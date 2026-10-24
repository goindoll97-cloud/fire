# -*- coding: utf-8 -*-
"""Pure-component property table for every molecular species the v2 model uses.

Species = every non-UVCB CAS on the 270 MSDS + every surrogate component in
v2/data/surrogates.csv. For each species we record SMILES, MW and the liquid
vapour pressure at 25 and 35 C with the correlation that produced it.

Vapour-pressure source priority (all from the open `chemicals`/`thermo` library,
Bell et al., which tabulates the original literature correlations):
    WAGNER_MCGARRY > WAGNER_POLING > ANTOINE_POLING > DIPPR_PERRY_8E >
    ANTOINE_EXTENDED_POLING > ANTOINE_WEBBOOK > LANDOLT > HEOS_FIT
If no correlation exists, a literature anchor point (LITERATURE_ANCHORS) or the
normal boiling point is combined with Clausius-Clapeyron using a tabulated
enthalpy of vaporisation (CC_* methods, flagged as estimates).
Every species also gets a QC column: relative spread between all available
correlations at 25 C.
"""
from __future__ import annotations

import csv
import math
import warnings
from pathlib import Path

warnings.filterwarnings('ignore')

from chemicals import Tb as tb_lookup
from chemicals.identifiers import search_chemical
from rdkit import Chem
from thermo import EnthalpyVaporization, VaporPressure

V2 = Path(__file__).resolve().parents[1]
BASE = V2.parent
COMP = BASE / 'data' / 'source_rebuild_20261023' / 'MSDS_구성성분_1454.csv'
SURR = V2 / 'data' / 'surrogates.csv'
OUT = V2 / 'data' / 'pure_properties.csv'

R = 8.314462618
T25, T35 = 298.15, 308.15
PRIORITY = ['WAGNER_MCGARRY', 'WAGNER_POLING', 'ANTOINE_POLING', 'DIPPR_PERRY_8E',
            'ANTOINE_EXTENDED_POLING', 'ANTOINE_WEBBOOK', 'LANDOLT', 'HEOS_FIT']

# UVCB / isomer-mixture CAS handled through surrogates (never as one molecule)
UVCB = {'1330-20-7', '25551-13-7', '64742-95-6', '64742-94-5', '64742-49-0', '31394-54-4',
        '64742-82-1', '64741-92-0', '64742-47-8', '8052-41-3'}
# Non-volatile constituents: kept in the liquid mass balance only as inert (see run_v2)
NONVOLATILE = {'26588-80-7': 'acrylic copolymer (polymer, non-volatile)'}

# Literature anchor points used only when no correlation is available.
# (P in kPa, T in K, citation)
LITERATURE_ANCHORS = {
    '108-65-6': (0.517, 298.15, 'MIDSIS-TROCS (REMPEC) PGMEA physical data: 517 Pa at 25 C'),
    '97-64-3': (0.16, 293.15, 'Merck/Sigma-Aldrich ethyl lactate product data: 1.6 hPa at 20 C'),
    '2110-78-3': (0.4608, 298.15, 'The Good Scents Company: 3.456 mmHg at 25 C (estimated value)'),
}
# Anchors that override a correlation found to be implausible in QC
# (the only LANDOLT correlation for dimethyl succinate extrapolates to 0.00037 kPa at 25 C,
#  ~50x below the supplier value).
OVERRIDE_ANCHORS = {
    '106-65-0': (0.02, 293.15, 'Merck/Sigma-Aldrich dimethyl succinate product data: 0.2 hPa at 20 C'),
}


def _hvap(cas, tb):
    h = EnthalpyVaporization(CASRN=cas)
    for m in ('DIPPR_PERRY_8E', 'CRC_HVAP_298', 'GHARAGHEIZI_HVAP_298', 'CRC_HVAP_TB', 'VDI_PPDS'):
        if m in h.all_methods:
            try:
                v = h.calculate(T25, m) if m not in ('CRC_HVAP_TB',) else h.calculate(tb, m)
                if v and v > 0:
                    return v, m
            except Exception:
                pass
    # Trouton-Kistiakowsky style estimate (Riedel-type constant 88 J/mol/K at Tb)
    return 88.0 * tb, 'TROUTON_88'


def species_list():
    cas = set()
    with COMP.open(encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            cas.add(r['CAS'])
    with SURR.open(encoding='utf-8') as f:
        surr = list(csv.DictReader(f))
    cas |= {r['component_CAS'] for r in surr}
    return sorted(c for c in cas if c not in UVCB and c not in NONVOLATILE)


def props(cas):
    m = search_chemical(cas)
    smi = m.smiles
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f'{cas}: invalid SMILES {smi}')
    smi = Chem.MolToSmiles(mol)
    vp = VaporPressure(CASRN=cas)
    avail = [x for x in PRIORITY if x in vp.all_methods]
    p25s = {}
    for meth in avail:
        try:
            v = vp.calculate(T25, meth)
            if v and v > 0 and math.isfinite(v):
                p25s[meth] = v / 1000
        except Exception:
            pass
    row = {'CAS': cas, 'name': m.common_name, 'SMILES': smi, 'MW_g_mol': round(m.MW, 4)}
    if p25s and cas not in OVERRIDE_ANCHORS:
        valid_both = [x for x in PRIORITY if x in p25s
                      and vp.test_method_validity(T25, x) and vp.test_method_validity(T35, x)]
        meth = valid_both[0] if valid_both else next(x for x in PRIORITY if x in p25s)
        p25 = p25s[meth]
        p35 = vp.calculate(T35, meth) / 1000
        valid = vp.test_method_validity(T25, meth)
        row.update(Psat25_kPa=p25, Psat35_kPa=p35, Psat_method=meth,
                   Psat_note='' if valid else 'extrapolated below correlation Tmin (sub-cooled liquid)',
                   n_correlations=len(p25s),
                   QC_spread_25C=(max(p25s.values()) / min(p25s.values()) - 1) if len(p25s) > 1 else '')
        return row
    tb = tb_lookup(cas)
    hv, hm = _hvap(cas, tb)
    anchors = {**LITERATURE_ANCHORS, **OVERRIDE_ANCHORS}
    if cas in anchors:
        p0, t0, cite = anchors[cas]
        meth = f'CC_ANCHOR+{hm}'
        note = f'estimate: Clausius-Clapeyron from {cite}'
    else:
        p0, t0, cite = 101.325, tb, f'normal boiling point {tb:.2f} K (chemicals.Tb)'
        meth = f'CC_TB+{hm}'
        note = f'estimate: Clausius-Clapeyron from {cite}'
    cc = lambda T: p0 * math.exp(-hv / R * (1 / T - 1 / t0))
    row.update(Psat25_kPa=cc(T25), Psat35_kPa=cc(T35), Psat_method=meth, Psat_note=note,
               n_correlations=len(p25s), QC_spread_25C='')
    return row


def main():
    rows = [props(c) for c in species_list()]
    cols = ['CAS', 'name', 'SMILES', 'MW_g_mol', 'Psat25_kPa', 'Psat35_kPa', 'Psat_method',
            'Psat_note', 'n_correlations', 'QC_spread_25C']
    with OUT.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    est = [r['CAS'] for r in rows if str(r['Psat_method']).startswith('CC_')]
    print(f'{len(rows)} species -> {OUT}; Clausius-Clapeyron estimates: {est}')


if __name__ == '__main__':
    main()
