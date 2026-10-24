# -*- coding: utf-8 -*-
"""Two-tier occupational exposure limit (OEL) assignment for every modelled species.

Tier A  = Korean MoEL notice <별표1> TWA only (parsed by parse_korean_oel.py).
          Exact CAS match, or an explicit isomer-group entry of the notice
          (xylene all isomers; trimethylbenzene mixed isomers; hexane other isomers).
Tier B  = Tier A + clearly labelled supplementary values:
          * foreign official limit (EU IOELV) where Korea has none, or
          * structural-analogue Korean value for surrogate components of UVCB solvents
            (reciprocal-calculation style; analogue named in the table).
Species with neither stay *uncovered*; their vapour share is reported as the OEL
coverage gap instead of being silently dropped.

UVCB substances that themselves carry an official Korean TWA (Stoddard solvent)
apply that TWA to all of their surrogate components (ORIGIN_OVERRIDES).
"""
from __future__ import annotations

import csv
from pathlib import Path

V2 = Path(__file__).resolve().parents[1]
OEL_TABLE = V2 / 'data' / 'korean_oel_table.csv'
PROPS = V2 / 'data' / 'pure_properties.csv'
OUT = V2 / 'data' / 'oel_assignment.csv'

# Isomer-group entries of the Korean notice that cover CAS numbers not printed in it.
ISOMER_GROUP = {
    '95-63-6': 603, '108-67-8': 603, '526-73-8': 603,   # 트리메틸 벤젠(혼합 이성체)
}

EU_2000_39 = 'EU IOELV, Commission Directive 2000/39/EC Annex (8 h TWA)'
TIER_B = {
    '108-65-6': (50.0, EU_2000_39 + ': 2-methoxy-1-methylethyl acetate 50 ppm'),
    '620-14-4': (25.0, 'Analogue: Korean TWA of trimethylbenzene (mixed isomers) 25 ppm (C9 alkylbenzene)'),
    '622-96-8': (25.0, 'Analogue: Korean TWA of trimethylbenzene (mixed isomers) 25 ppm (C9 alkylbenzene)'),
    '611-14-3': (25.0, 'Analogue: Korean TWA of trimethylbenzene (mixed isomers) 25 ppm (C9 alkylbenzene)'),
    '103-65-1': (25.0, 'Analogue: Korean TWA of trimethylbenzene (mixed isomers) 25 ppm (C9 alkylbenzene)'),
    '104-51-8': (25.0, 'Analogue: Korean TWA of trimethylbenzene (mixed isomers) 25 ppm (C10 alkylbenzene)'),
    '527-53-7': (25.0, 'Analogue: Korean TWA of trimethylbenzene (mixed isomers) 25 ppm (C10 alkylbenzene)'),
    '90-12-0': (10.0, 'Analogue: Korean TWA of naphthalene 10 ppm'),
    '591-76-4': (400.0, 'Analogue: Korean TWA of heptane 400 ppm (C7 isoparaffin)'),
    '589-34-4': (400.0, 'Analogue: Korean TWA of heptane 400 ppm (C7 isoparaffin)'),
    '108-08-7': (400.0, 'Analogue: Korean TWA of heptane 400 ppm (C7 isoparaffin)'),
    '124-18-5': (100.0, 'Analogue: Korean TWA of Stoddard solvent 100 ppm (C9-C13 paraffin)'),
    '1120-21-4': (100.0, 'Analogue: Korean TWA of Stoddard solvent 100 ppm (C9-C13 paraffin)'),
    '112-40-3': (100.0, 'Analogue: Korean TWA of Stoddard solvent 100 ppm (C9-C13 paraffin)'),
    '629-50-5': (100.0, 'Analogue: Korean TWA of Stoddard solvent 100 ppm (C9-C13 paraffin)'),
}
# UVCB carrying its own official Korean TWA: applies to all its surrogate components.
ORIGIN_OVERRIDES = {'8052-41-3': (100.0, 'Korean TWA 별표1 No.326 스토다드 용제 (applied to all surrogate components)')}


def korean_index():
    idx, rows = {}, {}
    with OEL_TABLE.open(encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            rows[int(r['일련번호'])] = r
            for c in filter(None, r['CAS'].split('|')):
                idx.setdefault(c, int(r['일련번호']))
    return idx, rows


def build():
    idx, rows = korean_index()
    out = []
    with PROPS.open(encoding='utf-8-sig') as f:
        species = list(csv.DictReader(f))
    for s in species:
        cas = s['CAS']
        no = idx.get(cas) or ISOMER_GROUP.get(cas)
        a_ppm, a_src = '', ''
        if no:
            r = rows[no]
            if r['TWA_ppm']:
                a_ppm = float(r['TWA_ppm'])
                a_src = f"Korean MoEL notice 별표1 No.{no} {r['국문']} ({r['영문']})"
                if cas in ISOMER_GROUP:
                    a_src += ' - isomer-group entry'
        if a_ppm != '':
            b_ppm, b_src = a_ppm, 'Tier A'
        elif cas in TIER_B:
            b_ppm, b_src = TIER_B[cas]
        else:
            b_ppm, b_src = '', ''
        out.append({'CAS': cas, 'name': s['name'], 'TierA_TWA_ppm': a_ppm, 'TierA_source': a_src,
                    'TierB_TWA_ppm': b_ppm, 'TierB_source': b_src})
    return out


def main():
    out = build()
    with OUT.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader(); w.writerows(out)
    na = sum(1 for r in out if r['TierA_TWA_ppm'] == '')
    nb = sum(1 for r in out if r['TierB_TWA_ppm'] == '')
    print(f'{len(out)} species; no Tier A: {na}; no Tier B: {nb} -> {OUT}')


if __name__ == '__main__':
    main()
