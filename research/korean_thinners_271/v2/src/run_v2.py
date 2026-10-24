# -*- coding: utf-8 -*-
"""v2 full-population run: 270 thinner MSDS x {Raoult, UNIFAC-Dortmund, HANNA} x {25, 35 C}.

What changed versus the v1 pipeline (see v2/README.md):
  * Isomer mixtures and UVCB solvents are expanded into documented surrogate components
    (v2/data/surrogates.csv) instead of blocking the whole product.
  * Pure-component properties and OELs are taken from sourced tables built by
    build_properties.py and build_oel.py (every value carries its source).
  * OEL coverage gap is computed explicitly (Tier A: Korean official; Tier B: + labelled
    supplementary values) instead of silently dropping species without a Korean TWA.
  * Modified UNIFAC (Dortmund) is added as the established non-ideal baseline.

Per product, N composition scenarios are sampled uniformly from the MSDS section-3
interval polytope (pairwise Gibbs, existing composition_sampling.py) and, inside each
UVCB, from its surrogate interval polytope. A deterministic 'central' scenario is added.
Outputs are equilibrium head-space *potentials*, not measured exposure.

Usage:
    python v2/src/run_v2.py --samples 100 --seed 20261024 --hanna-home /path/to/HANNA
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import random
import sys
import time
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np

warnings.filterwarnings('ignore')

V2 = Path(__file__).resolve().parents[1]
BASE = V2.parent
sys.path.insert(0, str(BASE / 'src'))
from composition_sampling import sample_bounded_simplex  # noqa: E402

DOCS = BASE / 'data' / 'source_rebuild_20261023' / 'MSDS_원본대조_270.csv'
COMP = BASE / 'data' / 'source_rebuild_20261023' / 'MSDS_구성성분_1454.csv'
USE = BASE / 'data' / '제품271건_용도재분류.csv'
PROPS = V2 / 'data' / 'pure_properties.csv'
OEL = V2 / 'data' / 'oel_assignment.csv'
SURR = V2 / 'data' / 'surrogates.csv'
RES = V2 / 'results'

P_ATM = 101.325
TEMPS = {'25C': 298.15, '35C': 308.15}
MODELS = ('Raoult', 'UNIFAC', 'HANNA')
NONVOLATILE = {'26588-80-7'}
STODDARD = '8052-41-3'
STODDARD_OEL = 100.0


def read(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def write(path, rows, cols=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = cols or list(rows[0])
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


def central(lo, hi):
    """Deterministic feasible point: lower bounds + remainder shared in proportion to range."""
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    rem, width = 100 - lo.sum(), (hi - lo).sum()
    return lo + (hi - lo) * (rem / width if width > 0 else 0)


def feasible_intervals(lo, hi):
    """Return intervals admitting a 100 % total; proportional repair for inconsistent MSDS."""
    lo, hi = list(map(float, lo)), list(map(float, hi))
    if sum(lo) <= 100 + 1e-9 and sum(hi) >= 100 - 1e-9:
        return lo, hi, ''
    k = 100.0 / ((sum(lo) + sum(hi)) / 2)
    lo, hi = [min(100, v * k) for v in lo], [min(100, v * k) for v in hi]
    return lo, hi, f'MSDS 범위합계 비정합(하한합 {sum(lo)/k:.1f}%, 상한합 {sum(hi)/k:.1f}%) -> 비례보정 k={k:.3f}'


def sample(lo, hi, n, rng):
    comps = [{'함량 하한(%)': a, '함량 상한(%)': b} for a, b in zip(lo, hi)]
    if len(lo) == 1:
        return np.full((n, 1), 100.0)
    return np.asarray(sample_bounded_simplex(comps, n, rng), float)


class Unifac:
    """Modified UNIFAC (Dortmund) via `thermo`, with explicit missing-parameter detection."""

    def __init__(self):
        from rdkit import Chem
        from thermo import unifac
        unifac.load_group_assignments_DDBST()
        self.Chem, self.u = Chem, unifac
        self.assign = unifac.DDBST_MODIFIED_UNIFAC_assignments

    def groups(self, smiles):
        key = self.Chem.MolToInchiKey(self.Chem.MolFromSmiles(smiles))
        return self.assign.get(key)

    def check(self, chemgroups):
        mains = sorted({self.u.DOUFSG[g].main_group_id for cg in chemgroups for g in cg})
        ip = self.u.DOUFIP2016
        missing = [(a, b) for i, a in enumerate(mains) for b in mains[i + 1:]
                   if not (a in ip and b in ip[a] and b in ip and a in ip[b])]
        return missing

    def model(self, chemgroups, T, xs):
        return self.u.UNIFAC.from_subgroups(T=T, xs=list(xs), chemgroups=chemgroups, version=1,
                                           interaction_data=self.u.DOUFIP2016,
                                           subgroups=self.u.DOUFSG)


class Hanna:
    def __init__(self, home):
        home = Path(home).resolve()
        if not (home / 'utils' / 'HANNA_predictor.py').is_file():
            raise FileNotFoundError(f'HANNA not found at {home}')
        sys.path.insert(0, str(home))
        prev = Path.cwd()
        os.chdir(home)
        try:
            from utils.HANNA_predictor import HANNA_Predictor
            self.m = HANNA_Predictor()
        finally:
            os.chdir(prev)
        self.home = home

    def gammas(self, smiles, x, T):
        prev = Path.cwd()
        os.chdir(self.home)
        try:
            lng, _ = self.m.predict(list(smiles), np.asarray(x, float), float(T))
        finally:
            os.chdir(prev)
        g = np.exp(np.asarray(lng, float))
        if g.shape != np.asarray(x).shape or not np.isfinite(g).all():
            raise ValueError('invalid HANNA output')
        return g


def load_inputs():
    props = {r['CAS']: r for r in read(PROPS)}
    oel = {r['CAS']: r for r in read(OEL)}
    surr = defaultdict(list)
    for r in read(SURR):
        surr[r['UVCB_CAS']].append(r)
    use = {r['제품ID']: r for r in read(USE)}
    parts = defaultdict(list)
    for r in read(COMP):
        parts[r['문서키']].append(r)
    return read(DOCS), parts, props, oel, surr, use


def expand(doc_parts, surr, n, rng):
    """Return (species list, origin list, mass array [n+1, n_species]) with row 0 = central."""
    lo, hi, note = feasible_intervals([p['함량하한_pct'] for p in doc_parts],
                                      [p['함량상한_pct'] for p in doc_parts])
    top = np.vstack([central(lo, hi), sample(lo, hi, n, rng)])  # product-level wt%
    cols, blocks = [], []
    for j, p in enumerate(doc_parts):
        cas = p['CAS']
        w = top[:, j:j + 1]
        if cas in NONVOLATILE:
            continue
        if cas in surr:
            sub = surr[cas]
            slo = [float(s['wt_lo']) for s in sub]
            shi = [float(s['wt_hi']) for s in sub]
            split = np.vstack([central(slo, shi), sample(slo, shi, n, rng)]) / 100.0
            for k, s in enumerate(sub):
                cols.append((s['component_CAS'], cas))
                blocks.append(w[:, 0] * split[:, k])
        else:
            cols.append((cas, cas))
            blocks.append(w[:, 0])
    mass = np.column_stack(blocks)
    return cols, mass, note


def run(args):
    t0 = time.time()
    docs, parts, props, oel, surr, use = load_inputs()
    uni = Unifac()
    hanna = Hanna(args.hanna_home) if args.hanna_home else None
    rng = random.Random(args.seed)
    scen_rows, species_rows, status_rows = [], [], []
    if args.limit:
        docs = docs[:args.limit]
    for di, d in enumerate(docs):
        key = d['문서키']
        cols, mass, note = expand(parts[key], surr, args.samples, rng)
        species = sorted({c for c, _ in cols})
        sidx = {c: i for i, c in enumerate(species)}
        # aggregate mass by species; remember Stoddard-origin share per species
        m_sp = np.zeros((mass.shape[0], len(species)))
        m_stod = np.zeros_like(m_sp)
        for (c, origin), col in zip(cols, mass.T):
            m_sp[:, sidx[c]] += col
            if origin == STODDARD:
                m_stod[:, sidx[c]] += col
        mw = np.array([float(props[c]['MW_g_mol']) for c in species])
        mol = m_sp / mw
        x = mol / mol.sum(axis=1, keepdims=True)
        stod_frac = np.divide(m_stod, m_sp, out=np.zeros_like(m_sp), where=m_sp > 0)
        smiles = [props[c]['SMILES'] for c in species]
        chemgroups = [uni.groups(s) for s in smiles]
        uni_issue = ''
        if any(g is None for g in chemgroups):
            uni_issue = 'UNIFAC 그룹할당 없음: ' + ','.join(c for c, g in zip(species, chemgroups) if g is None)
        else:
            miss = uni.check(chemgroups)
            if miss:
                uni_issue = 'UNIFAC 상호작용 파라미터 결측(주그룹쌍): ' + ';'.join(f'{a}-{b}' for a, b in miss)
        for tl, T in TEMPS.items():
            psat = np.array([float(props[c][f'Psat{tl[:2]}_kPa']) for c in species])
            gam = {'Raoult': np.ones_like(x)}
            if len(species) == 1:
                gam['UNIFAC'] = np.ones_like(x)
                gam['HANNA'] = np.ones_like(x)
            else:
                if not uni_issue:
                    um = uni.model(chemgroups, T, x[0])
                    gam['UNIFAC'] = np.array([um.to_T_xs(T, list(xi)).gammas() for xi in x])
                if hanna is not None:
                    gam['HANNA'] = hanna.gammas(smiles, x, T)
            for model, g in gam.items():
                p = x * g * psat
                ptot = p.sum(axis=1)
                out = {}
                for tier in ('A', 'B'):
                    oel_v = np.array([float(oel[c][f'Tier{tier}_TWA_ppm']) if oel[c][f'Tier{tier}_TWA_ppm'] else np.nan
                                      for c in species])
                    # reciprocal OEL with Stoddard-origin override
                    inv = np.where(np.isnan(oel_v), 0.0, 1.0 / np.where(np.isnan(oel_v), 1, oel_v))
                    inv_mat = (1 - stod_frac) * inv + stod_frac * (1 / STODDARD_OEL)
                    covered = (1 - stod_frac) * (~np.isnan(oel_v)) + stod_frac
                    terms = 1e6 * p / P_ATM * inv_mat
                    out[tier] = (terms.sum(axis=1), (p * covered).sum(axis=1) / ptot, terms)
                for i in range(x.shape[0]):
                    tA = out['A'][2][i]
                    scen_rows.append({
                        '문서키': key, '기존제품ID': d['기존제품ID'], 'scenario': 'central' if i == 0 else i,
                        'T': tl, 'model': model, 'Ptot_kPa': round(float(ptot[i]), 6),
                        'Csat_total_ppm': round(1e6 * float(ptot[i]) / P_ATM, 1),
                        'I_mix_TierA': round(float(out['A'][0][i]), 4),
                        'I_mix_TierB': round(float(out['B'][0][i]), 4),
                        'OELcov_TierA': round(float(out['A'][1][i]), 4),
                        'OELcov_TierB': round(float(out['B'][1][i]), 4),
                        'Hmax_CAS_TierA': species[int(np.argmax(tA))] if tA.max() > 0 else ''})
                    if i == 0:
                        for k, c in enumerate(species):
                            species_rows.append({
                                '문서키': key, 'T': tl, 'model': model, 'CAS': c, 'name': props[c]['name'],
                                'x_liq': round(float(x[0, k]), 6), 'gamma': round(float(g[0, k]), 5),
                                'p_kPa': round(float(p[0, k]), 7),
                                'y_vap': round(float(p[0, k] / ptot[0]), 6),
                                'term_TierA': round(float(out['A'][2][0, k]), 5),
                                'term_TierB': round(float(out['B'][2][0, k]), 5)})
        u = use.get(d['기존제품ID'], {})
        status_rows.append({
            '문서키': key, '기존제품ID': d['기존제품ID'], '제조사': d['제조사'], '제품명': d['제품명_파일명기반'],
            '용도군': u.get('1차 용도분류', ''), '세부분류': u.get('세부분류', ''),
            '계절변형': u.get('계절·건조 변형', ''), 'MSDS성분수': len(parts[key]),
            '모델화학종수': len(species), 'UVCB대리물질사용': any(c[1] != c[0] for c in cols),
            '고분자제외': any(p['CAS'] in NONVOLATILE for p in parts[key]),
            '조성비고': note, 'UNIFAC비고': uni_issue,
            'HANNA': 'computed' if hanna else 'not run',
            '조성지문': '|'.join(sorted(f"{p['CAS']}:{p['함량하한_pct']}-{p['함량상한_pct']}" for p in parts[key]))})
        if (di + 1) % 20 == 0:
            print(f'{di + 1}/{len(docs)} docs, {time.time() - t0:.0f}s', flush=True)
    RES.mkdir(parents=True, exist_ok=True)
    write(RES / 'scenarios_270.csv.gz', scen_rows)
    write(RES / 'species_central_270.csv.gz', species_rows)
    write(RES / 'products_270.csv', status_rows)
    meta = {'samples_per_product': args.samples, 'seed': args.seed, 'temperatures_K': TEMPS,
            'models': [m for m in MODELS if m != 'HANNA' or hanna], 'P_atm_kPa': P_ATM,
            'hanna_home': str(args.hanna_home) if args.hanna_home else None,
            'runtime_s': round(time.time() - t0, 1), 'products': len(docs)}
    (RES / 'run_meta.json').write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(meta, ensure_ascii=False))


def cli(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--samples', type=int, default=100)
    ap.add_argument('--seed', type=int, default=20261024)
    ap.add_argument('--hanna-home', default=os.environ.get('HANNA_HOME'))
    ap.add_argument('--limit', type=int, default=0, help='debug: first N documents only')
    run(ap.parse_args(argv))


if __name__ == '__main__':
    cli()
