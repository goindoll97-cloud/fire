# -*- coding: utf-8 -*-
"""화재 지표: 혼합물 인화점(Liaw 모델)과 포화증기 폭발하한 도달비(Le Chatelier)를 계산한다.

v2 노출 계산과 같은 입력(MSDS 범위 샘플, 대리물질, 활동도계수 Raoult/UNIFAC/HANNA)을 쓴다.
같은 시나리오를 재현하려면 run_v2.py와 같은 --samples, --seed를 쓴다 (난수 호출 순서가 같다).

지표
  Tfp_Liaw_C   Liaw 등(1997) 혼합물 인화점: Σ x_i γ_i Psat_i(T) / Psat_i(Tfp_i) = 1 을 푸는 T.
               γ_i는 25°C 값을 고정해 쓴다 (인화점 대부분이 0–45°C이므로 근사; 중앙조성에서 T-일관 계산과 비교).
  F25, F35     포화증기 화재지수 F = (Ptot/P) / LEL_mix,  LEL_mix = 1 / Σ(y_i / LFL_i) (Le Chatelier).
               F>=1은 혼합 기상 농도가 계산상 폭발하한(LEL)에 도달함을 뜻한다. 폭발상한(UEL) 미검토로 폭발범위 내부 여부를 확정하지 않는다.
  Ptot, y_i    노출 계산과 동일.
인화점은 폐쇄식(밀폐) 값을 가정한 성분별 값이다. 개방식 값이 섞이면 성분 값이 높게 나올 수 있다.

사용: python run_flammability.py --samples 100 --seed 20261024 --hanna-home <경로> [--limit N]
"""
from __future__ import annotations

import argparse
import json
import math
import random
import time
import warnings
from pathlib import Path

import numpy as np

import run_v2 as rv
from thermo import VaporPressure

warnings.filterwarnings('ignore')

V2 = rv.V2
FLAM = V2 / 'data' / 'pure_flammability.csv'
T_GRID = np.arange(233.15, 403.16, 2.5)          # -40 ~ 130 °C
R = 8.314462618


def psat_grid(props):
    """종별 Psat(T) [kPa] 격자. 상관식이 유효범위를 벗어나면 가장 가까운 유효 구간의 ln P - 1/T 직선으로 외삽."""
    out = {}
    for cas, r in props.items():
        meth = r['Psat_method']
        p25, p35 = float(r['Psat25_kPa']), float(r['Psat35_kPa'])
        if meth.startswith('CC_'):
            slope = math.log(p35 / p25) / (1 / 308.15 - 1 / 298.15)      # = -Hvap/R
            vals = np.array([p25 * math.exp(slope * (1 / T - 1 / 298.15)) for T in T_GRID])
        else:
            vp = VaporPressure(CASRN=cas)
            vals = np.full(T_GRID.shape, np.nan)
            for k, T in enumerate(T_GRID):
                try:
                    if vp.test_method_validity(T, meth):
                        v = vp.calculate(T, meth)
                        if v and v > 0 and math.isfinite(v):
                            vals[k] = v / 1000
                except Exception:
                    pass
            ok = np.isfinite(vals)
            if ok.sum() < 3:
                slope = math.log(p35 / p25) / (1 / 308.15 - 1 / 298.15)
                vals = np.array([p25 * math.exp(slope * (1 / T - 1 / 298.15)) for T in T_GRID])
            elif not ok.all():
                inv = 1 / T_GRID
                lp = np.log(vals[ok])
                idx = np.where(ok)[0]
                for side, sel in (('lo', idx[:4]), ('hi', idx[-4:])):
                    a, b = np.polyfit(inv[sel], np.log(vals[sel]), 1)
                    miss = ~ok & ((T_GRID < T_GRID[idx[0]]) if side == 'lo' else (T_GRID > T_GRID[idx[-1]]))
                    vals[miss] = np.exp(a * inv[miss] + b)
        out[cas] = vals
    return out


def lnp_at(grid, T):
    """격자(ln P vs 1/T) 선형보간"""
    return np.interp(1 / T, (1 / T_GRID)[::-1], np.log(grid)[::-1])


def liaw_flash(a, ratio_grid):
    """a: [S, n] = x_i γ_i, ratio_grid: [n, G] = Psat_i(T)/Psat_i(Tfp_i). 합이 1이 되는 T(K) 반환 [S]."""
    f = a @ ratio_grid            # [S, G]
    tf = np.full(a.shape[0], np.nan)
    for s in range(a.shape[0]):
        row = f[s]
        k = np.where(row >= 1)[0]
        if k.size == 0:
            continue
        j = k[0]
        if j == 0:
            tf[s] = T_GRID[0]
            continue
        # ln f 가 1/T에 거의 선형이므로 ln-ln 보간
        l0, l1 = math.log(row[j - 1]), math.log(row[j])
        w = (0 - l0) / (l1 - l0) if l1 != l0 else 0
        tf[s] = T_GRID[j - 1] + w * (T_GRID[j] - T_GRID[j - 1])
    return tf


def run(args):
    t0 = time.time()
    docs, parts, props, _, surr, _ = rv.load_inputs()
    flam = {r['CAS']: r for r in rv.read(FLAM)}
    grids = psat_grid(props)
    uni = rv.Unifac()
    hanna = rv.Hanna(args.hanna_home) if args.hanna_home else None
    rng = random.Random(args.seed)
    rows = []
    if args.limit:
        docs = docs[:args.limit]
    for di, d in enumerate(docs):
        key = d['문서키']
        cols, mass, note = rv.expand(parts[key], surr, args.samples, rng)
        species = sorted({c for c, _ in cols})
        sidx = {c: i for i, c in enumerate(species)}
        m_sp = np.zeros((mass.shape[0], len(species)))
        for (c, _o), col in zip(cols, mass.T):
            m_sp[:, sidx[c]] += col
        mw = np.array([float(props[c]['MW_g_mol']) for c in species])
        mol = m_sp / mw
        x = mol / mol.sum(axis=1, keepdims=True)
        smiles = [props[c]['SMILES'] for c in species]
        chemgroups = [uni.groups(s) for s in smiles]
        uni_ok = all(g is not None for g in chemgroups) and not uni.check(chemgroups)
        fp = np.array([float(flam[c]['FP_C']) + (args.fp_offset if flam[c].get('FP_grade', 'A') in ('B', 'C') else 0.0) + 273.15
                       for c in species])
        lfl = np.array([float(flam[c]['LFL_pct']) / 100 for c in species])
        est_any = any(flam[c]['estimated'] == 'True' for c in species)
        # Psat_i(T)/Psat_i(Tfp_i)
        ratio = np.array([np.exp(np.log(grids[c]) - lnp_at(grids[c], fp[k])) for k, c in enumerate(species)])
        gam25 = {'Raoult': np.ones_like(x)}
        gam35 = {'Raoult': np.ones_like(x)}
        if len(species) == 1:
            for g in (gam25, gam35):
                g['UNIFAC'] = np.ones_like(x)
                g['HANNA'] = np.ones_like(x)
        else:
            if uni_ok:
                for tl, T, g in (('25C', 298.15, gam25), ('35C', 308.15, gam35)):
                    um = uni.model(chemgroups, T, x[0])
                    g['UNIFAC'] = np.array([um.to_T_xs(T, list(xi)).gammas() for xi in x])
            if hanna is not None:
                gam25['HANNA'] = hanna.gammas(smiles, x, 298.15)
                gam35['HANNA'] = hanna.gammas(smiles, x, 308.15)
        for model in gam25:
            a = x * gam25[model]
            tf = liaw_flash(a, ratio)
            res = {}
            for tl, T, g in (('25C', 298.15, gam25), ('35C', 308.15, gam35)):
                psat = np.array([float(props[c][f'Psat{tl[:2]}_kPa']) for c in species])
                p = x * g[model] * psat
                ptot = p.sum(axis=1)
                y = p / ptot[:, None]
                lel_mix = 1 / (y / lfl).sum(axis=1)
                res[tl] = (ptot, lel_mix, (ptot / rv.P_ATM) / lel_mix)
            for i in range(x.shape[0]):
                rows.append({
                    '문서키': key, '기존제품ID': d['기존제품ID'], 'scenario': 'central' if i == 0 else i,
                    'model': model, 'Tfp_Liaw_C': '' if np.isnan(tf[i]) else round(tf[i] - 273.15, 2),
                    'F25': round(float(res['25C'][2][i]), 5), 'F35': round(float(res['35C'][2][i]), 5),
                    'LELmix25_pct': round(100 * float(res['25C'][1][i]), 4),
                    'Csat25_pct': round(100 * float(res['25C'][0][i] / rv.P_ATM), 4),
                    'FP추정성분포함': est_any})
        if (di + 1) % 20 == 0:
            print(f'{di + 1}/{len(docs)} docs, {time.time() - t0:.0f}s', flush=True)
    rv.write(rv.RES / f'flammability_scenarios_270{args.out_suffix}.csv.gz', rows)
    (rv.RES / 'run_flammability_meta.json').write_text(json.dumps(
        {'samples_per_product': args.samples, 'seed': args.seed, 'T_grid_C': [-40, 130, 2.5],
         'gamma_T_for_flash_K': 298.15, 'runtime_s': round(time.time() - t0, 1), 'products': len(docs)},
        ensure_ascii=False, indent=2), encoding='utf-8')
    print('done', round(time.time() - t0), 's')


def cli(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--samples', type=int, default=100)
    ap.add_argument('--seed', type=int, default=20261024)
    ap.add_argument('--hanna-home', default=None)
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--fp-offset', type=float, default=0.0, help='민감도: 출처 등급 B·C 성분 인화점에 더할 값(°C)')
    ap.add_argument('--out-suffix', default='')
    run(ap.parse_args(argv))


if __name__ == '__main__':
    cli()
