# -*- coding: utf-8 -*-
"""Paired 25 °C Raoult vs HANNA screen for eligible Korean thinner MSDS documents.

Use ONLY structurally identified, fully property-covered MSDS records. Composition
ranges are sampled as feasible hypothetical scenarios, not measured formulations.
Model accuracy, emission rate and worker exposure are NOT established by this tool.

Place beside run_raoult_271.py and hanna_gamma_adapter.py in
research/korean_thinners_271/src/, then run from the repository root:
  $env:HANNA_HOME = 'D:\\fire2\\HANNA'
  python research/korean_thinners_271/src/run_hanna_271.py --samples 20 --batch-size 8
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
import os
import random
import statistics
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from hanna_gamma_adapter import _get_model
from run_raoult_271 import feasible_compositions, quantile
from composition_sampling import sample_bounded_simplex

BASE = Path(__file__).resolve().parents[1]
DATA = BASE / 'data'
PHYS = Path(__file__).resolve().parents[3] / 'thinner_vle_reference_25C.csv'


def read_csv(path: Path):
    with path.open('r', encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def fingerprint(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def positive_property(row, key, cas):
    try:
        value = float(row[key])
    except (KeyError, ValueError, TypeError):
        raise ValueError(f'{cas}: {key} 누락/숫자 오류')
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f'{cas}: {key}는 양의 유한값이어야 합니다.')
    return value


def predict_ln_gamma_batch(model, smiles, fractions, temperature, batch_size):
    """Batch-compute ln gamma. No silent ideal-solution fallback."""
    results = []
    for start in range(0, len(fractions), batch_size):
        part = fractions[start:start + batch_size]
        ln_gamma, _gE = model.predict(smiles, part, temperature)
        ln_gamma = np.asarray(ln_gamma, dtype=float)
        if ln_gamma.shape != (len(part), len(smiles)):
            raise ValueError(f'HANNA prediction shape invalid: {ln_gamma.shape}')
        if not np.isfinite(ln_gamma).all():
            raise ValueError('HANNA returned NaN/Inf ln(gamma): fallback 금지')
        results.append(ln_gamma)
    return np.concatenate(results, axis=0)


def main():
    ap = argparse.ArgumentParser(description='Real MSDS eligible 25°C Raoult vs HANNA paired comparison')
    ap.add_argument('--samples', type=int, default=20, help='Feasible composition draws per eligible PDF (initial smoke: 20)')
    ap.add_argument('--batch-size', type=int, default=8, help='HANNA inference batch size (CPU: 8 or 16)')
    ap.add_argument('--seed', type=int, default=20261023)
    ap.add_argument('--sampler', choices=['legacy', 'gibbs'], default='legacy',
                    help='legacy conditional rejection or bounded-simplex Gibbs (exploratory)')
    ap.add_argument('--burn-in', type=int, default=300, help='Gibbs warmup pair updates')
    ap.add_argument('--thin', type=int, default=12, help='Pair updates per retained Gibbs draw')
    ap.add_argument('--output', type=Path, default=BASE / 'output')
    ap.add_argument('--all-eligible', action='store_true', help='Scan all 46 molecular candidates; missing properties excluded')
    args = ap.parse_args()
    if not 1 <= args.samples <= 100000:
        ap.error('--samples must be between 1 and 100000')
    if not 1 <= args.batch_size <= 1000:
        ap.error('--batch-size must be between 1 and 1000')
    if args.burn_in < 0 or args.thin < 1:
        ap.error('--burn-in must be >=0 and --thin >=1')

    classify_path = DATA / '제품271건_용도재분류.csv'
    components_path = DATA / '성분별함량_1462.csv'
    classifications = read_csv(classify_path)
    components_all = read_csv(components_path)
    props = {r['CAS'].strip(): r for r in read_csv(PHYS) if r.get('CAS', '').strip()}

    by_product = collections.defaultdict(list)
    for c in components_all:
        by_product[c['제품ID']].append(c)
    eligible = [r for r in classifications if r.get('분자모델 1차 후보') == '예']
    if not args.all_eligible:
        eligible = [r for r in eligible if r.get('기존물성 라울예비 대상') == '예']

    scenario_fields = ['제품ID', '제품명', '제품군', '시나리오ID', '총압력_Raoult_kPa', '총압력_HANNA_kPa', 'HANNA_대_Raoult_비', 'HANNA_압력차_kPa', '모델_상태']
    component_fields = ['제품ID', '시나리오ID', 'CAS', '화학물질명', '질량백분율', '액상몰분율', 'HANNA_ln_gamma', 'HANNA_gamma', 'Raoult_분압_kPa', 'HANNA_분압_kPa', 'Raoult_기상몰분율', 'HANNA_기상몰분율']
    summary_fields = ['제품ID', '제품명', '제품군', '시나리오수', 'Raoult_총압력_중앙값_kPa', 'HANNA_총압력_중앙값_kPa', 'HANNA_총압력_5백분위_kPa', 'HANNA_총압력_95백분위_kPa', '압력비_중앙값', '압력차_중앙값_kPa']
    excluded_fields = ['제품ID', '제품명', '제외사유']

    scenarios, details, summaries, excluded = [], [], [], []
    rng = random.Random(args.seed)
    model = None  # load only when at least one fully eligible case passes all checks

    for r in eligible:
        pid = r['제품ID']
        name = r.get('제품명(한글)', '')
        group = r.get('1차 용도분류', '')
        comps = by_product.get(pid, [])
        if len(comps) < 2:
            excluded.append({'제품ID': pid, '제품명': name, '제외사유': '성분 2개 미만/누락'})
            continue
        cas_list = [str(c.get('CAS 번호', '')).strip() for c in comps]
        if not all(cas_list) or len(set(cas_list)) != len(cas_list):
            excluded.append({'제품ID': pid, '제품명': name, '제외사유': 'CAS 누락 또는 중복; 수동 검토 필요'})
            continue
        missing = sorted(set(cas_list) - set(props))
        if missing:
            excluded.append({'제품ID': pid, '제품명': name, '제외사유': 'MW/VP/SMILES 물성표 없음: ' + ', '.join(missing)})
            continue
        try:
            smiles = []
            mw = []
            vp = []
            for cas in cas_list:
                pr = props[cas]
                smi = str(pr.get('SMILES', '')).strip()
                if not smi:
                    raise ValueError(f'{cas}: 확인된 순물질 SMILES 없음')
                smiles.append(smi)
                mw.append(positive_property(pr, 'MW_g_mol', cas))
                vp.append(positive_property(pr, 'VP_25C_kPa', cas))
            if args.sampler == 'legacy':
                weights = feasible_compositions(comps, args.samples, rng)
            else:
                weights = sample_bounded_simplex(comps, args.samples, rng, burn_in=args.burn_in, thin=args.thin)
        except (ValueError, RuntimeError, KeyError) as exc:
            excluded.append({'제품ID': pid, '제품명': name, '제외사유': str(exc)})
            continue

        weights = np.asarray(weights, dtype=float)
        moles = weights / np.asarray(mw)[None, :]
        x = moles / moles.sum(axis=1, keepdims=True)
        if not np.isfinite(x).all() or (x < -1e-10).any() or not np.allclose(x.sum(axis=1), 1, atol=1e-7, rtol=0):
            raise ValueError(f'{pid}: 생성된 몰분율 오류')
        x = np.maximum(x, 0)
        x = x / x.sum(axis=1, keepdims=True)

        if model is None:
            print('Loading HANNA ensemble (once)...')
            model = _get_model()
        # If any HANNA inference fails, stop and DO NOT mislabel Raoult as HANNA.
        ln_gamma = predict_ln_gamma_batch(model, smiles, x, 298.15, args.batch_size)
        with np.errstate(over='ignore', under='ignore', invalid='ignore'):
            gamma = np.exp(ln_gamma)
        if not np.isfinite(gamma).all() or (gamma <= 0).any():
            raise ValueError(f'{pid}: invalid HANNA gamma (no fallback)')

        raoult_partial = x * np.asarray(vp)[None, :]
        hanna_partial = raoult_partial * gamma
        r_total = raoult_partial.sum(axis=1)
        h_total = hanna_partial.sum(axis=1)
        if (r_total <= 0).any() or (h_total <= 0).any() or not np.isfinite(h_total).all():
            raise ValueError(f'{pid}: invalid total vapor pressure')
        ratio = h_total / r_total
        diff = h_total - r_total

        for i in range(len(weights)):
            scenarios.append({
                '제품ID': pid, '제품명': name, '제품군': group,
                '시나리오ID': i + 1,
                '총압력_Raoult_kPa': round(float(r_total[i]), 8),
                '총압력_HANNA_kPa': round(float(h_total[i]), 8),
                'HANNA_대_Raoult_비': round(float(ratio[i]), 8),
                'HANNA_압력차_kPa': round(float(diff[i]), 8),
                '모델_상태': 'HANNA_2026_실제모델_예측',
            })
            for j, c in enumerate(comps):
                details.append({
                    '제품ID': pid, '시나리오ID': i + 1, 'CAS': cas_list[j],
                    '화학물질명': c.get('화학물질명(한글)', ''),
                    '질량백분율': round(float(weights[i, j]), 8),
                    '액상몰분율': round(float(x[i, j]), 9),
                    'HANNA_ln_gamma': round(float(ln_gamma[i, j]), 9),
                    'HANNA_gamma': round(float(gamma[i, j]), 9),
                    'Raoult_분압_kPa': round(float(raoult_partial[i, j]), 8),
                    'HANNA_분압_kPa': round(float(hanna_partial[i, j]), 8),
                    'Raoult_기상몰분율': round(float(raoult_partial[i, j] / r_total[i]), 9),
                    'HANNA_기상몰분율': round(float(hanna_partial[i, j] / h_total[i]), 9),
                })
        summaries.append({
            '제품ID': pid, '제품명': name, '제품군': group,
            '시나리오수': len(weights),
            'Raoult_총압력_중앙값_kPa': round(float(statistics.median(r_total)), 8),
            'HANNA_총압력_중앙값_kPa': round(float(statistics.median(h_total)), 8),
            'HANNA_총압력_5백분위_kPa': round(float(quantile(h_total, 0.05)), 8),
            'HANNA_총압력_95백분위_kPa': round(float(quantile(h_total, 0.95)), 8),
            '압력비_중앙값': round(float(statistics.median(ratio)), 8),
            '압력차_중앙값_kPa': round(float(statistics.median(diff)), 8),
        })
        print(f'{pid}: {len(weights)} matched scenarios finished ({len(cas_list)} components)')

    out = args.output
    write_csv(out / 'hanna_라울_제품별_요약.csv', summaries, summary_fields)
    write_csv(out / 'hanna_라울_시나리오별_비교.csv', scenarios, scenario_fields)
    write_csv(out / 'hanna_라울_성분별_비교.csv', details, component_fields)
    write_csv(out / 'hanna_라울_계산제외.csv', excluded, excluded_fields)
    meta = {
        'script': 'run_hanna_271.py', 'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'model': '2026 multicomponent HANNA local ensemble',
        'model_home': os.environ.get('HANNA_HOME', '<adapter default>'),
        'temperature_K': 298.15, 'composition': 'hypothetical valid MSDS interval scenarios',
        'seed': args.seed, 'sampler': args.sampler, 'burn_in': args.burn_in, 'thin': args.thin,
        'samples_requested_per_document': args.samples,
        'batch_size': args.batch_size, 'all_eligible': args.all_eligible,
        'input_sha256': {
            'classification': fingerprint(classify_path),
            'components': fingerprint(components_path),
            'property_table': fingerprint(PHYS),
        },
        'documents_considered': len(eligible), 'documents_calculated': len(summaries),
        'documents_excluded': len(excluded),
        'limitations': [
            'Preliminary equilibrium screening only; not measured emission rate or workplace exposure.',
            'MSDS concentration intervals do not identify actual manufacturer formulations.',
            'Selected 5th/95th percentiles are distribution of generated scenarios, not statistical confidence intervals.',
            'Gibbs chain draws are autocorrelated; uniformity and mixing have not been proven for a finite run.',
            'Existing property table requires independent temperature, units and source verification.',
            'Unresolved complex substances and isomer mixtures are not replaced with arbitrary SMILES.',
            'HANNA predictive accuracy for this thinner application has not been externally validated.',
        ],
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / 'hanna_라울_실행정보.json').write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'완료: 시도 문서 {len(eligible)}, 계산 {len(summaries)}, 제외 {len(excluded)}, 시나리오 {len(scenarios)}')
    print('출력 폴더:', out.resolve())
    print('주의: 실측 증발률/작업장 노출/TWA 초과/AI 예측정확도 검증 결과가 아닙니다.')

if __name__ == '__main__':
    main()