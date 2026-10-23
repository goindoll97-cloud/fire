# -*- coding: utf-8 -*-
"""Reproducible sensitivity study for the 271-MSDS HANNA/Raoult research.

Only scenario-conditioned equilibrium predictions are assessed here. A paired
model comparison does not measure model accuracy, exposure or substitution safety.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
RESEARCH_ROOT = SCRIPT_DIR.parent
REPO_ROOT = RESEARCH_ROOT.parents[1]
CLASSIFICATION = RESEARCH_ROOT / 'data' / '제품271건_용도재분류.csv'
COMPONENTS = RESEARCH_ROOT / 'data' / '성분별함량_1462.csv'
PROPERTY = REPO_ROOT / 'thinner_vle_reference_25C.csv'
SCENARIO_FILE = 'hanna_라울_시나리오별_비교.csv'
COMPONENT_FILE = 'hanna_라울_성분별_비교.csv'
METADATA_FILE = 'hanna_라울_실행정보.json'


def read_csv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, cols):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def qtile(values, q):
    vals = sorted(values)
    if not vals:
        raise ValueError('Empty quantile series')
    t = (len(vals) - 1) * q
    i = int(t)
    return vals[i] + (vals[min(i+1, len(vals)-1)] - vals[i]) * (t-i)


def checked_number(v, name, *, strict_positive=False):
    try:
        x = float(v)
    except (TypeError, ValueError):
        raise ValueError(f'{name}: non-numeric value {v!r}')
    if not math.isfinite(x) or (strict_positive and x <= 0):
        raise ValueError(f'{name}: invalid number {v!r}')
    return x


def expected_sha():
    return {'classification': sha256(CLASSIFICATION), 'components': sha256(COMPONENTS),
            'property_table': sha256(PROPERTY)}


def run_one(out_dir, *, sampler, seed, samples, batch_size, burn_in, thin):
    cmd = [sys.executable, str(SCRIPT_DIR / 'run_hanna_271.py'),
           '--sampler', sampler, '--seed', str(seed), '--samples', str(samples),
           '--batch-size', str(batch_size), '--burn-in', str(burn_in),
           '--thin', str(thin), '--output', str(out_dir)]
    required = [out_dir / v for v in (SCENARIO_FILE, COMPONENT_FILE, METADATA_FILE)]
    if out_dir.exists() and any(out_dir.iterdir()):
        if not all(x.exists() for x in required):
            raise RuntimeError(f'Partial run directory (do not overwrite): {out_dir}')
        d = json.loads((out_dir / METADATA_FILE).read_text(encoding='utf-8-sig'))
        identity = (d.get('sampler') == sampler and d.get('seed') == seed and
                    d.get('samples_requested_per_document') == samples and
                    d.get('input_sha256') == expected_sha())
        if not identity:
            raise RuntimeError(f'Preexisting run has mismatched settings or source: {out_dir}')
        print(f'[SKIP] verified existing run metadata: {sampler}/{seed}', flush=True)
        return
    print(f'[START] {sampler} seed={seed} n={samples}', flush=True)
    subprocess.run(cmd, cwd=REPO_ROOT, check=True)
    if not all(x.exists() for x in required):
        raise RuntimeError(f'Missing results from {cmd}')


def validate_run(run_dir, meta, known_components, props, expected, top_tie_margin):
    if meta.get('input_sha256') != expected:
        raise ValueError(f'{run_dir}: underlying input SHA differs from current repository')
    if meta.get('model') != '2026 multicomponent HANNA local ensemble':
        raise ValueError(f'{run_dir}: invalid model identifier')
    if abs(checked_number(meta.get('temperature_K'), 'temperature_K') - 298.15) > 1e-6:
        raise ValueError(f'{run_dir}: only 25 C physical properties are supported')
    if meta.get('sampler') not in {'legacy', 'gibbs'}:
        raise ValueError(f'{run_dir}: missing valid sampler metadata')
    if not isinstance(meta.get('seed'), int):
        raise ValueError(f'{run_dir}: seed is missing or invalid')
    scenarios = read_csv(run_dir / SCENARIO_FILE)
    parts = read_csv(run_dir / COMPONENT_FILE)
    if len(scenarios) != meta.get('documents_calculated', -1) * meta.get('samples_requested_per_document', -2):
        raise ValueError(f'{run_dir}: unexpected number of scenarios')
    indexed_scenarios = {}
    indexed_parts = defaultdict(list)
    for s in scenarios:
        key = (s['제품ID'], s['시나리오ID'])
        if key in indexed_scenarios:
            raise ValueError(f'{run_dir}: duplicate scenario {key}')
        if s['모델_상태'] != 'HANNA_2026_실제모델_예측':
            raise ValueError(f'{run_dir}: fake/unknown HANNA output label')
        indexed_scenarios[key] = s
    for p in parts:
        key = (p['제품ID'], p['시나리오ID'])
        if key not in indexed_scenarios:
            raise ValueError(f'{run_dir}: unexpected component {key}')
        indexed_parts[key].append(p)
    if set(indexed_parts) != set(indexed_scenarios):
        raise ValueError(f'{run_dir}: missing component groups')
    results = defaultdict(lambda: {'ratios':[], 'diffs':[], 'raoult':[], 'hanna':[], 'nominal_flips':[],
                                   'robust_flips':[], 'ambiguous':[], 'changed_cas':Counter()})
    for key, sc in indexed_scenarios.items():
        pid, sid = key
        comp = indexed_parts[key]
        if pid not in known_components:
            raise ValueError(f'{run_dir}: product {pid} unknown')
        expected_cas = set(known_components[pid])
        if len(comp) != len(expected_cas) or {p['CAS'] for p in comp} != expected_cas:
            raise ValueError(f'{run_dir}: CAS coverage mismatch for {key}')
        total_r = checked_number(sc['총압력_Raoult_kPa'], 'Raoult', strict_positive=True)
        total_h = checked_number(sc['총압력_HANNA_kPa'], 'HANNA', strict_positive=True)
        ratio = checked_number(sc['HANNA_대_Raoult_비'], 'ratio', strict_positive=True)
        diff = checked_number(sc['HANNA_압력차_kPa'], 'diff')
        if not math.isclose(total_h/total_r, ratio, abs_tol=1e-6) or not math.isclose(total_h-total_r, diff, abs_tol=1e-5):
            raise ValueError(f'{run_dir}: inconsistent scenario pressures {key}')
        wr, xx, yr, yh, pp_r, pp_h = [], [], [], [], [], []
        for p in comp:
            cas = p['CAS']
            mass = checked_number(p['질량백분율'], 'mass')
            x = checked_number(p['액상몰분율'], 'mole fraction')
            y_r = checked_number(p['Raoult_기상몰분율'], 'Raoult y')
            y_h = checked_number(p['HANNA_기상몰분율'], 'HANNA y')
            gamma = checked_number(p['HANNA_gamma'], 'gamma', strict_positive=True)
            ln_gamma = checked_number(p['HANNA_ln_gamma'], 'ln_gamma')
            pr = checked_number(p['Raoult_분압_kPa'], 'Raoult partial')
            ph = checked_number(p['HANNA_분압_kPa'], 'HANNA partial')
            mw = checked_number(props[cas]['MW_g_mol'], 'MW', strict_positive=True)
            vp = checked_number(props[cas]['VP_25C_kPa'], 'P_sat', strict_positive=True)
            lo, hi = known_components[pid][cas]
            if mass < lo - 1e-5 or mass > hi + 1e-5 or x < -1e-7 or y_r < -1e-7 or y_h < -1e-7 or pr < -1e-6 or ph < -1e-6:
                raise ValueError(f'{run_dir}: composition/pressure bound violation {key}/{cas}')
            if not math.isclose(math.log(gamma), ln_gamma, abs_tol=2e-6):
                raise ValueError(f'{run_dir}: gamma transform mismatch {key}/{cas}')
            if not math.isclose(x * vp, pr, rel_tol=3e-6, abs_tol=1e-6):
                raise ValueError(f'{run_dir}: Raoult partial pressure mismatch {key}/{cas}')
            if not math.isclose(pr*gamma, ph, rel_tol=3e-6, abs_tol=2e-6):
                raise ValueError(f'{run_dir}: HANNA partial pressure mismatch {key}/{cas}')
            if not math.isclose(pr/total_r, y_r, abs_tol=3e-6) or not math.isclose(ph/total_h, y_h, abs_tol=3e-6):
                raise ValueError(f'{run_dir}: equilibrium vapor fraction mismatch {key}/{cas}')
            wr.append(mass); xx.append((x,mass,mw)); yr.append((y_r,cas)); yh.append((y_h,cas))
            pp_r.append(pr); pp_h.append(ph)
        if abs(sum(wr) - 100) > 2e-4 or abs(sum(x for x,_,_ in xx) - 1) > 2e-6:
            raise ValueError(f'{run_dir}: composition does not sum to 100% / 1 {key}')
        # Confirm mass-to-mole conversion using rounded output tolerances.
        mole_raw = [mass/mw for _,mass,mw in xx]
        mole_den = sum(mole_raw)
        if any(abs(x-(m/mole_den)) > 3e-6 for (x,_,_),m in zip(xx,mole_raw)):
            raise ValueError(f'{run_dir}: mass-to-mole mismatch {key}')
        if abs(sum(pp_r)-total_r) > 2e-5 or abs(sum(pp_h)-total_h) > 2e-5:
            raise ValueError(f'{run_dir}: sum of partial pressures mismatch {key}')
        if abs(sum(v for v,_ in yr)-1)>3e-6 or abs(sum(v for v,_ in yh)-1)>3e-6:
            raise ValueError(f'{run_dir}: sum of vapor fractions mismatch {key}')
        yr.sort(reverse=True); yh.sort(reverse=True)
        # Detect near-ties in either model, so a tiny numerical separation is not a robust "switch".
        margin_r = yr[0][0] - yr[1][0]
        margin_h = yh[0][0] - yh[1][0]
        ambiguous = min(margin_r, margin_h) < top_tie_margin
        flip = yr[0][1] != yh[0][1]
        item = results[pid]
        item['ratios'].append(ratio)
        item['diffs'].append(diff)
        item['raoult'].append(total_r)
        item['hanna'].append(total_h)
        item['nominal_flips'].append(int(flip))
        item['robust_flips'].append(int(flip and not ambiguous))
        item['ambiguous'].append(int(ambiguous))
        if flip:
            item['changed_cas'][(yr[0][1], yh[0][1])] += 1
    return results, len(scenarios), len(parts)


def summarize(run_root, out_dir, top_tie_margin, *, permit_incomplete=False):
    run_root, out_dir = Path(run_root), Path(out_dir)
    entries = sorted(run_root.glob('*/' + METADATA_FILE))
    if not entries:
        raise FileNotFoundError('No local HANNA runs in run-root; run paper_experiment_271.py run first')
    expected = expected_sha()
    props = {r['CAS']:r for r in read_csv(PROPERTY)}
    known_components = defaultdict(dict)
    for r in read_csv(COMPONENTS):
        pid, cas = r['제품ID'], r['CAS 번호']
        if cas in known_components[pid]:
            raise ValueError(f'Duplicate CAS in source MSDS product: {pid}/{cas}')
        known_components[pid][cas] = (float(r['함량 하한(%)']), float(r['함량 상한(%)']))
    classes = {r['제품ID']:r for r in read_csv(CLASSIFICATION)}
    per_run, scenarios_total, components_total = [], 0, 0
    seen = set()
    models = set()
    product_sets = set()
    for path in entries:
        meta = json.loads(path.read_text(encoding='utf-8-sig'))
        run_dir=path.parent
        identity=(meta.get('sampler'),meta.get('seed'))
        if identity in seen:
            raise ValueError('Duplicate sampler/seed experiment run')
        seen.add(identity)
        if int(meta.get('documents_excluded',-1)) and not permit_incomplete:
            raise ValueError(f'{run_dir}: excluded documents; use --permit-incomplete only if justified')
        validated, n_s, n_p = validate_run(run_dir, meta, known_components, props, expected, top_tie_margin)
        if not validated:
            raise ValueError(f'Empty run {run_dir}')
        product_sets.add(frozenset(validated))
        models.add(meta.get('model'))
        scenarios_total += n_s
        components_total += n_p
        for pid, a in sorted(validated.items()):
            ratios, diffs = a['ratios'], a['diffs']
            per_run.append({'표본추출':identity[0], '난수시드':identity[1], '제품ID':pid,
                            '제품명':classes[pid].get('제품명(한글)',''),
                            '제품군_잠정':classes[pid].get('1차 용도분류',''),
                            '시나리오수':len(ratios),
                            'Raoult_압력중앙값_kPa':statistics.median(a['raoult']),
                            'HANNA_압력중앙값_kPa':statistics.median(a['hanna']),
                            '압력비_중앙값':statistics.median(ratios),
                            '압력비_5백분위':qtile(ratios,.05),
                            '압력비_95백분위':qtile(ratios,.95),
                            '압력비_최솟값':min(ratios),'압력비_최댓값':max(ratios),
                            '압력비1초과_비율':sum(r>1 for r in ratios)/len(ratios),
                            '압력차_중앙값_kPa':statistics.median(diffs),
                            '최다증기성분_명목전환율':statistics.mean(a['nominal_flips']),
                            '최다증기성분_근접동률제외전환율':statistics.mean(a['robust_flips']),
                            '최다증기성분_근접동률비율':statistics.mean(a['ambiguous']),
                            '전환_방향_빈도':json.dumps({f'{k[0]}>{k[1]}':v for k,v in a['changed_cas'].items()},ensure_ascii=False)})
    if len(product_sets) != 1:
        raise ValueError('Runs covered different MSDS sets; pooling is invalid')
    by_doc = defaultdict(list)
    for r in per_run:
        by_doc[(r['표본추출'],r['제품ID'])].append(r)
    per_doc=[]
    for (sampler,pid),rs in sorted(by_doc.items()):
        m=[r['압력비_중앙값'] for r in rs]
        per_doc.append({'표본추출':sampler,'제품ID':pid,'제품명':rs[0]['제품명'],
                        '반복시드수':len(rs),'시드별_압력비중앙값_중앙값':statistics.median(m),
                        '시드별_압력비중앙값_최솟값':min(m), '시드별_압력비중앙값_최댓값':max(m),
                        '시드별_압력비중앙값_범위':max(m)-min(m),
                        '시드별_압력증가시나리오_비율_최솟값':min(r['압력비1초과_비율'] for r in rs),
                        '시드별_명목전환율_중앙값':statistics.median(r['최다증기성분_명목전환율'] for r in rs),
                        '시드별_근접동률제외전환율_중앙값':statistics.median(r['최다증기성분_근접동률제외전환율'] for r in rs)})
    compared=[]
    for pid in sorted({r['제품ID'] for r in per_doc}):
        vals={r['표본추출']:r['시드별_압력비중앙값_중앙값'] for r in per_doc if r['제품ID']==pid}
        if set(vals)=={'legacy','gibbs'}:
            compared.append({'제품ID':pid,'legacy_압력비중앙값':vals['legacy'],
                             'gibbs_압력비중앙값':vals['gibbs'],
                             '표본추출_변경에따른_압력비중앙값_차이':vals['gibbs']-vals['legacy']})
    out_dir.mkdir(parents=True,exist_ok=True)
    write_csv(out_dir/'paper_반복별_통계.csv',per_run, list(per_run[0]))
    write_csv(out_dir/'paper_제품별_안정성.csv',per_doc,list(per_doc[0]))
    write_csv(out_dir/'paper_표본추출_민감도.csv',compared, ['제품ID','legacy_압력비중앙값','gibbs_압력비중앙값','표본추출_변경에따른_압력비중앙값_차이'])
    products_covered = sorted(next(iter(product_sets)))
    runs_by_sampler = Counter(r['표본추출'] for r in per_run if r['제품ID']==products_covered[0])
    seeds_per_sampler=dict(runs_by_sampler)
    audit=read_csv(run_root/'baseline_audit'/'271_제품별_분석가능성_감사.csv') if (run_root/'baseline_audit'/'271_제품별_분석가능성_감사.csv').exists() else []
    external_status='NOT VALIDATED: no independent observed VLE/headspace data incorporated'
    metadata={'created_utc':datetime.now(timezone.utc).isoformat(),
              'status':'model-to-model comparison only; no empirical accuracy/worker exposure/substitution evidence',
              'input_sha256':expected,'run_folders':len(entries),'scenarios_checked':scenarios_total,
              'component_records_checked':components_total,'products_covered':products_covered,
              'seeds_per_sampler':seeds_per_sampler,
              'top_tie_margin_in_vapor_mole_fraction':top_tie_margin,
              'external_VLE_validation':external_status,
              'same_use_TWA_evidence':'not assessed: needs documented and verified OEL/compatibility',
              'optimization_evidence':'not assessed: requires experimentally justified process constraints',
              'generated_composition_percentiles':'descriptive sampling intervals; NOT population confidence intervals',
              'gibbs':'approximate constrained-space sampling only; finite-chain convergence unproven',
              'MSDS_audit_rows_available':len(audit)}
    (out_dir/'paper_실험재현정보.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')
    lines=['# 국내 신나 HANNA–Raoult 모형비교 민감도 분석', '',
           '## 결론의 범위','',
           f'- 검증된 계산: {len(entries)}개 실행, {scenarios_total:,}개 가상 조성 시나리오, {components_total:,}개 성분행.',
           f'- 문서 {len({r["제품ID"] for r in per_run})}건의 Raoult와 HANNA의 **동일 조성** 비교이며 독립 실측 검증이 아니다.',
           '- 실제 제조 배합, 작업장 노출농도, 대체제품 안전성, 최적 배합의 실험적 타당성을 증명하지 않는다.',
           '- 표본추출로 얻은 분위수와 시드 간 차이는 **설계된 가상 조성 분포**에 관한 기술통계다. 통계적 신뢰구간 또는 확률적 제품 순위가 아니다.',
           '- 같은 MSDS 명칭/개정본이 중복될 수 있어 6건을 독립 제품 6종이라고 부를 수 없다.',
           '', '## 안정성 요약','',
           '| 추출법 | 제품ID | 시드 수 | 압력비 중앙값의 시드 중앙값 | 시드 최솟값–최댓값 | 근접 동률 제외 우세성분 전환율(중앙) |',
           '|---|---|---:|---:|---|---:|']
    for row in per_doc:
        lines.append(f'| {row["표본추출"]} | {row["제품ID"]} | {row["반복시드수"]} | {row["시드별_압력비중앙값_중앙값"]:.4f} | {row["시드별_압력비중앙값_최솟값"]:.4f}–{row["시드별_압력비중앙값_최댓값"]:.4f} | {row["시드별_근접동률제외전환율_중앙값"]:.1%} |')
    lines+=['','## 논문 준비 여부', '',
            '- 모형 비교와 조성 민감도: 탐색적 결과로 보고 가능하되 실측 검증 없이 정확도 향상을 주장할 수 없음.',
            '- 독립 VLE 또는 headspace 참값 기반 MAE/RMSE 및 도료계 외삽 검증: **미완료**.',
            '- 한국 현행 TWA 및 독성영향 가산성/공정 호환성 검증: **미완료**.',
            '- 조성 상관, 제품 중복/개정본, 원물성 정확도, 샘플링 혼합 진단: 추가 연구 필요.',
            '- 따라서 현재 자료만으로 저위험 신나 구매 순위나 최적 혼합비를 논문 결론으로 주장해서는 안 됨.','']
    (out_dir/'paper_논문가치_판정.md').write_text('\n'.join(lines),encoding='utf-8')
    print(f'[COMPLETE] runs={len(entries)} scenarios={scenarios_total} components={components_total}; report={out_dir}')
    return metadata


def main(argv=None):
    parser = argparse.ArgumentParser(description='Paper-grade sensitivity workflow (no implied empirical validation)')
    sub = parser.add_subparsers(dest='action',required=True)
    r = sub.add_parser('run', help='Run multiple real local HANNA simulations and generate descriptive reports')
    r.add_argument('--root', type=Path, default=RESEARCH_ROOT/'output'/'paper_full')
    r.add_argument('--seeds', type=int, nargs='+', default=[20261023,20261024,20261025])
    r.add_argument('--samplers', nargs='+', choices=['legacy','gibbs'], default=['legacy','gibbs'])
    r.add_argument('--samples',type=int,default=200)
    r.add_argument('--batch-size',type=int,default=8)
    r.add_argument('--burn-in',type=int,default=300)
    r.add_argument('--thin',type=int,default=12)
    r.add_argument('--top-tie-margin',type=float,default=0.01,
                   help='A switch within 1 pct-point of either model top share counts as ambiguous')
    s = sub.add_parser('summarize', help='Analyze existing experiments; no HANNA model required')
    s.add_argument('--root',type=Path,default=RESEARCH_ROOT/'output'/'paper_full')
    s.add_argument('--top-tie-margin',type=float,default=0.01)
    s.add_argument('--permit-incomplete',action='store_true')
    args=parser.parse_args(argv)
    if not 0 <= args.top_tie_margin < 1:
        parser.error('top-tie-margin must be in [0,1)')
    if args.action=='run':
        if args.samples < 20 or args.batch_size < 1 or args.burn_in < 0 or args.thin < 1:
            parser.error('samples>=20, batch-size>0, burn-in>=0, thin>=1 required')
        if len(set(args.seeds))!=len(args.seeds) or len(set(args.samplers))!=len(args.samplers):
            parser.error('duplicate seed/sampler setting')
        for sampler in args.samplers:
            for seed in args.seeds:
                run_one(args.root/'runs'/f'{sampler}_seed_{seed}',sampler=sampler,seed=seed,
                        samples=args.samples,batch_size=args.batch_size,burn_in=args.burn_in,thin=args.thin)
        run_root=args.root/'runs'
        # The sensitivity results stand on their own; no automatic OEL approval or optimization.
        summarize(run_root,args.root/'reports',args.top_tie_margin)
    else:
        # Accept --root that is either the experiment base or the runs folder.
        run_root=args.root/'runs' if (args.root/'runs').exists() else args.root
        reports=args.root/'reports' if (args.root/'runs').exists() else args.root.parent/'reports'
        summarize(run_root,reports,args.top_tie_margin,permit_incomplete=args.permit_incomplete)

if __name__=='__main__':
    main()