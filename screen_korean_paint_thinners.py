# -*- coding: utf-8 -*-
"""Korean MSDS paint-thinner vapor screening / product shortlisting (25 C).

Does NOT determine workplace exposure, legal MSDS omissions, or safe substitution.
Domestic products and Korean legal OELs only. No NIOSH OEL fallback.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
DEFAULT_PRODUCTS = ROOT / 'domestic_product_inputs' / 'products.csv'
DEFAULT_OEL = ROOT / 'law_data' / '행정규칙' / '화학물질 및 물리적 인자의 노출기준.json'
DEFAULT_PROPERTIES = ROOT / 'thinner_vle_reference_25C.csv'
COLUMNS = ('product_id','product_name','manufacturer','market_country','product_evidence_url',
           'msds_revision','source_type','use_group','compatibility_group',
           'compatibility_verified','compatibility_evidence_url','CAS','component_name',
           'wt_pct_min','wt_pct_max','msds_disclosure')
OUT_COLUMNS = ('product_id','product_name','model','status','reason',
               'max_component_priority','top_component_CAS','top_component',
               'sum_component_priority_descriptive_only','total_eq_vapor_pressure_kPa',
               'use_group','compatibility_group','compatibility_verified','source_type')
CAS_PATTERN = re.compile(r'\[(\d{2,7}-\d{2}-\d)\]')


def read_csv(path):
    with path.open('r', encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(fields), extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def parse_korean_twa(path):
    """Read only Annex 1 key 000100: TWA ppm, never STEL or mg/m3.

    Duplicate CAS with conflicting limits is set to ambiguous and excluded.
    """
    data = json.loads(path.read_text(encoding='utf-8-sig'))['AdmRulService']
    annexes = [v for v in data['별표']['별표단위'] if str(v.get('별표키')) == '000100'
               and '노출기준' in v.get('별표제목','')]
    if len(annexes) != 1:
        raise ValueError('국내 노출기준 별표 1 (별표키 000100) 유일하게 확인되지 않음')
    def flatten(x):
        if isinstance(x, list):
            for y in x:
                yield from flatten(y)
        elif isinstance(x, str):
            yield x
    found = defaultdict(list)
    for line in flatten(annexes[0]['별표내용']):
        if '│' not in line or '[' not in line:
            continue
        parts = line.split('│')
        if len(parts) < 10 or not parts[1].strip().isdigit():
            continue
        cas_list = CAS_PATTERN.findall(line)
        try:
            ppm = float(parts[5].strip())
            if ppm <= 0 or not math.isfinite(ppm):
                continue
        except (ValueError, IndexError):
            continue
        for cas in cas_list:
            found[cas].append((ppm, parts[2].strip(), parts[3].strip()))
    result = {}
    for cas, vals in found.items():
        limits = {x[0] for x in vals}
        if len(limits) == 1:
            result[cas] = {'TWA_ppm': vals[0][0], 'name_kr': vals[0][1], 'name_en': vals[0][2],
                           'match_status':'MATCHED'}
        else:
            result[cas] = {'TWA_ppm': None, 'name_kr': '', 'name_en': '',
                           'match_status':'CONFLICTING_DOMESTIC_TWA'}
    return result, data['행정규칙기본정보']


def load_physchem(path):
    props = {}
    for row in read_csv(path):
        cas = row.get('CAS','').strip()
        if cas in props:
            raise ValueError('물성 CAS 중복: ' + cas)
        mw, vp = float(row['MW_g_mol']), float(row['VP_25C_kPa'])
        if not all(math.isfinite(x) and x > 0 for x in (mw, vp)) or not row['SMILES'].strip():
            raise ValueError('물성 누락: ' + cas)
        props[cas] = {'MW_g_mol':mw, 'VP_25C_kPa':vp, 'SMILES':row['SMILES'].strip(),
                      'chemical':row.get('chemical','').strip()}
    return props


def check_product(rows, props, oel):
    p = rows[0]
    def invalid(why):
        return why, []
    if any(x.get('product_id') != p.get('product_id') for x in rows):
        return invalid('PRODUCT_ID_INCONSISTENT')
    metadata = [x for x in COLUMNS[:11]]
    if any(any(row.get(col,'').strip() != p.get(col,'').strip() for col in metadata) for row in rows):
        return invalid('PRODUCT_METADATA_INCONSISTENT')
    if p.get('market_country','').strip().upper() != 'KR':
        return invalid('NOT_VERIFIED_KOREAN_MARKET')
    if not p.get('product_evidence_url','').strip() or not p.get('msds_revision','').strip():
        return invalid('MISSING_DOMESTIC_PRODUCT_OR_MSDS_PROVENANCE')
    if p.get('source_type') not in ('manufacturer_msds','supplier_msds','measured_verified_product'):
        return invalid('SOURCE_TYPE_NOT_VERIFIED')
    if not p.get('use_group','').strip():
        return invalid('MISSING_USE_GROUP')
    cas_seen, comps = set(), []
    for row in rows:
        cas = row.get('CAS','').strip()
        if not cas or cas in cas_seen:
            return invalid('EMPTY_OR_DUPLICATE_CAS')
        cas_seen.add(cas)
        try:
            lo = float(row.get('wt_pct_min',''))
            hi = float(row.get('wt_pct_max',''))
        except ValueError:
            return invalid('UNKNOWN_CONCENTRATION')
        if not all(math.isfinite(x) for x in (lo, hi)) or lo < 0 or hi > 100 or lo > hi:
            return invalid('INVALID_CONCENTRATION')
        if abs(lo-hi) > 1e-9:
            return invalid('CONCENTRATION_RANGE_NOT_RANKED')
        if cas not in props:
            return invalid('MISSING_PHYSICOCHEMICAL_DATA:' + cas)
        if cas not in oel or oel[cas]['match_status'] != 'MATCHED':
            return invalid('MISSING_OR_AMBIGUOUS_KOREAN_TWA:' + cas)
        comps.append((cas, lo))
    if len(comps) < 2 or abs(sum(x[1] for x in comps)-100) > 0.1:
        return invalid('INCOMPLETE_COMPOSITION_MASS_BALANCE')
    return '', comps


def screen_composition(comps, props, oel, gammas):
    cas = [x[0] for x in comps]
    gammas = np.asarray(gammas, dtype=float).reshape(-1)
    if gammas.shape != (len(cas),) or not np.all(np.isfinite(gammas)) or np.any(gammas <= 0):
        raise ValueError('INVALID_ACTIVITY_COEFFICIENTS')
    mole = np.asarray([wt / props[c]['MW_g_mol'] for c,wt in comps], float)
    xs = mole / mole.sum()
    partials = xs * gammas * [props[c]['VP_25C_kPa'] for c in cas]
    if not np.all(np.isfinite(partials)) or np.any(partials <= 0):
        raise ValueError('INVALID_PARTIAL_PRESSURES')
    total_p = float(sum(partials))
    if total_p >= 101.325:
        raise ValueError('TOTAL_EQUILIBRIUM_PRESSURE_GE_ATMOSPHERIC: model assumption requires review')
    priority = [float(x/101.325*1e6/oel[c]['TWA_ppm']) for c,x in zip(cas,partials)]
    # Component ratio of predicted saturated partial pressure (ppm) to Korean TWA.
    # NOT occupational air concentration, legal exceedance, or cumulative toxicological risk.
    details = []
    for i,(c,wt) in enumerate(comps):
        details.append(dict(CAS=c, chemical=props[c]['chemical'],liquid_wt_pct=wt,
                            liquid_mol_frac=float(xs[i]), gamma=float(gammas[i]),
                            predicted_partial_pressure_kPa=float(partials[i]),
                            predicted_vapor_mol_pct=float(100*partials[i]/total_p),
                            KOREAN_TWA_ppm=oel[c]['TWA_ppm'],
                            equilibrium_screening_priority=priority[i],
                            priority_interpretation='SCREENING_ONLY_NOT_WORKPLACE_EXPOSURE'))
    k = int(np.argmax(priority))
    summary = dict(max_component_priority=priority[k], top_component_CAS=cas[k],
                   top_component=props[cas[k]]['chemical'],
                   sum_component_priority_descriptive_only=float(sum(priority)),
                   total_eq_vapor_pressure_kPa=total_p)
    return summary, details


def run(input_csv, output_dir, model, hanna_repo, oel_path, physchem_path):
    from compare_thinner_raoult_unifac_hanna import GammaModels
    records = read_csv(input_csv)
    # An empty validated template is allowed; do not pretend products were analyzed.
    with Path(input_csv).open('r', encoding='utf-8-sig', newline='') as f:
        input_fields = set(csv.DictReader(f).fieldnames or [])
    missing = set(COLUMNS) - input_fields
    if missing:
        raise ValueError('제품 입력 컬럼 누락: '+', '.join(sorted(missing)))
    oel, legal_meta = parse_korean_twa(oel_path)
    props = load_physchem(physchem_path)
    batches = defaultdict(list)
    for r in records:
        batches[r['product_id'].strip()].append(r)
    model_ids = list(dict.fromkeys(model))
    predictor = GammaModels(Path(hanna_repo))
    audit, detail, summary, crosswalk = [], [], [], []
    for pid, rows in sorted(batches.items()):
        product = rows[0]
        reason, comps = check_product(rows, props, oel)
        for row in rows:
            c = row['CAS'].strip()
            d = oel.get(c, {})
            crosswalk.append(dict(product_id=pid,CAS=c,component_name=row.get('component_name',''),
                                  Korean_TWA_ppm=d.get('TWA_ppm',''),
                                  Korean_TWA_status=d.get('match_status','NOT_FOUND'),
                                  OEL_notice=legal_meta.get('발령번호',''),
                                  OEL_effective_date=legal_meta.get('시행일자','')))
        audit.append(dict(product_id=pid,product_name=product['product_name'],
                          market_country=product['market_country'], source_type=product['source_type'],
                          status='ELIGIBLE_FOR_MODEL' if not reason else 'EXCLUDED',reason=reason,
                          n_components=len(rows)))
        for mid in model_ids:
            base = dict(product_id=pid,product_name=product['product_name'],model=mid,
                        use_group=product['use_group'],compatibility_group=product['compatibility_group'],
                        compatibility_verified=product['compatibility_verified'],
                        source_type=product['source_type'])
            if reason:
                summary.append(dict(base,status='EXCLUDED',reason=reason))
                continue
            try:
                xs = np.array([w/props[c]['MW_g_mol'] for c,w in comps])
                xs /= xs.sum()
                gammas = predictor.predict(mid,[c for c,_ in comps],xs,298.15,props)
                smry, det = screen_composition(comps, props, oel, gammas)
                summary.append(dict(base,status='OK',reason='',**smry))
                detail.extend(dict(product_id=pid,product_name=product['product_name'],
                                   model=mid,**d) for d in det)
            except Exception as exc:
                summary.append(dict(base,status='MODEL_FAILED',
                                    reason=f'{type(exc).__name__}: {exc}'))
    candidates = []
    # Only verified, exact same application group and manufacturer-supported compatibility evidence.
    viable = [x for x in summary if x['model']=='hanna' and x['status']=='OK']
    by_group = defaultdict(list)
    for x in viable:
        p = batches[x['product_id']][0]
        if (p['compatibility_verified'].strip().lower() in ('yes','true','1')
                and p['compatibility_group'].strip() and p['compatibility_evidence_url'].strip()):
            by_group[(x['use_group'],x['compatibility_group'])].append(x)
    for group, members in by_group.items():
        if len(members) < 2:
            continue
        members.sort(key=lambda x:x['max_component_priority'])
        for rank,x in enumerate(members,1):
            candidates.append(dict(product_id=x['product_id'],product_name=x['product_name'],
                                   use_group=group[0],compatibility_group=group[1],
                                   rank_in_verified_group=rank,
                                   max_component_priority=x['max_component_priority'],
                                   selection_status='SCREENING_CANDIDATE_NOT_CERTIFIED_SUBSTITUTE',
                                   compatibility_evidence_url=batches[x['product_id']][0]['compatibility_evidence_url']))
    out = Path(output_dir)
    write_csv(out/'01_domestic_data_audit.csv',audit,('product_id','product_name','market_country','source_type','status','reason','n_components'))
    write_csv(out/'02_korean_twa_crosswalk.csv',crosswalk,('product_id','CAS','component_name','Korean_TWA_ppm','Korean_TWA_status','OEL_notice','OEL_effective_date'))
    write_csv(out/'03_component_vapor_screen.csv',detail,('product_id','product_name','model','CAS','chemical','liquid_wt_pct','liquid_mol_frac','gamma','predicted_partial_pressure_kPa','predicted_vapor_mol_pct','KOREAN_TWA_ppm','equilibrium_screening_priority','priority_interpretation'))
    write_csv(out/'04_product_screen_summary.csv',summary,OUT_COLUMNS)
    write_csv(out/'05_eligible_selection_candidates.csv',candidates,('product_id','product_name','use_group','compatibility_group','rank_in_verified_group','max_component_priority','selection_status','compatibility_evidence_url'))
    metadata = {'jurisdiction':'KR','product_data_requirement':'VERIFIED_KOREAN_MARKET','oel_source_notice':legal_meta.get('발령번호'),
                'oel_effective_date':legal_meta.get('시행일자'),'msds_notice':'2026-26',
                'temperature_K':298.15,'atmospheric_pressure_kPa':101.325,
                'physics_reference':str(physchem_path),'models_requested':model_ids,
                'us_niosh_oel_used':False,'price_required':False,
                'workplace_exposure_estimated':False,'regulatory_omission_determined':False,
                'caveat':'Equilibrium vapor screen only; MSDS ranges/unknown ingredients not modeled. Separate validation required.'}
    (out/'06_run_metadata.json').write_text(json.dumps(metadata,indent=2,ensure_ascii=False),encoding='utf-8')
    return {'total':len(batches),'eligible':sum(r['status']=='ELIGIBLE_FOR_MODEL' for r in audit),
            'model_ok':sum(r['status']=='OK' for r in summary),'candidates':len(candidates),'output_dir':str(out)}


def cli():
    p = argparse.ArgumentParser(description='Domestic Korean MSDS + Korean TWA + AI equilibrium vapor screening')
    p.add_argument('--input',type=Path,default=DEFAULT_PRODUCTS)
    p.add_argument('--output-dir',type=Path,default=ROOT/'domestic_screening_output')
    p.add_argument('--oel-json',type=Path,default=DEFAULT_OEL)
    p.add_argument('--reference',type=Path,default=DEFAULT_PROPERTIES)
    p.add_argument('--hanna-repo',type=Path,default=Path(r'D:\fire\HANNA'))
    p.add_argument('--models',nargs='+',choices=['hanna','unifac','raoult'],default=['hanna','unifac','raoult'])
    a=p.parse_args()
    print(json.dumps(run(a.input,a.output_dir,a.models,a.hanna_repo,a.oel_json,a.reference),ensure_ascii=False,indent=2))

if __name__=='__main__':
    cli()