# -*- coding: utf-8 -*-
"""Independent isothermal bubble-pressure validation against approved VLE data.

No built-in measurements. Must supply traced, manually verified experimental data
and a separately vetted 25C pure-component vapor pressure table. Else refuse.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import statistics
from pathlib import Path

from paper_experiment_271 import read_csv, write_csv, PROPERTY, checked_number

HEADERS = ['reference_id','source_url_or_doi','verification_status','temperature_K',
           'CAS_JSON','liquid_x_JSON','observed_total_pressure_kPa']
PROP_HEADERS = ['CAS','temperature_K','VP_25C_kPa','reference_url_or_doi','verification_status']


def parse_list(v, name):
    try:
        a = json.loads(v)
    except (TypeError,ValueError) as exc:
        raise ValueError(f'{name}: valid JSON array required') from exc
    if not isinstance(a,list):
        raise ValueError(f'{name}: JSON list required')
    return a


def validate_rows(measurements, evidence, base_props):
    evidence_map = {}
    for r in evidence:
        cas=r['CAS'].strip()
        if not cas or cas in evidence_map:
            raise ValueError('Evidence CAS empty or duplicated')
        if r['verification_status'].strip() != '승인':
            continue
        if not r['reference_url_or_doi'].strip() or not r['reference_url_or_doi'].startswith(('https://','http://','10.')):
            raise ValueError(f'{cas}: property provenance URL / DOI missing')
        t=checked_number(r['temperature_K'],'property temperature',strict_positive=True)
        vp=checked_number(r['VP_25C_kPa'],'property VP',strict_positive=True)
        if abs(t-298.15)>1e-4:
            raise ValueError(f'{cas}: property evidence must verify 298.15 K')
        if cas not in base_props or not math.isclose(vp,checked_number(base_props[cas]['VP_25C_kPa'],'existing VP'),rel_tol=1e-5):
            raise ValueError(f'{cas}: property evidence value disagrees with local property table')
        evidence_map[cas]=r
    if not measurements:
        raise ValueError('No independent VLE measurements supplied')
    cleaned=[]
    seen=set()
    for r in measurements:
        rid=r['reference_id'].strip()
        if not rid or rid in seen:
            raise ValueError('Missing or duplicate measurement row id')
        seen.add(rid)
        src=r['source_url_or_doi'].strip()
        if r['verification_status'].strip() != '승인' or not src.startswith(('https://','http://','10.')):
            raise ValueError(f'{rid}: original external VLE source and manual approval required')
        if abs(checked_number(r['temperature_K'],'measurement temperature')-298.15)>1e-4:
            raise ValueError(f'{rid}: only 25C models and saturation references supported')
        cas=[str(c).strip() for c in parse_list(r['CAS_JSON'],'CAS_JSON')]
        x=[checked_number(v,'mole fraction') for v in parse_list(r['liquid_x_JSON'],'liquid_x_JSON')]
        if len(cas)<2 or len(cas)!=len(x) or len(set(cas)) != len(cas):
            raise ValueError(f'{rid}: component order or identity invalid')
        if any(v < 0 or v>1 for v in x) or abs(sum(x)-1)>1e-7:
            raise ValueError(f'{rid}: liquid mole fractions invalid')
        if any(c not in evidence_map for c in cas):
            raise ValueError(f'{rid}: some CAS lack approved property evidence')
        if any(not base_props.get(c,{}).get('SMILES','').strip() for c in cas):
            raise ValueError(f'{rid}: missing SMILES in the research property table')
        obs=checked_number(r['observed_total_pressure_kPa'],'observed VLE pressure',strict_positive=True)
        cleaned.append((rid,src,cas,x,obs))
    return cleaned


def main(argv=None):
    p=argparse.ArgumentParser(description='Gated independent VLE accuracy evaluation')
    p.add_argument('--write-templates',type=Path)
    p.add_argument('--measurements',type=Path)
    p.add_argument('--property-evidence',type=Path)
    p.add_argument('--output',type=Path)
    args=p.parse_args(argv)
    if args.write_templates:
        args.write_templates.mkdir(parents=True,exist_ok=True)
        write_csv(args.write_templates/'독립VLE_실측입력_서식.csv',[],HEADERS)
        write_csv(args.write_templates/'순물질_포화증기압_출처승인_서식.csv',[],PROP_HEADERS)
        print('Created EMPTY verification templates; no VLE truth invented.')
        return
    if not (args.measurements and args.property_evidence and args.output):
        p.error('--measurements, --property-evidence, and --output required')
    base_props={r['CAS'].strip():r for r in read_csv(PROPERTY)}
    cleaned=validate_rows(read_csv(args.measurements),read_csv(args.property_evidence),base_props)
    from hanna_gamma_adapter import predict_gammas
    rows=[]
    for rid,source,cas,x,observed in cleaned:
        vp=[checked_number(base_props[c]['VP_25C_kPa'],'P_sat',strict_positive=True) for c in cas]
        smiles=[base_props[c]['SMILES'].strip() for c in cas]
        gamma=predict_gammas(smiles,x,298.15)
        if len(gamma)!=len(cas) or not all(math.isfinite(v) and v>0 for v in gamma):
            raise ValueError(f'{rid}: HANNA gamma invalid; no Raoult substitution')
        raoult=sum(a*b for a,b in zip(x,vp))
        hanna=sum(a*b*c for a,b,c in zip(x,vp,gamma))
        rows.append({'reference_id':rid,'source_url_or_doi':source,
                     'CAS_JSON':json.dumps(cas,ensure_ascii=False),
                     'n_components':len(cas),'observed_kPa':observed,
                     'Raoult_kPa':raoult,'HANNA_kPa':hanna,
                     'Raoult_error_kPa':raoult-observed,'HANNA_error_kPa':hanna-observed})
    errs_r=[r['Raoult_error_kPa'] for r in rows]
    errs_h=[r['HANNA_error_kPa'] for r in rows]
    metrics={'n_observations':len(rows),'n_independent_source_labels':len({r['source_url_or_doi'] for r in rows}),
             'temperature_K':298.15,
             'Raoult_MAE_kPa':statistics.mean(abs(x) for x in errs_r),
             'HANNA_MAE_kPa':statistics.mean(abs(x) for x in errs_h),
             'Raoult_RMSE_kPa':math.sqrt(statistics.mean(x*x for x in errs_r)),
             'HANNA_RMSE_kPa':math.sqrt(statistics.mean(x*x for x in errs_h)),
             'Raoult_bias_kPa':statistics.mean(errs_r),
             'HANNA_bias_kPa':statistics.mean(errs_h),
             'interpretation':'Source labels are not necessarily independent experiments; verify DOI, applicability, data quality, and within-source clustering before claims'}
    args.output.mkdir(parents=True,exist_ok=True)
    write_csv(args.output/'independent_vle_prediction_errors.csv',rows,list(rows[0]))
    (args.output/'independent_vle_metrics.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=2),encoding='utf8')
    print(f"External 25C bubble-pressure validation executed on {len(rows)} independently provided observations. Output: {args.output}")

if __name__=='__main__':
    main()