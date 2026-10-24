# -*- coding: utf-8 -*-
"""Domestic OEL sensitivity of the existing counterfactual thinner simulations.

Inputs are the existing local law_data/*.json and 03_component_predictions.csv.
The script DOES NOT determine whether an MSDS component can lawfully be omitted.
The resulting weighted fraction is not measured exposure or additive toxicity.
"""
from __future__ import annotations
import argparse
import json
import math
import re
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
LAW_DIR = ROOT / 'law_data' / '행정규칙'
OEL_JSON = LAW_DIR / '화학물질 및 물리적 인자의 노출기준.json'
MSDS_JSON = LAW_DIR / '화학물질의 분류·표시 및 물질안전보건자료에 관한 기준.json'
DEFAULT_INPUT = ROOT / 'chemexpo_paint_thinner_output' / 'multichem_screening_output' / '03_component_predictions.csv'
KEYS = ['product_id', 'target_CAS', 'target_wt_pct', 'model_id']


def regulatory_document(path: Path):
    if not path.is_file():
        raise FileNotFoundError(f'Required locally downloaded law JSON is missing: {path}')
    doc = json.loads(path.read_text(encoding='utf-8-sig'))['AdmRulService']
    return doc, doc['행정규칙기본정보']


def find_annex(doc: dict, annex_number: str) -> dict:
    """Select the *main* annex rather than an amended/deleted sub-annex.

    Korean regulatory JSON can contain e.g. two items numbered 0001:
    main Annex 1 with 별표키=000100 and deleted Annex 1-2 with 000102.
    Never silently choose the first duplicate without checking its identity.
    """
    entries = doc['별표']['별표단위']
    number = str(annex_number).zfill(4)
    numbered = [z for z in entries if str(z.get('별표번호', '')).zfill(4) == number]
    canonical_key = number + '00'
    keyed = [z for z in numbered if str(z.get('별표키', '')).zfill(6) == canonical_key]
    if keyed:
        selected = keyed
    else:
        # Compatibility with small fixtures and older JSON exports without keys.
        # Deleted annex items cannot be used as quantitative regulatory tables.
        selected = [z for z in numbered
                    if '삭제' not in str(z.get('별표제목', ''))
                    and isinstance(z.get('별표내용'), list)]
    if len(selected) != 1:
        brief = [(z.get('별표키'), z.get('별표제목')) for z in numbered]
        raise ValueError(
            f'Could not identify unique main annex {number}; candidates={brief}. '
            'Review regulatory JSON annex keys before calculation.'
        )
    return selected[0]


def annex_lines(annex: dict) -> list[str]:
    blocks = annex.get('별표내용')
    if not isinstance(blocks, list):
        raise ValueError(f"Annex {annex.get('별표키')} has no structured table lines")
    return [str(line) for block in blocks
            for line in (block if isinstance(block, list) else [block])]


def parse_twa_row(lines: list[str], cas: str):
    """Read numeric Korean TWA, verifying any synonymous CAS listings agree.

    Example in 2020-48: CAS 108-10-1 appears as Methyl isobutyl ketone
    and Hexone. Accept synonyms ONLY if all reported exposure-limit cells match.
    """
    needle = re.compile(r'\[' + re.escape(cas) + r'\]')
    matches = [line for line in lines if needle.search(line) and line.startswith('│')]
    if not matches:
        raise ValueError(f'CAS {cas} is absent from Korean Annex 1 (manual review required)')
    if len(set(matches)) != len(matches):
        raise ValueError(f'CAS {cas} has duplicate identical OEL table lines')
    parsed = []
    for line in matches:
        cells = [c.strip() for c in line.strip().strip('│').split('│')]
        # Identifier, Korean name, English name, chemical formula,
        # TWA ppm, TWA mg/m3, STEL ppm, STEL mg/m3, remarks.
        if len(cells) < 9:
            raise ValueError(f'Could not safely read TWA cell for CAS {cas}: {cells}')
        try:
            twa = float(cells[4])
        except ValueError as exc:
            raise ValueError(f'CAS {cas} has no unambiguous numeric TWA ppm: {cells[4]!r}') from exc
        if not math.isfinite(twa) or twa <= 0:
            raise ValueError(f'CAS {cas} invalid numeric TWA {twa}')
        try:
            source_order = int(cells[0])
        except ValueError:
            source_order = 10**9
        parsed.append((source_order, cells, line, twa))
    # Do not conflate two different numeric regulations that happen to share
    # a CAS number: require full identity of the four OEL columns.
    limits = {tuple(cells[4:8]) for _, cells, _, _ in parsed}
    if len(limits) != 1:
        raise ValueError(f'Conflicting Korean exposure-limit entries for CAS {cas}: {sorted(limits)}')
    parsed.sort(key=lambda item: (item[0], item[1][2]))
    _, cells, primary_line, twa = parsed[0]
    return {
        'CAS': cas,
        'KR_TWA_ppm': twa,
        'KR_chemical_name': cells[1],
        'KR_name_en': cells[2],
        'KR_OEL_listing_count': len(parsed),
        'KR_OEL_alias_names': ' | '.join(x[1][2] for x in parsed),
        'KR_source_table_row': primary_line,
        'KR_source_all_table_rows': '\n'.join(x[2] for x in parsed),
    }


def verify_msds_law(doc: dict):
    articles = doc['조문내용']
    full_text = json.dumps(articles, ensure_ascii=False)
    if '제11조(작성원칙)' not in full_text or '별표 6' not in full_text:
        raise ValueError('MSDS Article 11 could not be verified; review the legal source')
    annex6 = find_annex(doc, '0006')
    lines = annex_lines(annex6)
    # Verifies the annex actually contains the expected legal categories; this
    # does NOT classify individual CAS or compute legal disclosure obligations.
    joined = '\n'.join(lines)
    if not (re.search(r'발암성[^\n]*0\.1\s*%', joined)
            and re.search(r'급성 독성[^\n]*1\s*%', joined)):
        raise ValueError('MSDS annex 6 hazard cutoffs differ; manual review required')
    return lines


def build_crosswalk(df: pd.DataFrame, oel_doc: dict):
    lines = annex_lines(find_annex(oel_doc, '0001'))
    return pd.DataFrame([parse_twa_row(lines, cas) for cas in sorted(df['CAS'].unique())])


def apply_oel(df: pd.DataFrame, crosswalk: pd.DataFrame) -> pd.DataFrame:
    required = set(KEYS + ['CAS', 'vapor_mol_pct', 'OEL_weighted_share_pct', 'liquid_wt_pct'])
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f'Existing component result lacks expected columns: {sorted(missing)}')
    if df.duplicated(KEYS + ['CAS']).any():
        raise ValueError('Duplicate scenario-CAS entries would double count a compound')
    out = df.merge(crosswalk.drop(columns=['KR_source_table_row']), on='CAS', how='left', validate='many_to_one')
    if out['KR_TWA_ppm'].isna().any():
        raise ValueError(f'Missing Korean TWA CAS: {out.loc[out.KR_TWA_ppm.isna(), "CAS"].unique()}')
    if (out.vapor_mol_pct < 0).any() or out.vapor_mol_pct.isna().any():
        raise ValueError('Invalid source vapor mole fractions')
    group = out.groupby(KEYS, dropna=False)
    totals = group['vapor_mol_pct'].transform('sum')
    if not np.allclose(totals, 100., atol=1e-4):
        raise ValueError('Source vapor composition must sum to 100% for each scenario')
    out['KR_OEL_weighted_term'] = out['vapor_mol_pct'] / out['KR_TWA_ppm']
    sums = out.groupby(KEYS, dropna=False)['KR_OEL_weighted_term'].transform('sum')
    if (sums <= 0).any():
        raise ValueError('Invalid sum of OEL weighted terms')
    out['KR_OEL_weighted_share_pct'] = 100.0 * out['KR_OEL_weighted_term'] / sums
    out['KR_screening_rank'] = out.groupby(KEYS, dropna=False)['KR_OEL_weighted_share_pct'].rank(
        method='min', ascending=False).astype(int)
    if not np.allclose(out.groupby(KEYS)['KR_OEL_weighted_share_pct'].sum(), 100.0, atol=1e-6):
        raise AssertionError('KR weighted share total not 100%')
    out['delta_KR_minus_NIOSH_pct_points'] = out.KR_OEL_weighted_share_pct - out.OEL_weighted_share_pct
    out['regulatory_interpretation'] = 'screening_only_not_legal_omission_not_exposure_not_additive_toxicity'
    return out


def run(input_csv=DEFAULT_INPUT, oel_json=OEL_JSON, msds_json=MSDS_JSON, output_dir=None):
    df = pd.read_csv(input_csv, encoding='utf-8-sig', dtype={'product_id': str, 'target_CAS': str, 'CAS': str})
    oel_doc, oel_meta = regulatory_document(Path(oel_json))
    msds_doc, msds_meta = regulatory_document(Path(msds_json))
    annex6 = verify_msds_law(msds_doc)
    crosswalk = build_crosswalk(df, oel_doc)
    detailed = apply_oel(df, crosswalk)
    targets = detailed.loc[detailed.CAS == detailed.target_CAS].copy()
    if len(targets) != detailed[KEYS].drop_duplicates().shape[0]:
        raise AssertionError('Expected one target per scenario')
    out = Path(output_dir) if output_dir else Path(input_csv).parent / 'korean_regulatory_output'
    out.mkdir(parents=True, exist_ok=True)
    crosswalk.to_csv(out / '01_korean_oel_crosswalk.csv', index=False, encoding='utf-8-sig')
    detailed.to_csv(out / '02_korean_oel_component_comparison.csv', index=False, encoding='utf-8-sig')
    cols = KEYS + ['CAS', 'chemical', 'vapor_mol_pct', 'vapor_wt_pct',
                    'OEL_weighted_share_pct', 'KR_OEL_weighted_share_pct',
                    'screening_rank', 'KR_screening_rank', 'delta_KR_minus_NIOSH_pct_points']
    targets[cols].to_csv(out / '03_korean_oel_target_summary.csv', index=False, encoding='utf-8-sig')
    (out / '04_msds_annex6_verbatim.txt').write_text('\n'.join(annex6) + '\n', encoding='utf-8')
    metadata = {
        'OEL_notice': oel_meta.get('발령번호'),
        'OEL_effective_date': oel_meta.get('시행일자'),
        'MSDS_notice': msds_meta.get('발령번호'),
        'MSDS_effective_date': msds_meta.get('시행일자'),
        'OEL_metric': 'Korean time weighted average (TWA) ppm only',
        'OEL_source': str(oel_json), 'MSDS_source': str(msds_json),
        'rule_basis': '2026-26 Article 11(9) and Annex 6, plus enforcement rule Annex 18',
        'legal_omission_assessment': 'NOT_PERFORMED; per-CAS verified GHS categories and exceptions not established',
        'NIOSH_and_KR_comparison': 'Sensitivity of the same modeled vapor phase to the choice of weighting benchmark',
        'actual_exposure': False,
        'validity': 'Counterfactual mixtures; comparison index only, not aggregate mixture toxicity or compliance test',
        'successful_scenarios': len(targets),
        'complete_cas_count': len(crosswalk),
    }
    (out / '05_run_metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'Processed {len(detailed)} rows across {len(targets)} scenarios, {len(crosswalk)} CAS.')
    print('Korean notice references:',metadata['OEL_notice'],metadata['MSDS_notice'])
    print('OUTPUT', out)
    return crosswalk, targets


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Reweight existing thinner model output using local Korean regulation JSON')
    p.add_argument('--input', type=Path, default=DEFAULT_INPUT)
    p.add_argument('--oel-json', type=Path, default=OEL_JSON)
    p.add_argument('--msds-json', type=Path, default=MSDS_JSON)
    p.add_argument('--output-dir', type=Path)
    a = p.parse_args()
    run(a.input, a.oel_json, a.msds_json, a.output_dir)
