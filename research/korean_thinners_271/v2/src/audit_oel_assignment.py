# -*- coding: utf-8 -*-
"""CAS-to-OEL integrity check against the parsed CURRENT-LISTED 2020 notice.

This is an INTERNAL mapping audit (not an independent legal or experimental
verification of every Annex 1 row). The source notice must be checked against
the Ministry of Government Legislation before journal submission.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import build_oel as build

V2 = Path(__file__).resolve().parents[1]
OUT = V2 / 'results' / 'core'


def audit():
    rebuilt = {x['CAS']: x for x in build.build()}
    with build.OUT.open(encoding='utf-8-sig', newline='') as f:
        recorded = {x['CAS']: x for x in csv.DictReader(f)}
    if set(recorded) != set(rebuilt):
        raise ValueError('CAS rows differ between saved and rebuilt OEL tables')
    out = []
    for cas in sorted(recorded):
        saved = recorded[cas]
        new = rebuilt[cas]
        for prefix in ('TierA', 'TierB'):
            k = prefix + '_TWA_ppm'
            # CSV serializes floats; preserve missing data as missing, not zero.
            v1 = saved[k].strip()
            v2 = '' if new[k] == '' else str(float(new[k]))
            if v1 != v2:
                raise ValueError(f'OEL value differs for CAS {cas}: {prefix}: {v1!r} vs {v2!r}')
        if not saved['TierA_TWA_ppm']:
            type_ = 'foreign_or_analogue_TierB' if saved['TierB_TWA_ppm'] else 'uncovered'
        elif 'isomer-group' in saved['TierA_source']:
            type_ = 'official_isomer_group'
        else:
            type_ = 'official_exact_CAS'
        out.append({
            'CAS': cas, 'name': saved['name'], 'OEL_basis': type_,
            'TierA_TWA_ppm': saved['TierA_TWA_ppm'],
            'TierB_TWA_ppm': saved['TierB_TWA_ppm'],
            'TierA_source': saved['TierA_source'],
            'TierB_source': saved['TierB_source'],
            'source_status': '2020-48 in official current rule list; rows not independently re-extracted',
        })
    return out


def main():
    out = audit()
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / 'S4_oel_cas_mapping_audit.csv').open('w', newline='', encoding='utf-8-sig') as f:
        wr = csv.DictWriter(f, fieldnames=list(out[0]))
        wr.writeheader()
        wr.writerows(out)
    counts = {x: sum(r['OEL_basis'] == x for r in out) for x in sorted(
        set(r['OEL_basis'] for r in out))}
    print('NONEXPERIMENTAL_OEL_SOURCE_AUDIT ' + json.dumps({
        'total_CAS': len(out), 'assignment_matches_rebuilt_table': True,
        'status': 'Internal mapping check only, NOT independent revalidation of 2020 Annex',
        'counts': counts}, ensure_ascii=False))


if __name__ == '__main__':
    main()
