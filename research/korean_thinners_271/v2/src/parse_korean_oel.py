# -*- coding: utf-8 -*-
"""Parse 별표1 (화학물질의 노출기준) of the Korean MoEL notice JSON into a flat table.

Input : law_data/행정규칙/화학물질 및 물리적 인자의 노출기준.json (법제처 Open API dump in repo)
Output: v2/data/korean_oel_table.csv  (일련번호, 국문, 영문, TWA_ppm, TWA_mgm3, STEL_ppm, STEL_mgm3, CAS list, 비고)

The table is a box-drawing text table; one record spans several lines and ends at a
'├' or '└' separator line. Cells of continuation lines are concatenated.
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / 'law_data' / '행정규칙' / '화학물질 및 물리적 인자의 노출기준.json'
OUT = Path(__file__).resolve().parents[1] / 'data' / 'korean_oel_table.csv'
CAS_RE = re.compile(r'\b(\d{2,7}-\d{2}-\d)\b')


def _num(s: str):
    s = s.strip().replace(',', '')
    try:
        return float(s)
    except ValueError:
        return None


def parse(path: Path = SRC):
    d = json.loads(path.read_text(encoding='utf-8'))
    units = d['AdmRulService']['별표']['별표단위']
    unit = units if isinstance(units, list) else [units]
    lines = []
    for u in unit:
        if '화학물질의 노출기준' in u.get('별표제목', ''):
            for block in u['별표내용']:
                lines.extend(block if isinstance(block, list) else [block])
    records, cur = [], None
    for ln in lines:
        if not ln.startswith('│'):
            if cur and ln[:1] in '├└┝':
                records.append(cur)
                cur = None
            continue
        cells = [c for c in ln.strip().strip('│').split('│')]
        if len(cells) < 9:
            continue
        if cells[0].strip().isdigit():
            if cur:
                records.append(cur)
            cur = [c.strip() for c in cells[:9]]
        elif cur:
            for i in range(min(9, len(cells))):
                if cells[i].strip():
                    sep = '' if i in (1, 2) else ' '
                    cur[i] = (cur[i] + sep + cells[i].strip()).strip()
    if cur:
        records.append(cur)
    rows = []
    for r in records:
        no, ko, en, formula, twa_ppm, twa_mg, stel_ppm, stel_mg, note = r
        rows.append({'일련번호': int(no), '국문': ko, '영문': en, '화학식': formula,
                     'TWA_ppm': _num(twa_ppm), 'TWA_mgm3': _num(twa_mg),
                     'STEL_ppm': _num(stel_ppm), 'STEL_mgm3': _num(stel_mg),
                     'CAS': '|'.join(dict.fromkeys(CAS_RE.findall(note))), '비고': note})
    return rows


def main():
    rows = parse()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(f'{len(rows)} OEL records -> {OUT}')


if __name__ == '__main__':
    main()
