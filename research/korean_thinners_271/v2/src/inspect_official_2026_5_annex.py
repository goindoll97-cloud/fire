# -*- coding: utf-8 -*-
"""Fetch and audit the REAL NICS 2026-5 designation annex as an XLSX XML ZIP.

This is independent reference-document retrieval, not data-frame heuristics.
It records exact CAS-containing source rows for human/legal verification.
DOES NOT automatically infer concentration thresholds, effective-date exemptions,
or legal nonapplicability from absence of an exact CAS match.

Source NICS official bulletin (2026-03-27, effective 2026-07-01):
https://nics.mcee.go.kr/boardView.do?boardId=256&boardMasterId=4
Annex download as linked by the NICS bulletin:
https://nics.mcee.go.kr/common/userCmmnFileDown.do?fileId=9321&fileSeq=1
"""
from __future__ import annotations
import csv
import hashlib
import json
import re
import zipfile
from collections import defaultdict
from pathlib import Path
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET

import analyze_v2 as a2

SOURCE_URL = ('https://nics.mcee.go.kr/common/userCmmnFileDown.do'
              '?fileId=9321&fileSeq=1')
OUT = a2.RES / 'regulatory_gap' / 'official_2026_5'
COMPOSITION = (a2.V2.parent/'data'/'source_rebuild_20261023'
               /'MSDS_구성성분_1454.csv')
NS = {'x':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
CAS_RE = re.compile(r'(?<!\d)\d{2,7}-\d{2}-\d(?!\d)')


def excel_col(s):
    n=0
    for ch in s:
        if not ch.isalpha(): break
        n=n*26+ord(ch.upper())-ord('A')+1
    return n-1


def xml_xlsx_rows(raw: bytes):
    if raw[:2]!=b'PK':
        raise ValueError('Not a valid OOXML workbook ZIP')
    out=[]
    from io import BytesIO
    with zipfile.ZipFile(BytesIO(raw)) as z:
        shared=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            with z.open('xl/sharedStrings.xml') as f:
                root=ET.parse(f).getroot()
            for si in root.findall('x:si',NS):
                shared.append(''.join(t.text or '' for t in si.findall('.//x:t',NS)))
        sheets=sorted(n for n in z.namelist()
                      if re.match(r'xl/worksheets/sheet\d+\.xml$',n))
        if not sheets: raise ValueError('No worksheet data in official workbook')
        for sheet in sheets:
            with z.open(sheet) as f:
                root=ET.parse(f).getroot()
            count=0
            for row in root.findall('.//x:sheetData/x:row',NS):
                cells={}
                for c in row.findall('x:c',NS):
                    col=excel_col(c.get('r','A1'))
                    v=c.find('x:v',NS)
                    typ=c.get('t','')
                    inline=c.find('x:is',NS)
                    if typ=='s' and v is not None:
                        val=shared[int(v.text)]
                    elif inline is not None:
                        val=''.join(t.text or '' for t in inline.findall('.//x:t',NS))
                    else:
                        val=v.text if v is not None else ''
                    if val:cells[col]=val.strip()
                if cells:
                    items=[cells.get(i,'') for i in range(max(cells)+1)]
                    out.append({'sheet':sheet,'row':row.get('r',''),
                                'cells':items})
                    count+=1
            if count==0:raise ValueError(f'Empty official sheet {sheet}')
    return out


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    req=Request(SOURCE_URL,headers={'User-Agent':'Mozilla/5.0 (research source audit)'})
    with urlopen(req,timeout=35) as resp:
        raw=resp.read(5_000_000)
    if len(raw)<5000 or len(raw)>5_000_000:
        raise ValueError('Unexpected annex size; abort without approving source')
    rows=xml_xlsx_rows(raw)
    sha=hashlib.sha256(raw).hexdigest()
    annex=OUT/'NICS_2026_5_full_list_downloaded.xlsx'
    annex.write_bytes(raw)
    with COMPOSITION.open(encoding='utf-8-sig',newline='') as f:
        source_cas={x['CAS'] for x in csv.DictReader(f)}
    matched=[]
    for obj in rows:
        body=' | '.join(obj['cells'])
        found=set(CAS_RE.findall(body))&source_cas
        for cas in sorted(found):
            matched.append({
                'source_CAS':cas,'worksheet':obj['sheet'],
                'excel_row':obj['row'],'official_annex_row_raw':body,
                'law_verification_status':'source_CAS_matched_row_unverified_legal_threshold',
                'source_url':SOURCE_URL,'source_sha256':sha,
                'assumption_warning':'Exact CAS textual match only; aliases/group entries require independent review'
            })
    names=list(matched[0]) if matched else [
        'source_CAS','worksheet','excel_row','official_annex_row_raw',
        'law_verification_status','source_url','source_sha256','assumption_warning']
    with (OUT/'L3_annex_original_CAS_candidate_rows.csv').open(
        'w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=names)
        w.writeheader();w.writerows(matched)
    matched_cas={x['source_CAS'] for x in matched}
    ledger=[{
        'source_CAS':cas,
        'found_exact_CAS_in_2026_5_annex':cas in matched_cas,
        'legal_nonapplicability_inferred':False,
        'review_status':'source_text_row_candidate_only' if cas in matched_cas
                        else 'no_exact_CAS_match_not_an_exemption',
    } for cas in sorted(source_cas)]
    with (OUT/'L4_all_54_CAS_annex_text_search.csv').open(
        'w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(ledger[0]))
        w.writeheader();w.writerows(ledger)
    meta={
        'NICS_annex_source_url':SOURCE_URL,
        'NICS_2026_5_issued':'2026-03-27',
        'NICS_2026_5_effective':'2026-07-01',
        'file_sha256':sha,'xlsx_size_bytes':len(raw),
        'sheets':len({o['sheet'] for o in rows}),
        'nonblank_excel_rows':len(rows),'source_CAS_total':len(source_cas),
        'exact_CAS_matched':len(matched_cas),
        'matching_text_rows':len(matched),
        'matched_CAS_list':sorted(matched_cas),
        'must_not_infer_nonregulation_from_absence':True,
        'legal_thresholds_independently_approved':0
    }
    (OUT/'annex_ingest_meta.json').write_text(
        json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
    print('OFFICIAL_2026_5_ANNEX_INSPECT '+
          json.dumps(meta,ensure_ascii=False))


if __name__=='__main__':main()
