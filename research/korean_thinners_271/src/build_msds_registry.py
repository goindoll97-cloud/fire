# -*- coding: utf-8 -*-
"""Evidence-traceable MSDS research registry (NOT automatic chemistry verification).

Usage (local Windows or Linux):
 python build_msds_registry.py --noru-zip PATH --samhwa-zip PATH --output PATH

Prerequisite: pip install pymupdf. Optional --legacy-index matches SHA256 to existing 271 snapshot.
This ingests the two user-supplied ZIP files only, without publishing their PDF bytes.
All extracted records carry source hashes, page and review status. A parsed composition
range does not mean the ingredient list is an exhaustive formulation.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import json
import re
import sqlite3
import statistics
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import fitz

CAS_RE = re.compile(r'(?<!\d)(\d{2,7}-\d{2}-\d)(?!\d)')
SECTION_PATTERNS = {
    'NORU': re.compile(r'(?mi)^\s*3\s*[.．]\s*(?:구성\s*성분|구성성분|화학물질의 구성)'),
    'SAMHWA': re.compile(r'(?mi)^\s*3\s*[.．]\s*(?:CHEMICAL\s*COMPOSITION|COMPOSITION)'),
}
NEXT_SECTION = re.compile(r'(?mi)^\s*4\s*[.．]\s*(?:응급|HEALTH|FIRST|EMERGENCY|대응)')
DATE_RE = re.compile(r'(?<!\d)(?:19|20)\d{2}[-./]\d{1,2}[-./]\d{1,2}(?!\d)')
# Chemical identity groups: never silently represent petroleum/UVCB/isomer mixtures
# as a single structure; all other CAS are only *candidates* pending human identity QA.
MIXTURE_CAS = {
    '8052-41-3': 'Stoddard solvent; UVCB petroleum mixture (ATSDR)',
    '26588-80-7': 'acrylic copolymer; no single molecule (NITE J-CHECK)',
    '1330-20-7': 'xylene mixed isomers',
    '25551-13-7': 'trimethylbenzene isomer mixture',
    '31394-54-4': 'isoheptane isomer mixture',
    '64742-95-6': 'light aromatic solvent naphtha (UVCB)',
    '64742-94-5': 'heavy aromatic solvent naphtha (UVCB)',
    '64742-82-1': 'hydrodesulfurized heavy naphtha (UVCB)',
    '64742-47-8': 'hydrotreated light distillates (UVCB)',
    '64742-49-0': 'hydrotreated light naphtha (UVCB)',
    '64741-92-0': 'solvent-refined heavy naphtha (UVCB)',
    '64742-48-9': 'hydrotreated heavy naphtha (UVCB)',
    '64741-65-7': 'heavy alkylate naphtha (UVCB)',
}
PILOT_CAS = {'71-43-2', '67-64-1', '141-78-6', '108-88-3', '67-56-1', '123-86-4', '108-10-1'}


def csv_rows(path):
    with Path(path).open('r', encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, cols):
    with Path(path).open('w', encoding='utf-8-sig', newline='') as f:
        writer=csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)


def valid_cas(cas):
    try:
        body, digit = cas.rsplit('-', 1)
        code = body.replace('-', '')
        return sum((i + 1) * int(v) for i, v in enumerate(reversed(code))) % 10 == int(digit)
    except (ValueError, AttributeError):
        return False


def parse_interval(tail):
    """Only parse immediately following CAS; no imputed or 'typical' compositions."""
    s = re.sub(r'\s+', ' ', tail[:165]).strip()
    # Noru PDFs list the Korean inventory identifier immediately after CAS.
    s = re.sub(r'^/\s*(?:KE|KECL|KECI)\s*-\s*\d+\s*', '', s, flags=re.I)
    s = s.replace('～', '~').replace('∼', '~').replace('〜', '~').replace('–', '~')
    # A normalized pair of numbers with explicit tilde and recognized range words.
    m = re.match(r'(?i)^(?:More\s+than\s+|At\s+least\s+|over\s+)?'
                 r'(\d+(?:\.\d+)?)\s*(?:이상|초과|초과\s*~|이하)?\s*~\s*'
                 r'(?:(?:Under|Less\s+than|Below|No\s+more\s+than)\s*)?'
                 r'(\d+(?:\.\d+)?)\s*(?:%|이하|미만|이상|초과)?', s)
    if not m:
        # sometimes newline leads punctuation mixed with English: "More than 0.1 ~ Under 1%"
        return '', '', 'REVIEW_CONCENTRATION', s[:115]
    lo=float(m.group(1)); hi=float(m.group(2))
    if not (0 <= lo <= hi <= 100):
        return '', '', 'INVALID_RANGE', s[:115]
    # Both bounds may be open; retain raw expression and never equate endpoints
    # with guaranteed included manufacturer percentages.
    return lo, hi, 'RANGE_EXTRACTED_UNVERIFIED', s[:115]


def identify_uses(name):
    groups = {
      '세척·세정': ['세척','세정','세척용','cleaning','cleaner'],
      '우레탄': ['우레탄','urethane','pu ', 'p_u'],
      '에폭시': ['에폭시','epoxy'],
      '락카': ['락카','lacquer'],
      '에나멜': ['에나멜','enamel'],
      '자동차보수': ['자동차','차량','automotive','refinish'],
      '건조속도': ['속건','지건','초속건','초지건','retarder'],
      '계절': ['동절','하절','summer','winter'],
    }
    lower=name.casefold()
    matches=[g for g,keys in groups.items() if any(k.casefold() in lower for k in keys)]
    uses=[x for x in matches if x not in ('건조속도','계절')]
    return '|'.join(uses) if uses else '용도미확인', '|'.join(matches)


def load_legacy(path):
    if not path or not Path(path).is_file(): return {}
    rows=csv_rows(path)
    return {r['파일 SHA256'].strip():r for r in rows if r.get('파일 SHA256','').strip()}


def index_pdfs(archives, outdir, legacy=None):
    docs=[];ingredients=[];evidence=[]
    for brand, archive in archives.items():
        with zipfile.ZipFile(archive) as z:
            entries=sorted((x for x in z.infolist() if x.filename.lower().endswith('.pdf')),key=lambda x:x.filename.casefold())
            for entry in entries:
                raw=z.read(entry);sha=hashlib.sha256(raw).hexdigest()
                doc=fitz.open(stream=raw,filetype='pdf')
                # Section 3 is found on first pages for this particular corpus; never silently
                # process scans or missing tables as empty products.
                pages=[doc[i].get_text() for i in range(min(5,len(doc)))]
                text='\n'.join(pages)
                pattern=SECTION_PATTERNS[brand]
                section_starts=list(pattern.finditer(text))
                sec=section_starts[0] if section_starts else None
                # Some Samhwa PDFs reprint a '3.CHEMICAL COMPOSITION' continuation
                # after Section 4 material on the next page. Extract every explicit
                # Section 3 occurrence but never sweep Section 8 exposure/OEL values.
                section_parts=[]
                for si,msection in enumerate(section_starts):
                    cap=(section_starts[si+1].start() if si+1<len(section_starts) else len(text))
                    block=text[msection.end():cap]
                    end=NEXT_SECTION.search(block)
                    if end:block=block[:end.start()]
                    section_parts.append(block[:20000])
                secbody='\n'.join(section_parts)
                msds=re.search(r'AA\d{5}-\d{7,12}',text)
                product_filename=Path(entry.filename).name.rsplit('.',1)[0]
                if '_' in product_filename:product_filename=product_filename.split('_',1)[1]
                product_filename=re.sub(r'_GHS영문$', '', product_filename)
                direct, tags = identify_uses(product_filename)
                first=text[:min(len(text),4000)]
                issue_dates=DATE_RE.findall(first)
                dkey='D'+sha[:20]
                legacy_item=legacy.get(sha,{}) if legacy else {}
                docstatus='SECTION3_REVIEW_PENDING' if sec else 'SECTION3_NOT_FOUND'
                if len(doc) == 0 or sum(len(x) for x in pages)<200:
                    docstatus='TEXT_EXTRACTION_FAILED'
                matches=list(CAS_RE.finditer(secbody))
                caslist=[]
                for i,m in enumerate(matches):
                    cas=m.group(1)
                    if cas in caslist:continue
                    caslist.append(cas)
                    later=matches[i+1].start() if i+1<len(matches) else len(secbody)
                    tail=secbody[m.end():min(later,m.end()+165)]
                    lo,hi,flag,rawrange=parse_interval(tail)
                    before=secbody[max(0,m.start()-220):m.start()]
                    # evidence local only: raw text excerpt with page/adjacent synonyms
                    snippet=secbody[max(0,m.start()-170):min(len(secbody),m.end()+160)]
                    page_number=0
                    # best effort page location based on CAS detection per page; fallback unlocated
                    for pno,ptext in enumerate(pages,1):
                        if cas in ptext and ('3.' in ptext or pno>1):
                            page_number=pno;break
                    rec={
                     '문서키':dkey,'기존제품ID':legacy_item.get('제품ID',''),
                     '원본SHA256':sha,'제조사':brand,'CAS':cas,'CAS검증':('형식유효' if valid_cas(cas) else '체크디지트불일치'),
                     '함량하한_pct':lo,'함량상한_pct':hi,'함량추출상태':flag,
                     'CAS분류':('복합물/이성질체혼합물' if cas in MIXTURE_CAS else '순물질후보_미검증'),
                     '분자구조검증':'미검증','포화증기압검증':'미검증','국내TWA검증':'미검증',
                     '원문쪽수_잠정':page_number,
                     '연구판정':'복합물구조동정필요' if cas in MIXTURE_CAS else '수동동정/물성검증필요',
                    }
                    ingredients.append(rec)
                    evidence.append({'문서키':dkey,'CAS':cas,'상대위치':m.start(),
                        '원문쪽수_잠정':page_number,'범위원문':rawrange,'원문문맥':snippet.replace('\x00',''),
                        '추출확인상태':'미검수'})
                mix=[c for c in caslist if c in MIXTURE_CAS]
                intervalrows=[r for r in ingredients if r['문서키']==dkey]
                complete=bool(intervalrows) and all(r['함량추출상태']=='RANGE_EXTRACTED_UNVERIFIED' for r in intervalrows)
                lowtotal=sum(float(r['함량하한_pct']) for r in intervalrows if r['함량하한_pct']!='')
                hitotal=sum(float(r['함량상한_pct']) for r in intervalrows if r['함량상한_pct']!='')
                feasible=complete and lowtotal<=100+1e-7 and hitotal>=100-1e-7
                if not matches:docstatus='SECTION3_NO_VALID_CAS'
                if not complete and matches:docstatus='INTERVAL_MANUAL_REVIEW'
                docs.append({
                   '문서키':dkey,'기존제품ID':legacy_item.get('제품ID',''),
                   '제조사':brand,'원본PDF파일명':entry.filename.split('/')[-1],
                   '원본SHA256':sha,'PDF쪽수':len(doc),'제품명_파일명기반':product_filename,
                   '용도태그_잠정':direct,'전체태그_잠정':tags,'제품군검증':'미검증',
                   '원문MSDS번호_잠정':msds.group() if msds else '',
                   '1~5쪽날짜_참고': '|'.join(dict.fromkeys(issue_dates[:4])),
                   '제3항찾음':'예' if sec else '아니오', '추출CAS종수':len(caslist),
                   '추출하한합계_pct':round(lowtotal,5),'추출상한합계_pct':round(hitotal,5),
                   '함량합계100가능_표기성분만':'예' if feasible else '확인필요',
                   '복합물성분종수':len(mix),
                   '기존7종물성만으로_잠정구성충족':'예' if feasible and len(caslist)>=2 and set(caslist).issubset(PILOT_CAS) else '아니오',
                   '구조미동정복합물있음':'예' if mix else '아니오',
                   '검수상태':docstatus,
                   '한계':'MSDS_공개함량범위≠제조사정확배합/미공개성분존재가능',
                })
    return docs,ingredients,evidence


def make_cas_master(documents, components):
    groups=collections.defaultdict(list)
    for c in components:groups[c['CAS']].append(c)
    result=[]
    for cas,recs in sorted(groups.items(),key=lambda kv:(-len(kv[1]),kv[0])):
        result.append({'CAS':cas,'포함PDF수':len(set(r['문서키'] for r in recs)),
           '유효형식':recs[0]['CAS검증'],
           '성분분류':recs[0]['CAS분류'],
           '복합물또는이성질체':MIXTURE_CAS.get(cas,''),
           '기존파일럿7종포함':'예' if cas in PILOT_CAS else '아니오',
           '원문범위추출완료문서수':len({r['문서키'] for r in recs if r['함량추출상태']=='RANGE_EXTRACTED_UNVERIFIED'}),
           'SMILES_인증':'미검증', 'MW_인증':'미검증', 'VP_온도출처_인증':'미검증',
           '현재국내TWA_인증':'미검증','수집해야할원문근거':'물성측정/공신력있는DB/공식국내노출기준',
           '모델입력승인':'미승인'})
    return result


def add_sqlite(path, datasets):
    with sqlite3.connect(path) as con:
        con.execute('PRAGMA foreign_keys=ON')
        for name, (rows, fields) in datasets.items():
            con.execute(f'DROP TABLE IF EXISTS "{name}"')
            col_sql=','.join(f'"{f}" TEXT' for f in fields)
            con.execute(f'CREATE TABLE "{name}" ({col_sql})')
            if rows:
                binds=','.join('?' for _ in fields)
                con.executemany(f'INSERT INTO "{name}" VALUES ({binds})', [[str(r.get(f,'')) for f in fields] for r in rows])
        con.execute('CREATE INDEX IF NOT EXISTS idx_components_cas ON 성분근거(CAS)')
        con.execute('CREATE INDEX IF NOT EXISTS idx_components_doc ON 성분근거(문서키)')


def main(argv=None):
    ap=argparse.ArgumentParser()
    ap.add_argument('--noru-zip',type=Path,required=True)
    ap.add_argument('--samhwa-zip',type=Path,required=True)
    ap.add_argument('--legacy-index',type=Path,help='Optional existing research/.../data/제품목록_271.csv')
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args(argv)
    args.output.mkdir(parents=True,exist_ok=True)
    default_legacy=Path(__file__).resolve().parents[1]/'data'/'제품목록_271.csv'
    legacy_path=args.legacy_index or (default_legacy if default_legacy.is_file() else None)
    legacy=load_legacy(legacy_path)
    docs,components,evidence=index_pdfs({'NORU':args.noru_zip,'SAMHWA':args.samhwa_zip},args.output,legacy)
    master=make_cas_master(docs,components)
    unresolved=[r for r in docs if r['검수상태']!='SECTION3_REVIEW_PENDING' or r['복합물성분종수']]
    ds={'문서원본':docs,'성분근거':components,'CAS기본마스터':master}
    for n,rr in ds.items():write_csv(args.output/(n+'.csv'),rr,list(rr[0]) if rr else [])
    write_csv(args.output/'검수우선문서.csv',unresolved,list(docs[0]))
    if legacy_path:
        known={d['원본SHA256']:d for d in docs}
        reconciliation=[]
        for r in csv_rows(legacy_path):
            sha=r.get('파일 SHA256','').strip()
            current=known.get(sha)
            reconciliation.append({
                '기존제품ID':r.get('제품ID',''),
                '기존원본파일명':r.get('원본 PDF 파일명',''),
                '원본SHA256':sha,
                '이번ZIP에포함':'예' if current else '아니오',
                '이번추출제3항CAS수':current['추출CAS종수'] if current else '',
                '기존성분기록수':r.get('성분 기록수',''),
                '성분수일치':'예' if current and str(current['추출CAS종수'])==str(r.get('성분 기록수','')) else '확인필요',
                '검수대상':'ZIP미포함_기존원본확인' if not current else ('정상대조' if str(current['추출CAS종수'])==str(r.get('성분 기록수','')) else '원문성분수차이'),
            })
        write_csv(args.output/'기존271_vs_원본ZIP_재대조.csv',reconciliation,list(reconciliation[0]))

    write_csv(args.output/'원문근거_로컬전용.csv',evidence,list(evidence[0]))
    add_sqlite(args.output/'국내신나_연구DB.sqlite',{n:(rr,list(rr[0])) for n,rr in ds.items()})
    cc=collections.Counter(d['제조사'] for d in docs)
    obs_count=sum(1 for d in docs if d['검수상태']=='SECTION3_REVIEW_PENDING')
    feasible=sum(1 for d in docs if d['함량합계100가능_표기성분만']=='예')
    candidate=sum(1 for d in docs if d['기존7종물성만으로_잠정구성충족']=='예')
    mixtures=sum(1 for d in docs if d['복합물성분종수']>0)
    info={
      'source':'supplied ZIP PDFs (original bytes not uploaded to Github)',
      'timestamp_utc':datetime.now(timezone.utc).isoformat(),
      'noru_zip_sha256':hashlib.sha256(args.noru_zip.read_bytes()).hexdigest(),
      'samhwa_zip_sha256':hashlib.sha256(args.samhwa_zip.read_bytes()).hexdigest(),
      'pdf_count':len(docs),'brand_counts':dict(cc),'distinct_cas_regex_valid':len(master),
      'pdf_section3_with_cas':sum(d['추출CAS종수']>0 for d in docs),
      'component_cas_document_records':len(components),
      'pdf_all_numeric_bounds_parsed':obs_count,
      'pdf_sum100_compatible_reported_components':feasible,
      'pdf_with_named_isomer_or_uvcb_mixture':mixtures,
      'pdf_with_only_pilot7_reference_cas_and_feasible_intervals':candidate,
      'manual_review_needed_docs':len(unresolved),
      'legacy_matched_exact_sha256':sum(bool(d['기존제품ID']) for d in docs),
      'limitations':['Section 3 values extracted from PDFs; no chemical/properties human approval claimed',
          'Private evidence file contains source snippets and should not be published automatically',
          'Valid MSDS composition intervals need not reflect actual formulation',
          'Some undisclosed ingredients may be missing from MSDS; no exposure rate derived',
          'CAS regex result is not complete ontology of chemical entities',
          'No imputation for UVCB, mixed-isomers, unknown structures or Korean OELs',
          'Only reference 7 CAS have prior pilot physical properties; all source validations pending']}
    (args.output/'원본대조_요약.json').write_text(json.dumps(info,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(info,ensure_ascii=False,indent=2))
    return info
if __name__=='__main__':main()