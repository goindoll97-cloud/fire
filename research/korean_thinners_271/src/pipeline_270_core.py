# -*- coding: utf-8 -*-
"""Evidence-gated Korean thinner MSDS research pipeline (270 original PDFs).

None of the stages asserts that a CAS, physical property, exposure limit,
actual product compatibility, or measured risk is verified without human evidence.
Run numbered 00--09 wrappers from the repository root; generated files are ignored.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import random
import re
import shutil
import statistics
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[3]
RESEARCH = ROOT / 'research' / 'korean_thinners_271'
SOURCE = RESEARCH / 'data' / 'source_rebuild_20261023'
WORK = RESEARCH / 'output' / 'full270'
DOCFILE = 'MSDS_원본대조_270.csv'
COMPFILE = 'MSDS_구성성분_1454.csv'
MASTERFILE = 'CAS_54종_출처상태.csv'
ATM = 101.325
PRODUCT_FIELDS = [
    '문서키', '기존제품ID', '제조사', '제품명', '용도ID', '비교군ID',
    'MSDS제1항용도검증', 'MSDS제3항함량검증', '전체모델성분확정',
    '중복개정판검토', '독립제품검증', '공정수지호환성검증', 'TDS성능제약검증', '근거문서'
]
PROPERTY_FIELDS = [
    'CAS', '복합물또는이성질체', 'SMILES', '분자량_g_mol', '25C_포화증기압_kPa',
    '온도_K', '구조공식출처', '물성출처와측정조건', '구조검증상태',
    '물성검증상태', 'HANNA입력승인'
]
TWA_FIELDS = [
    'CAS', '한국_TWA_ppm', '기준시점', '고시명번호', '공식출처',
    '공식원문대조', '검증상태'
]
PRED_SCEN = [
    '문서키', '기존제품ID', '제품명', '시나리오ID', '온도_K', 'Raoult_총압_kPa',
    'HANNA_총압_kPa', '압력비_HANNA_Raoult', '모델명', 'seed', 'sampler'
]
PRED_PART = [
    '문서키', '시나리오ID', 'CAS', '질량_pct', '액상몰분율', 'HANNA_gamma',
    'Raoult_분압_kPa', 'HANNA_분압_kPa', 'Raoult_기상몰분율', 'HANNA_기상몰분율'
]


def read(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def write(path, rows, columns):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def jsonfile(path, obj):
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding='utf-8')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def csv_index(path, key):
    ans = {}
    for row in read(path):
        name = row.get(key, '').strip()
        if not name or name in ans:
            raise ValueError(f'{path}: 빈 값 또는 중복 {key}: {name!r}')
        ans[name] = row
    return ans


def ensure_source():
    for name in (DOCFILE, COMPFILE, MASTERFILE):
        if not (SOURCE / name).is_file():
            raise FileNotFoundError(f'원본 DB 누락: {SOURCE/name} (git pull origin main)')
    docs, parts, master = read(SOURCE/DOCFILE), read(SOURCE/COMPFILE), read(SOURCE/MASTERFILE)
    if len(docs) != 270 or len(parts) != 1454 or len(master) != 54:
        raise ValueError(f'원본 자료 무결성 오류: {len(docs)} PDF, {len(parts)} 성분, {len(master)} CAS (기대 270/1454/54)')
    dids = {d['문서키'] for d in docs}
    if len(dids) != 270 or len({d['원본SHA256'] for d in docs}) != 270:
        raise ValueError('문서키/원본 SHA-256 중복 또는 누락')
    if any(c['문서키'] not in dids for c in parts):
        raise ValueError('원본 PDF에 연결되지 않는 성분기록 존재')
    if {c['CAS'] for c in parts} != {m['CAS'] for m in master}:
        raise ValueError('성분 CAS와 CAS 마스터 불일치')
    if any(d['기존제품ID'] == 'NR007' for d in docs):
        raise ValueError('비신나 NR007이 연구 DB에 포함됨')
    return docs, parts, master


def load_env():
    """Load only the root .env; existing shell environment takes precedence.
    Never print secrets or construct authenticated requests on a failure path.
    """
    p = ROOT/'.env'
    if not p.is_file(): return
    for line in p.read_text(encoding='utf-8-sig').splitlines():
        t = line.strip()
        if not t or t.startswith('#'): continue
        if t.startswith('export '): t=t[7:].strip()
        if '=' not in t: continue
        k, v = t.split('=', 1)
        k = k.strip()
        if not re.fullmatch(r'[A-Z][A-Z0-9_]*', k): continue
        v = v.strip()
        if len(v)>1 and v[0] == v[-1] and v[0] in "\"'": v=v[1:-1]
        os.environ.setdefault(k, v)


def require_previous(stage):
    if stage != '00' and not (WORK/'01_source'/'source_manifest.json').exists() and stage != '01':
        raise RuntimeError('먼저 01 원본 대조를 실행하세요.')


def stage00(a):
    template = ROOT/'.env.example'
    local = ROOT/'.env'
    if not local.exists():
        if not template.exists(): raise FileNotFoundError('.env.example is missing')
        # Explicitly local: ignored by Git, empty keys; never overwrite existing .env
        shutil.copyfile(template, local)
        try: local.chmod(0o600)
        except OSError: pass
        print('로컬 .env 생성됨 (빈 값만 포함). 필요한 값을 로컬에서 입력하세요.')
    load_env()
    summary = {'local_env_exists':local.exists(),
               'LAW_OC_configured':bool(os.environ.get('LAW_OC')),
               'KOSHA_MSDS_SERVICE_KEY_configured':bool(os.environ.get('KOSHA_MSDS_SERVICE_KEY')),
               'HANNA_HOME_configured':bool(os.environ.get('HANNA_HOME')),
               'security':'Values never printed, .env ignored, no automatic source verification'}
    jsonfile(WORK/'00_environment'/'status.json',summary)
    print(json.dumps(summary,ensure_ascii=False,indent=2))


def stage01(a):
    docs, parts, master = ensure_source()
    # Cross-check supplied raw PDF ZIP file bytes against committed per-file SHA-256.
    args = [a.noru_zip, a.samhwa_zip]
    if any(args) and not all(args):
        raise ValueError('--noru-zip and --samhwa-zip must be provided together')
    zips = []
    if all(args):
        seen=[]
        for archive, brand in ((a.noru_zip, 'NORU'),(a.samhwa_zip,'SAMHWA')):
            if not Path(archive).is_file(): raise FileNotFoundError(archive)
            with zipfile.ZipFile(archive) as z:
                pdfs=[n for n in z.namelist() if n.lower().endswith('.pdf') and not n.startswith('__MACOSX')]
                for entry in pdfs:
                    raw = z.read(entry)
                    seen.append((brand,hashlib.sha256(raw).hexdigest()))
            zips.append({'brand':brand,'archive_sha256':sha(archive),'pdf_count':len(pdfs)})
        truth={(d['제조사'],d['원본SHA256']) for d in docs}
        if set(seen) != truth or len(seen)!=270:
            raise ValueError('ZIP 원본 SHA와 정제 270건 레지스트리가 불일치. PDF 수/변경 내역 확인 필요')
        # Optional full raw PDF re-extraction via the earlier tested Section 3 parser.
        # Keep extracted evidence local; never commit it to GitHub.
        from build_msds_registry import index_pdfs, load_legacy
        legacy_path=RESEARCH/'data'/'제품목록_271.csv'
        legacy=load_legacy(legacy_path if legacy_path.is_file() else None)
        raw_docs,raw_parts,raw_evidence=index_pdfs({'NORU':a.noru_zip,'SAMHWA':a.samhwa_zip},WORK/'01_source',legacy)
        mapping={d['문서키']:d['원본SHA256'] for d in raw_docs}
        if len(raw_docs)!=270 or len(raw_parts)!=1454:
            raise ValueError('원본 PDF 재추출 수가 GitHub 대조 DB와 다릅니다. 수동 비교 필요')
        def normalize_comp(rows, keyfunc, casfunc, lofunc, hifunc):
            result={}
            for r in rows:
                key=(keyfunc(r),casfunc(r))
                if key in result:raise ValueError('하나의 PDF에 CAS 중복 존재: '+str(key))
                result[key]=(float(lofunc(r)),float(hifunc(r)))
            return result
        committed_doc={d['문서키']:d['원본SHA256'] for d in docs}
        source_cmp=normalize_comp(parts,lambda r:committed_doc[r['문서키']],lambda r:r['CAS'],
                                 lambda r:r['함량하한_pct'],lambda r:r['함량상한_pct'])
        raw_cmp=normalize_comp(raw_parts,lambda r:mapping[r['문서키']],lambda r:r['CAS'],
                              lambda r:r['함량하한_pct'],lambda r:r['함량상한_pct'])
        if source_cmp.keys()!=raw_cmp.keys() or any(abs(source_cmp[k][j]-raw_cmp[k][j])>1e-6
                  for k in source_cmp for j in (0,1)):
            raise ValueError('원본 PDF 재추출과 GitHub 구성성분/함량 DB 불일치: 원문 검수 필요')
        private=WORK/'01_source'/'local_pdf_evidence'
        private.mkdir(parents=True,exist_ok=True)
        write(private/'원문근거_로컬전용.csv',raw_evidence,list(raw_evidence[0]))
        jsonfile(private/'reextraction_audit.json',{'pdf':len(raw_docs),'compositions':len(raw_parts),
                 'sha_and_interval_match':'YES, AUTOMATED ONLY; not human-approved'})
    folder=WORK/'01_source';folder.mkdir(parents=True,exist_ok=True)
    out={
        'document_count':len(docs),'component_rows':len(parts),'cas_count':len(master),
        'source_state':'SOURCE_PDF_SHA_VERIFIED' if all(args) else 'GITHUB_REBUILD_SNAPSHOT_ONLY',
        'zip_verification':zips, 'source_sha256':{x:sha(SOURCE/x) for x in (DOCFILE,COMPFILE,MASTERFILE)},
        'pdf_not_shipped':True, 'extracted_content_review':'NOT HUMAN APPROVED',
        'created_at_utc':datetime.now(timezone.utc).isoformat(),
    }
    jsonfile(folder/'source_manifest.json',out)
    print('01:',out['source_state'],len(docs),'documents',len(parts),'components',len(master),'CAS')


def stage02(a):
    docs, parts, master=ensure_source()
    dst=WORK/'02_review';dst.mkdir(parents=True,exist_ok=True)
    reviews=[]
    for d in docs:
        reviews.append({'문서키':d['문서키'],'기존제품ID':d['기존제품ID'],'제조사':d['제조사'],
                        '제품명':d['제품명_파일명기반'],'용도ID':d['용도태그_잠정'],
                        '비교군ID':'',
                        'MSDS제1항용도검증':'미검증','MSDS제3항함량검증':'미검증',
                        '전체모델성분확정':'미검증','중복개정판검토':'미검증',
                        '독립제품검증':'미검증','공정수지호환성검증':'미검증',
                        'TDS성능제약검증':'미검증','근거문서':''})
    props=[]
    twa=[]
    for row in master:
        cas=row['CAS']
        props.append({'CAS':cas,'복합물또는이성질체':'예' if row['복합물또는이성질체'].strip() else '아니오',
                      'SMILES':'','분자량_g_mol':'','25C_포화증기압_kPa':'',
                      '온도_K':'298.15','구조공식출처':'','물성출처와측정조건':'',
                      '구조검증상태':'미검증','물성검증상태':'미검증','HANNA입력승인':'미승인'})
        twa.append({'CAS':cas,'한국_TWA_ppm':'','기준시점':'','고시명번호':'',
                    '공식출처':'','공식원문대조':'미검증','검증상태':'미검증'})
    for name,rows,fields in (
        ('제품_원문TDS_검증.csv',reviews,PRODUCT_FIELDS),
        ('CAS_물성_검증.csv',props,PROPERTY_FIELDS),
        ('CAS_한국TWA_검증.csv',twa,TWA_FIELDS)):
        path=dst/name
        if not path.exists(): write(path,rows,fields)
        else:
            current=read(path)
            if {r['CAS' if name.startswith('CAS') else '문서키'] for r in current} != {r['CAS' if name.startswith('CAS') else '문서키'] for r in rows}:
                raise ValueError(f'기존 검수 파일의 대상 키가 원본과 다릅니다: {path} (보존 위해 중단)')
    bydoc=defaultdict(list)
    for c in parts: bydoc[c['문서키']].append(c)
    statuses=[]
    for d in docs:
        cs=bydoc[d['문서키']]; mix=[c['CAS'] for c in cs if c['CAS분류']=='복합물/이성질체혼합물']
        why=[]
        if d['함량합계100가능_표기성분만']!='예':why.append('함량합계100_불가능')
        if mix:why.append('UVCB_혼합이성질체:'+','.join(mix))
        if any(c['함량추출상태']!='RANGE_EXTRACTED_UNVERIFIED' for c in cs):why.append('성분함량자동추출검토')
        statuses.append({'문서키':d['문서키'],'기존제품ID':d['기존제품ID'],
                         '용도_잠정':d['용도태그_잠정'], '성분기록수':len(cs),
                         '모델링_원자료_보류사유':' ; '.join(why) if why else '원문/물성승인대기'})
    write(dst/'270건_성분및용도_QA.csv',statuses,list(statuses[0]))
    print('02: 270 product review rows + 54 property/TWA rows (all unapproved by default)')


def http_xml(url, params):
    # URL may include secrets; never log request URL or exception containing query string.
    query=urllib.parse.urlencode(params)
    request=urllib.request.Request(url+'?'+query,headers={'User-Agent':'thinner-msds-research/1.0'})
    try:
        with urllib.request.urlopen(request,timeout=20) as response:
            raw=response.read(3_000_001)
    except Exception as exc:
        raise RuntimeError(f'기관 API 요청 실패: {type(exc).__name__} (인증값 미노출)') from None
    if len(raw)>3_000_000:raise ValueError('기관 API 응답이 제한 크기를 초과했습니다')
    try: return ET.fromstring(raw)
    except ET.ParseError: raise RuntimeError('기관 API가 XML을 반환하지 않았습니다. 인증·응답형식 확인 필요') from None


def stage03(a):
    load_env();out=WORK/'03_official';out.mkdir(parents=True,exist_ok=True)
    lawkey=os.environ.get('LAW_OC','').strip()
    koskey=os.environ.get('KOSHA_MSDS_SERVICE_KEY','').strip()
    # data.go.kr offers both URL-encoded and plain keys. urlencode() below
    # must encode exactly once (as proven by the working single-CAS probe).
    if re.search(r'%[0-9A-Fa-f]{2}', koskey):
        koskey=urllib.parse.unquote(koskey)
    rows=[]
    if lawkey and not a.offline:
        for query in ('화학물질 및 물리적 인자의 노출기준','산업안전보건법'):
            target='admrul' if '노출기준' in query else 'eflaw'
            el=http_xml('https://www.law.go.kr/DRF/lawSearch.do',
                        {'OC':lawkey,'target':target,'query':query,'type':'XML','display':10,'page':1})
            # Record XML source excerpt as evidence only, never infer per-CAS TWA.
            reply=ET.tostring(el,encoding='utf-8')
            safe_name='oel_admrul' if target=='admrul' else 'occupational_law'
            (out/(safe_name+'.xml')).write_bytes(reply)
            rows.append({'source':safe_name,'status':'RAW_XML_DOWNLOADED_REVIEW_REQUIRED',
                         'xml_sha256':hashlib.sha256(reply).hexdigest(),
                         'note':'Legal metadata not a per-CAS TWA table; manually review current annex'})
    else:rows.append({'source':'LAW_OC','status':'OFFLINE' if a.offline else 'KEY_NOT_SET','xml_sha256':'','note':'No law evidence approval'})
    # Official data.go.kr catalogue 15157612 specifies KOSHA msdschem1/getChemList001.
    # Its query parameter must be supplied from that service's current API specification.
    if koskey and not a.offline and a.kosha_query_param:
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]{1,40}',a.kosha_query_param):
            raise ValueError('KOSHA 공식 API 조회 매개변수명 형식 오류')
        for c in read(SOURCE/MASTERFILE)[:min(a.kosha_limit,54)]:
            cas=c['CAS']
            xml=http_xml('https://apis.data.go.kr/B552468/msdschem1/getChemList001',
                {'serviceKey':koskey,a.kosha_query_param:cas,'pageNo':1,'numOfRows':5})
            raw=ET.tostring(xml,encoding='utf-8')
            (out/f'kosha_{cas}.xml').write_bytes(raw)
            rows.append({'source':'KOSHA_MSDS_'+cas,'status':'RAW_REFERENCE_REVIEW_REQUIRED',
                         'xml_sha256':hashlib.sha256(raw).hexdigest(),
                         'note':'Reference-only chemical MSDS; manufacturer original and current statutory TWA not certified'})
    else:
        rows.append({'source':'KOSHA_MSDS_SERVICE_KEY',
                     'status':'OFFLINE' if a.offline else (
                        'QUERY_PARAMETER_REVIEW_REQUIRED' if koskey else 'KEY_NOT_SET'),
                     'xml_sha256':'',
                     'note':'Official listing endpoint identified, specific search parameter needs provider spec; never infer TWA'})
    write(out/'기관API_증거수집현황.csv',rows,['source','status','xml_sha256','note'])
    print('03:',', '.join(f"{r['source']}={r['status']}" for r in rows))


def finite_pos(val, message):
    try: x=float(val)
    except (ValueError, TypeError): raise ValueError(message) from None
    if not math.isfinite(x) or x<=0:raise ValueError(message)
    return x


def approved_properties():
    path=WORK/'02_review'/'CAS_물성_검증.csv'
    if not path.exists():raise RuntimeError('02 CAS_물성_검증.csv 파일이 없습니다')
    accepted={}
    for r in read(path):
        if r['HANNA입력승인']!='승인':continue
        cas=r['CAS']
        if r['복합물또는이성질체']!='아니오':raise ValueError(f'{cas} 혼합물을 단일 SMILES로 입력 금지')
        if r['구조검증상태']!='승인' or r['물성검증상태']!='승인':
            raise ValueError(f'{cas} 근거 없는 물성 승인')
        if not (r['SMILES'].strip() and r['구조공식출처'].strip() and r['물성출처와측정조건'].strip()):
            raise ValueError(f'{cas} SMILES/공식 구조·물성출처 누락')
        if abs(float(r['온도_K'])-298.15)>0.001:raise ValueError(f'{cas} 포화증기압 온도 불일치')
        r['_mw']=finite_pos(r['분자량_g_mol'],f'{cas} 분자량 오류')
        r['_vp']=finite_pos(r['25C_포화증기압_kPa'],f'{cas} 포화증기압 오류')
        accepted[cas]=r
    return accepted


def approved_twa():
    path=WORK/'02_review'/'CAS_한국TWA_검증.csv'
    if not path.exists():raise RuntimeError('02 CAS_한국TWA_검증.csv 파일이 없습니다')
    accepted={}
    for r in read(path):
        if r['검증상태']!='승인':continue
        cas=r['CAS']
        if r['공식원문대조']!='승인' or not all(r[x].strip() for x in ('기준시점','고시명번호','공식출처')):
            raise ValueError(f'{cas} 국내 TWA 원문증거 미확인')
        accepted[cas]=finite_pos(r['한국_TWA_ppm'],f'{cas} 국내 TWA 값 오류')
    return accepted


def reviews_index():
    return csv_index(WORK/'02_review'/'제품_원문TDS_검증.csv','문서키')


def stage04(a):
    docs, parts, master=ensure_source()
    props=approved_properties();twas=approved_twa();rev=reviews_index()
    bydoc=defaultdict(list)
    for p in parts:bydoc[p['문서키']].append(p)
    status=[]
    for doc in docs:
        key=doc['문서키']; cs=bydoc[key]; r=rev[key];why=[]
        if len(cs)<2:why.append('성분수 부족')
        if doc['함량합계100가능_표기성분만']!='예':why.append('MSDS 함량 100% 범위 부적합')
        if r['MSDS제3항함량검증']!='승인':why.append('MSDS 제3항 원문 미승인')
        if r['전체모델성분확정']!='승인':why.append('전체 모델 성분 구성 미확정')
        if any(c['CAS분류']=='복합물/이성질체혼합물' for c in cs):why.append('복합물/이성질체 분자동정 불가')
        missing=sorted({c['CAS'] for c in cs}-set(props))
        if missing:why.append('검증 물성 누락: '+','.join(missing))
        twa_missing=sorted({c['CAS'] for c in cs}-set(twas))
        status.append({'문서키':key,'기존제품ID':doc['기존제품ID'],
                       'HANNA계산허용':'예' if not why else '아니오',
                       'TWA전성분승인':'예' if not twa_missing else '아니오',
                       '미승인TWA':'|'.join(twa_missing),'계산보류사유':'; '.join(why)})
    out=WORK/'04_eligibility';write(out/'270건_모델입력_판정.csv',status,list(status[0]))
    jsonfile(out/'approved_provenance.json',{'properties_approved':len(props),'twa_approved':len(twas),
       'model_eligible':sum(x['HANNA계산허용']=='예' for x in status),
       'input_sha256':{'docs':sha(SOURCE/DOCFILE),'components':sha(SOURCE/COMPFILE),
                       'property_review':sha(WORK/'02_review'/'CAS_물성_검증.csv'),
                       'twa_review':sha(WORK/'02_review'/'CAS_한국TWA_검증.csv'),
                       'product_review':sha(WORK/'02_review'/'제품_원문TDS_검증.csv')}})
    print('04: model eligible',sum(x['HANNA계산허용']=='예' for x in status),'/270',
          'approved properties',len(props),'approved TWA',len(twas))


def stage05(a):
    # Only approved, entire identified mixture properties. No fallback to Raoult.
    docs, parts,_=ensure_source();props=approved_properties()
    approved_meta=json.loads((WORK/'04_eligibility'/'approved_provenance.json').read_text(encoding='utf-8'))
    expected={'docs':sha(SOURCE/DOCFILE),'components':sha(SOURCE/COMPFILE),
              'property_review':sha(WORK/'02_review'/'CAS_물성_검증.csv'),
              'twa_review':sha(WORK/'02_review'/'CAS_한국TWA_검증.csv'),
              'product_review':sha(WORK/'02_review'/'제품_원문TDS_검증.csv')}
    if approved_meta['input_sha256']!=expected:
        raise ValueError('04 모델 적용성 평가 이후 검증 자료가 변경됐습니다. 04부터 재실행하세요')
    allowed={r['문서키'] for r in read(WORK/'04_eligibility'/'270건_모델입력_판정.csv') if r['HANNA계산허용']=='예'}
    if not allowed:
        out=WORK/'05_predictions'
        write(out/'HANNA_Raoult_시나리오.csv',[],PRED_SCEN)
        write(out/'HANNA_Raoult_성분.csv',[],PRED_PART)
        write(out/'계산제외.csv',[{'문서키':d['문서키'],'사유':'검증된 전성분 물성/원문 조성 대기'} for d in docs],['문서키','사유'])
        jsonfile(out/'prediction_provenance.json',{'scenarios':0,'documents':0,'status':'WAITING_APPROVAL',
           'source_sha256':{n:sha(SOURCE/n) for n in (DOCFILE,COMPFILE)},
           'property_review_sha256':sha(WORK/'02_review'/'CAS_물성_검증.csv'),
           'product_review_sha256':sha(WORK/'02_review'/'제품_원문TDS_검증.csv')})
        print('05: WAITING_APPROVAL (270 document screening, 0 HANNA predictions)');return
    import numpy as np
    from hanna_gamma_adapter import _get_model
    from run_raoult_271 import feasible_compositions
    from composition_sampling import sample_bounded_simplex
    bydoc=defaultdict(list)
    for p in parts:bydoc[p['문서키']].append(p)
    model=None; scen=[];det=[];excluded=[];rng=random.Random(a.seed)
    for d in docs:
        key=d['문서키']
        if key not in allowed:continue
        cs=bydoc[key]
        smiles=[props[c['CAS']]['SMILES'] for c in cs]
        mw=np.asarray([props[c['CAS']]['_mw'] for c in cs],dtype=float)
        vp=np.asarray([props[c['CAS']]['_vp'] for c in cs],dtype=float)
        intervals=[{'함량 하한(%)':c['함량하한_pct'],'함량 상한(%)':c['함량상한_pct']} for c in cs]
        try:
            weights=(feasible_compositions(intervals,a.samples,rng) if a.sampler=='legacy' else
                     sample_bounded_simplex(intervals,a.samples,rng,burn_in=a.burn_in,thin=a.thin))
        except (RuntimeError,ValueError) as e:
            excluded.append({'문서키':key,'사유':'샘플링 오류: '+str(e)});continue
        weights=np.asarray(weights,dtype=float)
        if (weights<0).any() or not np.allclose(weights.sum(axis=1),100,atol=1e-5):
            raise ValueError(f'{key}: 가상 조성 합계 오류')
        moles=weights/mw[None,:]; x=moles/moles.sum(axis=1,keepdims=True)
        if model is None: model=_get_model()
        all_lng=[]
        for off in range(0,len(x),a.batch_size):
            vals,_=model.predict(smiles,x[off:off+a.batch_size],298.15)
            vals=np.asarray(vals,dtype=float)
            if vals.shape!=(len(x[off:off+a.batch_size]),len(cs)) or not np.isfinite(vals).all():
                raise ValueError(f'{key}: HANNA ln(gamma) 형식/비유한값 오류; fallback 금지')
            all_lng.append(vals)
        gamma=np.exp(np.concatenate(all_lng,axis=0))
        if not np.isfinite(gamma).all() or (gamma<=0).any():raise ValueError(f'{key}: HANNA gamma 오류')
        pr=x*vp[None,:];ph=pr*gamma;tr=pr.sum(axis=1);th=ph.sum(axis=1)
        if (th<=0).any() or (tr<=0).any():raise ValueError(f'{key} 총압력 오류')
        for i in range(len(x)):
            scen.append({'문서키':key,'기존제품ID':d['기존제품ID'],
                         '제품명':d['제품명_파일명기반'],'시나리오ID':i+1,'온도_K':298.15,
                         'Raoult_총압_kPa':float(tr[i]),'HANNA_총압_kPa':float(th[i]),
                         '압력비_HANNA_Raoult':float(th[i]/tr[i]),'모델명':'HANNA 실제 로컬 앙상블',
                         'seed':a.seed,'sampler':a.sampler})
            for j,c in enumerate(cs):
                det.append({'문서키':key,'시나리오ID':i+1,'CAS':c['CAS'],
                            '질량_pct':float(weights[i,j]),'액상몰분율':float(x[i,j]),
                            'HANNA_gamma':float(gamma[i,j]),'Raoult_분압_kPa':float(pr[i,j]),
                            'HANNA_분압_kPa':float(ph[i,j]),'Raoult_기상몰분율':float(pr[i,j]/tr[i]),
                            'HANNA_기상몰분율':float(ph[i,j]/th[i])})
        print(f'05: {d["기존제품ID"] or key} -> {len(x)} real HANNA scenarios')
    out=WORK/'05_predictions';write(out/'HANNA_Raoult_시나리오.csv',scen,PRED_SCEN)
    write(out/'HANNA_Raoult_성분.csv',det,PRED_PART)
    write(out/'계산제외.csv',excluded,['문서키','사유'])
    jsonfile(out/'prediction_provenance.json',{'model':'actual external HANNA 2026 ensemble',
         'temperature_K':298.15,'seed':a.seed,'sampler':a.sampler,'samples_per_document':a.samples,
         'source_sha256':{n:sha(SOURCE/n) for n in (DOCFILE,COMPFILE)},
         'property_review_sha256':sha(WORK/'02_review'/'CAS_물성_검증.csv'),
         'product_review_sha256':sha(WORK/'02_review'/'제품_원문TDS_검증.csv'),
         'scenarios':len(scen),'component_rows':len(det),'documents':len({s['문서키'] for s in scen})})
    print('05 completed:',len(scen),'scenarios; other documents remain unranked')


def stage06(a):
    p=WORK/'05_predictions';out=WORK/'06_hazard'
    scen=read(p/'HANNA_Raoult_시나리오.csv') if (p/'HANNA_Raoult_시나리오.csv').exists() else []
    details=read(p/'HANNA_Raoult_성분.csv') if (p/'HANNA_Raoult_성분.csv').exists() else []
    if scen and not (p/'prediction_provenance.json').exists():raise ValueError('HANNA 실행정보가 없는 예측결과 거부')
    if scen:
        manifest=json.loads((p/'prediction_provenance.json').read_text(encoding='utf-8'))
        if manifest.get('source_sha256')!={n:sha(SOURCE/n) for n in (DOCFILE,COMPFILE)}:
            raise ValueError('MSDS 원본 DB SHA가 계산 시점과 불일치')
        if manifest.get('property_review_sha256')!=sha(WORK/'02_review'/'CAS_물성_검증.csv'):
            raise ValueError('물성승인 기록이 계산 후 변경됨: 05 재실행 필요')
        if manifest.get('product_review_sha256')!=sha(WORK/'02_review'/'제품_원문TDS_검증.csv'):
            raise ValueError('제품 원문승인 기록이 계산 후 변경됨: 05 재실행 필요')
    twa=approved_twa()
    bykey=defaultdict(list)
    for r in details:bykey[(r['문서키'],r['시나리오ID'])].append(r)
    scored=[];pending=[]
    for s in scen:
        key=(s['문서키'],s['시나리오ID']);cs=bykey[key]
        if not cs: raise ValueError('성분이 없는 시나리오 '+str(key))
        p_h=sum(float(x['HANNA_분압_kPa']) for x in cs)
        p_r=sum(float(x['Raoult_분압_kPa']) for x in cs)
        if not math.isclose(p_h,float(s['HANNA_총압_kPa']),rel_tol=1e-6,abs_tol=1e-6) or not math.isclose(p_r,float(s['Raoult_총압_kPa']),rel_tol=1e-6,abs_tol=1e-6):
            raise ValueError('성분 분압 합계/총압력 불일치 '+str(key))
        missing=sorted(set(x['CAS'] for x in cs)-set(twa))
        if missing:
            pending.append({'문서키':key[0],'시나리오ID':key[1],'보류사유':'TWA 미승인 '+','.join(missing)});continue
        hr=[(1e6*float(x['HANNA_분압_kPa'])/(ATM*twa[x['CAS']]),x['CAS']) for x in cs]
        rr=[(1e6*float(x['Raoult_분압_kPa'])/(ATM*twa[x['CAS']]),x['CAS']) for x in cs]
        hmax,hcas=max(hr);rmax,rcas=max(rr)
        scored.append({'문서키':key[0],'시나리오ID':key[1],
                       'HANNA_최대성분정규화평형분압':hmax,'HANNA_최대기여CAS':hcas,
                       'Raoult_최대성분정규화평형분압':rmax,'Raoult_최대기여CAS':rcas,
                       '주의':'노출농도·TWA 초과판정 아님, 평형 기반 스크리닝'})
    write(out/'성분별_잠재지표_시나리오.csv',scored,
          ['문서키','시나리오ID','HANNA_최대성분정규화평형분압','HANNA_최대기여CAS',
           'Raoult_최대성분정규화평형분압','Raoult_최대기여CAS','주의'])
    write(out/'보류.csv',pending,['문서키','시나리오ID','보류사유'])
    jsonfile(out/'hazard_provenance.json',{'verified_twa_review_sha256':sha(WORK/'02_review'/'CAS_한국TWA_검증.csv'),
        'predictions_sha256':sha(p/'HANNA_Raoult_시나리오.csv') if scen else '',
        'scenarios_scored':len(scored),'scenarios_pending':len(pending),
        'interpretation':'Equilibrium partial-pressure based potential; NOT breathing-zone exposure'})
    print('06: scores',len(scored),'TWA pending',len(pending),'NO worker exposure claim')


def stage07(a):
    out=WORK/'07_same_use'
    src=WORK/'06_hazard'/'성분별_잠재지표_시나리오.csv'
    scores=read(src) if src.exists() else []
    if scores:
        mp=json.loads((WORK/'05_predictions'/'prediction_provenance.json').read_text(encoding='utf-8'))
        if mp.get('product_review_sha256')!=sha(WORK/'02_review'/'제품_원문TDS_검증.csv'):
            raise ValueError('제품 검증내역 변경 후 HANNA 계산 재수행 필요')
        hp=json.loads((WORK/'06_hazard'/'hazard_provenance.json').read_text(encoding='utf-8'))
        if hp['verified_twa_review_sha256']!=sha(WORK/'02_review'/'CAS_한국TWA_검증.csv'):
            raise ValueError('TWA 승인정보 변경됨: 06 재실행 필요')
    docs,_,_=ensure_source(); dmap={d['문서키']:d for d in docs}
    rev=reviews_index(); bydoc=defaultdict(list)
    for x in scores:bydoc[x['문서키']].append(float(x['HANNA_최대성분정규화평형분압']))
    groups=defaultdict(list);blocked=[]
    checks=('MSDS제1항용도검증','MSDS제3항함량검증','전체모델성분확정',
            '중복개정판검토','독립제품검증','공정수지호환성검증','TDS성능제약검증')
    for key,d in dmap.items():
        r=rev[key];why=[]
        if len(bydoc[key])<20:why.append('충분한 가상조성+전 성분 TWA 승인+HANNA 결과 없음 (20개 이상 필요)')
        for k in checks:
            if r.get(k)!='승인':why.append(k+' 미승인')
        if not r['비교군ID'].strip():why.append('같은 수지·공정 비교군 미확정')
        if not r['근거문서'].strip():why.append('원본/TDS 근거 미기록')
        if d['함량합계100가능_표기성분만']!='예':why.append('원문 조성오류')
        if why:
            blocked.append({'문서키':key,'제품명':d['제품명_파일명기반'],'이유':' ; '.join(why)})
        else:
            v=bydoc[key]
            groups[(r['비교군ID'],r['용도ID'])].append({'문서키':key,'기존제품ID':d['기존제품ID'],
                  '제품명':d['제품명_파일명기반'], '비교군ID':r['비교군ID'], '용도ID':r['용도ID'],
                  '시나리오수':len(v),'중앙값':statistics.median(v),'최솟값':min(v),'최댓값':max(v)})
    rank=[];pairs=[]
    for group,items in sorted(groups.items()):
        if len(items)<2:
            for item in items:blocked.append({'문서키':item['문서키'],'제품명':item['제품명'],'이유':'검증된 같은 용도/공정 내 비교제품 2건 미만'})
            continue
        items.sort(key=lambda x:(x['중앙값'],x['문서키']))
        for i,item in enumerate(items,1):
            rank.append(dict(item, 참고순위=i, 해석='평형분압 기반 잠재지표, 저위험·동등기능 인증 아님'))
        for i,lo in enumerate(items):
            for hi in items[i+1:]:
                pairs.append({'비교군ID':group[0],'용도ID':group[1],
                  '잠재지표낮은_문서키':lo['문서키'],'잠재지표높은_문서키':hi['문서키'],
                  '중앙값비_높음나눔낮음':hi['중앙값']/lo['중앙값'],
                  '가상조성범위분리':'예' if lo['최댓값']<hi['최솟값'] else '아니오',
                  '주의':'함량구간 가상 시나리오만 비교; 실제 대체제품 승인 아님'})
    write(out/'같은용도_낮은잠재지표_후보.csv',rank,
          ['문서키','기존제품ID','제품명','비교군ID','용도ID','시나리오수','중앙값','최솟값','최댓값','참고순위','해석'])
    write(out/'같은용도_제품쌍_비교.csv',pairs,
          ['비교군ID','용도ID','잠재지표낮은_문서키','잠재지표높은_문서키','중앙값비_높음나눔낮음','가상조성범위분리','주의'])
    write(out/'비교보류_270건.csv',blocked,['문서키','제품명','이유'])
    print('07: approved same-use candidates',len(rank),'pairs',len(pairs),'waiting',len(blocked))


def stage08(a):
    """Prepare an UNAPPROVED optimization template even in a partial local checkout.

    Do not silently approve properties/constraints or claim an optimization result.
    The actual optimizer module is required for any later manual optimization.
    """
    out=WORK/'08_optimization';out.mkdir(parents=True,exist_ok=True)
    requirements=['성분구조검증','동일공정호환성검증','물성출처검증',
                  '한국TWA최신성검증','기능제약실측검증']
    optimizer_available=True
    try:
        from optimize_verified_blend import template
        blank=template()
    except ModuleNotFoundError as exc:
        # Only tolerate a missing optimizer file. Do NOT suppress a missing
        # dependency *inside* an available optimizer or any other import error.
        if exc.name!='optimize_verified_blend':
            raise
        optimizer_available=False
        blank={
          '대상제품군':'예: 특정 도료 및 공정에서 호환성 검증된 그룹',
          '근거문서':'TDS/실측/법령 출처',
          '검증':{key:False for key in requirements},
          '온도_K':298.15,
          '목표_총평형증기압범위_kPa':[0.0,101.325],
          '성분':[{'CAS':'REPLACE','SMILES':'REPLACE','분자량_g_mol':None,
                  '포화증기압_kPa':None,'한국_TWA_ppm':None,
                  '최소질량비_pct':None,'최대질량비_pct':None}],
        }
        print('08: WARNING optimize_verified_blend.py is missing locally; '
              'writing UNAPPROVED template only. Manual optimization is DISABLED.')
    # Preserve manually edited, unapproved templates instead of overwriting them.
    cfg=out/'검증전_Pareto_제약서식.json'
    if not cfg.exists():jsonfile(cfg,blank)
    pairs=read(WORK/'07_same_use'/'같은용도_제품쌍_비교.csv')
    jsonfile(out/'최적화_진행조건.json',{
      'same_use_pairs_evidence':len(pairs),
      'auto_optimization_performed':False,
      'optimizer_module_available':optimizer_available,
      'optimizer_missing_status':'DISABLED_MISSING_LOCAL_MODULE' if not optimizer_available else '',
      'conditions_required':requirements,
      'action':('Fill constraints with measured evidence and run optimize_verified_blend.py explicitly'
                if optimizer_available else
                'Restore research/korean_thinners_271/src/optimize_verified_blend.py '
                'from the GitHub main branch before any manual optimization'),
      'not_proved':'Optimal real-world mixture, product substitution, regulatory exposure'})
    print('08: Pareto template exists, unapproved; optimization not executed')


def stage09(a):
    docs,parts,master=ensure_source();out=WORK/'09_report';out.mkdir(parents=True,exist_ok=True)
    coverage=read(WORK/'04_eligibility'/'270건_모델입력_판정.csv')
    preds=read(WORK/'05_predictions'/'HANNA_Raoult_시나리오.csv') if (WORK/'05_predictions'/'HANNA_Raoult_시나리오.csv').exists() else []
    risk=read(WORK/'06_hazard'/'성분별_잠재지표_시나리오.csv')
    ranks=read(WORK/'07_same_use'/'같은용도_낮은잠재지표_후보.csv')
    pairs=read(WORK/'07_same_use'/'같은용도_제품쌍_비교.csv')
    summary={'genuine_thinner_pdfs':len(docs),'ingredient_records':len(parts),'unique_cas':len(master),
     'model_eligible_approved_docs':sum(x['HANNA계산허용']=='예' for x in coverage),
     'hanna_computed_docs':len({r['문서키'] for r in preds}), 'hanna_scenarios':len(preds),
     'verified_twa_headspace_scenarios':len(risk),'same_use_ranked_documents':len(ranks),
     'same_use_pairs':len(pairs),'real_product_risk_validated':False,
     'external_VLE_accuracy_validated':False,
     'interpretation':'Automated source-linked screening. Regulatory, measured and functional claims require independent evidence.'}
    jsonfile(out/'논문_연구진행_요약.json',summary)
    text='\n'.join(['# 국내 신나 원본 270건 통합 논문 파이프라인 점검','',
        '| 지표 | 값 |','|---|---:|']+[f'| {k} | {v} |' for k,v in summary.items() if isinstance(v,int)])
    text+='\n\n## 논문 해석\n\n- 270개 PDF는 독립 시판 제품 270종을 의미하지 않습니다.\n- 한 번의 MSDS 구간 표본추출은 실제 제조 배합/제품 간 확률분포가 아닙니다.\n- HANNA-Raoult 차이는 독립 실측 VLE 정확도가 아닙니다.\n- TWA 정규화 평형분압은 작업장 실제 노출/법정 기준 위반을 뜻하지 않습니다.\n- 같은 용도 내 순위는 수지·공정·TDS·구성성분·TWA 근거가 모두 승인된 경우에만 참고용으로 산출됩니다.\n- 모델 미적용 복합물·이성질체·미공개 조성을 안전한 것으로 처리하지 않습니다.\n- 최적 배합비는 기능 제약이 독립 검증되기 전에는 주장하지 않습니다.\n'
    (out/'논문용_진행보고서.md').write_text(text,encoding='utf-8')
    print('09: reports generated ->',out)

STAGES={'00':stage00,'01':stage01,'02':stage02,'03':stage03,
        '04':stage04,'05':stage05,'06':stage06,'07':stage07,'08':stage08,'09':stage09}


def cli(stage=None, argv=None):
    ap=argparse.ArgumentParser(description='Source-linked 270 MSDS research pipeline, 00..09')
    if stage is None:ap.add_argument('stage',choices=list(STAGES)+['all'])
    ap.add_argument('--noru-zip',type=Path,help='Optional original Noru ZIP, needs --samhwa-zip')
    ap.add_argument('--samhwa-zip',type=Path,help='Optional original Samhwa ZIP')
    ap.add_argument('--offline',action='store_true',help='Skip official API network retrieval')
    ap.add_argument('--kosha-query-param',help='Verified list-search parameter name in data.go.kr service 15157612 guide')
    ap.add_argument('--kosha-limit',type=int,default=3,help='Max CAS KOSHA evidence requests, 1-54')
    ap.add_argument('--samples',type=int,default=200)
    ap.add_argument('--batch-size',type=int,default=8)
    ap.add_argument('--seed',type=int,default=20261023)
    ap.add_argument('--sampler',choices=['legacy','gibbs'],default='gibbs')
    ap.add_argument('--burn-in',type=int,default=300)
    ap.add_argument('--thin',type=int,default=12)
    a=ap.parse_args(argv)
    if a.samples<1 or a.batch_size<1 or a.thin<1 or a.burn_in<0 or not 1<=a.kosha_limit<=54:ap.error('Invalid samples or KOSHA request limit')
    selected=stage or a.stage
    for s in (list(STAGES) if selected=='all' else [selected]):
        require_previous(s)
        print(f'\n[STAGE {s}]')
        STAGES[s](a)

if __name__=='__main__':cli()