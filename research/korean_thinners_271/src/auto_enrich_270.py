# -*- coding: utf-8 -*-
"""Evidence-aware 270 MSDS enrichment. No fictitious CAS, TWA, VLE or approval.

PubChem public PUG REST supplies candidate structures/MW, PUG VIEW annotations are
reference-only. KOSHA MSDS is NOT the current Korean statutory exposure standard.
Law official document + annex snippets are not accepted as numeric TWA without
human confirmation. Source API keys never enter output or logs.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import time
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError
import xml.etree.ElementTree as ET
import zipfile

BASE=Path(__file__).resolve().parents[1]
ROOT=BASE.parents[1]
WORK=BASE/'output'/'full270'
SOURCEDIR=BASE/'data'/'source_rebuild_20261023'
PROPS='CAS_물성_검증.csv'
TWA='CAS_한국TWA_검증.csv'
PRODUCT='제품_원문TDS_검증.csv'
CAS_PATTERN=re.compile(r'(?<![0-9])\d{2,7}-\d{2}-\d(?![0-9])')
# Known non-single-molecular entries: ATSDR Stoddard solvent; NITE acrylic polymer.
NON_SINGLE_CAS={'8052-41-3','26588-80-7'}


def read(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(p,rs,fields):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rs)
def dump(p,obj):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
def sha(raw):return hashlib.sha256(raw).hexdigest()
def load_env():
    # Delegate to the repository's existing non-overwriting loader.
    from pipeline_270_core import load_env as safe_load_env
    safe_load_env()
def get_bytes(url, params=None, timeout=20, fetch=None):
    if fetch is not None:return fetch(url, params or {})
    if params:url+='?'+urlencode(params)
    req=Request(url,headers={'User-Agent':'Research chemical properties project (academic; PubChem public API)','Accept':'application/json, application/xml;q=0.8'})
    data=None
    for attempt in range(4):
        if 'pubchem.ncbi.nlm.nih.gov' in url:time.sleep(0.26)  # PUG service: <=5 requests/sec
        try:
            with urlopen(req,timeout=timeout) as response:
                data=response.read(3_000_001)
            break
        except HTTPError as exc:
            if exc.code in (429,500,502,503,504) and attempt<3:
                time.sleep(0.6*(2**attempt))
                continue
            raise RuntimeError('기관 요청 실패: HTTP_'+str(exc.code)) from None
        except Exception as exc:
            # Secret-bearing full URL/exception must NEVER be logged.
            raise RuntimeError('기관 요청 실패: '+type(exc).__name__) from None
    if len(data)>3_000_000:raise RuntimeError('응답 크기 제한 초과')
    return data


def cid_lookup(cas, fetch=None):
    base='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+quote(cas,safe='')
    try:
        cids=json.loads(get_bytes(base+'/cids/JSON',fetch=fetch)).get('IdentifierList',{}).get('CID',[])
    except (RuntimeError,ValueError,KeyError):return {'CAS':cas,'상태':'조회실패'}
    if len(cids)!=1:return {'CAS':cas,'상태':'CID_없음' if not cids else '여러_CID_동정필요','CID_목록':'|'.join(map(str,cids[:20]))}
    cid=int(cids[0]); prefix=f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}'
    try:
        ps=json.loads(get_bytes(prefix+'/property/SMILES,MolecularWeight,InChIKey/JSON',fetch=fetch))['PropertyTable']['Properties']
        syn=json.loads(get_bytes(prefix+'/synonyms/JSON',fetch=fetch))['InformationList']['Information'][0].get('Synonym',[])
        if len(ps)!=1:return {'CAS':cas,'CID':cid,'상태':'구조_여러결과_검토'}
    except (RuntimeError,ValueError,KeyError,IndexError,TypeError):return {'CAS':cas,'CID':cid,'상태':'구조_조회실패'}
    if cas not in set(s.strip() for s in syn):return {'CAS':cas,'CID':cid,'상태':'CAS_동의어_미확인'}
    p=ps[0]
    smi=p.get('SMILES') or p.get('CanonicalSMILES') or p.get('ConnectivitySMILES')
    mw=p.get('MolecularWeight')
    if not smi or not mw:return {'CAS':cas,'CID':cid,'상태':'구조또는분자량없음'}
    try:
        if float(mw)<=0:raise ValueError()
    except (ValueError,TypeError):return {'CAS':cas,'CID':cid,'상태':'분자량오류'}
    return {'CAS':cas,'CID':cid,'상태':'CAS_동의어_일치_구조후보',
            'SMILES_후보':str(smi),'분자량_후보_g_mol':str(mw),
            'InChIKey':p.get('InChIKey',''),
            'PubChem_URL':f'https://pubchem.ncbi.nlm.nih.gov/compound/{cid}'}


def flatten_annotated_strings(node):
    """Restrict to actual vapor-pressure heading; retain exact source quotation."""
    hits=[]
    def walk(n,active=False):
        if isinstance(n,dict):
            active = active or ('vapor pressure' in str(n.get('TOCHeading','')).casefold())
            if active and isinstance(n.get('StringWithMarkup'),list):
                for x in n['StringWithMarkup']:
                    if isinstance(x,dict) and x.get('String'):
                        hits.append(str(x['String'])[:450])
            for k,v in n.items():
                if k not in ('StringWithMarkup',):walk(v,active)
        elif isinstance(n,list):
            for x in n:walk(x,active)
    walk(node)
    return list(dict.fromkeys(hits))[:30]


def pubchem_vp_reference(cid,fetch=None):
    url=f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON'
    try:
        data=json.loads(get_bytes(url,{'heading':'Vapor Pressure'},fetch=fetch))
        strings=flatten_annotated_strings(data)
        return {'문구':' | '.join(strings[:6]),'주석문구수':len(strings),'URL':url+'?heading=Vapor%20Pressure',
                '상태':'VP_원문조건검토' if strings else 'VP_주석없음'}
    except (RuntimeError,ValueError,KeyError,TypeError):
        return {'문구':'','주석문구수':0,'URL':url,'상태':'VP_조회실패'}


def direct_25c_vp_candidates(annotations):
    """Conservative candidate-only conversion; not a claim of experimental accuracy.

    Requires an explicit 25 C (or 298.15 K) on the SAME short statement as
    a recognizable pressure quantity. Never extrapolate values at 20 C.
    """
    factors={'mmhg':0.13332236842105263,'torr':0.13332236842105263,
             'kpa':1.0,'pa':0.001,'hpa':0.1,'bar':100.0,'atm':101.325}
    pat=re.compile(r'(?i)(?<![\w.+\-])(\d+(?:\.\d+)?(?:e[+\-]?\d+|[x×]\s*10\s*\^?\s*[+\-]?\d+)?)\s*'
                   r'(mm\s*hg|torr|kpa|hpa|pa|bar|atm)\b')
    tpat=re.compile(r'(?i)\b25\s*(?:°\s*|deg\s*|degrees\s*)?c(?:elsius)?\b|'
                    r'\b298\.15\s*k\b')
    out=[]
    for ann in annotations:
        for statement in re.split(r'[;\n]',ann):
            if not tpat.search(statement):continue
            matches=list(pat.finditer(statement))
            if len(matches)!=1:continue  # avoid multiple units/numbers in one phrase
            m=matches[0];raw_unit=m.group(2).replace(' ','').lower()
            try:
                num=m.group(1).replace('×','x')
                if re.search(r'(?i)x\s*10',num):
                    base,exp=re.split(r'(?i)x\s*10\s*\^?\s*',num,maxsplit=1)
                    magnitude=float(base)*(10**int(exp))
                else:magnitude=float(num)
                val=magnitude*factors[raw_unit]
            except (ValueError,KeyError,OverflowError):continue
            if not (0 < val < 1000):continue
            out.append((val,statement.strip()[:240]))
    # Multiple inconsistent reports should remain unresolved.
    if not out:return None, 'VP_25C_직접수치없음'
    vals=[x[0] for x in out]
    if max(vals)-min(vals)>max(0.02*min(vals),1e-5):
        return None, 'VP_25C_출처별수치불일치'
    return (sum(vals)/len(vals),out[0][1]),'VP_25C_직접조건후보_검증필요'


def property_enrich(*,fetch=None,limit=54,sleep=0.35,work=WORK,apply=True):
    review=Path(work)/'02_review'
    path=review/PROPS;rows=read(path)
    out=Path(work)/'02a_pubchem';out.mkdir(parents=True,exist_ok=True)
    candidates=[];raw=[];repairs=[]
    for entry in rows[:limit]:
        cas=entry['CAS']
        if cas in NON_SINGLE_CAS:
            if entry['HANNA입력승인']=='승인':
                raise ValueError(f'{cas}: known polymer/UVCB cannot be approved as one molecule')
            entry['복합물또는이성질체']='예'
        mixture=entry['복합물또는이성질체']!='아니오'
        source=entry.get('물성출처와측정조건','')
        if (entry['HANNA입력승인']!='승인'
                and entry.get('25C_포화증기압_kPa','').strip()
                and ' ; original: ' in source
                and '자동후보/출처 원문 미승인' in source):
            raw_phrase=source.split(' ; original: ',1)[1].split(' ; 자동후보',1)[0]
            parsed,_status=direct_25c_vp_candidates([raw_phrase])
            try:
                stored=float(entry['25C_포화증기압_kPa'])
                mismatch=parsed is not None and abs(stored-parsed[0])>max(parsed[0]*0.02,1e-5)
            except ValueError:
                mismatch=True
            if mismatch:
                repairs.append({'CAS':cas,'기존값_kPa':entry['25C_포화증기압_kPa'],
                    '원문':raw_phrase,'파서재판독_kPa':format(parsed[0],'.8g') if parsed else '',
                    '상태':'미승인 자동후보 정정·재검토 필요'})
                entry['25C_포화증기압_kPa']=''
                entry['물성출처와측정조건']=''
        if mixture:
            cand={'CAS':cas,'상태':'UVCB_혼합이성질체_단일구조차단','CID':'','SMILES_후보':'','분자량_후보_g_mol':'','InChIKey':'','PubChem_URL':'','VP_후보문구':'','VP_상태':'차단'}
        else:
            c=cid_lookup(cas,fetch=fetch)
            cand={'CAS':cas,**c,'VP_후보문구':'','VP_상태':''}
            if c.get('상태')=='CAS_동의어_일치_구조후보':
                vp=pubchem_vp_reference(c['CID'],fetch=fetch)
                cand.update({'VP_후보문구':vp['문구'],'VP_상태':vp['상태'], 'VP_출처':vp['URL']})
                direct, direct_status=direct_25c_vp_candidates(vp['문구'].split(' | '))
                cand['VP_25C_직접후보_kPa']=format(direct[0],'.8g') if direct else ''
                cand['VP_25C_자동판독상태']=direct_status
                raw.append({'CAS':cas,'PubChem_CID':c['CID'],'vapor_pressure_annotations':vp['문구'],
                            'direct_25C_candidate':direct[0] if direct else None,
                            'review_note':'온도 및 출처·실측/추정 구분 검증 전 승인 금지'})
                if apply and entry.get('HANNA입력승인')!='승인':
                    # Populate only empty provisional fields, never modify existing approvals.
                    if not entry['SMILES'].strip():entry['SMILES']=str(c['SMILES_후보'])
                    if not entry['분자량_g_mol'].strip():entry['분자량_g_mol']=str(c['분자량_후보_g_mol'])
                    if not entry['구조공식출처'].strip():entry['구조공식출처']=str(c['PubChem_URL'])+' (CAS synonym matched; unverified)'
                    # A DIRECT 25C annotation may be provisionally prefilled; it
                    # stays UNAPPROVED and will not unlock HANNA or a hazard ranking.
                    if direct and not entry['25C_포화증기압_kPa'].strip():
                        entry['25C_포화증기압_kPa']=format(direct[0],'.8g')
                        if not entry['물성출처와측정조건'].strip():
                            entry['물성출처와측정조건']=(vp['URL']+' ; original: '+direct[1]
                                                       +' ; 자동후보/출처 원문 미승인')

            if sleep and fetch is None:time.sleep(sleep)
        candidates.append(cand)
    fields=['CAS','상태','CID','CID_목록','SMILES_후보','분자량_후보_g_mol','InChIKey','PubChem_URL','VP_상태','VP_후보문구','VP_출처','VP_25C_직접후보_kPa','VP_25C_자동판독상태']
    write(out/'PubChem_54CAS_구조분자량_자동후보.csv',candidates,fields)
    write(out/'자동VP후보_오류정정감사.csv',repairs,
          ['CAS','기존값_kPa','원문','파서재판독_kPa','상태'])
    dump(out/'PubChem_VP_원문후보_검수.json',raw)
    if apply:write(path,rows,list(rows[0]))
    summary={'processed':len(candidates),'structure_candidates':sum(c.get('상태')=='CAS_동의어_일치_구조후보' for c in candidates),
             'mixed_blocked':sum(c.get('상태')=='UVCB_혼합이성질체_단일구조차단' for c in candidates),
             'VP_25C_direct_candidates':sum(bool(c.get('VP_25C_직접후보_kPa')) for c in candidates),
             'VP_25C_values_automatically_approved':0,'property_approvals_created':0,
             'VP_candidates_corrected':len(repairs)}
    dump(out/'summary.json',summary)
    print('02a:',summary)
    return summary


def safe_xml_cas_snippets(xml_data, cas, window=190):
    try:
        root=ET.fromstring(xml_data)
        s=' '.join(t.strip() for t in root.itertext() if t and t.strip())
    except ET.ParseError:
        s=xml_data.decode('utf-8','replace')
    results=[]
    for m in re.finditer(re.escape(cas),s):
        results.append(s[max(0,m.start()-window):min(len(s),m.end()+window)])
    return list(dict.fromkeys(results))[:8]


def official_twa_evidence(*,fetch=None,limit=54,sleep=0.25,work=WORK,offline=False):
    load_env();out=Path(work)/'02b_twa_evidence';out.mkdir(parents=True,exist_ok=True)
    rows=read(Path(work)/'02_review'/TWA)
    key=os.environ.get('KOSHA_MSDS_SERVICE_KEY','').strip()
    law=os.environ.get('LAW_OC','').strip()
    evidence=[];law_evidence=[]
    if not offline and law:
        # Use the exact regulation title. The response may be a SEARCH LIST, not annex text!
        try:
            raw=get_bytes('https://www.law.go.kr/DRF/lawSearch.do',{
                'OC':law,'target':'admrul','query':'화학물질 및 물리적 인자의 노출기준',
                'nw':1,'type':'XML','display':30,'page':1},fetch=fetch)
            (out/'법제처_행정규칙_검색응답.xml').write_bytes(raw.replace(law.encode(),b'[REDACTED]'))
            law_evidence.append({'source':'LAW_SEARCH','status':'RAW_SEARCH_ONLY_ANNEX_NOT_VALIDATED','sha256':sha(raw),'message':'별표별 원문/시행일을 검증하지 않았음'})
            snippets={cas:safe_xml_cas_snippets(raw,cas) for cas in (r['CAS'] for r in rows)}
            for cas,ss in snippets.items():
                if ss:law_evidence.append({'source':'LAW_SEARCH_'+cas,'status':'CAS_TEXT_CANDIDATE_UNVERIFIED','sha256':sha(raw),'message':' | '.join(ss)[:500]})
            # Retrieve exactly matching present-day regulation details and annex lists;
            # no guessing about which regulation number or annex is the current one.
            try:
                root=ET.fromstring(raw)
                candidates=[]
                for child in root.iter():
                    title=(child.findtext('행정규칙명') or '').strip()
                    if title.replace(' ','')=='화학물질및물리적인자의노출기준':
                        serial=(child.findtext('행정규칙일련번호') or '').strip()
                        if serial.isdigit():candidates.append(serial)
                candidates=list(dict.fromkeys(candidates))
                if len(candidates)==1:
                    detail=get_bytes('https://www.law.go.kr/DRF/lawService.do',{
                        'OC':law,'target':'admrul','ID':candidates[0],'type':'XML'},fetch=fetch)
                    clean=detail.replace(law.encode(),b'[REDACTED]')
                    (out/'법제처_현행노출기준_행정규칙원문.xml').write_bytes(clean)
                    law_evidence.append({'source':'LAW_CURRENT_DETAIL','status':'FULL_TEXT_AND_ANNEX_REVIEW_REQUIRED',
                        'sha256':sha(clean),'message':'행정규칙일련번호='+candidates[0]+'; 별표의 시행일/CAS/ppm 검토 필요'})
                    for ent in rows:
                        snippets=safe_xml_cas_snippets(detail,ent['CAS'])
                        if snippets:
                            law_evidence.append({'source':'LAW_DETAIL_'+ent['CAS'],
                              'status':'CAS_TEXT_CANDIDATE_UNVERIFIED','sha256':sha(clean),
                              'message':' | '.join(snippets)[:500]})
                else:
                    law_evidence.append({'source':'LAW_CURRENT_DETAIL','status':'EXACT_REGULATION_ID_NOT_UNIQUE',
                        'sha256':'','message':str(len(candidates))+' current exact matches; manual match required'})
                # Annex listing is a separate official endpoint, NOT a normalized TWA table.
                annex=get_bytes('https://www.law.go.kr/DRF/lawSearch.do',{
                    'OC':law,'target':'admbyl','search':2,
                    'query':'화학물질 및 물리적 인자의 노출기준','type':'XML','display':100,'page':1},fetch=fetch)
                annex_clean=annex.replace(law.encode(),b'[REDACTED]')
                (out/'법제처_별표서식_목록.xml').write_bytes(annex_clean)
                law_evidence.append({'source':'LAW_ANNEX_LIST','status':'ANNEX_LINKS_NEED_CONTENT_AND_EFFECTIVE_DATE',
                  'sha256':sha(annex_clean),'message':'별표 목록 응답 확보, 실제 별표 표 데이터 자동 검증 전'})
            except (RuntimeError,ET.ParseError,ValueError) as exc:
                law_evidence.append({'source':'LAW_DETAIL_OR_ANNEX','status':'REQUEST_OR_PARSE_FAILED',
                  'sha256':'','message':type(exc).__name__+' (key hidden)'})
        except RuntimeError as exc:law_evidence.append({'source':'LAW_SEARCH','status':'REQUEST_FAILED','sha256':'','message':str(exc)})
    else:law_evidence.append({'source':'LAW_SEARCH','status':'OFFLINE_OR_KEY_MISSING','sha256':'','message':'공식 고시 원문 확인 필요'})
    # KOSHA current catalogue CAS search: searchCnd=1 and searchWrd=CAS.
    for entry in rows[:limit]:
        cas=entry['CAS']
        state='OFFLINE_OR_KEY_MISSING';raw=b'';snippets=[]
        if not offline and key:
            try:
                # Two provider-published paths exist across catalogue revisions.
                # Try the current path first; fallback only on request failure.
                params={'serviceKey':key,'searchWrd':cas,'searchCnd':'1',
                        'numOfRows':10,'pageNo':1}
                endpoints=('https://apis.data.go.kr/B552468/msdschem1/getChemList001',
                           'https://apis.data.go.kr/B552468/msdschem/getChemList')
                for endpoint in endpoints:
                    try:
                        raw=get_bytes(endpoint,params,fetch=fetch)
                        try:
                            root=ET.fromstring(raw)
                            code=(root.findtext('.//resultCode') or root.findtext('.//returnReasonCode') or '').strip()
                            message=(root.findtext('.//resultMsg') or root.findtext('.//returnAuthMsg') or '').strip().upper()
                            bad=(code not in ('','00','0','000','NORMAL_SERVICE') or
                                 any(z in message for z in
                                    ('ERROR','INVALID','UNREGISTERED','SERVICE_KEY','AUTHORIZATION')))
                            if bad:
                                if endpoint==endpoints[0]:continue
                                raise RuntimeError('PROVIDER_XML_ERROR_'+(code or 'MESSAGE'))
                        except ET.ParseError:
                            if endpoint==endpoints[0]:continue
                            raise RuntimeError('PROVIDER_XML_INVALID') from None
                        break
                    except RuntimeError:
                        if endpoint == endpoints[-1]:raise
                snippets=safe_xml_cas_snippets(raw,cas)
                state='KOSHA_REFERENCE_FOUND' if snippets else 'NO_EXACT_CAS_IN_RESPONSE'
                (out/f'KOSHA_MSDS_{cas}.xml').write_bytes(raw.replace(key.encode(),b'[REDACTED]'))
            except RuntimeError as exc:state='REQUEST_FAILED_'+str(exc).replace(' ','_')[:45]
        evidence.append({'CAS':cas,'KOSHA_상태':state,'KOSHA_조회원문_해시':sha(raw.replace(key.encode(),b'[REDACTED]')) if raw and key else (sha(raw) if raw else ''),
                         'KOSHA_CAS_일치원문발췌':' | '.join(snippets[:3])[:700],
                         '공식_국내_TWA_ppm_후보':'','법정고시_별표_확인상태':'미확인',
                         '국내TWA_승인상태':'미검증',
                         'KOSHA_API_근거':'https://www.data.go.kr/data/15157612/openapi.do',
                         '주의':'KOSHA 참고 MSDS는 대한민국 법정 TWA 표와 동일하지 않음'})
        if sleep and fetch is None and key and not offline:time.sleep(sleep)
    write(out/'54CAS_KOSHA_공식TWA_검토큐.csv',evidence,list(evidence[0]))
    write(out/'법제처_공식고시_원문검토상태.csv',law_evidence,['source','status','sha256','message'])
    from collections import Counter
    summary={'KOSHA_status_counts':dict(Counter(x['KOSHA_상태'] for x in evidence)),
             'LAW_status_counts':dict(Counter(x['status'] for x in law_evidence)),
             'CAS_considered':len(evidence),'KOSHA_reference_candidates':sum(x['KOSHA_상태']=='KOSHA_REFERENCE_FOUND' for x in evidence),
             'approved_Korean_TWA_added':0,'law_annex_validated':False,
             'notice':'TWA 수치 자동 승인 불가: 현행 별표·CAS·단위 원문 검증 필요'}
    dump(out/'summary.json',summary)
    print('02b:',summary)
    return summary


def product_evidence(noru_zip,samhwa_zip,*,work=WORK,limit=None):
    """Extract section-1 use evidence and flag TDS/compatibility still missing."""
    try:import fitz
    except ImportError:raise RuntimeError('python -m pip install pymupdf 필요') from None
    review=Path(work)/'02_review';docs=read(SOURCEDIR/'MSDS_원본대조_270.csv')
    mapping={x['원본SHA256']:x for x in docs};current={r['문서키']:r for r in read(review/PRODUCT)}
    out=Path(work)/'02c_product_evidence';out.mkdir(parents=True,exist_ok=True)
    queues=[]
    # The Section-1 text is extracted ONLY from source PDFs and includes no invented TDS content.
    for brand,archive in [('NORU',noru_zip),('SAMHWA',samhwa_zip)]:
        with zipfile.ZipFile(archive) as z:
            for info in z.infolist():
                if not info.filename.lower().endswith('.pdf'):continue
                raw=z.read(info);digest=sha(raw);d=mapping.get(digest)
                if not d:raise ValueError('등록된 원본 SHA와 불일치: '+info.filename)
                if d['제조사']!=brand:raise ValueError('제조사 매칭 실패: '+info.filename)
                pdf=fitz.open(stream=raw,filetype='pdf')
                s='\n'.join(pdf[i].get_text() for i in range(min(2,len(pdf))))
                # Raw limited line evidence; an MSDS title or use is not verified TDS.
                lines=[x.strip() for x in s.splitlines() if x.strip()]
                keywords=['제품의 권고 용도','사용상의 제한','recommended use','specific end use',
                          '사용 용도','권고 용도','usage','application','thinner','신나']
                hints=[x for x in lines if any(k in x.casefold() for k in keywords)]
                section_1_end=next((i for i,x in enumerate(lines) if re.search(
                    r'^\s*2\s*[.．]\s*(?:유해|hazard)',x,re.I)),min(120,len(lines)))
                section_1=lines[:section_1_end]
                anchor=next((i for i,x in enumerate(section_1) if any(k in x.casefold() for k in
                  ['권고 용도','recommended use','product use','사용 용도','권고용도'])),None)
                use_block=' | '.join(section_1[anchor:anchor+15])[:1500] if anchor is not None else ''
                queues.append({'문서키':d['문서키'],'기존제품ID':d['기존제품ID'],'제품명':d['제품명_파일명기반'],
                 '원본SHA256':digest,'제1항_자동발췌':' | '.join(hints[:10])[:1100],
                 '제1항_권고용도_주변문맥':use_block,
                 '기존_용도ID':current[d['문서키']]['용도ID'],
                 '동일용도_비교군_검수필요':'예','TDS_공정호환성_미확인':'예',
                 '검증상태':'AUTO_EXTRACTED_NOT_APPROVED'})
                pdf.close()
                if limit and len(queues)>=limit:break
        if limit and len(queues)>=limit:break
    write(out/'제품270건_제1항_용도자동근거.csv',queues,list(queues[0]))
    report={'processed':len(queues),'manual_TDS_approvals_added':0,'same_use_products_ranked':0,
           'source':'SHA256 matched original MSDS first 1–2 pages',
           'limitations':'TDS·수지 호환성·제품 개정/중복과 실제 비교군은 별도 검토'}
    dump(out/'summary.json',report)
    print('02c:',report)
    return report


def report_queue(*,work=WORK):
    review=Path(work)/'02_review';out=Path(work)/'02d_review_queue';out.mkdir(parents=True,exist_ok=True)
    props=read(review/PROPS);twa=read(review/TWA);products=read(review/PRODUCT)
    cands=read(Path(work)/'02a_pubchem'/'PubChem_54CAS_구조분자량_자동후보.csv') if (Path(work)/'02a_pubchem'/'PubChem_54CAS_구조분자량_자동후보.csv').exists() else []
    ix={r['CAS']:r for r in cands}
    items=[]
    for r in props:
        candidate=ix.get(r['CAS'],{})
        reasons=[]
        if r['복합물또는이성질체']=='예':reasons.append('복합물/이성질체 동정 불가')
        if not r['SMILES'] or not r['분자량_g_mol']:reasons.append('SMILES/MW 없음')
        if not r['25C_포화증기압_kPa']:reasons.append('25C 포화증기압 미확보')
        if r['HANNA입력승인']!='승인':reasons.append('물성 원문 검증 미승인')
        items.append({'범주':'물성','키':r['CAS'],'상태':'보류' if reasons else '입력승인',
                      '검수우선항목':' ; '.join(reasons),'자동근거':candidate.get('PubChem_URL','')})
    for r in twa:
        items.append({'범주':'한국TWA','키':r['CAS'],'상태':'보류' if r['검증상태']!='승인' else '승인',
                      '검수우선항목':'법제처 현행 별표의 CAS/단위(ppm)/고시 시행일 대조' if r['검증상태']!='승인' else '',
                      '자동근거':r['공식출처']})
    for r in products:
        flags=['MSDS제1항용도검증','MSDS제3항함량검증','전체모델성분확정','중복개정판검토','독립제품검증','공정수지호환성검증','TDS성능제약검증']
        missing=[k for k in flags if r.get(k)!='승인']
        if not r.get('비교군ID'):missing.append('비교군ID')
        items.append({'범주':'제품','키':r['문서키'],'상태':'보류' if missing else '승인',
                      '검수우선항목':' ; '.join(missing),'자동근거':r.get('근거문서','')})
    write(out/'논문연구_자동입력후_검수대상.csv',items,list(items[0]))
    stats={'property_rows':len(props),'structure_candidates_in_review':sum(bool(r['SMILES']) for r in props),
           'VP_25C_present':sum(bool(r['25C_포화증기압_kPa']) for r in props),
           'HANNA_approved_CAS':sum(r['HANNA입력승인']=='승인' for r in props),
           'Korean_TWA_approved_CAS':sum(r['검증상태']=='승인' for r in twa),
           'same_use_product_evidence_approved':sum(r['공정수지호환성검증']=='승인' for r in products),
           'review_queue_count':len(items)}
    dump(out/'summary.json',stats);print('02d:',stats);return stats


def main(argv=None):
    parser=argparse.ArgumentParser(description='270 MSDS automated evidence collection; NO automatic approvals')
    parser.add_argument('step',choices=['pubchem','twa','products','queue'])
    parser.add_argument('--limit',type=int,default=None,help='Use 1..54 for CAS pilot')
    parser.add_argument('--offline',action='store_true')
    parser.add_argument('--no-apply',action='store_true',help='Do not populate blank SMILES/MW in review input')
    parser.add_argument('--noru-zip',type=Path)
    parser.add_argument('--samhwa-zip',type=Path)
    a=parser.parse_args(argv)
    if a.limit is not None and a.limit<1:parser.error('--limit must be >=1')
    if a.step=='pubchem':property_enrich(limit=min(a.limit or 54,54),apply=not a.no_apply)
    elif a.step=='twa':official_twa_evidence(limit=min(a.limit or 54,54),offline=a.offline)
    elif a.step=='products':
        if not a.noru_zip or not a.samhwa_zip:parser.error('--noru-zip and --samhwa-zip are required for source evidence')
        product_evidence(a.noru_zip,a.samhwa_zip,limit=a.limit)
    else:report_queue()
if __name__=='__main__':main()