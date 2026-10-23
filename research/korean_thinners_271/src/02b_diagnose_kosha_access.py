# -*- coding: utf-8 -*-
"""One-request KOSHA MSDS OpenAPI diagnostic. NEVER print or save secrets/URLs.

Official data.go.kr service: 15157612, B552468/msdschem1/getChemList001.
This is a reference MSDS API, not an authoritative Korean statutory TWA table.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import re
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

from pipeline_270_core import load_env

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'output'/'full270'/'02b_twa_evidence'
ENDPOINT='https://apis.data.go.kr/B552468/msdschem1/getChemList001'
KNOWN_CODES={
    '01':'GW_APPLICATION_ERROR',
    '04':'HTTP_REQUEST_OR_PROVIDER_ERROR',
    '05':'SERVICE_TIMEOUT',
    '10':'INVALID_REQUEST_PARAMETER',
    '12':'API_SERVICE_NOT_AVAILABLE',
    '20':'ACCESS_DENIED_OR_MISSING_KEY',
    '22':'DAILY_QUOTA_EXCEEDED',
    '23':'REQUEST_RATE_LIMIT',
    '29':'IP_BLOCKED',
    '30':'SERVICE_KEY_NOT_REGISTERED',
    '31':'SERVICE_KEY_EXPIRED',
}
def response_metadata(raw):
    """Extract only whitelisted error codes; do not echo free-form server messages."""
    typ='XML'
    try:
        root=ET.fromstring(raw)
    except ET.ParseError:
        text=raw[:160].lstrip().lower()
        return {'body_format':'HTML' if b'<html' in text or b'<!doctype html' in text else 'OTHER',
                'provider_code':'NOT_XML','provider_status':'UNKNOWN'}
    codes=[]
    for elem in root.iter():
        if elem.tag.split('}')[-1] in ('resultCode','returnReasonCode','errCode','errorCode'):
            v=(elem.text or '').strip().upper()
            if re.fullmatch(r'[A-Z_0-9-]{1,40}',v):
                codes.append(v)
    code=next((v for v in codes if v not in ('00','0','000','NORMAL_SERVICE')),None)
    if code is None:
        code=codes[0] if codes else 'NONE'
    return {'body_format':typ,'provider_code':code,
            'provider_status':KNOWN_CODES.get(code,'SUCCESS_OR_UNCLASSIFIED' if code in ('00','0','000','NORMAL_SERVICE','NONE') else 'UNRECOGNIZED_PROVIDER_ERROR')}

def normalize_key(raw):
    """If a *percent-encoded* portal key was placed in .env, decode once.

    urllib.parse.urlencode will then encode exactly once when building the query.
    """
    s=raw.strip()
    encoded=bool(re.search(r'%[0-9A-Fa-f]{2}',s))
    return (urllib.parse.unquote(s) if encoded else s),encoded

def probe(*,cas='71-43-2',fetch=None):
    load_env()
    key=os.environ.get('KOSHA_MSDS_SERVICE_KEY','').strip()
    result={
        'checked_utc':datetime.now(timezone.utc).isoformat(),
        'service':'data.go.kr/15157612; KOSHA MSDS reference',
        'cas_probe':cas,
        'key_configured':bool(key),
        'key_was_percent_encoded':False,
        'http_status':None,
        'body_format':'',
        'provider_code':'',
        'provider_status':'',
        'diagnosis':'',
        'korean_TWA_validated':False,
    }
    if not key:
        result['diagnosis']='KOSHA_MSDS_SERVICE_KEY_NOT_CONFIGURED'
        return result
    key,was_encoded=normalize_key(key)
    result['key_was_percent_encoded']=was_encoded
    params={'serviceKey':key,'searchWrd':cas,'searchCnd':'1','numOfRows':'1','pageNo':'1'}
    try:
        if fetch is not None:
            status,body=fetch(ENDPOINT,params)
        else:
            url=ENDPOINT+'?'+urllib.parse.urlencode(params)
            req=urllib.request.Request(url,headers={'User-Agent':'ThinnerResearch-KOSHA-Access-Diagnostic/1.0'})
            try:
                with urllib.request.urlopen(req,timeout=18) as response:
                    status=response.getcode()
                    body=response.read(8192)
            except urllib.error.HTTPError as err:
                status=err.code
                body=err.read(8192)
        result['http_status']=status
        result.update(response_metadata(body))
        if status==403:
            result['diagnosis']='HTTP_403_CHECK_APPLICATION_APPROVAL_KEY_AND_GATEWAY'
        elif status==401:
            result['diagnosis']='HTTP_401_AUTHORIZATION_REJECTED'
        elif status!=200:
            result['diagnosis']='HTTP_NON_200_CHECK_PROVIDER_CODE'
        elif result['provider_code'] in ('20','30','31','12','10','22','23','29'):
            result['diagnosis']='HTTP_200_BUT_PROVIDER_ERROR'
        else:
            result['diagnosis']='HTTP_200_CHECK_DATA_ROWS_MANUALLY'
    except Exception as exc:
        # Some Python URL exceptions may embed the secret-bearing URL.
        # Store ONLY the exception's type, never str(exc).
        result['diagnosis']='NETWORK_FAILURE_'+type(exc).__name__
    return result

def main(argv=None):
    p=argparse.ArgumentParser(description='Diagnose exactly one KOSHA request without showing credentials')
    p.add_argument('--cas',default='71-43-2',help='One CAS to probe, default benzene')
    p.add_argument('--out',type=Path,default=DEST/'KOSHA_단일요청_접근진단.json')
    args=p.parse_args(argv)
    if not re.fullmatch(r'\d{2,7}-\d{2}-\d',args.cas):
        p.error('Invalid CAS syntax')
    report=probe(cas=args.cas)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))
    print('Keys and full request URLs are deliberately omitted. A KOSHA MSDS '
          'lookup is not official statutory Korean TWA evidence.')
if __name__=='__main__':main()
