# -*- coding: utf-8 -*-
"""Audit Korean SDS composition ranges; no inferred formulation."""
from __future__ import annotations
import argparse
import csv
import hashlib
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DEFAULT_PRODUCTS = ROOT/'domestic_product_inputs'/'korean_msds_products_7.csv'
DEFAULT_COMPONENTS = ROOT/'domestic_product_inputs'/'korean_msds_components_39.csv'

def load(path):
    with Path(path).open('r',encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))

def audit(products, components, pdf_dir=None):
    if not products or not components:
        raise ValueError('제품/성분 데이터 없음')
    ids=[p['product_id'] for p in products]
    if len(set(ids))!=len(ids):
        raise ValueError('제품 ID 중복')
    product_map={p['product_id']:p for p in products}
    groups=defaultdict(list)
    for c in components:
        if c['product_id'] not in product_map:
            raise ValueError('등록되지 않은 제품의 성분')
        groups[c['product_id']].append(c)
    results=[]
    for p in products:
        pid=p['product_id']; rows=groups[pid]
        if len(rows)!=int(p['components_reported']) or len({c['CAS'] for c in rows})!=len(rows):
            raise ValueError('제품별 성분 수 또는 CAS 중복 오류: '+pid)
        lo=hi=0.0
        for c in rows:
            if any(c[k]!=p[k] for k in ('product_id','source_file','source_pdf_sha256','msds_revision')):
                raise ValueError('원본 문서 출처 불일치: '+pid)
            a,b=float(c['wt_pct_min']),float(c['wt_pct_max'])
            if not all(map(math.isfinite,(a,b))) or not 0<=a<=b<=100:
                raise ValueError('함량 범위 오류: '+pid)
            if c['min_inclusive'] not in ('yes','no') or c['max_inclusive'] not in ('yes','no'):
                raise ValueError('하한/상한 경계 값 누락')
            lo+=a;hi+=b
        if abs(lo-float(p['declared_min_sum_wt_pct']))>1e-8 or abs(hi-float(p['declared_max_sum_wt_pct']))>1e-8:
            raise ValueError('문서별 함량 합계 불일치: '+pid)
        possible=lo<=100+1e-9<=hi+1e-9
        if p['bounds_allow_100']!=('yes' if possible else 'no'):
            raise ValueError('합계 100% 가능 여부 오류: '+pid)
        if p['all_components_exact']!='no':
            raise ValueError('범위형 MSDS를 실제 확정 배합비로 잘못 표시함: '+pid)
        if (any(c['chemical_identity']=='mixed_isomers' for c in rows)) != (p['contains_mixed_isomers']=='yes'):
            raise ValueError('혼합 자일렌 이성질체 표시 오류: '+pid)
        if (any(c['chemical_identity']=='UVCB' for c in rows)) != (p['contains_uvcb']=='yes'):
            raise ValueError('석유계 UVCB 표시 오류: '+pid)
        if pdf_dir is not None:
            source=Path(pdf_dir)/p['source_file']
            if hashlib.sha256(source.read_bytes()).hexdigest()!=p['source_pdf_sha256']:
                raise ValueError('원본 PDF SHA256 불일치: '+source.name)
        results.append((pid,len(rows),lo,hi,possible))
    return results

def main():
    parser=argparse.ArgumentParser(description='국내 MSDS 원문 함량 데이터 무결성 검증')
    parser.add_argument('--products',type=Path,default=DEFAULT_PRODUCTS)
    parser.add_argument('--components',type=Path,default=DEFAULT_COMPONENTS)
    parser.add_argument('--pdf-dir',type=Path,default=None)
    args=parser.parse_args()
    rows=audit(load(args.products),load(args.components),args.pdf_dir)
    for pid,n,lo,hi,possible in rows:
        print(f'{pid}: n={n}, sum-bounds=[{lo:.1f}, {hi:.1f}], contains-100={possible}')
    print(f'PASS: {len(rows)} products, {sum(x[1] for x in rows)} component rows, HANNA not executed')

if __name__=='__main__':
    main()
