# -*- coding: utf-8 -*-
"""Full-population, evidence-gated analysis for three thinner-paper questions.

RQ1: Same-use product contrasts in component saturated-vapor potential and
     mixture headspace index (sum of component partial-pressure/OEL ratios).
RQ2: Raoult vs REAL HANNA product rank reversals and component contributions.
RQ3: MSDS constrained low/middle/high composition anchors and sampled robustness.

All 270 MSDS documents enter the eligibility census. No approval, OEL or gamma
is ever imputed. These are *equilibrium potential* scores, not workplace dose.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import statistics
import sys

BASE = Path(__file__).resolve().parents[1]
DATA = BASE / 'data' / 'source_rebuild_20261023'
WORK = BASE / 'output' / 'full270'
OUT = WORK / '10_13_paper_questions'
ATM = 101.325
REQUIRED_REVIEW = ('MSDS제1항용도검증', 'MSDS제3항함량검증',
                   '전체모델성분확정', '중복개정판검토', '독립제품검증',
                   '공정수지호환성검증', 'TDS성능제약검증')
MIXED_NAME = '복합물/이성질체혼합물'


def read(path):
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f'필수 입력 누락: {p}')
    with p.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def write(path, rows, fields):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def dump(path, obj):
    p=Path(path); p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')


def number(v):
    try:
        x=float(v)
        return x if math.isfinite(x) and x>0 else None
    except (ValueError,TypeError,OverflowError):
        return None


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def index_unique(rows,key):
    d={}
    for r in rows:
        k=r.get(key,'').strip()
        if not k or k in d:
            raise ValueError(f'중복 또는 빈 {key}: {k}')
        d[k]=r
    return d


def load_inputs(work=WORK,data=DATA):
    docs=read(data/'MSDS_원본대조_270.csv')
    parts=read(data/'MSDS_구성성분_1454.csv')
    products=read(work/'02_review'/'제품_원문TDS_검증.csv')
    props=read(work/'02_review'/'CAS_물성_검증.csv')
    twa=read(work/'02_review'/'CAS_한국TWA_검증.csv')
    if (len(docs),len(parts))!=(270,1454):
        raise ValueError('전체 270 MSDS/1454 성분 모집단이 아니므로 중단')
    d=index_unique(docs,'문서키');r=index_unique(products,'문서키')
    p=index_unique(props,'CAS');t=index_unique(twa,'CAS')
    if set(d)!=set(r):
        raise ValueError('제품 검토 CSV 문서키와 270건 모집단이 불일치')
    for x in parts:
        if x['문서키'] not in d:
            raise ValueError('원본 문서에 존재하지 않는 CAS 성분 연결')
    return dict(docs=docs,parts=parts,product=r,property=p,twa=t,
                docmap=d,work=Path(work),data=Path(data))


def approved_twa(t):
    v=number(t.get('한국_TWA_ppm'))
    return (v if t.get('검증상태')=='승인' and
            t.get('공식원문대조')=='승인' and
            all(t.get(k,'').strip() for k in ('기준시점','고시명번호','공식출처')) else None)


def approved_property(p):
    vp=number(p.get('25C_포화증기압_kPa')); mw=number(p.get('분자량_g_mol'))
    return (vp,mw,p['SMILES']) if (
        vp and mw and p.get('복합물또는이성질체')=='아니오' and
        p.get('HANNA입력승인')=='승인' and p.get('구조검증상태')=='승인' and
        p.get('물성검증상태')=='승인' and p.get('SMILES','').strip() and
        p.get('구조공식출처','').strip() and p.get('물성출처와측정조건','').strip() and
        abs(float(p.get('온도_K') or 0)-298.15)<0.001) else None


def provenance_check(ctx):
    """Reject stale or synthetic stage05 outputs without source-linked manifest."""
    stage=ctx['work']/'05_predictions'
    scenarios=read(stage/'HANNA_Raoult_시나리오.csv')
    parts=read(stage/'HANNA_Raoult_성분.csv')
    if scenarios:
        meta_path=stage/'prediction_provenance.json'
        if not meta_path.exists():
            raise ValueError('예측결과 provenance 없음: 05 재실행 필요')
        meta=json.loads(meta_path.read_text(encoding='utf-8'))
        checks={
          'property_review_sha256':ctx['work']/'02_review'/'CAS_물성_검증.csv',
          'product_review_sha256':ctx['work']/'02_review'/'제품_원문TDS_검증.csv',
        }
        for k,path in checks.items():
            if meta.get(k)!=sha(path):
                raise ValueError(f'{k} 변경: 05 재실행 필요')
        sources=meta.get('source_sha256') or {}
        for name in ('MSDS_원본대조_270.csv','MSDS_구성성분_1454.csv'):
            if sources.get(name)!=sha(ctx['data']/name):
                raise ValueError('원본 DB SHA 변경: 05 재실행 필요')
        if meta.get('model')!='actual external HANNA 2026 ensemble':
            raise ValueError('실제 HANNA 모델 provenance 미확인')
    return scenarios,parts


def empirical_scores(ctx):
    """Only scores with full product & CAS evidence; paired Raoult/HANNA draws."""
    scen,parts=provenance_check(ctx)
    per_scen={
        (r['문서키'],str(r['시나리오ID'])):r for r in scen}
    if len(per_scen)!=len(scen):
        raise ValueError('중복된 예측 시나리오')
    perparts=defaultdict(list)
    for row in parts:
        perparts[(row['문서키'],str(row['시나리오ID']))].append(row)
    res=defaultdict(list);delta=defaultdict(lambda:defaultdict(list)); comp=[]
    tp={cas:approved_twa(r) for cas,r in ctx['twa'].items()}
    vp={cas:approved_property(r) for cas,r in ctx['property'].items()}
    for key,scenario in per_scen.items():
        if key[0] not in ctx['docmap']:
            raise ValueError('270모집단 외 예측 시나리오 입력 금지')
        ss=perparts.get(key,[])
        if not ss:
            raise ValueError('성분별 분압 없는 예측 시나리오')
        names=[s['CAS'] for s in ss]
        if len(names)!=len(set(names)):
            raise ValueError('한 시나리오 내 중복 CAS')
        if not all(tp.get(c) and vp.get(c) for c in names):
            continue  # No numeric OEL substitution, even for a subset of components.
        ra_h=[]; ha_h=[]; cas_contrib={}
        for p in ss:
            cas=p['CAS'];r=number(p['Raoult_분압_kPa']);h=number(p['HANNA_분압_kPa'])
            if not r or not h:raise ValueError('시나리오 성분 분압 0/음수/비유한')
            kr=1e6*r/(ATM*tp[cas]);kh=1e6*h/(ATM*tp[cas]);
            ra_h.append(kr);ha_h.append(kh)
            cas_contrib[cas]=(kr,kh)
            comp.append(dict(문서키=key[0],시나리오ID=key[1],CAS=cas,
                 순물질_포화증기위험비=1e6*vp[cas][0]/(ATM*tp[cas]),
                 Raoult_성분평형기여=kr,HANNA_성분평형기여=kh,
                 HANNA_차이_기여=kh-kr))
        # Completeness and pressure consistency: do NOT silently drop components.
        if len(ss)!=len([p for p in ctx['parts'] if p['문서키']==key[0]]):
            raise ValueError('05 예측 성분 수와 원문 성분 수 불일치')
        for label,li,col in [('Raoult',ra_h,'Raoult_총압_kPa'),('HANNA',ha_h,'HANNA_총압_kPa')]:
            pressure=sum(float(p[f'{label}_분압_kPa']) for p in ss)
            total=number(scenario.get(col))
            if total is None or not math.isclose(pressure,total,rel_tol=1e-5,abs_tol=1e-5):
                raise ValueError(f'{label} 성분 분압/총압 불일치')
        res[key[0]].append(dict(문서키=key[0],시나리오ID=key[1],
            Raoult_혼합물지수=sum(ra_h),HANNA_혼합물지수=sum(ha_h),
            Raoult_최대성분위험비=max(ra_h),HANNA_최대성분위험비=max(ha_h),
            HANNA_대_Raoult_혼합물지수비=sum(ha_h)/sum(ra_h)))
        for cas,(r,h) in cas_contrib.items():
            delta[key[0]][cas].append(h-r)
    return res,delta,comp


def census(ctx,score):
    bydoc=defaultdict(list)
    for p in ctx['parts']:bydoc[p['문서키']].append(p)
    rows=[];groups=defaultdict(list)
    for doc in ctx['docs']:
        key=doc['문서키']; review=ctx['product'][key]; parts=bydoc[key]; reasons=[]
        if doc.get('함량합계100가능_표기성분만')!='예':
            reasons.append('원문 함량구간으로 100% 조성 불가능')
        if len(parts)<2:reasons.append('다성분 비교 대상 아님')
        if any(p.get('CAS분류')==MIXED_NAME for p in parts):
            reasons.append('복합물/혼합이성질체/고분자 포함')
        if any(not approved_property(ctx['property'].get(p['CAS'],{})) for p in parts):
            reasons.append('전 성분 구조/25C 증기압 승인 부족')
        if any(not approved_twa(ctx['twa'].get(p['CAS'],{})) for p in parts):
            reasons.append('전 성분 한국 법정 TWA 승인 부족')
        for item in REQUIRED_REVIEW:
            if review.get(item)!='승인':reasons.append(item+' 미승인')
        if not review.get('용도ID','').strip() or not review.get('비교군ID','').strip():
            reasons.append('같은 도료/수지/공정 비교군 미검증')
        if not review.get('근거문서','').strip():
            reasons.append('MSDS/TDS 근거문서 누락')
        series=score.get(key,[])
        if len(series)<20:reasons.append('실제 HANNA 및 TWA 완비 시나리오 20개 미만')
        group=(review.get('용도ID','').strip(),review.get('비교군ID','').strip())
        record=dict(문서키=key,기존제품ID=doc['기존제품ID'],제조사=doc['제조사'],
            제품명=doc['제품명_파일명기반'],용도ID=group[0],비교군ID=group[1],
            원본성분수=len(parts),완전평형시나리오수=len(series),
            제품단위검증='승인' if not reasons else '보류',보류사유='; '.join(reasons))
        rows.append(record)
        if not reasons:groups[group].append(record)
    for group,rr in groups.items():
        count=Counter(r['기존제품ID'] for r in rr)
        if any(v>1 or not k for k,v in count.items()):
            for r in rr:r['제품단위검증']='보류';r['보류사유']='같은 비교군 독립 제품 ID 중복'
            groups[group]=[]
        elif len(rr)<2:
            for r in rr:r['제품단위검증']='보류';r['보류사유']='같은 공정 내 승인된 비교제품 2개 미만'
            groups[group]=[]
    return rows,{k:v for k,v in groups.items() if len(v)>=2}


def summarize_products(groups,scores):
    output=[]
    for group,eligible in sorted(groups.items()):
        for r in eligible:
            v=scores[r['문서키']]
            output.append(dict(문서키=r['문서키'],기존제품ID=r['기존제품ID'],
              제조사=r['제조사'],제품명=r['제품명'],용도ID=group[0],비교군ID=group[1],
              시나리오수=len(v),
              Raoult_혼합물지수_중앙값=statistics.median(s['Raoult_혼합물지수'] for s in v),
              HANNA_혼합물지수_중앙값=statistics.median(s['HANNA_혼합물지수'] for s in v),
              Raoult_혼합물지수_최솟값=min(s['Raoult_혼합물지수'] for s in v),
              Raoult_혼합물지수_최댓값=max(s['Raoult_혼합물지수'] for s in v),
              HANNA_혼합물지수_최솟값=min(s['HANNA_혼합물지수'] for s in v),
              HANNA_혼합물지수_최댓값=max(s['HANNA_혼합물지수'] for s in v),
              Raoult_최대성분위험비_중앙값=statistics.median(s['Raoult_최대성분위험비'] for s in v),
              HANNA_최대성분위험비_중앙값=statistics.median(s['HANNA_최대성분위험비'] for s in v)))
    return output


def model_order(rows,model):
    key=f'{model}_혼합물지수_중앙값'
    return sorted(rows,key=lambda r:(r[key],r['문서키']))


def compare_pairs(products,delta):
    bygroup=defaultdict(list)
    for p in products:bygroup[(p['용도ID'],p['비교군ID'])].append(p)
    rq1=[];rq2=[];gsummary=[]
    for (use,gid),rows in sorted(bygroup.items()):
        ra=model_order(rows,'Raoult'); ha=model_order(rows,'HANNA')
        ranks={model:{p['문서키']:i+1 for i,p in enumerate(order)}
               for model,order in [('Raoult',ra),('HANNA',ha)]}
        def relative_grade(position):
            # Within-group exploration only; not a statutory or toxicological class.
            percentile=(position-1)/(len(rows)-1)
            return ('상대낮음' if percentile<=1/3 else
                    '상대중간' if percentile<2/3 else '상대높음')
        for p in rows:
            p['Raoult_순위']=ranks['Raoult'][p['문서키']]
            p['HANNA_순위']=ranks['HANNA'][p['문서키']]
            p['순위변화']=p['HANNA_순위']-p['Raoult_순위']
            p['Raoult_상대등급']=relative_grade(p['Raoult_순위'])
            p['HANNA_상대등급']=relative_grade(p['HANNA_순위'])
            p['상대등급변화']='YES' if p['Raoult_상대등급']!=p['HANNA_상대등급'] else 'NO'
        gsummary.append(dict(용도ID=use,비교군ID=gid,검증제품수=len(rows),
            Raoult_최고대최저_배수=ra[-1]['Raoult_혼합물지수_중앙값']/ra[0]['Raoult_혼합물지수_중앙값'],
            HANNA_최고대최저_배수=ha[-1]['HANNA_혼합물지수_중앙값']/ha[0]['HANNA_혼합물지수_중앙값'],
            Raoult_지표낮은제품=ra[0]['문서키'],Raoult_지표높은제품=ra[-1]['문서키'],
            HANNA_지표낮은제품=ha[0]['문서키'],HANNA_지표높은제품=ha[-1]['문서키']))
        for x,y in itertools.combinations(sorted(rows,key=lambda r:r['문서키']),2):
            rd=x['Raoult_혼합물지수_중앙값']-y['Raoult_혼합물지수_중앙값']
            hd=x['HANNA_혼합물지수_중앙값']-y['HANNA_혼합물지수_중앙값']
            # Meaningful sign changes only, no rank flipping claimed for ties.
            def sign(diff,a,b):
                tol=1e-9*max(a,b,1.0)
                return 0 if abs(diff)<=tol else (1 if diff>0 else -1)
            rs=sign(rd,x['Raoult_혼합물지수_중앙값'],y['Raoult_혼합물지수_중앙값'])
            hs=sign(hd,x['HANNA_혼합물지수_중앙값'],y['HANNA_혼합물지수_중앙값'])
            flip='YES' if rs*hs==-1 else ('TIE_UNRESOLVED' if not rs or not hs else 'NO')
            cx=delta.get(x['문서키'],{});cy=delta.get(y['문서키'],{})
            comps=[]
            for cas in sorted(set(cx)|set(cy)):
                ax=statistics.mean(cx[cas]) if cx.get(cas) else 0.0
                ay=statistics.mean(cy[cas]) if cy.get(cas) else 0.0
                comps.append((cas,ay-ax))
            driver=max(comps,key=lambda z:abs(z[1])) if comps else ('',0.0)
            rq1.append(dict(용도ID=use,비교군ID=gid,문서키_A=x['문서키'],문서키_B=y['문서키'],
               Raoult_A대B_중앙값비=x['Raoult_혼합물지수_중앙값']/y['Raoult_혼합물지수_중앙값'],
               HANNA_A대B_중앙값비=x['HANNA_혼합물지수_중앙값']/y['HANNA_혼합물지수_중앙값'],
               해석='평형기반 잠재지표의 탐색적 제품쌍 비교; 실제 노출/안전 인증 아님'))
            rq2.append(dict(용도ID=use,비교군ID=gid,문서키_A=x['문서키'],문서키_B=y['문서키'],
               Raoult_A순위=ranks['Raoult'][x['문서키']],Raoult_B순위=ranks['Raoult'][y['문서키']],
               HANNA_A순위=ranks['HANNA'][x['문서키']],HANNA_B순위=ranks['HANNA'][y['문서키']],
               순위역전=flip,변화설명_최대기여CAS=driver[0],
               기여차이_B빼기A=driver[1],
               설명='시나리오 평균 성분별 HANNA-Raoult 차이에 대한 기여, 인과 확정 아님'))
    return rq1,rq2,gsummary


def feasible_projection(lower,upper,target=100.0,fix=None):
    """Project interval midpoints on bounded mass simplex (convex projection)."""
    n=len(lower)
    if len(upper)!=n or n<2:raise ValueError('2개 이상 성분 필요')
    if any(l<0 or u<l or u>100 for l,u in zip(lower,upper)):
        raise ValueError('조성 하한/상한 오류')
    lows=list(lower);ups=list(upper)
    if fix:
        k,value=fix
        if value<lows[k]-1e-8 or value>ups[k]+1e-8:raise ValueError('유효 범위 외 고정값')
        lows[k]=ups[k]=value
    if sum(lows)>target+1e-8 or sum(ups)<target-1e-8:
        raise ValueError('지정된 조성값으로 100% 배합 불가')
    mids=[(l+u)/2 for l,u in zip(lows,ups)]
    left=-200.0;right=200.0
    for _ in range(100):
        lam=(left+right)/2
        vals=[min(u,max(l,m+lam)) for l,u,m in zip(lows,ups,mids)]
        if sum(vals)<target:left=lam
        else:right=lam
    vals=[min(u,max(l,m+(left+right)/2)) for l,u,m in zip(lows,ups,mids)]
    if not math.isclose(sum(vals),target,abs_tol=1e-6):
        raise ValueError('프로젝션 합계 100% 실패')
    return vals


def composition_anchors(rows):
    """One central point plus each constituent's feasible LOW/HIGH endpoint.

    'Low/mid/high' is NOT simultaneous every-CAS low/high (usually impossible).
    """
    if len(rows)<2:raise ValueError('다성분 필요')
    lower=[float(r['함량하한_pct']) for r in rows]
    upper=[float(r['함량상한_pct']) for r in rows]
    ret=[('CENTRAL','MID','','',feasible_projection(lower,upper))]
    for i,part in enumerate(rows):
        lo=max(lower[i],100-sum(upper[j] for j in range(len(rows)) if j!=i))
        hi=min(upper[i],100-sum(lower[j] for j in range(len(rows)) if j!=i))
        if lo>hi+1e-7:raise ValueError('성분별 실현가능 구간 없음')
        for level,value in [('LOW',lo),('HIGH',hi)]:
            ret.append((part['CAS'],level,lo,hi,feasible_projection(lower,upper,fix=(i,value))))
    return ret


def anchors_hanna(ctx, eligible_keys,*,batch_size=16):
    """Real external HANNA is mandatory; no gamma=1 substitute for HANNA."""
    if not eligible_keys:return []
    import numpy as np
    from hanna_gamma_adapter import _get_model
    model=_get_model()
    bydoc=defaultdict(list)
    for r in ctx['parts']:bydoc[r['문서키']].append(r)
    outs=[]
    for key in sorted(eligible_keys):
        rows=bydoc[key];score=[]
        props=[approved_property(ctx['property'][r['CAS']]) for r in rows]
        twas=[approved_twa(ctx['twa'][r['CAS']]) for r in rows]
        if not all(props) or not all(twas):continue
        anchors=composition_anchors(rows)
        masses=np.asarray([a[4] for a in anchors]);mw=np.array([p[1] for p in props])
        mole=masses/mw[None,:];x=mole/mole.sum(axis=1,keepdims=True)
        ln=[]
        for offset in range(0,len(x),batch_size):
            result,_=model.predict([p[2] for p in props],x[offset:offset+batch_size],298.15)
            result=np.asarray(result,dtype=float)
            if result.shape!=(len(x[offset:offset+batch_size]),len(rows)) or not np.isfinite(result).all():
                raise ValueError('실제 HANNA 앵커 모델 출력 오류')
            ln.extend(result.tolist())
        gamma=np.exp(np.asarray(ln))
        if not np.isfinite(gamma).all() or (gamma<=0).any():
            raise ValueError('실제 HANNA 감마 유효성 오류')
        vp=np.asarray([p[0] for p in props]);twa=np.asarray(twas)
        rp=x*vp[None,:];hp=rp*gamma
        rs=(1e6*rp/(ATM*twa[None,:])).sum(axis=1)
        hs=(1e6*hp/(ATM*twa[None,:])).sum(axis=1)
        for i,anchor in enumerate(anchors):
            focus,level,lo,hi,mass=anchor
            outs.append(dict(문서키=key,표본종류='MSDS_제약_확정적_하중상',
              대상CAS=focus,수준=level,실현가능하한_pct=lo,실현가능상한_pct=hi,
              실제배합가정_pct='|'.join(f'{r["CAS"]}:{w:.6f}' for r,w in zip(rows,mass)),
              Raoult_혼합물지수=float(rs[i]),HANNA_혼합물지수=float(hs[i]),
              경고='각 CAS의 극단을 개별 고정하고 나머지를 100%로 보정; 동시 전성분 하한·상한 아님'))
    return outs


def uncertainty(products,scores,anchors):
    byprod=defaultdict(list)
    for s in anchors:byprod[s['문서키']].append(s)
    series={r['문서키']:r for r in products}
    env=[];bounds={}
    for key,p in sorted(series.items()):
        draws=scores.get(key,[]);aa=byprod.get(key,[])
        if not draws:continue
        if not aa:
            raise ValueError('실제 HANNA 하한·중간·상한 앵커 누락; 12 재실행 필요')
        v={}
        for model in ('Raoult','HANNA'):
            samples=[s[f'{model}_혼합물지수'] for s in draws]
            allv=samples+[s[f'{model}_혼합물지수'] for s in aa]
            v[model]=(min(allv),max(allv))
        mid=next((a for a in aa if a['수준']=='MID'),None)
        env.append(dict(문서키=key,용도ID=p['용도ID'],비교군ID=p['비교군ID'],
           실측_제품배합_알수없음='예',랜덤시나리오수=len(draws),
           조성앵커수=len(aa),Raoult_최저관찰지표=v['Raoult'][0],
           Raoult_중앙조성지표=mid['Raoult_혼합물지수'] if mid else '',
           Raoult_최고관찰지표=v['Raoult'][1],HANNA_최저관찰지표=v['HANNA'][0],
           HANNA_중앙조성지표=mid['HANNA_혼합물지수'] if mid else '',
           HANNA_최고관찰지표=v['HANNA'][1],
           한계='표본+성분별 실현가능 경계 앵커에서 관찰된 범위이며 수학적 전역 보증 아님'))
        bounds[key]=v
    bygroup=defaultdict(list)
    for row in products:bygroup[(row['용도ID'],row['비교군ID'])].append(row)
    pairs=[]
    for group,rs in sorted(bygroup.items()):
        for a,b in itertools.combinations(sorted(rs,key=lambda x:x['문서키']),2):
            ab=bounds[a['문서키']];bb=bounds[b['문서키']]
            classification={}
            for model in ('Raoult','HANNA'):
                alo,ahi=ab[model];blo,bhi=bb[model]
                classification[model]=('A_LOWER_B' if ahi < blo else
                                       'B_LOWER_A' if bhi < alo else 'OVERLAP_OR_TIE')
            agreed=(classification['Raoult']==classification['HANNA'] and
                    classification['Raoult']!='OVERLAP_OR_TIE')
            pairs.append(dict(용도ID=group[0],비교군ID=group[1],
                문서키_A=a['문서키'],문서키_B=b['문서키'],
                Raoult_시나리오범위판정=classification['Raoult'],
                HANNA_시나리오범위판정=classification['HANNA'],
                두모델_같은방향_비중첩='YES' if agreed else 'NO_OR_UNRESOLVED',
                해석='관찰된 표본/앵커 범위 분리일 뿐, 전 조성구간 순위 유지의 수학적 증명 아님'))
    return env,pairs


CENSUS_FIELDS=['문서키','기존제품ID','제조사','제품명','용도ID','비교군ID',
               '원본성분수','완전평형시나리오수','제품단위검증','보류사유']
PRODUCT_FIELDS=['문서키','기존제품ID','제조사','제품명','용도ID','비교군ID',
                '시나리오수','Raoult_혼합물지수_중앙값','HANNA_혼합물지수_중앙값',
                'Raoult_혼합물지수_최솟값','Raoult_혼합물지수_최댓값',
                'HANNA_혼합물지수_최솟값','HANNA_혼합물지수_최댓값',
                'Raoult_최대성분위험비_중앙값','HANNA_최대성분위험비_중앙값',
                'Raoult_순위','HANNA_순위','순위변화',
                'Raoult_상대등급','HANNA_상대등급','상대등급변화']
COMP_FIELDS=['문서키','시나리오ID','CAS','순물질_포화증기위험비',
             'Raoult_성분평형기여','HANNA_성분평형기여','HANNA_차이_기여']
GROUP_FIELDS=['용도ID','비교군ID','검증제품수','Raoult_최고대최저_배수',
              'HANNA_최고대최저_배수','Raoult_지표낮은제품','Raoult_지표높은제품',
              'HANNA_지표낮은제품','HANNA_지표높은제품']
RQ1_PAIR_FIELDS=['용도ID','비교군ID','문서키_A','문서키_B','Raoult_A대B_중앙값비',
                'HANNA_A대B_중앙값비','해석']
RQ2_FIELDS=['용도ID','비교군ID','문서키_A','문서키_B','Raoult_A순위',
            'Raoult_B순위','HANNA_A순위','HANNA_B순위','순위역전',
            '변화설명_최대기여CAS','기여차이_B빼기A','설명']
ANCHOR_FIELDS=['문서키','표본종류','대상CAS','수준','실현가능하한_pct',
               '실현가능상한_pct','실제배합가정_pct','Raoult_혼합물지수',
               'HANNA_혼합물지수','경고']
ENV_FIELDS=['문서키','용도ID','비교군ID','실측_제품배합_알수없음',
            '랜덤시나리오수','조성앵커수','Raoult_최저관찰지표',
            'Raoult_중앙조성지표','Raoult_최고관찰지표',
            'HANNA_최저관찰지표','HANNA_중앙조성지표',
            'HANNA_최고관찰지표','한계']
ROBUST_FIELDS=['용도ID','비교군ID','문서키_A','문서키_B',
               'Raoult_시나리오범위판정','HANNA_시나리오범위판정',
               '두모델_같은방향_비중첩','해석']


def analysis(work=WORK,data=DATA,output=None,*,with_anchors=True):
    """Produce 270-census and RQ1–RQ3 outputs; no approved data => explicit HOLD."""
    global OUT
    out=Path(output) if output is not None else Path(work)/'10_13_paper_questions'
    ctx=load_inputs(work,data)
    scores,delta,components=empirical_scores(ctx)
    counted,groups=census(ctx,scores)
    products=summarize_products(groups,scores)
    rq1pairs,rq2pairs,gsummary=compare_pairs(products,delta)
    if with_anchors:
        eligible_keys={r['문서키'] for r in products}
        anchors=anchors_hanna(ctx,eligible_keys)
    else:
        anchors=[]
    if products and not with_anchors:
        # Allow synthetic/test coverage without local HANNA weights, never pass as RQ3.
        bounds=[];rankstable=[]
    else:
        bounds,rankstable=uncertainty(products,scores,anchors)
    write(out/'00_전체270건_모델적용성및보류.csv',counted,CENSUS_FIELDS)
    write(out/'RQ1_같은용도_제품별평형지표.csv',products,PRODUCT_FIELDS)
    write(out/'RQ1_같은용도_최고최저배수.csv',gsummary,GROUP_FIELDS)
    write(out/'RQ1_같은용도_제품쌍비.csv',rq1pairs,RQ1_PAIR_FIELDS)
    write(out/'RQ1_RQ2_성분별증기위험기여.csv',components,COMP_FIELDS)
    write(out/'RQ2_Raoult_HANNA_제품쌍순위역전.csv',rq2pairs,RQ2_FIELDS)
    write(out/'RQ3_MSDS_실현가능_하중상_앵커.csv',anchors,ANCHOR_FIELDS)
    write(out/'RQ3_제품별_시나리오범위.csv',bounds,ENV_FIELDS)
    write(out/'RQ3_제품쌍_관찰범위_분리.csv',rankstable,ROBUST_FIELDS)
    flips=sum(r['순위역전']=='YES' for r in rq2pairs)
    counts=Counter(r['보류사유'].split('; ')[0] for r in counted if r['보류사유'])
    result=dict(source_documents=270,source_component_records=1454,
        unique_cas=len(ctx['property']),docs_with_all_approved_vle_and_twa_scenarios=len(scores),
        documents_in_same_use_verified_comparisons=len(products),
        verified_same_use_groups=len(gsummary),ranked_pairs_RQ1=len(rq1pairs),
        reversed_rank_pairs_RQ2=flips,
        relative_grade_changed_products_RQ2=sum(p['상대등급변화']=='YES' for p in products),
        question3_anchor_scenarios=len(anchors),
        empirically_range_separated_pairs_RQ3=sum(r['두모델_같은방향_비중첩']=='YES' for r in rankstable),
        question1_status='ANSWERABLE' if gsummary else 'ON_HOLD_NO_VERIFIED_SAME_USE_COMPARISON',
        question2_status='ANSWERABLE' if rq2pairs else 'ON_HOLD_NO_VERIFIED_PAIRS',
        question3_status=('EXPLORATORY_ENVELOPES_NOT_GLOBAL_PROOF' if rankstable else
                          'ON_HOLD_NO_VERIFIED_PAIRS_OR_ANCHORS'),
        top_blockers=dict(counts.most_common(7)),
        scientific_limitations=[
          'MSDS 공개 함량구간 기반 모델 시나리오는 실제 제조 배합의 통계적 표본이 아님',
          '순물질 포화증기위험비는 성분별 p_sat/P_atm/OEL, 혼합물 지수는 평형분압/OEL의 합계',
          '모든 지표는 밀폐 평형 상부공간 잠재치이며 작업장 노출/기준 초과/안전 인증이 아님',
          '같은 용도의 기능·공정 호환성은 MSDS/TDS 원문 및 제품 독립성 증거로 검증해야 함',
          '하/중/상 상대등급은 같은 비교군 내부 순위에만 기반하며 법적/독성학적 등급이 아님',
          '한계 조성의 모든 가능한 조합에 대한 수학적 보증은 아님',
          '승인된 한국 TWA가 없는 성분은 임의의 해외 OEL로 치환하지 않음',
        ])
    dump(out/'RQ_논문실험_성립여부_판정.json',result)
    md=['# 신나 270건: 논문 연구질문 1–3 분석 현황','',
        '전체 원본 270개를 모집단으로 유지하며 제외 문서는 보류 사유가 CSV에 남습니다.','',
        '## RQ1 · 같은 용도 제품 간 증기위험비와 혼합물 평형지표','',
        f'- 상태: **{result["question1_status"]}**; 검증된 비교군 {len(gsummary)}개; 비교 제품 {len(products)}개',
        '- 혼합물 지표 = Σ[1,000,000 × 성분 평형분압(kPa)/(101.325 × 한국 TWA(ppm))].',
        '- 순물질 증기위험비 = 1,000,000 × 해당 성분 25°C 포화증기압(kPa)/(101.325 × 한국 TWA(ppm)).',
        '- 혼합물 지표의 합계와 최대 성분 지표는 구별합니다.','',
        '## RQ2 · Raoult 대비 실제 HANNA 순위·상대등급 변동','',
        f'- 상태: **{result["question2_status"]}**; 제품쌍 역전 {flips}건 (검증 비교 가능 쌍 {len(rq2pairs)}건).',
        '- 순위는 동일 용도/수지/공정 비교군 내 제품별 시나리오 지표 중앙값에 기초.',
        '- 같은 비교군 내부 순위 백분위를 3개 상대등급으로 나눈 탐색적 등급변화만 제시. 법정 유해성 등급/노출기준 초과판정이 아님.','',
        '## RQ3 · MSDS 함량구간 하한·중간·상한 불확실성','',
        f'- 상태: **{result["question3_status"]}**; 앵커 {len(anchors)}개; 관찰범위 분리 제품쌍 {len(rankstable)}개 중 {result["empirically_range_separated_pairs_RQ3"]}개.',
        '- 각 CAS를 실현가능한 하한 또는 상한에 놓고 나머지 성분은 총 100%가 되도록 조정. 중앙조성은 제약 투영으로 구함.',
        '- 표본 및 앵커의 범위 분리는 모든 가능한 실제 배합에 대한 수학적 전역 보증은 아님.','',
        '## 논문 진행 결론','',
        ('- 세 질문을 평가할 동일 용도 검증제품쌍이 아직 없어 **질문은 미해결**입니다. '
         '원본 270개라는 표본 숫자만으로 결과가 있다고 주장할 수 없습니다.' if not rq1pairs else
         '- 검증된 제품쌍에 대해 세 연구질문 산출물을 생성했습니다. 모델 예측 정확도는 별도 실측 비교가 필요합니다.'),
        '- 다음 필수조건: 독립 제품·용도·TDS 호환성, 전성분 구조/25°C 물성, 현행 한국 TWA 확인, 실제 HANNA 가동 및 조성 앵커.',
        '', '## 주요 보류 원인', '']
    md += [f'- {reason}: {n}건' for reason,n in counts.most_common(8)]
    (out/'RQ_논문실험_성립여부_판정.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
    print('PAPER RQ ANALYSIS:',json.dumps({k:result[k] for k in (
          'source_documents','documents_in_same_use_verified_comparisons',
          'verified_same_use_groups','ranked_pairs_RQ1',
          'reversed_rank_pairs_RQ2','question3_anchor_scenarios',
          'question1_status','question2_status','question3_status')},ensure_ascii=False))
    return result


def main(argv=None):
    parser=argparse.ArgumentParser(description='Full 270 MSDS paper RQ1–RQ3; evidence-gated, actual HANNA only')
    parser.add_argument('--work',type=Path,default=WORK)
    parser.add_argument('--data',type=Path,default=DATA)
    parser.add_argument('--output',type=Path,default=None)
    args=parser.parse_args(argv)
    return analysis(args.work,args.data,args.output)

if __name__=='__main__':main()