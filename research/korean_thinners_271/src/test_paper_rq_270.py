# -*- coding: utf-8 -*-
"""Synthetic regression: no external HANNA checkpoint, CAS approvals or internet."""
from collections import defaultdict
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).parent))
import paper_rq_270 as rq

class ThreeQuestionTests(unittest.TestCase):
    def test_bounded_simplex_anchors_are_feasible(self):
        rows=[{'CAS':'A','함량하한_pct':'30','함량상한_pct':'70'},
              {'CAS':'B','함량하한_pct':'20','함량상한_pct':'60'},
              {'CAS':'C','함량하한_pct':'10','함량상한_pct':'30'}]
        anchors=rq.composition_anchors(rows)
        self.assertEqual(len(anchors),7)
        self.assertEqual(anchors[0][1],'MID')
        for target,level,lo,hi,m in anchors:
            self.assertAlmostEqual(sum(m),100,places=6)
            for v,r in zip(m,rows):
                self.assertGreaterEqual(v,float(r['함량하한_pct'])-1e-7)
                self.assertLessEqual(v,float(r['함량상한_pct'])+1e-7)
            if level in ('LOW','HIGH'):
                self.assertAlmostEqual(m[[r['CAS'] for r in rows].index(target)],
                                       lo if level=='LOW' else hi,places=6)
        with self.assertRaises(ValueError):
            rq.composition_anchors([{'CAS':'A','함량하한_pct':'75','함량상한_pct':'85'},
                                    {'CAS':'B','함량하한_pct':'50','함량상한_pct':'70'}])

    def test_rq1_ratios_rq2_rank_and_relative_tier_flip(self):
        def p(key,r,h):
            return dict(문서키=key,기존제품ID=key,용도ID='우레탄',비교군ID='동일공정',
                        Raoult_혼합물지수_중앙값=r,HANNA_혼합물지수_중앙값=h)
        products=[p('A',1,4),p('B',2,1),p('C',3,2)]
        d={'A':{'CAS1':[3]},'B':{'CAS1':[-1]},'C':{'CAS1':[-1]}}
        pairs,flips,summary=rq.compare_pairs(products,d)
        self.assertEqual(len(pairs),3)
        self.assertEqual(summary[0]['Raoult_최고대최저_배수'],3)
        self.assertEqual(summary[0]['HANNA_최고대최저_배수'],4)
        self.assertEqual(sum(x['순위역전']=='YES' for x in flips),2)
        self.assertEqual(next(x for x in flips if x['문서키_A']=='A' and x['문서키_B']=='B')['변화설명_최대기여CAS'],'CAS1')
        self.assertEqual(next(x for x in products if x['문서키']=='A')['상대등급변화'],'YES')

    def test_rq3_observed_disjoint_not_claimed_global_guarantee(self):
        products=[{'문서키':c,'용도ID':'우레탄','비교군ID':'same'} for c in ('A','B')]
        scores={'A':[{'Raoult_혼합물지수':1.0,'HANNA_혼합물지수':1.1}],
                'B':[{'Raoult_혼합물지수':5.0,'HANNA_혼합물지수':6.0}]}
        anchors=[]
        for key,r,h in [('A',1.1,1.2),('B',5.1,6.2)]:
            anchors.append({'문서키':key,'수준':'MID','Raoult_혼합물지수':r,'HANNA_혼합물지수':h})
        env,pairs=rq.uncertainty(products,scores,anchors)
        self.assertEqual(pairs[0]['두모델_같은방향_비중첩'],'YES')
        self.assertIn('수학적',pairs[0]['해석'])
        scores['B'][0]['HANNA_혼합물지수']=1.15
        env,pairs=rq.uncertainty(products,scores,anchors)
        self.assertEqual(pairs[0]['두모델_같은방향_비중첩'],'NO_OR_UNRESOLVED')

    def test_rq1_rq2_full_census_with_real_model_provenance_fixture(self):
        # Source-linked synthetic 270-document universe; no external HANNA
        # model is called, but results must carry actual model provenance.
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);data=root/'data';work=root/'work';data.mkdir()
            review=work/'02_review';pred=work/'05_predictions'
            review.mkdir(parents=True);pred.mkdir(parents=True)
            docs=[dict(문서키=f'D{i:03}',기존제품ID=f'X{i:03}',
                제조사='SYNTH',제품명_파일명기반=f'paper{i}',
                함량합계100가능_표기성분만='예') for i in range(270)]
            parts=[]
            for i,d in enumerate(docs):
                for cas in (('A','B') if i<2 else ('C',)):
                    parts.append(dict(문서키=d['문서키'],CAS=cas,
                        CAS분류='순물질후보_미검증',함량하한_pct='30',함량상한_pct='70'))
            parts += [dict(문서키='D002',CAS='C',CAS분류='순물질후보_미검증',
                           함량하한_pct='0',함량상한_pct='100') for _ in range(1454-len(parts))]
            products=[]
            for i,d in enumerate(docs):
                rev=dict(문서키=d['문서키'],용도ID='우레탄' if i<2 else '',
                         비교군ID='same_resin' if i<2 else '',
                         근거문서='SYNTH_Fixture' if i<2 else '')
                for flag in rq.REQUIRED_REVIEW:rev[flag]='승인' if i<2 else '미검증'
                products.append(rev)
            props=[];twas=[]
            for cas in ('A','B','C'):
                props.append(dict(CAS=cas,분자량_g_mol='100',**{'25C_포화증기압_kPa':'8'},
                    온도_K='298.15',SMILES='CC',복합물또는이성질체='아니오',
                    구조검증상태='승인',물성검증상태='승인',HANNA입력승인='승인',
                    구조공식출처='SYNTH_Fixture',물성출처와측정조건='SYNTH_Fixture'))
                twas.append(dict(CAS=cas,한국_TWA_ppm='100',기준시점='2026-10-24',
                    고시명번호='SYNTH',공식출처='SYNTH',공식원문대조='승인',검증상태='승인'))
            rq.write(data/'MSDS_원본대조_270.csv',docs,list(docs[0]))
            rq.write(data/'MSDS_구성성분_1454.csv',parts,list(parts[0]))
            rq.write(review/'제품_원문TDS_검증.csv',products,list(products[0]))
            rq.write(review/'CAS_물성_검증.csv',props,list(props[0]))
            rq.write(review/'CAS_한국TWA_검증.csv',twas,list(twas[0]))
            ps=[];cs=[]
            for key in ('D000','D001'):
                for i in range(20):
                    ideal=(1,1) if key=='D000' else (1.5,1.5)
                    han=(1.9,1.9) if key=='D000' else (1.2,1.2)
                    ps.append(dict(문서키=key,시나리오ID=i+1,
                        Raoult_총압_kPa=sum(ideal),HANNA_총압_kPa=sum(han)))
                    for j,cas in enumerate(('A','B')):
                        cs.append(dict(문서키=key,시나리오ID=i+1,CAS=cas,
                            Raoult_분압_kPa=ideal[j],HANNA_분압_kPa=han[j]))
            rq.write(pred/'HANNA_Raoult_시나리오.csv',ps,list(ps[0]))
            rq.write(pred/'HANNA_Raoult_성분.csv',cs,list(cs[0]))
            rq.dump(pred/'prediction_provenance.json',dict(
                model='actual external HANNA 2026 ensemble',
                property_review_sha256=rq.sha(review/'CAS_물성_검증.csv'),
                product_review_sha256=rq.sha(review/'제품_원문TDS_검증.csv'),
                source_sha256={n:rq.sha(data/n) for n in
                  ('MSDS_원본대조_270.csv','MSDS_구성성분_1454.csv')}))
            result=rq.analysis(work,data,with_anchors=False)
            self.assertEqual(result['source_documents'],270)
            self.assertEqual(result['documents_in_same_use_verified_comparisons'],2)
            self.assertEqual(result['ranked_pairs_RQ1'],1)
            self.assertEqual(result['reversed_rank_pairs_RQ2'],1)
            self.assertEqual(result['relative_grade_changed_products_RQ2'],2)
            self.assertEqual(result['question3_status'],'ON_HOLD_NO_VERIFIED_PAIRS_OR_ANCHORS')
            product_rows=rq.read(work/'10_13_paper_questions'/'RQ1_같은용도_제품별평형지표.csv')
            self.assertEqual(len(product_rows),2)
            self.assertGreater(float(product_rows[0]['HANNA_혼합물지수_중앙값']),
                               float(product_rows[1]['HANNA_혼합물지수_중앙값']))
            self.assertLess(float(product_rows[0]['Raoult_혼합물지수_중앙값']),
                            float(product_rows[1]['Raoult_혼합물지수_중앙값']))

    def test_main_preserves_270_even_when_no_approved_scores(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); data=root/'data'; work=root/'work';data.mkdir();
            (work/'02_review').mkdir(parents=True);(work/'05_predictions').mkdir(parents=True)
            docs=[{'문서키':f'D{i:03}','기존제품ID':f'X{i:03}',
                   '제조사':'TEST','제품명_파일명기반':f'test{i}',
                   '함량합계100가능_표기성분만':'예'} for i in range(270)]
            parts=[{'문서키':d['문서키'],'CAS':'TESTCAS','CAS분류':'순물질후보_미검증'} for d in docs]
            parts += [{'문서키':'D000','CAS':'TESTCAS','CAS분류':'순물질후보_미검증'} for _ in range(1454-len(parts))]
            reviews=[{'문서키':d['문서키'],'용도ID':'검증전','비교군ID':'',
                      'MSDS제1항용도검증':'미검증','독립제품검증':'미검증'} for d in docs]
            props=[{'CAS':'TESTCAS','복합물또는이성질체':'아니오','HANNA입력승인':'미승인'}]
            twa=[{'CAS':'TESTCAS','검증상태':'미검증','공식원문대조':'미검증'}]
            rq.write(data/'MSDS_원본대조_270.csv',docs,list(docs[0]))
            rq.write(data/'MSDS_구성성분_1454.csv',parts,list(parts[0]))
            rq.write(work/'02_review'/'제품_원문TDS_검증.csv',reviews,list(reviews[0]))
            rq.write(work/'02_review'/'CAS_물성_검증.csv',props,list(props[0]))
            rq.write(work/'02_review'/'CAS_한국TWA_검증.csv',twa,list(twa[0]))
            rq.write(work/'05_predictions'/'HANNA_Raoult_시나리오.csv',[],['문서키','시나리오ID'])
            rq.write(work/'05_predictions'/'HANNA_Raoult_성분.csv',[],['문서키','시나리오ID'])
            result=rq.analysis(work,data)
            self.assertEqual(result['source_documents'],270)
            self.assertEqual(result['documents_in_same_use_verified_comparisons'],0)
            self.assertEqual(result['question1_status'],'ON_HOLD_NO_VERIFIED_SAME_USE_COMPARISON')
            out=work/'10_13_paper_questions'
            self.assertEqual(len(rq.read(out/'00_전체270건_모델적용성및보류.csv')),270)
            self.assertEqual(len(rq.read(out/'RQ2_Raoult_HANNA_제품쌍순위역전.csv')),0)
            self.assertEqual(len(rq.read(out/'RQ3_제품쌍_관찰범위_분리.csv')),0)

if __name__=='__main__':unittest.main()