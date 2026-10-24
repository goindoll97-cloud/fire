# -*- coding: utf-8 -*-
"""Regression, feasibility and nonlegal-gate tests for inverse design."""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

import optimize_msds_gap as opt
import compare_thinner_raoult_unifac_hanna as core

HERE=Path(__file__).resolve().parent


def fixture_law_files(folder: Path):
    """Reproduce the 7 actual OEL table rows from prior Korean crosswalk.

    This is a schema/logic fixture, not a real law-document integration test.
    Run the CLI with the original law_data/*.json to perform that test.
    """
    # Minimal lines reproducing CAS / TWA cells verified against the user's
    # Korean crosswalk. No external generated CSV is needed to run unit tests.
    official_rows = [
      ('71-43-2', '벤젠', 'Benzene', 0.5, '2.5'),
      ('67-56-1', '메틸 알코올', 'Methanol', 200, '250'),
      ('67-64-1', '아세톤', 'Acetone', 500, '750'),
      ('123-86-4', '노말-초산 부틸', 'n-Butyl acetate', 150, '200'),
      ('141-78-6', '초산 에틸', 'Ethyl acetate', 400, '-'),
      ('108-88-3', '톨루엔', 'Toluene', 50, '150'),
      ('108-10-1', '메틸 이소부틸 케톤', 'Methyl isobutyl ketone', 50, '75'),
      ('108-10-1', '헥손', 'Hexone', 50, '75'),
    ]
    oel_rows=[f'│{i} │{ko} │{en} │FORM │{twa} │- │{stel} │- │[{cas}] │'
              for i,(cas,ko,en,twa,stel) in enumerate(official_rows,1)]
    oel={'AdmRulService':{
      '행정규칙기본정보':{'발령번호':'2020-48','시행일자':'20200116'},
      '별표':{'별표단위':[
        {'별표번호':'0001','별표키':'000100','별표제목':'화학물질의 노출기준','별표내용':[oel_rows]},
        {'별표번호':'0001','별표키':'000102','별표제목':'삭제','별표내용':'<별표 1-2> <삭제>'}]}}}
    msds={'AdmRulService':{
      '행정규칙기본정보':{'발령번호':'2026-26','시행일자':'20260424'},
      '조문내용':{'x':'제11조(작성원칙) 제9항 별표 6'},
      '별표':{'별표단위':[{'별표번호':'0006','별표키':'000600',
        '별표제목':'건강 및 환경 유해성 분류에 대한 한계농도 기준',
        '별표내용':[['┃건강 │7. 발암성 │0.1%┃','┃건강 │1. 급성 독성 │1%┃']]}]}}}
    a=folder/'oel.json';b=folder/'msds.json'
    a.write_text(json.dumps(oel,ensure_ascii=False),encoding='utf-8')
    b.write_text(json.dumps(msds,ensure_ascii=False),encoding='utf-8')
    return a,b


class AlgorithmTests(unittest.TestCase):
    def setUp(self):
        self.base=[('108-10-1',47.),('108-88-3',53.)]
        self.p0=opt.normalized_base(self.base)
        self.policy=opt.BoundsPolicy()
        self.props=core.load_reference(opt.ROOT/'thinner_vle_reference_25C.csv')
        self.predictor=core.GammaModels(opt.ROOT/'HANNA_NOT_INSTALLED')
        self.twa={'108-10-1':50.,'108-88-3':50.,'71-43-2':0.5}

    def test_mixture_mass_balance_and_fixed_trace(self):
        p=self.p0.copy()
        comps=opt.mixture(self.base,'71-43-2',.099,p)
        self.assertAlmostEqual(sum(w for _,w in comps),100.,places=10)
        self.assertAlmostEqual(dict(comps)['71-43-2'],.099,places=12)
        self.assertTrue(opt.feasible(p,self.p0,self.policy))
        with self.assertRaisesRegex(ValueError,'already'):
            opt.mixture(self.base,'108-10-1',.099,p)

    def test_reallocation_limit_rejects_extreme_compositions(self):
        p=np.array([.20,.80])
        self.assertFalse(opt.feasible(p,self.p0,self.policy))
        self.assertFalse(opt.feasible(np.array([.30,.69]),self.p0,self.policy))

    def test_raoult_optimal_never_worse_and_constraints(self):
        p0,pnew,initial,best,audit=opt.optimize_one(self.base,'71-43-2',.099,'raoult',
            self.predictor,self.props,self.twa,self.policy,starts=2,maxiter=40,seed=7)
        self.assertTrue(opt.feasible(pnew,p0,self.policy))
        self.assertGreaterEqual(best['KR_share_pct']+1e-7,initial['KR_share_pct'])
        self.assertTrue(0. <= .5*np.abs(pnew-p0).sum()*100. <= 10.00001)
        self.assertGreaterEqual(len(audit),2)
        self.assertEqual(best['KR_TWA_ppm'],.5)

    def test_registry_unreviewed_does_not_get_automatic_limit(self):
        registry=opt.multi.load_registry(opt.ROOT/'multichem_screening_registry.csv')
        rows=['┃건강 │7. 발암성 │0.1%┃']
        v,kind=opt.trace_cutoff('71-43-2',registry,rows,None)
        self.assertEqual(v,.1)
        self.assertIn('PRELIMINARY',kind)
        with self.assertRaisesRegex(ValueError,'No verified numeric cutoff'):
            opt.trace_cutoff('67-56-1',registry,rows,None)
        v,kind=opt.trace_cutoff('67-56-1',registry,rows,1.)
        self.assertEqual(v,1.)
        self.assertIn('NOT_REGULATORY',kind)

    def test_missing_hanna_fails_rather_than_raoult_fallback(self):
        comps=opt.mixture(self.base,'71-43-2',.099,self.p0)
        with self.assertRaisesRegex(RuntimeError,'HANNA'):
            opt.evaluate(comps,'71-43-2','hanna',self.predictor,self.props,self.twa)

    def test_real_law_json_integration_if_present(self):
        if not opt.kr.OEL_JSON.is_file() or not opt.kr.MSDS_JSON.is_file():
            self.skipTest('Live Korean law JSON absent; sample-format fixture tested instead')
        with tempfile.TemporaryDirectory() as folder:
            args=opt.parse_args(['--output-dir',str(Path(folder)/'out'),
                                 '--optimize-model','raoult','--check-models','raoult',
                                 '--starts','1','--maxiter','30',
                                 '--doses','.099','--product-ids','1183746'])
            result=opt.run(args)
            self.assertEqual(result['successful_scenarios'],1)
            self.assertEqual(result['failed_scenarios'],0)

    def test_end_to_end_law_fixture_raoult(self):
        with tempfile.TemporaryDirectory() as folder:
            folder=Path(folder)
            oj,mj=fixture_law_files(folder)
            args=opt.parse_args([
              '--oel-json',str(oj),'--msds-json',str(mj),
              '--output-dir',str(folder/'out'),
              '--optimize-model','raoult','--check-models','raoult',
              '--starts','1','--maxiter','35','--doses','.099',
              '--product-ids','1183746'])
            result=opt.run(args)
            self.assertEqual(result['successful_scenarios'],1)
            self.assertEqual(result['failed_scenarios'],0)
            summary=pd.read_csv(folder/'out'/'02_optimization_summary.csv')
            self.assertAlmostEqual(summary.iloc[0].dose_wt_pct,.099)
            self.assertGreaterEqual(summary.iloc[0].optimized_KR_share_pct,
                                    summary.iloc[0].baseline_KR_share_pct-1e-7)
            report=pd.read_csv(folder/'out'/'04_cross_model_checks.csv')
            self.assertEqual(len(report),2)
            self.assertEqual(set(report.status),{'OK'})
            composition=pd.read_csv(folder/'out'/'03_baseline_vs_optimized_compositions.csv')
            self.assertAlmostEqual(composition.optimized_trace_wt_pct.sum(),100.,places=6)
            self.assertEqual(json.loads((folder/'out'/'06_run_metadata.json').read_text())['legal_omission_assessed'],False)

if __name__=='__main__':
    unittest.main()
