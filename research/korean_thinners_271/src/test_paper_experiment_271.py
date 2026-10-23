# -*- coding: utf-8 -*-
"""No model checkpoints, no external downloads; run on GitHub Actions stdlib."""
import csv
import json
import math
import random
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import composition_sampling as sampler
import paper_experiment_271 as paper
import validate_external_vle_271 as external


class CompositionTests(unittest.TestCase):
    def test_feasible_box_simplex_many_samples(self):
        bounds=[{'함량 하한(%)':'10','함량 상한(%)':'60'},
                {'함량 하한(%)':'20','함량 상한(%)':'70'},
                {'함량 하한(%)':'5','함량 상한(%)':'60'}]
        for seed in (13,42,20261023):
            samples=sampler.sample_bounded_simplex(bounds,180,random.Random(seed),burn_in=80,thin=8)
            self.assertEqual(len(samples),180)
            for row in samples:
                self.assertAlmostEqual(sum(row),100,places=7)
                for x,b in zip(row,bounds):
                    self.assertGreaterEqual(x,float(b['함량 하한(%)'])-1e-8)
                    self.assertLessEqual(x,float(b['함량 상한(%)'])+1e-8)
        self.assertNotEqual(samples[0],samples[1])

    def test_degenerate_and_invalid_bounds(self):
        yes=[{'함량 하한(%)':'30','함량 상한(%)':'30'},
             {'함량 하한(%)':'70','함량 상한(%)':'70'}]
        self.assertEqual(sampler.sample_bounded_simplex(yes,4,random.Random(1)),[[30.,70.]]*4)
        bad=[{'함량 하한(%)':'80','함량 상한(%)':'90'},
             {'함량 하한(%)':'30','함량 상한(%)':'50'}]
        with self.assertRaises(ValueError):
            sampler.sample_bounded_simplex(bad,2,random.Random(1))


class FullPaperQualityGates(unittest.TestCase):
    def setUp(self):
        self.t=tempfile.TemporaryDirectory()
        self.addCleanup(self.t.cleanup)
        root=Path(self.t.name)
        self.research=root/'research'/'korean_thinners_271'
        d=self.research/'data'
        d.mkdir(parents=True)
        self.props=root/'thinner_vle_reference_25C.csv'
        self.cls=d/'제품271건_용도재분류.csv'
        self.comp=d/'성분별함량_1462.csv'
        self.save(self.props,[{'CAS':'C1','MW_g_mol':100,'VP_25C_kPa':2,'SMILES':'CC'},
                              {'CAS':'C2','MW_g_mol':50,'VP_25C_kPa':4,'SMILES':'CO'}])
        self.save(self.cls,[{'제품ID':'P1','제품명(한글)':'FAKE SAMPLE','1차 용도분류':'테스트'}])
        self.save(self.comp,[{'제품ID':'P1','CAS 번호':'C1','함량 하한(%)':20,'함량 상한(%)':80},
                             {'제품ID':'P1','CAS 번호':'C2','함량 하한(%)':20,'함량 상한(%)':80}])
        self.patcher=patch.multiple(paper,RESEARCH_ROOT=self.research,REPO_ROOT=root,
                                   CLASSIFICATION=self.cls,COMPONENTS=self.comp,PROPERTY=self.props)
        self.patcher.start()
        self.addCleanup(self.patcher.stop)
        self.runroot=root/'output'/'runs'
        self.out=root/'output'/'reports'

    @staticmethod
    def save(path,rows):
        path.parent.mkdir(exist_ok=True,parents=True)
        with path.open('w',newline='',encoding='utf-8-sig') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

    def make_run(self,sampler='legacy',seed=1,bad=False):
        d=self.runroot/f'{sampler}_seed_{seed}'
        d.mkdir(exist_ok=True,parents=True)
        scenarios=[];components=[]
        for i,w in enumerate([20.,35.,50.,80.],1):
            v=[w,100-w]
            mol=[v[0]/100,v[1]/50]
            x=[m/sum(mol) for m in mol]
            vp=[2.,4.]
            gamma=[1.2,1.3]
            raoult=[a*b for a,b in zip(x,vp)]
            hanna=[a*b for a,b in zip(raoult,gamma)]
            rt=sum(raoult);ht=sum(hanna)
            scenarios.append({'제품ID':'P1','시나리오ID':i,'총압력_Raoult_kPa':round(rt,8),
                              '총압력_HANNA_kPa':round(ht,8),
                              'HANNA_대_Raoult_비':round(ht/rt,8),'HANNA_압력차_kPa':round(ht-rt,8),
                              '모델_상태':'HANNA_2026_실제모델_예측'})
            for cas,weight,xi,g,pr,ph in zip(['C1','C2'],v,x,gamma,raoult,hanna):
                components.append({'제품ID':'P1','시나리오ID':i,'CAS':cas,'질량백분율':round(weight,8),
                                   '액상몰분율':round(xi,9),'HANNA_ln_gamma':round(math.log(g),9),
                                   'HANNA_gamma':round(g,9),'Raoult_분압_kPa':round(pr,8),
                                   'HANNA_분압_kPa':round(ph if not(bad and i==1 and cas=='C2') else 9,8),
                                   'Raoult_기상몰분율':round(pr/rt,9),
                                   'HANNA_기상몰분율':round(ph/ht,9)})
        self.save(d/paper.SCENARIO_FILE,scenarios)
        self.save(d/paper.COMPONENT_FILE,components)
        meta={'input_sha256':paper.expected_sha(),'model':'2026 multicomponent HANNA local ensemble',
              'temperature_K':298.15,'seed':seed,'sampler':sampler,
              'documents_calculated':1,'documents_excluded':0,'samples_requested_per_document':4}
        (d/paper.METADATA_FILE).write_text(json.dumps(meta),encoding='utf8')
        return d

    def test_repeated_seed_and_sampler_summary_and_gates(self):
        for sampler in ('legacy','gibbs'):
            for seed in (1,2):self.make_run(sampler,seed)
        out=paper.summarize(self.runroot,self.out,top_tie_margin=.01)
        self.assertEqual(out['scenarios_checked'],16)
        self.assertEqual(out['component_records_checked'],32)
        self.assertEqual(out['seeds_per_sampler'],{'gibbs':2,'legacy':2})
        self.assertEqual(len(paper.read_csv(self.out/'paper_제품별_안정성.csv')),2)
        self.assertEqual(len(paper.read_csv(self.out/'paper_표본추출_민감도.csv')),1)
        self.assertIn('NOT VALIDATED',out['external_VLE_validation'])

    def test_detect_fake_hanna_pressure(self):
        self.make_run(bad=True)
        with self.assertRaisesRegex(ValueError,'partial pressure mismatch'):
            paper.summarize(self.runroot,self.out,.01)

    def test_reject_mismatched_provenance(self):
        d=self.make_run()
        info=d/paper.METADATA_FILE
        o=json.loads(info.read_text())
        o['input_sha256']['components']='broken'
        info.write_text(json.dumps(o))
        with self.assertRaisesRegex(ValueError,'SHA'):
            paper.summarize(self.runroot,self.out,.01)

    def test_external_vle_requires_manual_property_evidence(self):
        base={r['CAS']:r for r in paper.read_csv(self.props)}
        rows=[{'reference_id':'R1','source_url_or_doi':'https://example.org/real-doi',
               'verification_status':'승인','temperature_K':'298.15',
               'CAS_JSON':'["C1","C2"]','liquid_x_JSON':'[0.5,0.5]',
               'observed_total_pressure_kPa':'4'}]
        with self.assertRaisesRegex(ValueError,'property evidence'):
            external.validate_rows(rows,[],base)
        evidence=[{'CAS':'C1','temperature_K':'298.15','VP_25C_kPa':'2',
                   'reference_url_or_doi':'https://example.org/VP','verification_status':'승인'},
                  {'CAS':'C2','temperature_K':'298.15','VP_25C_kPa':'4',
                   'reference_url_or_doi':'https://example.org/VP','verification_status':'승인'}]
        self.assertEqual(len(external.validate_rows(rows,evidence,base)),1)
        rows[0]['verification_status']='미검증'
        with self.assertRaisesRegex(ValueError,'manual approval'):
            external.validate_rows(rows,evidence,base)

if __name__=='__main__':unittest.main()