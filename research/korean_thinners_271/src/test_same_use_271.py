# -*- coding: utf-8 -*-
"""Standard-library, model-free safety regressions for conditional comparisons."""
import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import screen_same_use_271 as app


class SameUseGuardTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        base = Path(self.tmp.name)
        self.data = base / 'research' / 'korean_thinners_271' / 'data'
        self.data.mkdir(parents=True)
        self.in_dir = base / 'results'
        self.in_dir.mkdir()
        self.out_dir = base / 'reports'
        self.props = base / 'thinner_vle_reference_25C.csv'
        self.props.write_text('CAS,MW_g_mol,VP_25C_kPa\nC1,100,1\nC2,80,2\n', encoding='utf-8')
        self.products = [
            {'제품ID': pid, '제품명(한글)': '가상 '+pid,
             '1차 용도분류': group, '분류 신뢰도': '명칭명시', '분류 근거': 'fixture',
             '제1항·제품명 불일치': '원문 불일치' if pid == 'C' else ''}
            for pid, group in [('A', '우레탄계'), ('B', '우레탄계'), ('C', '우레탄계')]
        ]
        self.components = [{'제품ID': pid, 'CAS 번호': cas} for pid in 'ABC' for cas in ['C1', 'C2']]
        self.write(self.data / '제품271건_용도재분류.csv', self.products)
        self.write(self.data / '성분별함량_1462.csv', self.components)
        scenarios, details = [], []
        for pid, factor in [('A', 1.), ('B', 2.), ('C', 3.)]:
            for sid in range(1, 21):
                scenarios.append({'제품ID': pid, '시나리오ID': sid,
                                  '총압력_Raoult_kPa': 1 * factor, '총압력_HANNA_kPa': 1.2 * factor,
                                  '모델_상태': 'HANNA_2026_실제모델_예측'})
                for cas, share in [('C1', .25), ('C2', .75)]:
                    details.append({'제품ID': pid, '시나리오ID': sid, 'CAS': cas,
                                    'HANNA_분압_kPa': 1.2 * factor * share})
        self.write(self.in_dir / 'hanna_라울_시나리오별_비교.csv', scenarios)
        self.write(self.in_dir / 'hanna_라울_성분별_비교.csv', details)
        meta = {'input_sha256': {'classification': app.sha(self.data / '제품271건_용도재분류.csv'),
                                 'components': app.sha(self.data / '성분별함량_1462.csv'),
                                 'property_table': app.sha(self.props)},
                'temperature_K': 298.15, 'model': '2026 multicomponent HANNA local ensemble'}
        (self.in_dir / 'hanna_라울_실행정보.json').write_text(json.dumps(meta), encoding='utf-8')
        self.review = base / 'review.csv'
        self.oel = base / 'oel.csv'
        self.write(self.review, [{'제품ID': pid, '비교군ID': 'urethane_verified',
                                  '용도검증': '승인', '공정호환성검증': '승인',
                                  '중복개정검토': '승인', '근거문서': '가상 fixture 검증'} for pid in 'ABC'])
        self.write(self.oel, [{'CAS': cas, '한국_TWA_ppm': twa, '기준시점': '2026-10-01',
                               '공식출처': 'fixture://not-an-official-source', '검증상태': '승인'}
                              for cas, twa in [('C1', '100'), ('C2', '50')]])
        self.patcher = patch.multiple(app, DATA=self.data, REPO_ROOT=base)
        self.patcher.start()
        self.addCleanup(self.patcher.stop)

    @staticmethod
    def write(path, rows):
        with path.open('w', newline='', encoding='utf-8-sig') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    def cmd(self, *args):
        app.main(['--input-dir', str(self.in_dir), '--output-dir', str(self.out_dir), *args])

    def read(self, name):
        return app.read_csv(self.out_dir / name)

    def test_no_manual_approval_means_no_hazard_ranks(self):
        self.cmd('--write-templates')
        self.assertEqual(len(self.read('271_제품별_분석가능성_감사.csv')), 3)
        self.assertEqual(len(self.read('제품군별_증기상_기술비교.csv')), 3)
        self.assertEqual(self.read('동일용도_저유해성_후보선별.csv'), [])
        self.assertEqual(len(self.read('국내TWA_공식검증입력_서식.csv')), 2)

    def test_incomplete_oel_blocks_all_ranking(self):
        self.write(self.oel, [{'CAS': 'C1', '한국_TWA_ppm': '100', '기준시점': '2026-10-01',
                               '공식출처': 'fixture', '검증상태': '승인'}])
        self.cmd('--review-csv', str(self.review), '--verified-oel-csv', str(self.oel))
        self.assertEqual(self.read('동일용도_저유해성_후보선별.csv'), [])
        self.assertIn('C2', self.read('271_제품별_분석가능성_감사.csv')[0]['차단사유'])

    def test_two_approved_product_compared_and_mismatch_excluded(self):
        self.cmd('--review-csv', str(self.review), '--verified-oel-csv', str(self.oel))
        haz = self.read('동일용도_저유해성_후보선별.csv')
        self.assertEqual(len(haz), 2)
        self.assertEqual([(r['제품ID'], r['참고순위']) for r in haz], [('A', '1'), ('B', '2')])
        pairs = self.read('동일용도_제품쌍_잠재지표비교.csv')
        self.assertEqual(len(pairs), 1)
        self.assertAlmostEqual(float(pairs[0]['중앙값_높음나눔낮음']), 2)
        self.assertEqual(self.read('271_제품별_분석가능성_감사.csv')[2]['평가상태'], '보류')

    def test_metadata_mismatch_blocks_rank(self):
        meta = self.in_dir / 'hanna_라울_실행정보.json'
        d = json.loads(meta.read_text(encoding='utf-8'))
        d['input_sha256']['components'] = 'mismatch'
        meta.write_text(json.dumps(d), encoding='utf-8')
        self.cmd('--review-csv', str(self.review), '--verified-oel-csv', str(self.oel))
        self.assertEqual(self.read('동일용도_저유해성_후보선별.csv'), [])

    def test_partial_pressure_inconsistent_rejected(self):
        path = self.in_dir / 'hanna_라울_성분별_비교.csv'
        rows = app.read_csv(path)
        rows[0]['HANNA_분압_kPa'] = '999'
        self.write(path, rows)
        with self.assertRaisesRegex(ValueError, '성분분압 합계 불일치'):
            self.cmd()

    def test_optimizer_requires_evidence_before_any_search(self):
        from optimize_verified_blend import template, load_config
        path = Path(self.tmp.name) / 'optim.json'
        path.write_text(json.dumps(template(), ensure_ascii=False), encoding='utf-8')
        with self.assertRaisesRegex(ValueError, '최적화 실행 보류'):
            load_config(path)


if __name__ == '__main__':
    unittest.main()