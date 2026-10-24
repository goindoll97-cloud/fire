# -*- coding: utf-8 -*-
"""성분별 인화점(FP)·폭발한계(LFL/UFL) 표를 만든다 -> data/pure_flammability.csv

우선순위:
 1) 문헌 앵커(ANCHORS): 라이브러리에 값이 없는 성분을 공급사 SDS·ICSC·OSHA 등에서 직접 확인한 값
 2) chemicals.safety 표 (IEC 60079-20-1, NFPA 497, Serat/DIPPR 2017, Wikidata 순)
 3) 추정값: LFL은 Crowl–Louvar(1990) 화학양론식, FP는 유사 물질 유추 (FP_source에 'ESTIMATE' 표기)
추정값을 쓴 성분은 estimated=True로 남겨 민감도 분석에서 따로 다룬다.
"""
import csv
from pathlib import Path

from chemicals import safety
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

V2 = Path(__file__).resolve().parent.parent

# (FP_C, FP_출처, LFL_%, UFL_%, 폭발한계 출처)  None = 해당 없음
# GESTIS 항목은 2026-10-24 사용자 PC의 내장 브라우저로 gestis-database.dguv.de에서 직접 확인한 폐쇄식 인화점이다.
# '<' 로 보고된 값은 그 상한을 사용한다(보수적으로 화재 위험이 낮게 잡힘).
ANCHORS = {
    '616-38-6': (18.0, 'ICSC 1080 (open cup 18°C; closed cup 17°C 이하로 보고되나 개방식 값 사용)', 4.2, 12.9, 'ICSC 1080'),
    '107-98-2': (36.1, 'OSHA Chemical Data 72 (97°F)', 1.6, 13.8, 'OSHA Chemical Data 72'),
    '2110-78-3': (45.0, 'TCI H0256 SDS', None, None, ''),
    '106-65-0': (94.0, 'GESTIS (DGUV) 492698, closed cup', 1.0, 8.5, 'GESTIS (DGUV) 492698'),
    '1119-40-0': (97.0, 'GESTIS (DGUV) 530305, closed cup', None, None, ''),
    '627-93-0': (107.0, 'GESTIS (DGUV) 032130, closed cup', 0.8, 8.1, 'GESTIS (DGUV) 032130'),
    '104-51-8': (58.0, 'NITE-CHRIP (GESTIS cited), closed cup', None, None, ''),
    '107-83-5': (-7.0, 'GESTIS (DGUV) 037070, closed cup, 보고값 < -7°C (상한으로 사용)', 1.2, None, 'GESTIS (DGUV) 037070'),
    '589-34-4': (-11.0, 'GESTIS (DGUV) 490240, closed cup, ca. -11°C', 1.0, 7.0, 'GESTIS (DGUV) 490240'),
    # 출처 미확인: 유사 물질 유추 (estimated)
    '108-08-7': (-20.0, 'GESTIS (DGUV) 496176, closed cup, 보고값 < -20°C (상한으로 사용)', 0.9, 6.9, 'GESTIS (DGUV) 496176'),
    '527-53-7': (63.0, 'GESTIS (DGUV) 570009 (1,2,3,5-tetramethylbenzene), closed cup', None, None, ''),
    '112-40-3': (80.0, 'GESTIS (DGUV) 013560, closed cup', 0.6, None, 'GESTIS (DGUV) 013560'),
}


def formula_atoms(smiles):
    from rdkit.Chem import rdMolDescriptors
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    atoms = {}
    for a in m.GetAtoms():
        atoms[a.GetSymbol()] = atoms.get(a.GetSymbol(), 0) + 1
    return atoms


def crowl_louvar_pct(a):
    """LFL(vol%) = 0.55 * Cst, Cst = 1/(1 + 4.76*(C + (H-X)/4 - O/2)) (Crowl & Louvar, 2nd ed.)"""
    x = sum(a.get(h, 0) for h in ('F', 'Cl', 'Br', 'I'))
    n_o2 = a.get('C', 0) + (a.get('H', 0) - x) / 4 - a.get('O', 0) / 2
    return 100 * 0.55 / (1 + 4.76 * n_o2)


def source_grade(fp_src, est):
    """A: 표준·공공 데이터(IEC, NFPA, DIPPR/Serat, ICSC, OSHA, GESTIS 인용 NITE). B: 공급사 SDS·위키 등 2차 출처. C: 추정."""
    if est or fp_src.startswith('ESTIMATE'):
        return 'C'
    if fp_src.startswith(('IEC', 'NFPA', 'Serat', 'ICSC', 'OSHA', 'NITE', 'GESTIS')):
        return 'A'
    return 'B'


def main():
    out = []
    with (V2 / 'data' / 'pure_properties.csv').open(encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        cas, smi = r['CAS'], r['SMILES']
        fp = fp_src = lfl = ufl = lfl_src = None
        est = False
        if cas in ANCHORS:
            fp, fp_src, lfl, ufl, lfl_src = ANCHORS[cas]
            est = fp_src.startswith('ESTIMATE')
        if fp is None:
            t = safety.T_flash(CASRN=cas)
            if t is not None:
                fp = t - 273.15
                ms = safety.T_flash_methods(cas)
                fp_src = ms[0] if ms else 'chemicals.safety'
        if lfl is None:
            l = safety.LFL(CASRN=cas)
            if l is not None:
                lfl = l * 100
                u = safety.UFL(CASRN=cas)
                ufl = u * 100 if u is not None else None
                lfl_src = 'chemicals.safety (NFPA/IEC/Serat)'
        if lfl is None:
            lfl = crowl_louvar_pct(formula_atoms(smi))
            lfl_src = 'ESTIMATE: Crowl-Louvar LFL=0.55*Cst (화학양론)'
            est = True
        out.append({'CAS': cas, 'name': r['name'], 'FP_C': '' if fp is None else round(fp, 1),
                    'FP_source': fp_src or '', 'LFL_pct': '' if lfl is None else round(lfl, 3),
                    'UFL_pct': '' if ufl is None else round(ufl, 2), 'LFL_source': lfl_src or '',
                    'estimated': est, 'FP_grade': source_grade(fp_src or '', fp_src is not None and fp_src.startswith('ESTIMATE'))})
    cols = list(out[0])
    with (V2 / 'data' / 'pure_flammability.csv').open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, cols)
        w.writeheader()
        w.writerows(out)
    miss = [o['name'] for o in out if o['FP_C'] == '' or o['LFL_pct'] == '']
    print(f'{len(out)} 화학종, 인화점/LFL 누락: {miss}')


if __name__ == '__main__':
    main()
