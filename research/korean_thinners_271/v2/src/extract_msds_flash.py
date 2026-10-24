# -*- coding: utf-8 -*-
"""MSDS 제9항(물리화학적 특성)의 인화점·폭발범위를 PDF 텍스트에서 추출한다.

사용: python extract_msds_flash.py <MSDS PDF 폴더 1> [폴더 2 ...]
출력: ../data/msds_flash.csv  (기존제품ID, 인화점_C, 인화점_원문, 하한_pct, 상한_pct, 폭발범위_원문)
원본 PDF는 pdftotext -layout으로 텍스트화한다. 인화점 값이 없거나 숫자가 아니면 빈 칸으로 두고 원문을 남긴다.
"""
import csv
import re
import subprocess
import sys
from pathlib import Path

V2 = Path(__file__).resolve().parent.parent
REPO = V2.parent
RECON = REPO / 'data' / 'source_rebuild_20261023' / 'MSDS_원본대조_270.csv'

NUM = r'(-?\d+(?:\.\d+)?)'
FP_LINE = re.compile(r'(?:인화점(?:\(℃\))?|Flash[ \t]*point)[ \t]*[:：]?[ \t]*([^\n]*)', re.I)


def pdf_text(p):
    return subprocess.run(['pdftotext', '-layout', str(p), '-'], capture_output=True,
                          text=True, errors='ignore').stdout


def flash(text):
    """첫 번째 숫자형 인화점 줄을 반환 (값, 원문). 제9항 이전의 소화 안내문은 숫자가 없어 걸러진다."""
    for m in FP_LINE.finditer(text):
        raw = m.group(1).strip()
        n = re.search(NUM, raw)
        if n and len(raw) < 40:
            return float(n.group(1)), raw
    return None, ''


def limits(text):
    m = re.search(r'폭발\s*범위의\s*상한\s*/\s*하한\s*[:：]?\s*([^\n]*)', text)
    if m:
        raw = m.group(1).strip()
        nums = re.findall(NUM, raw)
        if len(nums) >= 2:
            return float(nums[1]), float(nums[0]), raw      # 상한/하한 순서
        return None, None, raw
    m = re.search(r'Upper\s+expolosion\s+limit\s*/\s*Lower\s+explosion\s+limit\s*[:：]?\s*([^\n]*)', text, re.I)
    if m:
        raw = m.group(1).strip()
        nums = re.findall(NUM, raw)
        if len(nums) >= 2:
            return float(nums[1]), float(nums[0]), raw
        return None, None, raw
    return None, None, ''


def main(dirs):
    files = {}
    for d in dirs:
        for p in Path(d).iterdir():
            if p.suffix.lower() == '.pdf':
                files[p.stem.lower()] = p
    out = []
    with RECON.open(encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            stem = Path(r['원본PDF파일명']).stem.lower()
            p = files.get(stem)
            if p is None:   # 삼화 파일명은 코드 앞부분만 일치하기도 한다
                code = stem.split('_')[0]
                p = next((q for k, q in files.items() if k.split('_')[0] == code), None)
            if p is None:
                out.append({'기존제품ID': r['기존제품ID'], 'PDF': '', '인화점_C': '', '인화점_원문': 'PDF 없음'})
                continue
            t = pdf_text(p)
            fp, fraw = flash(t)
            up, lo, lraw = limits(t)
            out.append({'기존제품ID': r['기존제품ID'], 'PDF': p.name, '인화점_C': fp if fp is not None else '',
                        '인화점_원문': fraw, '하한_pct': lo if lo is not None else '',
                        '상한_pct': up if up is not None else '', '폭발범위_원문': lraw})
    dst = V2 / 'data' / 'msds_flash.csv'
    cols = ['기존제품ID', 'PDF', '인화점_C', '인화점_원문', '하한_pct', '상한_pct', '폭발범위_원문']
    with dst.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, cols)
        w.writeheader()
        w.writerows(out)
    n = sum(1 for o in out if o['인화점_C'] != '')
    print(f'{len(out)} 제품, 인화점 숫자 추출 {n}건 -> {dst}')


if __name__ == '__main__':
    main(sys.argv[1:])
