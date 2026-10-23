# -*- coding: utf-8 -*-
"""MSDS 271건 연구 CSV 구조 검증. 추출·분류 정확도의 원문 수작업 검증을 대신하지 않음."""
import csv
import collections
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

def read(name):
    with (DATA / name).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def main():
    p = read("제품목록_271.csv")
    c = read("성분별함량_1462.csv")
    r = read("제품271건_용도재분류.csv")
    g = read("제품군별_현황.csv")
    z = read("물성보완_62CAS.csv")
    assert len(p) == 271 and len(c) == 1462 and len(r) == 271
    assert len(g) == 14 and len(z) == 62
    product_ids = {row["제품ID"] for row in p}
    assert len(product_ids) == 271
    assert {row["제품ID"] for row in r} == product_ids
    assert all(row["제품ID"] in product_ids for row in c)
    assert len({row["CAS 번호"] for row in c}) == 62
    assert sum(int(row["MSDS 문서 건수"]) for row in g) == 271
    assert sum(row["분자모델 1차 후보"] == "예" for row in r) == 46
    assert sum(row["기존물성 라울예비 대상"] == "예" for row in r) == 6
    assert sum(row["MSDS 개정일"] == "미확인" for row in r) == 8
    counts = collections.Counter(x["제품ID"] for x in c)
    for row in p:
        pid = row["제품ID"]
        assert counts[pid] == int(row["성분 기록수"]), (pid, counts[pid])
    print("통과: 271 PDF 문서 / 1462 성분 / 62 CAS / 14 제품군 / HANNA 예비 46 / 라울 6")
    print("주의: 원문 추출오류, 혼합물 식별, 최신 국내 TWA 및 물성 검증은 별도입니다.")

if __name__ == "__main__":
    main()
