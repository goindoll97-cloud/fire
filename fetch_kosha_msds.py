r"""
KOSHA MSDS Open API 수집 스크립트 (D:\fire 에서 실행)

목적
1) 모델에 쓰는 용제 7종의 KOSHA MSDS에서 물리화학적 특성(9번 항목: 인화점, 폭발한계, 증기압 등)을 받아
   논문의 물성값 출처를 'KOSHA MSDS'로 통일한다.
2) '시너', '신나', '세척제', '희석제'로 검색해 혼합물(시판 제품) MSDS가 이 API에 있는지 확인한다.

.env (같은 폴더)
    KOSHA_MSDS_SERVICE_KEY=공공데이터포털_인증키
    * 공공데이터포털 '한국산업안전보건공단_물질안전보건자료 조회 서비스'(15157612) 활용신청 키 (planner와 같은 키)
    * 인코딩/디코딩 키 어느 쪽을 넣어도 됩니다.

실행
    cd /d D:\fire
    python fetch_kosha_msds.py

결과 (D:\fire\kosha_msds)
    raw/            응답 원문 XML (문제 생기면 이 폴더를 올려주세요)
    chem_ids.csv    CAS → chemId 매칭 결과
    section09.csv   물리화학적 특성 항목 (물질, 항목명, 값)
    product_search.csv  '시너' 등 검색 결과 목록
"""
import csv
import sys
import time
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path

import requests

BASE = Path(__file__).resolve().parent
OUT = BASE / "kosha_msds"
RAW = OUT / "raw"
# 2026년 공공데이터포털 서비스 이전 후 주소 (planner 저장소 engine/__init__.py와 동일)
API = "https://apis.data.go.kr/B552468/msdschem"
OP_LIST = "getChemList"
OP_DETAIL = "getChemDetail{:02d}"
# data.go.kr 게이트웨이가 python-requests 기본 User-Agent를 봇으로 막으므로 브라우저 UA 사용
HEADERS = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
                         "(KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36"}
SLEEP = 0.4

SUBSTANCES = {
    "108-88-3": "톨루엔",
    "1330-20-7": "크실렌",
    "100-41-4": "에틸벤젠",
    "78-93-3": "메틸에틸케톤",
    "141-78-6": "아세트산에틸",
    "67-64-1": "아세톤",
    "67-56-1": "메탄올",
    # 확장 후보 (문헌에서 시너 성분으로 자주 검출)
    "123-86-4": "아세트산부틸",
    "108-10-1": "메틸이소부틸케톤",
    "67-63-0": "이소프로필알코올",
    "110-54-3": "노말헥산",
    "71-43-2": "벤젠",
}
PRODUCT_WORDS = ["시너", "신나", "세척제", "희석제", "락카"]


# 윈도우는 파일명 끝의 마침표를 지우므로 "env."는 실제로 "env"로 저장됨 → 둘 다 찾음
ENV_NAMES = [".env", "env.", "env", ".env.txt", "env.txt"]


def find_env() -> Path:
    for n in ENV_NAMES:
        p = BASE / n
        if p.is_file():
            return p
    sys.exit(f"설정 파일이 없습니다. {BASE} 폴더에 다음 중 하나가 필요합니다: {', '.join(ENV_NAMES)}")


def read_env_value(name: str) -> str:
    env = find_env()
    for line in env.read_text(encoding="utf-8-sig").splitlines():
        line = line.strip()
        if "=" in line and line.split("=", 1)[0].strip() == name:
            val = line.split("=", 1)[1].strip().strip('"').strip("'")
            # 인증키는 '='나 '.'로 끝날 수 있어 마침표를 지우지 않음. 앞뒤 공백만 제거
            if val:
                return val
    sys.exit(f"{env.name}에서 {name} 값을 찾지 못했습니다. 예: {name}=값")


def read_key() -> str:
    # 인코딩 키가 들어와도 디코딩해서 사용
    return urllib.parse.unquote(read_env_value("KOSHA_MSDS_SERVICE_KEY"))


def call(op: str, params: dict, key: str, tag: str):
    q = {"serviceKey": key, **params}
    try:
        r = requests.get(f"{API}/{op}", params=q, timeout=30, headers=HEADERS)
    except requests.RequestException as e:
        print(f"  [오류] {op} {tag}: {e}")
        return None
    (RAW / f"{op}_{tag}.xml").write_text(r.text, encoding="utf-8")
    try:
        root = ET.fromstring(r.content)
    except ET.ParseError:
        print(f"  [오류] {op} {tag}: XML 아님 -> {r.text[:200]}")
        return None
    code = (root.findtext(".//resultCode") or root.findtext(".//returnReasonCode") or "").strip()
    if code and code not in ("0", "00", "0000") and not code.upper().startswith("NORMAL"):
        msg = root.findtext(".//resultMsg") or root.findtext(".//returnAuthMsg") or ""
        print(f"  [API 오류] {op} {tag}: {code} {msg}")
        return None
    return root


def items(root):
    out = []
    for it in root.iter():
        if it.tag.rsplit("}", 1)[-1].lower() == "item":
            row = {c.tag.rsplit("}", 1)[-1]: (c.text or "").strip() for c in it}
            if row:
                out.append(row)
    return out


def pick(row, *names):
    for n in names:
        if row.get(n):
            return row[n]
    return ""


def main():
    key = read_key()
    RAW.mkdir(parents=True, exist_ok=True)

    # 1) CAS → chemId
    print("[1] CAS로 chemId 찾기")
    ids = []
    for cas, ko in SUBSTANCES.items():
        root = call(OP_LIST, {"searchWrd": cas, "searchCnd": "1", "numOfRows": "10", "pageNo": "1"}, key, cas)
        if root is None:
            continue
        rows = [r for r in items(root) if pick(r, "casNo", "cas") == cas] or items(root)
        if not rows:
            print(f"  [없음] {ko} {cas}")
            continue
        r = rows[0]
        cid = pick(r, "chemId", "chemNo", "chem_id", "msdsId")
        ids.append({"CAS": cas, "물질": ko, "chemId": cid, "KOSHA명": pick(r, "chemNameKor", "chemNmKor", "chemName"),
                    "최종갱신": pick(r, "lastDate", "updtDt")})
        print(f"  {ko} {cas} -> chemId {cid}")
        time.sleep(SLEEP)
    with open(OUT / "chem_ids.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["CAS", "물질", "chemId", "KOSHA명", "최종갱신"])
        w.writeheader()
        w.writerows(ids)

    # 2) 9번 항목(물리화학적 특성)
    print("\n[2] 물리화학적 특성(9번 항목)")
    sec = []
    for r in ids:
        root = call(OP_DETAIL.format(9), {"chemId": r["chemId"], "chemNo": r["chemId"]}, key, r["CAS"])
        if root is None:
            continue
        for it in items(root):
            name = pick(it, "msdsItemNameKor", "itemNameKor", "itemNmKor", "msdsItemName", "itemName")
            val = pick(it, "itemDetail", "detail", "contents", "content", "value")
            sec.append({"CAS": r["CAS"], "물질": r["물질"], "항목": name, "값": val})
        print(f"  {r['물질']}: {sum(1 for s in sec if s['CAS'] == r['CAS'])}개 항목")
        time.sleep(SLEEP)
    with open(OUT / "section09.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["CAS", "물질", "항목", "값"])
        w.writeheader()
        w.writerows(sec)

    # 3) 혼합물 제품이 있는지 이름 검색
    print("\n[3] 제품명 검색 (혼합물 MSDS 존재 여부 확인)")
    prod = []
    for wd in PRODUCT_WORDS:
        root = call(OP_LIST, {"searchWrd": wd, "searchCnd": "0", "numOfRows": "100", "pageNo": "1"}, key, f"name_{wd}")
        if root is None:
            continue
        rows = items(root)
        total = root.findtext(".//totalCount")
        print(f"  '{wd}': 총 {total}건")
        for x in rows:
            prod.append({"검색어": wd, **x})
        time.sleep(SLEEP)
    cols = sorted({k for p in prod for k in p}) or ["검색어"]
    with open(OUT / "product_search.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(prod)

    print(f"\n완료 -> {OUT}\nkosha_msds 폴더를 zip으로 압축해서 올려주세요.")


if __name__ == "__main__":
    main()
