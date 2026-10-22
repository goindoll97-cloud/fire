r"""
국가법령정보 Open API 통합 수집 스크립트
(위험물안전관리법 · 화학물질관리법 · 산업안전보건법 + 관련 고시)

사용법 (Windows 명령 프롬프트):
    cd /d D:\fire
    pip install requests
    python fetch_all.py              # 이미 받은 법령 JSON은 다시 받지 않음
    python fetch_all.py --refresh    # 전부 새로 받기

.env 파일(같은 폴더)에 다음 한 줄이 있어야 합니다.
    LAW_OC=신청한_이메일ID

결과 (D:\fire\law_data):
    법령/<법령명>.json              법령 본문 원본
    별표/<법령명>_별표N.txt          별표 텍스트 (서식은 _서식N 으로 따로 저장)
    별표_목록.csv                    별표·서식 전체 목록
    행정규칙_검색결과.csv            고시 검색 결과 (추가로 받을 고시 고를 때 사용)
    행정규칙/<고시명>.json, .txt     선택한 고시 본문
"""
import csv
import json
import re
import sys
import time
from pathlib import Path

import requests

BASE = Path(__file__).resolve().parent
OUT = BASE / "law_data"
API_SEARCH = "https://www.law.go.kr/DRF/lawSearch.do"
API_SERVICE = "https://www.law.go.kr/DRF/lawService.do"
SLEEP = 0.5

# ── 1. 받을 법령 (정확한 법령명) ─────────────────────────────
LAWS = [
    "위험물안전관리법",
    "위험물안전관리법 시행령",
    "위험물안전관리법 시행규칙",
    "화학물질관리법",
    "화학물질관리법 시행령",
    "화학물질관리법 시행규칙",
    "산업안전보건법",
    "산업안전보건법 시행령",
    "산업안전보건법 시행규칙",
    "산업안전보건기준에 관한 규칙",
]

# ── 2. 고시 검색어 (결과는 목록만 저장) ─────────────────────
ADMRUL_QUERIES = ["유해화학물질", "인체등유해성물질", "위험물", "노출기준", "화학물질"]

# ── 3. 본문까지 받을 고시: (저장이름, 행정규칙일련번호) ──────
#    추가하려면 행정규칙_검색결과.csv의 '행정규칙일련번호'를 넣으면 됩니다.
ADMRULS = [
    ("위험물안전관리에 관한 세부기준", "2100000285200"),
    ("유해화학물질 소량 취급시설에 관한 고시", "2100000262632"),
    ("유해화학물질의 규정수량에 관한 규정_2026-4_현행", "2100000278502"),
    ("유해화학물질의 규정수량에 관한 규정_2026-20_시행예정", "2100000286916"),
    ("화학물질 통계조사에 관한 규정", "2100000266554"),
    ("화학물질 조사결과 및 정보공개제도 운영에 관한 규정", "2100000269068"),
    ("화학물질 및 물리적 인자의 노출기준", "2100000186058"),
    ("화학물질의 분류·표시 및 물질안전보건자료에 관한 기준", "2100000278320"),
]


# ── 공통 함수 ────────────────────────────────────────────────
def read_oc() -> str:
    env = BASE / ".env"
    if not env.exists():
        sys.exit(f".env 파일이 없습니다: {env}")
    for line in env.read_text(encoding="utf-8-sig").splitlines():
        line = line.strip()
        if line.startswith("LAW_OC"):
            val = line.split("=", 1)[1].strip().strip('"').strip("'")
            if val:
                return val
    sys.exit(".env에서 LAW_OC 값을 찾지 못했습니다. 예: LAW_OC=abcd1234")


def get_json(url: str, params: dict, label: str = "") -> dict | None:
    try:
        r = requests.get(url, params=params, timeout=30)
        r.raise_for_status()
    except requests.RequestException as e:
        print(f"  [오류] {label}: {e}")
        return None
    try:
        return r.json()
    except ValueError:
        print(
            f"  [오류] {label}: JSON 응답이 아닙니다. OC 값, API 신청 항목, 등록 IP를 확인하세요.\n"
            f"         응답 앞부분: {r.text[:200]}"
        )
        return None


def as_list(x):
    if x is None:
        return []
    return x if isinstance(x, list) else [x]


def flatten(v) -> str:
    """문자열/리스트/딕셔너리가 섞인 값을 텍스트로 변환"""
    if v is None:
        return ""
    if isinstance(v, str):
        return v
    if isinstance(v, list):
        return "\n".join(flatten(i) for i in v)
    if isinstance(v, dict):
        return "\n".join(flatten(i) for i in v.values())
    return str(v)


def safe(name: str) -> str:
    return re.sub(r'[\\/:*?"<>|]', "_", name)


def num(x):
    s = str(x or "").strip()
    return int(s) if s.isdigit() else s


def write_csv(path: Path, rows: list[dict], cols: list[str]):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


# ── 1단계: 법령 본문 ─────────────────────────────────────────
def find_law(oc: str, name: str):
    data = get_json(
        API_SEARCH,
        {"OC": oc, "target": "law", "type": "JSON", "query": name, "display": 100},
        label=f"검색 {name}",
    )
    if not data:
        return None
    items = as_list(data.get("LawSearch", {}).get("law"))
    exact = [i for i in items if i.get("법령명한글", "").replace(" ", "") == name.replace(" ", "")]
    if not exact:
        return None
    cur = [i for i in exact if i.get("현행연혁코드") == "현행"] or exact
    return cur[0]


def fetch_laws(oc: str, refresh: bool):
    d = OUT / "법령"
    d.mkdir(parents=True, exist_ok=True)
    print("\n[1] 법령 본문")
    for name in LAWS:
        path = d / f"{safe(name)}.json"
        if path.exists() and not refresh:
            print(f"  [있음] {name}")
            continue
        meta = find_law(oc, name)
        if not meta:
            print(f"  [못 찾음] {name}")
            continue
        mst = meta.get("법령일련번호")
        body = get_json(
            API_SERVICE, {"OC": oc, "target": "law", "type": "JSON", "MST": mst}, label=name
        )
        if body:
            path.write_text(json.dumps(body, ensure_ascii=False, indent=1), encoding="utf-8")
            print(f"  [받음] {name} (시행일 {meta.get('시행일자')})")
        time.sleep(SLEEP)


# ── 2단계: 별표 추출 (별표/서식 구분) ───────────────────────
def extract_annexes():
    src, dst = OUT / "법령", OUT / "별표"
    dst.mkdir(parents=True, exist_ok=True)
    for old in dst.glob("*.txt"):  # 예전 버전 파일 정리
        old.unlink()
    rows, kinds = [], {}
    for jf in sorted(src.glob("*.json")):
        name = jf.stem
        law = json.loads(jf.read_text(encoding="utf-8")).get("법령", {})
        info = law.get("기본정보", {}) or {}
        for i, u in enumerate(as_list((law.get("별표") or {}).get("별표단위"))):
            kind = (u.get("별표구분") or "").strip() or "구분없음"
            no, sub = num(u.get("별표번호")), num(u.get("별표가지번호"))
            tag = f"{kind}{no}" + (f"의{sub}" if isinstance(sub, int) and sub > 0 else "")
            title = (u.get("별표제목") or "").strip()
            text = flatten(u.get("별표내용"))
            fn = f"{safe(name)}_{tag}.txt"
            if (dst / fn).exists():
                fn = f"{safe(name)}_{tag}_{i}.txt"
            (dst / fn).write_text(f"{title}\n\n{text}", encoding="utf-8")
            kinds[kind] = kinds.get(kind, 0) + 1
            rows.append({
                "법령명": name,
                "시행일자": info.get("시행일자", ""),
                "구분": kind,
                "번호": tag,
                "제목": title,
                "텍스트길이": len(text),
                "파일": fn,
                "원본파일링크": u.get("별표서식파일링크", "") or "",
            })
    write_csv(OUT / "별표_목록.csv", rows,
              ["법령명", "시행일자", "구분", "번호", "제목", "텍스트길이", "파일", "원본파일링크"])
    empty = sum(1 for r in rows if r["텍스트길이"] < 200)
    print(f"\n[2] 별표 추출: {len(rows)}개 {kinds} (텍스트 거의 없음 {empty}개)")


# ── 3단계: 고시 검색 목록 ────────────────────────────────────
def search_admrul(oc: str):
    rows, seen = [], set()
    for q in ADMRUL_QUERIES:
        data = get_json(
            API_SEARCH,
            {"OC": oc, "target": "admrul", "type": "JSON", "query": q, "display": 100},
            label=f"고시 검색 {q}",
        )
        for i in as_list((data or {}).get("AdmRulSearch", {}).get("admrul")):
            sid = i.get("행정규칙일련번호")
            if sid in seen:
                continue
            seen.add(sid)
            rows.append({
                "검색어": q,
                "행정규칙명": i.get("행정규칙명"),
                "종류": i.get("행정규칙종류"),
                "발령번호": i.get("발령번호"),
                "시행일자": i.get("시행일자"),
                "소관부처": i.get("소관부처명"),
                "행정규칙ID": i.get("행정규칙ID"),
                "행정규칙일련번호": sid,
            })
        time.sleep(SLEEP)
    write_csv(OUT / "행정규칙_검색결과.csv", rows,
              ["검색어", "행정규칙명", "종류", "발령번호", "시행일자", "소관부처", "행정규칙ID", "행정규칙일련번호"])
    print(f"\n[3] 고시 검색: {len(rows)}건 (중복 제거)")


# ── 4단계: 선택한 고시 본문 ──────────────────────────────────
def fetch_admrul(oc: str, refresh: bool):
    d = OUT / "행정규칙"
    d.mkdir(parents=True, exist_ok=True)
    print("\n[4] 고시 본문")
    for name, sid in ADMRULS:
        jp = d / f"{safe(name)}.json"
        if jp.exists() and not refresh:
            print(f"  [있음] {name}")
            continue
        data = get_json(API_SERVICE, {"OC": oc, "target": "admrul", "type": "JSON", "ID": sid}, label=name)
        if not data:
            continue
        text = flatten(data)
        jp.write_text(json.dumps(data, ensure_ascii=False, indent=1), encoding="utf-8")
        (d / f"{safe(name)}.txt").write_text(text, encoding="utf-8")
        print(f"  [받음] {name}: {len(text):,}자")
        time.sleep(SLEEP)


def main():
    refresh = "--refresh" in sys.argv
    oc = read_oc()
    OUT.mkdir(exist_ok=True)
    fetch_laws(oc, refresh)
    extract_annexes()
    search_admrul(oc)
    fetch_admrul(oc, refresh)
    print(f"\n완료 -> {OUT}\nlaw_data 폴더를 zip으로 압축해서 올려주세요.")


if __name__ == "__main__":
    main()
