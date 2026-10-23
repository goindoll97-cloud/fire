# -*- coding: utf-8 -*-
"""Evidence-gated screening of Korean thinner documents by verified same-use group.

Only standard library is required. This is equilibrium *potential*, never a
workplace exposure or regulatory compliance calculation. No hazard ranking is
produced without explicit use/duplicate review AND complete verified Korean OELs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

RESEARCH_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[3]
DATA = RESEARCH_ROOT / "data"
ATM_KPA = 101.325
TRUE_VALUES = {"승인"}  # deliberate, no implicit approval by 'yes'/'true'
BLOCKED_CATEGORIES = {"명칭만으로 판별 불가", "검색 오탐·비신나"}

AUDIT_COLUMNS = ["제품ID", "제품명", "원분류", "분류신뢰도", "분류근거", "모델결과존재", "모델시나리오수",
                 "모델검증정보일치", "비교군ID", "동일용도검증", "공정호환성검증", "중복개정검토",
                 "국내TWA전성분검증", "평가상태", "차단사유"]
PHYSICAL_COLUMNS = ["제품ID", "제품명", "원분류", "시나리오수", "Raoult_중앙값_kPa", "HANNA_중앙값_kPa",
                    "HANNA_5백분위_kPa", "HANNA_95백분위_kPa", "압력비_중앙값", "비교해석"]
HAZARD_COLUMNS = ["비교군ID", "제품ID", "제품명", "시나리오수", "최대TWA정규화평형분압_중앙값",
                  "최대TWA정규화평형분압_최솟값", "최대TWA정규화평형분압_최댓값", "최대TWA정규화평형분압_5백분위",
                  "최대TWA정규화평형분압_95백분위", "참고순위", "해석상태"]
PAIR_COLUMNS = ["비교군ID", "지표_낮은제품ID", "지표_높은제품ID", "중앙값_높음나눔낮음", "시나리오범위_분리", "해석"]
SCENARIO_COLUMNS = ["제품ID", "비교군ID", "시나리오ID", "최대TWA정규화평형분압", "최대기여성분_CAS"]
OEL_HEADERS = ["CAS", "한국_TWA_ppm", "기준시점", "공식출처", "검증상태"]
REVIEW_HEADERS = ["제품ID", "비교군ID", "용도검증", "공정호환성검증", "중복개정검토", "근거문서"]


def read_csv(path, required=True):
    path = Path(path)
    if not path.exists():
        if required:
            raise FileNotFoundError(str(path))
        return []
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def save_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def num(value, name, positive=False):
    try:
        x = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name}: 숫자 오류 {value!r}") from exc
    if not math.isfinite(x) or (positive and x <= 0):
        raise ValueError(f"{name}: 유한한 양수 필요")
    return x


def qtile(v, q):
    a = sorted(v)
    t = (len(a) - 1) * q
    i = int(t)
    return a[i] + (a[min(i + 1, len(a) - 1)] - a[i]) * (t - i)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def report_provenance(input_dir, classification, component_data):
    path = input_dir / "hanna_라울_실행정보.json"
    if not path.is_file():
        return False, "실행정보 JSON 없음"
    try:
        d = json.loads(path.read_text(encoding="utf-8-sig"))
        hashes = d["input_sha256"]
        match = (hashes["classification"] == sha(classification)
                 and hashes["components"] == sha(component_data)
                 and hashes["property_table"] == sha(REPO_ROOT / "thinner_vle_reference_25C.csv")
                 and d["model"] == "2026 multicomponent HANNA local ensemble"
                 and abs(num(d["temperature_K"], "온도") - 298.15) < 1e-6)
        return match, "원본 SHA-256 또는 모델·온도 불일치" if not match else ""
    except (KeyError, ValueError, OSError, json.JSONDecodeError) as exc:
        return False, f"실행정보 검증 불가: {exc}"


def build_templates(target, products, components):
    review = [{"제품ID": p["제품ID"], "비교군ID": "", "용도검증": "미검증",
               "공정호환성검증": "미검증", "중복개정검토": "미검증", "근거문서": ""} for p in products]
    cases = sorted({c.get("CAS 번호", "").strip() for c in components if c.get("CAS 번호", "").strip()})
    oel = [{"CAS": cas, "한국_TWA_ppm": "", "기준시점": "", "공식출처": "", "검증상태": "미검증"} for cas in cases]
    save_csv(target / "동일용도_검증입력_서식.csv", review, REVIEW_HEADERS)
    save_csv(target / "국내TWA_공식검증입력_서식.csv", oel, OEL_HEADERS)


def load_verified_oel(path):
    verified = {}
    if path is None:
        return verified
    for r in read_csv(path):
        cas = r.get("CAS", "").strip()
        if not cas:
            continue
        if cas in verified:
            raise ValueError(f"검증 TWA CAS 중복: {cas}")
        if r.get("검증상태", "").strip() != "승인":
            continue
        # Human-reviewed dates and official links; never infer 'current' from a URL alone.
        if not r.get("공식출처", "").strip() or not r.get("기준시점", "").strip():
            continue
        verified[cas] = num(r.get("한국_TWA_ppm"), f"{cas} TWA", positive=True)
    return verified


def load_review(path, products):
    if path is None:
        return {}
    ids = {p["제품ID"] for p in products}
    reviews = {}
    for r in read_csv(path):
        pid = r.get("제품ID", "").strip()
        if pid not in ids:
            raise ValueError(f"검증 파일의 알 수 없는 제품ID: {pid}")
        if pid in reviews:
            raise ValueError(f"검증 파일 제품ID 중복: {pid}")
        reviews[pid] = r
    return reviews


def validate_prediction_rows(scenarios, details):
    by_sc = {}
    by_detail = defaultdict(list)
    for r in scenarios:
        key = (r.get("제품ID", "").strip(), str(r.get("시나리오ID", "")).strip())
        if not key[0] or key in by_sc:
            raise ValueError(f"시나리오 제품ID/키 중복 또는 누락: {key}")
        if r.get("모델_상태") != "HANNA_2026_실제모델_예측":
            raise ValueError(f"HANNA 실제모델 출력으로 표기되지 않음: {key}")
        pressure = num(r["총압력_HANNA_kPa"], f"{key} 총압력", positive=True)
        num(r["총압력_Raoult_kPa"], "Raoult 압력", positive=True)
        by_sc[key] = (r, pressure)
    for r in details:
        key = (r.get("제품ID", "").strip(), str(r.get("시나리오ID", "")).strip())
        if key not in by_sc:
            raise ValueError(f"시나리오가 없는 성분행: {key}")
        by_detail[key].append(r)
    for key, (_, total) in by_sc.items():
        rows = by_detail[key]
        if len(rows) < 2 or len({r.get("CAS", "").strip() for r in rows}) != len(rows):
            raise ValueError(f"성분 누락·CAS 중복: {key}")
        pressures = [num(r["HANNA_분압_kPa"], f"{key} 성분분압") for r in rows]
        if any(p < 0 for p in pressures) or not math.isclose(sum(pressures), total, abs_tol=1e-5, rel_tol=1e-6):
            raise ValueError(f"성분분압 합계 불일치: {key}")
    return by_sc, by_detail


def main(argv=None):
    ap = argparse.ArgumentParser(description="동일 용도 신나의 증기상 잠재성: 확인된 제품만 조건부 비교")
    ap.add_argument("--input-dir", type=Path, default=RESEARCH_ROOT / "output")
    ap.add_argument("--output-dir", type=Path, default=RESEARCH_ROOT / "output")
    ap.add_argument("--review-csv", type=Path, help="명시적으로 승인된 용도/공정/중복 검토 CSV")
    ap.add_argument("--verified-oel-csv", type=Path, help="현행 한국 TWA를 수작업 검증한 승인 CSV")
    ap.add_argument("--write-templates", action="store_true", help="271개 검토행과 62개 CAS TWA 빈 서식 출력")
    ap.add_argument("--min-scenarios", type=int, default=20)
    args = ap.parse_args(argv)
    if args.min_scenarios < 1:
        ap.error("--min-scenarios must be >= 1")
    classification = DATA / "제품271건_용도재분류.csv"
    component_data = DATA / "성분별함량_1462.csv"
    products = read_csv(classification)
    component_rows = read_csv(component_data)
    if args.write_templates:
        build_templates(args.output_dir, products, component_rows)
        print("미승인 검토 서식 2개 생성 (승인 상태는 자동으로 채우지 않음)")
    scenarios = read_csv(args.input_dir / "hanna_라울_시나리오별_비교.csv", required=False)
    details = read_csv(args.input_dir / "hanna_라울_성분별_비교.csv", required=False)
    by_sc, by_detail = validate_prediction_rows(scenarios, details) if scenarios or details else ({}, {})
    valid_meta, meta_issue = report_provenance(args.input_dir, classification, component_data)
    rev = load_review(args.review_csv, products)
    oel = load_verified_oel(args.verified_oel_csv)
    by_product_sc = defaultdict(list)
    for key in by_sc:
        by_product_sc[key[0]].append(key)
    product_details_cas = defaultdict(set)
    for key, rows in by_detail.items():
        product_details_cas[key[0]].update(r.get("CAS", "").strip() for r in rows)
    classification_ids = {p["제품ID"] for p in products}
    unknown = set(by_product_sc) - classification_ids
    if unknown:
        raise ValueError(f"분류표에 없는 계산 제품ID: {sorted(unknown)}")

    audit, physical, scores, scenario_scores = [], [], [], []
    for p in products:
        pid = p["제품ID"]
        name = p.get("제품명(한글)", "")
        original_group = p.get("1차 용도분류", "")
        keys = sorted(by_product_sc[pid], key=lambda z: int(z[1]))
        group_review = rev.get(pid, {})
        group = group_review.get("비교군ID", "").strip()
        reasons = []
        if not keys:
            reasons.append("HANNA 계산결과 없음")
        elif len(keys) < args.min_scenarios:
            reasons.append(f"시나리오 부족({len(keys)}<{args.min_scenarios})")
        if not valid_meta:
            reasons.append(meta_issue)
        if original_group in BLOCKED_CATEGORIES:
            reasons.append("용도 미확정 또는 비대상")
        if p.get("제1항·제품명 불일치", "").strip():
            reasons.append("MSDS 제1항/제품명 용도 불일치")
        if not group:
            reasons.append("비교군 검증 없음")
        for field, description in [("용도검증", "동일용도"), ("공정호환성검증", "공정호환성"), ("중복개정검토", "중복/개정")]:
            if group_review.get(field, "").strip() not in TRUE_VALUES:
                reasons.append(f"{description} 미승인")
        if not group_review.get("근거문서", "").strip():
            reasons.append("근거문서 없음")
        missing_oel = sorted(product_details_cas[pid] - set(oel))
        if keys and missing_oel:
            reasons.append("공식검증 TWA 누락: " + ",".join(missing_oel))
        if keys:
            rp = [num(by_sc[k][0]["총압력_Raoult_kPa"], "Raoult", positive=True) for k in keys]
            hp = [by_sc[k][1] for k in keys]
            physical.append({"제품ID": pid, "제품명": name, "원분류": original_group, "시나리오수": len(keys),
                             "Raoult_중앙값_kPa": statistics.median(rp), "HANNA_중앙값_kPa": statistics.median(hp),
                             "HANNA_5백분위_kPa": qtile(hp, 0.05), "HANNA_95백분위_kPa": qtile(hp, 0.95),
                             "압력비_중앙값": statistics.median([h/r for h, r in zip(hp, rp)]),
                             "비교해석": "압력만 기술·저유해성 순위 아님"})
        approved = not reasons
        audit.append({"제품ID": pid, "제품명": name, "원분류": original_group,
                      "분류신뢰도": p.get("분류 신뢰도", ""), "분류근거": p.get("분류 근거", ""),
                      "모델결과존재": "예" if keys else "아니오", "모델시나리오수": len(keys),
                      "모델검증정보일치": "예" if valid_meta else "아니오",
                      "비교군ID": group, "동일용도검증": group_review.get("용도검증", "미검증"),
                      "공정호환성검증": group_review.get("공정호환성검증", "미검증"),
                      "중복개정검토": group_review.get("중복개정검토", "미검증"),
                      "국내TWA전성분검증": "예" if keys and not missing_oel else "아니오",
                      "평가상태": "검증완료_비교군수확인전" if approved else "보류",
                      "차단사유": "; ".join(reasons)})
        if not approved:
            continue
        values = []
        for key in keys:
            parts = by_detail[key]
            normalized = [(1e6 * num(row["HANNA_분압_kPa"], "분압") / (ATM_KPA * oel[row["CAS"].strip()]), row["CAS"].strip()) for row in parts]
            val, cas_max = max(normalized)
            values.append(val)
            scenario_scores.append({"제품ID": pid, "비교군ID": group, "시나리오ID": key[1],
                                    "최대TWA정규화평형분압": val, "최대기여성분_CAS": cas_max})
        scores.append({"비교군ID": group, "제품ID": pid, "제품명": name, "시나리오수": len(values),
                       "최대TWA정규화평형분압_중앙값": statistics.median(values),
                       "최대TWA정규화평형분압_최솟값": min(values),
                       "최대TWA정규화평형분압_최댓값": max(values),
                       "최대TWA정규화평형분압_5백분위": qtile(values, .05),
                       "최대TWA정규화평형분압_95백분위": qtile(values, .95),
                       "참고순위": "", "해석상태": "비교군 수 부족"})

    by_group = defaultdict(list)
    for row in scores:
        by_group[row["비교군ID"]].append(row)
    pairs = []
    for group, rows in by_group.items():
        if len(rows) < 2:
            continue
        rows.sort(key=lambda x: (x["최대TWA정규화평형분압_중앙값"], x["제품ID"]))
        for rank, row in enumerate(rows, 1):
            row["참고순위"] = rank
            row["해석상태"] = "상대 잠재지표 참고순위·저위험 인증 아님"
        for i, a in enumerate(rows):
            for b in rows[i + 1:]:
                sep = a["최대TWA정규화평형분압_최댓값"] < b["최대TWA정규화평형분압_최솟값"]
                pairs.append({"비교군ID": group, "지표_낮은제품ID": a["제품ID"],
                              "지표_높은제품ID": b["제품ID"],
                              "중앙값_높음나눔낮음": b["최대TWA정규화평형분압_중앙값"] / a["최대TWA정규화평형분압_중앙값"],
                              "시나리오범위_분리": "예" if sep else "아니오",
                              "해석": "샘플링 범위에서만 분리·실측 검증 필요" if sep else "샘플링 범위 중첩·서열 불확실"})
    for a in audit:
        if a["평가상태"] == "검증완료_비교군수확인전":
            if len(by_group.get(a["비교군ID"], [])) >= 2:
                a["평가상태"] = "조건부_비교가능"
            else:
                a["평가상태"] = "보류"
                a["차단사유"] = "같은 검증 비교군의 제품 2개 이상 필요"
    out = args.output_dir
    save_csv(out / "271_제품별_분석가능성_감사.csv", audit, AUDIT_COLUMNS)
    save_csv(out / "제품군별_증기상_기술비교.csv", physical, PHYSICAL_COLUMNS)
    save_csv(out / "동일용도_저유해성_후보선별.csv", scores, HAZARD_COLUMNS)
    save_csv(out / "동일용도_제품쌍_잠재지표비교.csv", pairs, PAIR_COLUMNS)
    save_csv(out / "동일용도_시나리오_잠재지표.csv", scenario_scores, SCENARIO_COLUMNS)
    print(f"문서 {len(products)}건 감사; HANNA 계산 {len(physical)}건; 조건부 평가 {len(scores)}건; "
          f"동일 비교군 내 제품쌍 {len(pairs)}쌍")
    print("검증입력 없으면 비교순위 생성하지 않음. 실제 작업장 노출·법규위반·대체제품 승인 아님.")


if __name__ == "__main__":
    main()