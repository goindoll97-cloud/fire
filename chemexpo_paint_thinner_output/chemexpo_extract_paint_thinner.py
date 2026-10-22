# -*- coding: utf-8 -*-
"""EPA ChemExpo 대용량 Composition Data에서 도료 희석제 성분자료만 추출.

외부 라이브러리 불필요 (Python 3.9+). CSV/TSV/TXT/ZIP/GZ 입력 지원.

사용 (Windows PowerShell / CMD):
    python chemexpo_extract_paint_thinner.py "ChemExpo_bulk_composition_chemicals.csv"

또는 아래 INPUT_FILE 경로를 수정한 뒤:
    python chemexpo_extract_paint_thinner.py

출력: 원본 파일과 같은 폴더의 chemexpo_paint_thinner_output/ 하위 CSV 5개.

중요: ChemExpo Reported는 원천 문서에 보고된 값입니다. 반드시 실험 실측값은 아닙니다.
      공개 Composition 자료로 SDS 미기재 미량성분의 부재/존재를 증명할 수 없습니다.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import re
import sys
import zipfile
from collections import Counter, defaultdict
from contextlib import contextmanager
from pathlib import Path

# >>> 파일 경로만 바꾸면 실행 가능. 또는 명령줄에 파일 경로를 입력하세요.
INPUT_FILE = r"ChemExpo_bulk_composition_chemicals.csv"

# PUC 검색은 paint thinner, lacquer thinner 등 관련 제품 유형/제품군을 대상으로 함.
PUC_REGEX = re.compile(
    r"\b(?:paint|lacquer|enamel|coating|varnish|urethane|epoxy|automotive)\b"
    r".{0,75}\b(?:thinners?|reducers?)\b"
    r"|\b(?:thinners?|reducers?)\b.{0,75}"
    r"\b(?:paint|lacquer|enamel|coating|varnish|urethane|epoxy|automotive)\b",
    re.IGNORECASE,
)
# 제품명만 보고 찾은 자료는 도료 희석제 여부가 확실하지 않으므로 별도 보관.
NAME_REGEX = re.compile(r"\b(?:thinners?|reducers?)\b", re.IGNORECASE)

RAW_COLUMNS = [
    "Data Source", "Data Document ID", "Data Document Title", "Data Document Subtitle",
    "Document Date", "Product Name", "PUC Kind", "PUC General Category",
    "PUC Product Family", "PUC Product Type", "PUC Classification Method",
    "Raw Chemical Name", "Raw CAS", "DTXSID", "Curated Chemical Name",
    "Curated CAS", "Provisional DTXSID", "Raw Min Comp", "Raw Max Comp",
    "Raw Central Comp", "Unit Type", "Lower Weight Fraction",
    "Upper Weight Fraction", "Central Weight Fraction", "Weight Fraction Type",
]
REQUIRED_COLUMNS = [
    "Product Name", "PUC Product Type", "Curated Chemical Name",
    "Curated CAS", "Lower Weight Fraction", "Upper Weight Fraction",
    "Central Weight Fraction", "Weight Fraction Type",
]
ADDED_COLUMNS = [
    "extract_group", "extract_reason", "product_key", "component_CAS_used",
    "CAS_source", "weight_fraction_method", "wf_for_screening", "wt_percent_for_screening",
    "reported_flag", "bounds_valid", "row_warning",
]


def normalize(s: object) -> str:
    return re.sub(r"\s+", " ", str(s or "").replace("\ufeff", "").strip()).casefold()


def safe_float(s: object):
    value = str(s or "").strip().replace(",", "")
    if not value or value.casefold() in ("nan", "null", "none", "na"):
        return None
    try:
        out = float(value)
        return out if 0 <= out <= 1 else None
    except (ValueError, OverflowError):
        return None


def fmt(n):
    return "" if n is None else f"{n:.10g}"


def get_source_encoding(prefix: bytes):
    if prefix.startswith((b"\xff\xfe", b"\xfe\xff")):
        return "utf-16"
    if prefix.startswith(b"\xef\xbb\xbf"):
        return "utf-8-sig"
    return "utf-8-sig"


@contextmanager
def open_bulk_file(path: Path):
    """1GB 파일/ZIP 압축을 풀지 않고 스트리밍으로 읽기."""
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as z:
            items = [i for i in z.infolist()
                     if not i.is_dir() and Path(i.filename).suffix.lower() in (".csv", ".tsv", ".txt")]
            if not items:
                raise ValueError("ZIP 안에 CSV/TSV/TXT 파일을 찾을 수 없습니다.")
            member = max(items, key=lambda x: x.file_size)
            print(f"[ZIP 내부파일] {member.filename}", flush=True)
            with z.open(member) as raw:
                with io.BufferedReader(raw) as buf:
                    encoding = get_source_encoding(buf.peek(4)[:4])
                    with io.TextIOWrapper(buf, encoding=encoding, errors="replace", newline="") as f:
                        yield f
    elif path.suffix.lower() == ".gz":
        with gzip.open(path, "rb") as raw:
            with io.BufferedReader(raw) as buf:
                encoding = get_source_encoding(buf.peek(4)[:4])
                with io.TextIOWrapper(buf, encoding=encoding, errors="replace", newline="") as f:
                    yield f
    else:
        with open(path, "rb") as raw:
            with io.BufferedReader(raw) as buf:
                encoding = get_source_encoding(buf.peek(4)[:4])
                with io.TextIOWrapper(buf, encoding=encoding, errors="replace", newline="") as f:
                    yield f


def choose_delimiter(header_line: str):
    return max([",", "\t", ";"], key=header_line.count)


def identify_product(row):
    family = row.get("PUC Product Family", "")
    ptype = row.get("PUC Product Type", "")
    family_type = f"{family} {ptype}"
    if PUC_REGEX.search(family_type):
        return "PUC_matched", "PUC 제품군/유형에 paint/lacquer 등과 thinner/reducer 일치"
    # 단독으로 더 구체적 분류명일 수도 있음.
    if normalize(ptype) in ("paint thinner", "paint thinners", "lacquer thinner", "lacquer thinners"):
        return "PUC_matched", "PUC 제품유형 정확히 일치"
    name = row.get("Product Name", "")
    if NAME_REGEX.search(name):
        return "Name_candidate", "제품명에 thinner/reducer 포함; 도료용 제품인지 수동 검토"
    return "", ""


def add_derived(row, group, reason, product_key):
    """누락값을 임의로 채우지 않음. 원본 정제 질량분율만 계산에 사용."""
    out = dict(row)
    lo = safe_float(row.get("Lower Weight Fraction"))
    hi = safe_float(row.get("Upper Weight Fraction"))
    ce = safe_float(row.get("Central Weight Fraction"))
    wtype = normalize(row.get("Weight Fraction Type"))
    reported = wtype.startswith("reported")
    cas_c = str(row.get("Curated CAS", "") or "").strip()
    cas_r = str(row.get("Raw CAS", "") or "").strip()
    cas = cas_c or cas_r
    cas_source = "curated" if cas_c else ("raw_unverified" if cas_r else "missing")
    warnings = []
    bounds_valid = not (lo is not None and hi is not None and lo > hi)
    if not bounds_valid:
        warnings.append("lower_greater_than_upper")
    if ce is not None and lo is not None and ce < lo - 1e-12:
        warnings.append("central_below_lower")
    if ce is not None and hi is not None and ce > hi + 1e-12:
        warnings.append("central_above_upper")
    if not cas:
        warnings.append("CAS_missing")
    if not reported:
        warnings.append("not_reported_data")
    wf = None
    method = "missing_or_range_only"
    if bounds_valid:
        if ce is not None:
            wf = ce
            method = "central_weight_fraction"
        elif lo is not None and hi is not None and abs(lo - hi) <= 1e-12:
            wf = lo
            method = "equal_lower_upper"
    if not bounds_valid or any(x in warnings for x in ["central_below_lower", "central_above_upper"]):
        wf = None
        method = "invalid_bounds"
    out.update({
        "extract_group": group,
        "extract_reason": reason,
        "product_key": product_key,
        "component_CAS_used": cas,
        "CAS_source": cas_source,
        "weight_fraction_method": method,
        "wf_for_screening": fmt(wf),
        "wt_percent_for_screening": fmt(wf * 100) if wf is not None else "",
        "reported_flag": "1" if reported else "0",
        "bounds_valid": "1" if bounds_valid else "0",
        "row_warning": ";".join(warnings),
    })
    return out


def product_id(row):
    # 서로 다른 SDS 문서의 같은 제품명을 합쳐 함량을 더하는 실수를 방지.
    return " || ".join([
        str(row.get("Data Source", "") or "").strip(),
        str(row.get("Data Document ID", "") or "").strip(),
        str(row.get("Data Document Title", "") or "").strip(),
        str(row.get("Product Name", "") or "").strip(),
    ])


def summarize_product(records):
    first = records[0]
    n = len(records)
    valid = [safe_float(r["wf_for_screening"]) for r in records]
    known = [x for x in valid if x is not None]
    have_both = []
    for r in records:
        lo = safe_float(r.get("Lower Weight Fraction"))
        hi = safe_float(r.get("Upper Weight Fraction"))
        have_both.append((lo, hi) if lo is not None and hi is not None and lo <= hi else None)
    wt_types = Counter(normalize(r.get("Weight Fraction Type")) or "unknown" for r in records)
    cas_values = [str(r["component_CAS_used"]).strip() for r in records if r["component_CAS_used"]]
    n_dup_cas = sum(cnt - 1 for cnt in Counter(cas_values).values() if cnt > 1)
    fully_known_reported = all(r["reported_flag"] == "1" and valid[i] is not None
                               and r["CAS_source"] == "curated" and not r["row_warning"].replace("CAS_missing", "").strip(";")
                               for i, r in enumerate(records))
    known_sum = sum(known) if known else None
    plausible_sum = known_sum is not None and 0.95 <= known_sum <= 1.05
    candidate = (n >= 2 and fully_known_reported and n_dup_cas == 0 and plausible_sum)
    if candidate:
        status = "candidate_manual_review_required"
    elif len(known) < 2:
        status = "not_ready_less_than_2_quantified_chemicals"
    elif len(known) < n:
        status = "not_ready_missing_exact_fractions"
    elif any(r["reported_flag"] != "1" for r in records):
        status = "not_ready_predicted_or_unknown_fractions"
    elif n_dup_cas:
        status = "not_ready_duplicate_chemical_rows"
    elif not plausible_sum:
        status = "not_ready_fraction_total_not_near_100pct"
    else:
        status = "not_ready_incomplete_identifier_or_bounds"
    result = {
        "product_key": first["product_key"],
        "extract_group": first["extract_group"],
        "Data Source": first.get("Data Source", ""),
        "Data Document ID": first.get("Data Document ID", ""),
        "Data Document Title": first.get("Data Document Title", ""),
        "Document Date": first.get("Document Date", ""),
        "Product Name": first.get("Product Name", ""),
        "PUC Product Family": first.get("PUC Product Family", ""),
        "PUC Product Type": first.get("PUC Product Type", ""),
        "n_chemical_rows": n,
        "n_distinct_CAS": len(set(cas_values)),
        "n_duplicate_CAS_rows": n_dup_cas,
        "n_rows_with_numeric_point_fraction": len(known),
        "n_rows_with_range_bounds": sum(x is not None for x in have_both),
        "n_rows_reported": sum(r["reported_flag"] == "1" for r in records),
        "n_rows_predicted": sum("predicted" in normalize(r.get("Weight Fraction Type")) for r in records),
        "n_rows_unknown_type": wt_types.get("unknown", 0),
        "sum_of_available_point_fraction": fmt(known_sum),
        "sum_of_available_point_wt_pct": fmt(known_sum * 100) if known_sum is not None else "",
        "sum_of_lower_bounds_if_all_present": fmt(sum(x[0] for x in have_both)) if all(x is not None for x in have_both) else "",
        "sum_of_upper_bounds_if_all_present": fmt(sum(x[1] for x in have_both)) if all(x is not None for x in have_both) else "",
        "model_input_candidate": "YES_manual_review" if candidate else "NO",
        "quality_status": status,
        "note": "EPA reported does not imply laboratory-measured; check source SDS, missing trace chemicals, and document duplicates",
    }
    return result, candidate


def write_csv_file(path, headers, rows):
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=headers, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def run(input_path: Path, output_dir: Path):
    if not input_path.is_file():
        raise FileNotFoundError(f"원본 파일이 없습니다: {input_path}")
    output_dir.mkdir(parents=True, exist_ok=True)
    try:
        csv.field_size_limit(sys.maxsize)
    except OverflowError:
        csv.field_size_limit(2**31 - 1)
    puc_records = []
    name_records = []
    read_n = 0
    counts = Counter()
    with open_bulk_file(input_path) as f:
        first_header_line = f.readline()
        if not first_header_line:
            raise ValueError("입력 파일이 비어 있습니다.")
        delimiter = choose_delimiter(first_header_line)
        # 한 줄만 읽어보고 파일 처음으로 복귀. 압축 파일도 seek(0)을 지원함.
        try:
            f.seek(0)
        except (OSError, io.UnsupportedOperation):
            # seek가 불가능한 스트림의 경우 첫 줄을 다시 포함하는 iterator 사용
            import itertools
            reader = csv.DictReader(itertools.chain([first_header_line], f), delimiter=delimiter)
        else:
            reader = csv.DictReader(f, delimiter=delimiter)
        fieldnames = [(x or "").strip().lstrip("\ufeff") for x in (reader.fieldnames or [])]
        reader.fieldnames = fieldnames
        missing = [c for c in REQUIRED_COLUMNS if c not in fieldnames]
        if missing:
            raise ValueError(
                f"필수 열 누락: {missing}\n인식한 열: {fieldnames}\n"
                "파일이 PUC 분류표인지 Composition 자료인지 확인하세요."
            )
        print(f"[읽기 시작] 구분자={repr(delimiter)} / {len(fieldnames)}개 열 / 파일={input_path}", flush=True)
        for row in reader:
            read_n += 1
            if None in row:
                counts["malformed_extra_columns"] += 1
            group, reason = identify_product(row)
            if group:
                clean = add_derived(row, group, reason, product_id(row))
                if group == "PUC_matched":
                    puc_records.append(clean)
                else:
                    name_records.append(clean)
                counts[group] += 1
            if read_n % 100000 == 0:
                print(f"읽은 행 {read_n:,} / PUC 일치 {counts['PUC_matched']:,} / 제품명 후보 {counts['Name_candidate']:,}", flush=True)

    cols = list(dict.fromkeys(fieldnames + ADDED_COLUMNS))
    write_csv_file(output_dir / "01_paint_thinner_PUC_matched.csv", cols, puc_records)
    write_csv_file(output_dir / "02_thinner_name_only_review.csv", cols, name_records)

    products = defaultdict(list)
    for rec in puc_records:
        products[rec["product_key"]].append(rec)
    summaries = []
    model_ready = []
    for pk, recs in products.items():
        s, candidate = summarize_product(recs)
        summaries.append(s)
        if candidate:
            for rec in recs:
                # 검토용 후보만 추출: 확정된 완전 조성 실측자료가 아님.
                model_ready.append({
                    "product_key": rec["product_key"],
                    "Product Name": rec["Product Name"],
                    "Data Document ID": rec.get("Data Document ID", ""),
                    "Data Source": rec.get("Data Source", ""),
                    "chemical_name": rec.get("Curated Chemical Name", "") or rec.get("Raw Chemical Name", ""),
                    "CAS": rec["component_CAS_used"],
                    "DTXSID": rec.get("DTXSID", ""),
                    "weight_fraction": rec["wf_for_screening"],
                    "wt_percent": rec["wt_percent_for_screening"],
                    "Weight Fraction Type": rec.get("Weight Fraction Type", ""),
                    "data_caveat": "composition candidate only; independent source verification required",
                })
    summaries.sort(key=lambda x: (x["model_input_candidate"] != "YES_manual_review", x["Product Name"]))
    summary_cols = list(summaries[0]) if summaries else [
        "product_key", "extract_group", "Data Source", "Data Document ID", "Product Name",
        "n_chemical_rows", "model_input_candidate", "quality_status"
    ]
    write_csv_file(output_dir / "03_product_quality_summary.csv", summary_cols, summaries)
    model_headers = [
        "product_key", "Product Name", "Data Document ID", "Data Source", "chemical_name",
        "CAS", "DTXSID", "weight_fraction", "wt_percent", "Weight Fraction Type", "data_caveat"
    ]
    write_csv_file(output_dir / "04_model_input_candidates_REVIEW.csv", model_headers, model_ready)
    audit = [
        ("input_file", str(input_path)),
        ("total_rows_read", read_n),
        ("PUC_matched_rows", len(puc_records)),
        ("name_only_candidate_rows", len(name_records)),
        ("PUC_matched_products_by_document", len(products)),
        ("candidate_products_manual_review", sum(s["model_input_candidate"] == "YES_manual_review" for s in summaries)),
        ("malformed_rows_with_extra_columns", counts["malformed_extra_columns"]),
        ("caution_1", "Reported means in-source-reported, not necessarily laboratory measured."),
        ("caution_2", "Composition data may omit trace chemicals; cannot test legal omission directly."),
        ("caution_3", "Fraction totals near 1 do not establish complete true composition."),
        ("caution_4", "PUC matched is a category screen, not confirmation of product identity."),
    ]
    write_csv_file(output_dir / "05_readme_and_audit.csv", ["metric", "value"],
                   [{"metric": k, "value": v} for k, v in audit])
    print("\n[완료]", flush=True)
    for k, v in audit[:6]:
        print(f"  {k}: {v}")
    print(f"\n생성된 출력 폴더: {output_dir.resolve()}")
    for p in sorted(output_dir.glob("*.csv")):
        print(f"  - {p.name} ({p.stat().st_size:,} bytes)")


def main():
    parser = argparse.ArgumentParser(description="EPA ChemExpo 대용량 데이터에서 도료 희석제 성분 추출")
    parser.add_argument("input", nargs="?", default=INPUT_FILE, help="원본 .csv/.tsv/.txt/.zip/.gz 경로")
    parser.add_argument("--out", default=None, help="출력 폴더 (생략 시 원본 폴더 하위) ")
    args = parser.parse_args()
    path = Path(args.input).expanduser()
    output = Path(args.out).expanduser() if args.out else path.parent / "chemexpo_paint_thinner_output"
    run(path, output)


if __name__ == "__main__":
    main()
