# -*- coding: utf-8 -*-
"""Guarded blend scenario search. Real HANNA coefficients require a verified external model.

No suitability / substitution claim without experimentally backed process constraints.
"""
import argparse
import csv
import importlib
import json
import math
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = ("성분구조검증", "동일공정호환성검증", "물성출처검증", "한국TWA최신성검증", "기능제약실측검증")
ATM_KPA = 101.325

def template():
    return {
        "대상제품군": "예: 특정 도료 및 공정에서 호환성 검증된 그룹",
        "근거문서": "TDS/실측/법령 출처",
        "검증": {item: False for item in REQUIRED},
        "온도_K": 298.15,
        "목표_총평형증기압범위_kPa": [0.0, 101.325],
        "성분": [
            {
                "CAS": "REPLACE", "SMILES": "REPLACE", "분자량_g_mol": None,
                "포화증기압_kPa": None, "한국_TWA_ppm": None,
                "최소질량비_pct": None, "최대질량비_pct": None,
            }
        ],
    }

def load_config(path):
    d = json.loads(path.read_text(encoding="utf-8-sig"))
    flags = d.get("검증", {})
    failed = [key for key in REQUIRED if flags.get(key) is not True]
    if failed:
        raise ValueError("최적화 실행 보류: 검증되지 않은 조건: " + ", ".join(failed))
    comps = d.get("성분") or []
    if len(comps) < 2:
        raise ValueError("2종 이상 성분이 필요함")
    cas = [x.get("CAS", "") for x in comps]
    if not all(cas) or len(set(cas)) != len(cas):
        raise ValueError("CAS 누락/중복")
    for x in comps:
        for k in ("분자량_g_mol", "포화증기압_kPa", "한국_TWA_ppm"):
            if not isinstance(x.get(k), (float, int)) or x[k] <= 0:
                raise ValueError("검증된 양의 물성/TWA 누락: " + str((x["CAS"], k)))
        for k in ("최소질량비_pct", "최대질량비_pct"):
            if not isinstance(x.get(k), (float, int)):
                raise ValueError("조성 제약 누락: " + str((x["CAS"], k)))
        if not 0 <= x["최소질량비_pct"] <= x["최대질량비_pct"] <= 100:
            raise ValueError("조성 범위 오류: " + x["CAS"])
        if not x.get("SMILES") or x["SMILES"] == "REPLACE":
            raise ValueError("분자구조(SMILES) 미확인: " + x["CAS"])
    if sum(x["최소질량비_pct"] for x in comps) > 100 or sum(x["최대질량비_pct"] for x in comps) < 100:
        raise ValueError("조성합계 100% 불가능")
    band = d["목표_총평형증기압범위_kPa"]
    if len(band) != 2 or not 0 <= band[0] < band[1]:
        raise ValueError("기능 검증된 총증기압 허용범위 필요")
    if abs(float(d["온도_K"]) - 298.15) > 0.001:
        raise ValueError("현재 스크립트는 물성표가 25°C인 경우만 지원")
    return d, comps, band

def draw(comps, rng):
    ids = list(range(len(comps)))
    rng.shuffle(ids)
    masses = [0.0] * len(comps)
    for i in ids[:-1]:
        c = comps[i]
        masses[i] = rng.uniform(c["최소질량비_pct"], c["최대질량비_pct"])
    masses[ids[-1]] = 100 - sum(masses)
    j = ids[-1]
    if comps[j]["최소질량비_pct"] <= masses[j] <= comps[j]["최대질량비_pct"]:
        return masses
    return None

def gamma_predictor(spec):
    if not spec:
        return lambda smiles, x, temp: [1.0] * len(x), "Raoult (gamma=1, AI 아님)"
    if ":" not in spec:
        raise ValueError("--gamma-predictor는 module:function 형식")
    module, func = spec.split(":", 1)
    predictor = getattr(importlib.import_module(module), func)
    return predictor, "외부 실검증된 gamma 예측기: " + spec

def pareto(rows):
    """Return nondominated points minimizing both objective columns."""
    result = []
    for row in sorted(rows, key=lambda z: (z["총평형증기압_kPa"], z["최대TWA정규화분압"])):
        if not result or row["최대TWA정규화분압"] < result[-1]["최대TWA정규화분압"] - 1e-12:
            result.append(row)
    return result

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-template", type=Path)
    ap.add_argument("--config", type=Path)
    ap.add_argument("--gamma-predictor", help="External VERIFIED 2026 HANNA adapter module:function")
    ap.add_argument("--n", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20261023)
    ap.add_argument("--out", type=Path, default=ROOT / "output" / "배합탐색_Pareto.csv")
    args = ap.parse_args()
    if args.write_template:
        args.write_template.parent.mkdir(parents=True, exist_ok=True)
        args.write_template.write_text(json.dumps(template(), ensure_ascii=False, indent=2), encoding="utf-8")
        print("미검증 초기 템플릿 생성:", args.write_template)
        return
    if not args.config:
        ap.error("--config 또는 --write-template 지정")
    if not 1 <= args.n <= 100000:
        ap.error("--n must be 1..100000")
    d, comps, band = load_config(args.config)
    fn, model = gamma_predictor(args.gamma_predictor)
    rng = random.Random(args.seed)
    results = []
    tries = 0
    while len(results) < args.n and tries < args.n * 2000:
        tries += 1
        masses = draw(comps, rng)
        if masses is None:
            continue
        mol = [w / c["분자량_g_mol"] for c, w in zip(comps, masses)]
        mol_sum = sum(mol)
        x = [a / mol_sum for a in mol]
        gammas = list(fn([c["SMILES"] for c in comps], x, 298.15))
        if len(gammas) != len(x) or not all(math.isfinite(g) and g > 0 for g in gammas):
            raise ValueError("HANNA predictor invalid gamma values: refuse fallback")
        partial = [xi * g * c["포화증기압_kPa"] for xi, g, c in zip(x, gammas, comps)]
        total = sum(partial)
        if not band[0] <= total <= band[1]:
            continue
        max_index = max(1e6 * pi / (ATM_KPA * c["한국_TWA_ppm"]) for pi, c in zip(partial, comps))
        results.append({
            "총평형증기압_kPa": round(total, 8),
            "최대TWA정규화분압": round(max_index, 8),
            "질량비_pct_JSON": json.dumps({c["CAS"]: round(w, 7) for c, w in zip(comps, masses)},
                                        ensure_ascii=False),
            "모델": model,
            "주의": "검증된 가상배합 탐색; 증발속도·작업장농도·실제품 권고 아님",
        })
    if not results:
        raise RuntimeError("No feasible scenario satisfying model/functional band; no optimum claimed")
    answers = pareto(results)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(answers[0]))
        writer.writeheader()
        writer.writerows(answers)
    print(f"모델: {model}; feasible {len(results)}, Pareto {len(answers)}; 저장 {args.out}")
    print("TWA 정규화 평형분압은 실제 노출초과를 나타내지 않습니다.")

if __name__ == "__main__":
    main()
