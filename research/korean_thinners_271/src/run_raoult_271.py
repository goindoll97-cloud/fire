# -*- coding: utf-8 -*-
"""25°C Raoult equilibrium *potential* from feasible MSDS composition ranges.

Not measured evaporation rate, occupational exposure or an AI/HANNA result.
Conditionally sampled feasible scenarios are NOT confidence intervals.
"""
import argparse
import collections
import csv
import random
import statistics
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
DATA = BASE / "data"
PHYS = Path(__file__).resolve().parents[3] / "thinner_vle_reference_25C.csv"

def read(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def output_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

def quantile(values, q):
    a = sorted(values)
    k = (len(a) - 1) * q
    low = int(k)
    high = min(len(a) - 1, low + 1)
    return a[low] + (a[high] - a[low]) * (k - low)

def feasible_compositions(comps, count, rng, attempts_multiplier=1500):
    """Produce valid scenarios, no unsupported assumption about true SDS weights."""
    lows = [float(r["함량 하한(%)"]) for r in comps]
    highs = [float(r["함량 상한(%)"]) for r in comps]
    if any(lo < 0 or hi > 100 or lo > hi for lo, hi in zip(lows, highs)):
        raise ValueError("Invalid MSDS composition bound")
    if sum(lows) > 100.000001 or sum(highs) < 99.999999:
        raise ValueError("Composition intervals cannot sum to 100%")
    if len(comps) < 2:
        raise ValueError("At least two identified components required")
    feasible = []
    attempts = 0
    while len(feasible) < count and attempts < count * attempts_multiplier:
        attempts += 1
        permutation = list(range(len(comps)))
        rng.shuffle(permutation)
        fractions = [0.0] * len(comps)
        for j in permutation[:-1]:
            fractions[j] = rng.uniform(lows[j], highs[j])
        last = permutation[-1]
        fractions[last] = 100.0 - sum(fractions)
        if lows[last] - 1e-8 <= fractions[last] <= highs[last] + 1e-8:
            feasible.append(fractions)
    if len(feasible) != count:
        raise RuntimeError("Insufficient feasible samples: acceptance too low")
    return feasible

def main():
    parser = argparse.ArgumentParser(description="Raoult equilibrium scenario audit at 25 degrees C")
    parser.add_argument("--samples", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261023)
    parser.add_argument("--output", type=Path, default=BASE / "output")
    parser.add_argument("--all-eligible", action="store_true",
                        help="Use all 46 molecular candidates, skipping missing-property rows (default only 6)")
    args = parser.parse_args()
    if not 1 <= args.samples <= 100000:
        parser.error("--samples must be between 1 and 100000")

    classifications = read(DATA / "제품271건_용도재분류.csv")
    all_components = read(DATA / "성분별함량_1462.csv")
    props = {r["CAS"]: r for r in read(PHYS)}
    by_id = collections.defaultdict(list)
    for r in all_components:
        by_id[r["제품ID"]].append(r)

    eligible = [r for r in classifications
                if r["분자모델 1차 후보"] == "예"]
    if not args.all_eligible:
        eligible = [r for r in eligible if r["기존물성 라울예비 대상"] == "예"]
    rng = random.Random(args.seed)
    scenario_rows = []
    component_rows = []
    summary_rows = []
    excluded_rows = []
    for r in eligible:
        pid = r["제품ID"]
        components = by_id[pid]
        cas_list = [c["CAS 번호"] for c in components]
        missing = sorted(set(cas_list) - set(props))
        if missing:
            excluded_rows.append({"제품ID": pid, "제외사유": "25C VP/MW 미확보: " + ", ".join(missing)})
            continue
        try:
            samples = feasible_compositions(components, args.samples, rng)
        except (ValueError, RuntimeError) as e:
            excluded_rows.append({"제품ID": pid, "제외사유": str(e)})
            continue

        pressures = []
        for scenario_id, weights in enumerate(samples, 1):
            mw = [float(props[c]["MW_g_mol"]) for c in cas_list]
            vp = [float(props[c]["VP_25C_kPa"]) for c in cas_list]
            if not all(x > 0 for x in mw + vp):
                raise ValueError("Non-positive property in " + pid)
            amounts = [w / m for w, m in zip(weights, mw)]
            mol_sum = sum(amounts)
            mole_fractions = [a / mol_sum for a in amounts]
            partial = [x * p for x, p in zip(mole_fractions, vp)]
            total = sum(partial)
            pressures.append(total)
            liquid_top = components[max(range(len(weights)), key=weights.__getitem__)]["CAS 번호"]
            vapor_top = components[max(range(len(partial)), key=partial.__getitem__)]["CAS 번호"]
            scenario_rows.append({
                "제품ID": pid, "제품군": r["1차 용도분류"], "시나리오ID": scenario_id,
                "총평형증기압_kPa": round(total, 7),
                "액상최다성분_CAS": liquid_top, "증기상최다성분_CAS": vapor_top,
                "두우세성분_상이": int(liquid_top != vapor_top),
                "모델": "Raoult 25C (gamma=1)",
            })
            for c, weight, xi, pi in zip(components, weights, mole_fractions, partial):
                component_rows.append({
                    "제품ID": pid, "시나리오ID": scenario_id, "CAS": c["CAS 번호"],
                    "화학물질명": c["화학물질명(한글)"], "MSDS질량백분율": round(weight, 7),
                    "액상몰분율": round(xi, 9), "평형분압_kPa": round(pi, 7),
                    "평형기상몰분율": round(pi / total, 9), "모델": "Raoult 25C",
                })
        summary_rows.append({
            "제품ID": pid, "제품명": r["제품명(한글)"], "제품군": r["1차 용도분류"],
            "시나리오수": len(pressures), "총증기압_5백분위_kPa": round(quantile(pressures, .05), 5),
            "총증기압_중앙값_kPa": round(statistics.median(pressures), 5),
            "총증기압_95백분위_kPa": round(quantile(pressures, .95), 5),
            "모델": "Raoult 25C (no gamma or rate)",
        })
    out = args.output
    output_csv(out / "raoult_제품별_요약.csv", summary_rows, [
        "제품ID", "제품명", "제품군", "시나리오수", "총증기압_5백분위_kPa",
        "총증기압_중앙값_kPa", "총증기압_95백분위_kPa", "모델"])
    output_csv(out / "raoult_시나리오.csv", scenario_rows, [
        "제품ID", "제품군", "시나리오ID", "총평형증기압_kPa", "액상최다성분_CAS",
        "증기상최다성분_CAS", "두우세성분_상이", "모델"])
    output_csv(out / "raoult_성분별.csv", component_rows, [
        "제품ID", "시나리오ID", "CAS", "화학물질명", "MSDS질량백분율",
        "액상몰분율", "평형분압_kPa", "평형기상몰분율", "모델"])
    output_csv(out / "raoult_계산제외.csv", excluded_rows, ["제품ID", "제외사유"])
    print(f"완료: Raoult 제품 {len(summary_rows)}건, 시나리오 {len(scenario_rows)}개, 제외 {len(excluded_rows)}건")
    print("주의: 통계적 신뢰구간·실제 시간당 증발률·HANNA 예측·실측노출 결과 아님.")

if __name__ == "__main__":
    main()
