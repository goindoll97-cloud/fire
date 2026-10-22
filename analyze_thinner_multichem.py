# -*- coding: utf-8 -*-
"""Multi-chemical counterfactual screening for ChemExpo paint thinners.

Reuses the real Raoult/UNIFAC/HANNA gamma implementations from
compare_thinner_raoult_unifac_hanna.py. Vapor pressure is fixed at 25 C.
This is NOT an exposure assessment or a finding about legal SDS omission.

Default input (relative to this file):
  chemexpo_paint_thinner_output/04_model_input_candidates_REVIEW.csv
Default output:
  chemexpo_paint_thinner_output/multichem_screening_output/
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

import compare_thinner_raoult_unifac_hanna as core

ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "chemexpo_paint_thinner_output" / "04_model_input_candidates_REVIEW.csv"
DEFAULT_REGISTRY = ROOT / "multichem_screening_registry.csv"
DEFAULT_OUT = ROOT / "chemexpo_paint_thinner_output" / "multichem_screening_output"
DEFAULT_CANDIDATES = ("71-43-2", "67-56-1", "67-64-1")  # benzene, methanol, acetone
DEFAULT_DOSES = (0.039, 0.099, 0.100, 0.990, 1.000)


def arguments(argv=None):
    p = argparse.ArgumentParser(description="Multi-chemical hypothetical vapor screening at 25C")
    p.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    p.add_argument("--reference", type=Path, default=ROOT / "thinner_vle_reference_25C.csv")
    p.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    p.add_argument("--hanna-repo", type=Path, default=Path(r"D:\fire\HANNA"))
    p.add_argument("--product-ids", nargs="+", default=list(core.RESEARCH_PRODUCT_IDS))
    p.add_argument("--candidates", nargs="+", default=list(DEFAULT_CANDIDATES), help="Trace candidate CAS")
    p.add_argument("--doses", nargs="+", type=float, default=list(DEFAULT_DOSES), help="Hypothetical wt%%")
    p.add_argument("--models", nargs="+", choices=core.MODEL_IDS,
                   default=["raoult", "unifac", "hanna"])
    p.add_argument("--total-air-ppm", type=float, default=100.0,
                   help="Hypothetical reference only, never actual personal exposure")
    p.add_argument("--strict-models", action="store_true", help="Nonzero exit if a selected model fails")
    return p.parse_args(argv)


def load_registry(path: Path) -> dict:
    table = pd.read_csv(path, dtype={"CAS": str}, encoding="utf-8-sig", keep_default_na=False)
    mandatory = {"CAS", "chemical_label", "screening_cutoff_wt_pct",
                 "screening_cutoff_basis", "legal_review_status"}
    if not mandatory.issubset(table.columns):
        raise ValueError(f"Registry is missing columns: {sorted(mandatory-set(table.columns))}")
    if table.CAS.duplicated().any():
        raise ValueError("Duplicated CAS in screening registry")
    return {str(row.CAS).strip(): row.to_dict() for _, row in table.iterrows()}


def make_trace_mixture(base: list[tuple[str, float]], target: str, dose: float):
    """Reduce all baseline mass fractions proportionally and insert new CAS.

    Refuse to silently replace an existing baseline chemical or create duplicates.
    """
    if not math.isfinite(dose) or not 0 < dose < 100:
        raise ValueError(f"Trace dose must satisfy 0 < wt%% < 100, got {dose!r}")
    if any(cas == target for cas, _ in base):
        raise ValueError(f"ALREADY_PRESENT: {target} occurs in the real EPA input")
    if not np.isclose(sum(w for _, w in base), 100.0, atol=0.05):
        raise ValueError("Base mass fractions do not sum to 100 wt%")
    return [(cas, wt * (100-dose)/100) for cas, wt in base] + [(target, float(dose))]


def cutoff_relation(dose: float, entry: dict) -> str:
    raw = str(entry.get("screening_cutoff_wt_pct", "")).strip()
    if not raw:
        return "NO_REVIEWED_NUMERIC_CUTOFF"
    cutoff = float(raw)
    if not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("Invalid screening cutoff")
    # NOTE: deliberately no legal omission label. Other rules need a professional review.
    if math.isclose(dose, cutoff, abs_tol=1e-10):
        return "AT_NUMERIC_CUTOFF_NOT_LEGAL_DETERMINATION"
    return ("BELOW_NUMERIC_CUTOFF_NOT_LEGAL_DETERMINATION" if dose < cutoff
            else "ABOVE_NUMERIC_CUTOFF_NOT_LEGAL_DETERMINATION")


def classify_error(exc: Exception) -> str:
    s = str(exc).lower()
    if isinstance(exc, (ImportError, ModuleNotFoundError)):
        return "UNAVAILABLE"
    if any(t in s for t in ("missing", "absent", "not found", "unavailable", "not installed")):
        return "UNAVAILABLE"
    return "FAILED"


def compute_one(components, props, predictor, model_id, total_air_ppm):
    cass = [cas for cas, _ in components]
    xs = core.liquid_mole_fractions(components, props)
    gammas = predictor.predict(model_id, cass, xs, 298.15, props)
    partial, ys, terms, shares = core.evaluate_model(components, props, gammas, xs, total_air_ppm)
    mws = np.asarray([props[cas]["MW_g_mol"] for cas in cass], dtype=float)
    vap_wt = (ys * mws) / float(np.dot(ys, mws))
    if not np.isclose(ys.sum(), 1.0, atol=1e-10):
        raise AssertionError("Vapor mole fraction mass balance failed")
    if not np.isclose(vap_wt.sum(), 1.0, atol=1e-10):
        raise AssertionError("Vapor mass fraction mass balance failed")
    if not np.isclose(shares.sum(), 1.0, atol=1e-10):
        raise AssertionError("Screening shares mass balance failed")
    return xs, np.asarray(gammas), partial, ys, vap_wt, terms, shares


def run(args):
    if not 0 < args.total_air_ppm < 1e6:
        raise ValueError("total-air-ppm must satisfy 0 < ppm < 1000000")
    doses = list(dict.fromkeys(float(d) for d in args.doses))
    if any(not math.isfinite(d) or d <= 0 or d >= 100 for d in doses):
        raise ValueError("All injected doses must be finite and satisfy 0 < wt% < 100")
    candidates = list(dict.fromkeys(args.candidates))
    models = list(dict.fromkeys(args.models))
    if "raoult" not in models:
        models.insert(0, "raoult")  # explicit reference, not fallback
    props = core.load_reference(args.reference)
    registry = load_registry(args.registry)
    for target in candidates:
        if target not in props:
            raise ValueError(f"No documented property row for candidate CAS {target}")
        if target not in registry:
            raise ValueError(f"Missing legal-screening review entry for {target}")
    products = core.load_product_compositions(args.input, [str(p) for p in args.product_ids])
    predictor = core.GammaModels(args.hanna_repo)
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    coverage, status, details, target_rows = [], [], [], []
    for pid, base in products.items():
        for target in candidates:
            exists = any(c == target for c, _ in base)
            coverage.append(dict(product_id=pid, target_CAS=target,
                                 chemical=registry[target]["chemical_label"],
                                 result="SKIP_ALREADY_IN_EPA_BASELINE" if exists else "ELIGIBLE_COUNTERFACTUAL",
                                 original_component_count=len(base),
                                 explanation="Do not duplicate or overwrite EPA-reported ingredients" if exists else
                                 "Addition is hypothetical, never analytically detected in this product"))
            if exists:
                continue
            for dose in doses:
                comps = make_trace_mixture(base, target, dose)
                cass = [c for c, _ in comps]
                relation = cutoff_relation(dose, registry[target])
                for model in models:
                    t0 = time.perf_counter()
                    st = dict(product_id=pid, target_CAS=target, target_wt_pct=dose,
                              model_id=model, status="OK", error_message="",
                              elapsed_seconds=None, base_data="EPA_ChemExpo_reported_SDS",
                              target_data="hypothetical_counterfactual")
                    try:
                        xs, gamma, pp, ys, vap_wt, terms, shares = compute_one(
                            comps, props, predictor, model, args.total_air_ppm)
                        liquid_rank = pd.Series([-w for _, w in comps]).rank(method="min").astype(int)
                        index_rank = pd.Series(-shares).rank(method="min").astype(int)
                        for i, (cas, wt) in enumerate(comps):
                            row = dict(product_id=pid, target_CAS=target, target_wt_pct=dose,
                                       model_id=model, CAS=cas, chemical=props[cas]["chemical"],
                                       composition_source="COUNTERFACTUAL_INSERTION" if cas == target else "EPA_REPORTED_RESCALED",
                                       liquid_wt_pct=wt, liquid_mol_pct=float(xs[i]*100),
                                       gamma=float(gamma[i]), vapor_pressure_kPa=props[cas]["VP_25C_kPa"],
                                       partial_pressure_kPa=float(pp[i]),
                                       vapor_mol_pct=float(ys[i]*100), vapor_wt_pct=float(vap_wt[i]*100),
                                       vapor_wt_enrichment_fold=float(vap_wt[i]*100/wt),
                                       OEL_NIOSH_ppm=props[cas]["OEL_NIOSH_ppm"],
                                       hypothetical_component_air_ppm=float(ys[i]*args.total_air_ppm),
                                       OEL_weighted_term=float(terms[i]), OEL_weighted_share_pct=float(shares[i]*100),
                                       liquid_rank=int(liquid_rank[i]), screening_rank=int(index_rank[i]),
                                       rank_change=int(liquid_rank[i]-index_rank[i]),
                                       interpretation="comparative_screening_not_exposure_or_mixture_toxicity")
                            details.append(row)
                            if cas == target:
                                target_rows.append(dict(product_id=pid, target_CAS=target,
                                        chemical=registry[target]["chemical_label"], target_wt_pct=dose,
                                        model_id=model, gamma=float(gamma[i]),
                                        vapor_mol_pct=float(ys[i]*100), vapor_wt_pct=float(vap_wt[i]*100),
                                        vapor_wt_enrichment_fold=float(vap_wt[i]*100/wt),
                                        OEL_NIOSH_ppm=props[cas]["OEL_NIOSH_ppm"],
                                        OEL_weighted_share_pct=float(shares[i]*100),
                                        liquid_rank=int(liquid_rank[i]), screening_rank=int(index_rank[i]),
                                        rank_change=int(liquid_rank[i]-index_rank[i]),
                                        screening_cutoff_wt_pct=registry[target]["screening_cutoff_wt_pct"],
                                        cutoff_relation=relation,
                                        legal_review_status=registry[target]["legal_review_status"],
                                        legal_omission_assessment="NOT_PERFORMED",
                                        real_product_detection="NOT_DETERMINED_COUNTERFACTUAL"))
                    except Exception as exc:
                        st["status"] = classify_error(exc)
                        st["error_message"] = f"{type(exc).__name__}: {str(exc)[:450]}"
                    st["elapsed_seconds"] = round(time.perf_counter()-t0, 4)
                    status.append(st)
                    print(f"{pid}: {target} {dose:.3f} wt%, {model}: {st['status']}", flush=True)

    cv = pd.DataFrame(coverage)
    st = pd.DataFrame(status)
    det = pd.DataFrame(details)
    tar = pd.DataFrame(target_rows)
    core.csv_save(cv, output / "01_scenario_coverage.csv")
    core.csv_save(st, output / "02_model_status.csv")
    if not det.empty:
        core.csv_save(det, output / "03_component_predictions.csv")
    if not tar.empty:
        core.csv_save(tar, output / "04_trace_screening.csv")
        key = ["product_id", "target_CAS", "target_wt_pct"]
        grouping = tar.groupby(key, as_index=False).agg(
            models_success=("model_id", "nunique"),
            vapor_wt_min=("vapor_wt_pct", "min"), vapor_wt_max=("vapor_wt_pct", "max"),
            share_min=("OEL_weighted_share_pct", "min"),
            share_max=("OEL_weighted_share_pct", "max"),
            rank_min=("screening_rank", "min"),rank_max=("screening_rank", "max"))
        grouping["vapor_wt_range_pct_points"] = grouping.vapor_wt_max-grouping.vapor_wt_min
        grouping["share_range_pct_points"] = grouping.share_max-grouping.share_min
        core.csv_save(grouping, output / "05_model_agreement.csv")
        core.csv_save(tar[["product_id", "target_CAS", "target_wt_pct", "model_id",
                           "liquid_rank", "screening_rank", "rank_change", "OEL_weighted_share_pct"]],
                      output / "06_rank_discordance.csv")
    metadata = dict(date_policy="2026-10-22; separate full legal verification required",
                    input=str(args.input), input_status="EPA-reported, not independently measured",
                    property_reference=str(args.reference), candidate_registry=str(args.registry),
                    model_ids=models, target_cas=candidates, doses_wt_pct=doses,
                    temperature_C=25, vp_temperatures="25C only, not temperature extrapolation",
                    reference_total_air_ppm=args.total_air_ppm, reference_total_air_is_measured=False,
                    exposure_limit_source="US NIOSH REL, not Korean OEL",
                    quantified_worker_exposure=False, toxicological_additivity_verified=False,
                    tested_regulatory_omission=False, original_HANNA_validation_done=False,
                    successfully_computed_runs=int((st.status == "OK").sum()),
                    unavailable_or_failed_runs=int((st.status != "OK").sum()),
                    eligible_product_chemical_pairs=int((cv.result=="ELIGIBLE_COUNTERFACTUAL").sum()),
                    interpretation="hypothetical screening; NOT law compliance, hazard proof, or personal exposure")
    (output / "07_analysis_metadata.json").write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding="utf-8")
    print("\nOUTPUT", output.resolve(), "successful", metadata["successfully_computed_runs"],
          "failed", metadata["unavailable_or_failed_runs"])
    if args.strict_models and (st.status != "OK").any():
        raise RuntimeError("At least one requested model failed; see 02_model_status.csv")
    return metadata


if __name__ == "__main__":
    try:
        run(arguments())
    except Exception as error:
        print(f"\n[ERROR] {type(error).__name__}: {error}", file=sys.stderr)
        raise SystemExit(1)
