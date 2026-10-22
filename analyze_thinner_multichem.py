# -*- coding: utf-8 -*-
"""Counterfactual MULTI-CHEMICAL screening of EPA reported paint thinners.
25 C ONLY. Not observed benzene/methanol/acetone in these EPA products.
Uses the installed Raoult/UNIFAC/HANNA backends from the previous script.
Results are comparative OEL-weighted indices, NOT measured exposure or law compliance.
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
BASE_INPUT = ROOT / "chemexpo_paint_thinner_output" / "04_model_input_candidates_REVIEW.csv"
BASE_REGISTRY = ROOT / "multichem_screening_registry.csv"
BASE_OUTPUT = ROOT / "chemexpo_paint_thinner_output" / "multichem_screening_output"
CANDIDATES = ("71-43-2", "67-56-1", "67-64-1")
DOSES = (0.039, 0.099, 0.100, 0.990, 1.000)


def arguments(argv=None):
    p = argparse.ArgumentParser(description="Multi-chemical 25C vapor screening")
    p.add_argument("--input", type=Path, default=BASE_INPUT)
    p.add_argument("--reference", type=Path, default=ROOT / "thinner_vle_reference_25C.csv")
    p.add_argument("--registry", type=Path, default=BASE_REGISTRY)
    p.add_argument("--output-dir", type=Path, default=BASE_OUTPUT)
    p.add_argument("--hanna-repo", type=Path, default=Path(r"D:\fire\HANNA"))
    p.add_argument("--product-ids", nargs="+", default=list(core.RESEARCH_PRODUCT_IDS))
    p.add_argument("--candidates", nargs="+", default=list(CANDIDATES))
    p.add_argument("--doses", nargs="+", type=float, default=list(DOSES))
    p.add_argument("--models", nargs="+", choices=core.MODEL_IDS, default=["raoult", "unifac", "hanna"])
    p.add_argument("--total-air-ppm", type=float, default=100.)
    p.add_argument("--strict-models", action="store_true")
    return p.parse_args(argv)


def registry_load(path):
    data = pd.read_csv(path, dtype={"CAS": str}, encoding="utf-8-sig", keep_default_na=False)
    required = {"CAS", "chemical_label", "screening_cutoff_wt_pct",
                "screening_cutoff_basis", "legal_review_status"}
    if not required.issubset(data.columns):
        raise ValueError(f"Missing registry columns: {required-set(data.columns)}")
    if data.CAS.duplicated().any():
        raise ValueError("Duplicated registry CAS")
    return {str(row.CAS).strip(): row.to_dict() for _, row in data.iterrows()}


def make_trace_mixture(base, cas, dose):
    if not math.isfinite(dose) or not 0 < dose < 100:
        raise ValueError("Injection dose must satisfy 0 < wt% < 100")
    if any(c == cas for c, _ in base):
        raise ValueError(f"ALREADY_PRESENT: {cas}")
    if not np.isclose(sum(w for _, w in base), 100., atol=.05):
        raise ValueError("EPA baseline fractions do not sum to 100")
    return [(c, w*(100.-dose)/100.) for c, w in base]+[(cas, float(dose))]


def cutoff_relation(dose, entry):
    value = str(entry.get("screening_cutoff_wt_pct", "")).strip()
    if not value:
        return "NO_REVIEWED_NUMERIC_CUTOFF"
    cutoff = float(value)
    if not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("Invalid registry cutoff")
    if math.isclose(dose, cutoff, abs_tol=1e-10):
        return "AT_NUMERIC_CUTOFF_NOT_LEGAL_DETERMINATION"
    return ("BELOW_NUMERIC_CUTOFF_NOT_LEGAL_DETERMINATION" if dose < cutoff
            else "ABOVE_NUMERIC_CUTOFF_NOT_LEGAL_DETERMINATION")


def compute_one(components, props, engine, model, total_air):
    cas = [c for c, _ in components]
    xs = core.liquid_mole_fractions(components, props)
    gamma = np.asarray(engine.predict(model, cas, xs, 298.15, props), dtype=float)
    part, ys, terms, shares = core.evaluate_model(components, props, gamma, xs, total_air)
    mws = np.asarray([props[c]["MW_g_mol"] for c in cas])
    vap_wt = ys*mws/float(np.dot(ys, mws))
    for vector in (xs, ys, vap_wt, shares):
        if not np.isclose(np.sum(vector), 1., atol=1e-10):
            raise AssertionError("Mass/mole/screening share balance failure")
    return xs, gamma, part, ys, vap_wt, terms, shares


def fail_status(err):
    s = str(err).lower()
    if isinstance(err, (ImportError, ModuleNotFoundError)) or any(
            t in s for t in ("missing", "absent", "not found", "unavailable", "not installed")):
        return "UNAVAILABLE"
    return "FAILED"


def run(args):
    if not math.isfinite(args.total_air_ppm) or not 0 < args.total_air_ppm < 1e6:
        raise ValueError("Invalid assumed total-air-ppm")
    doses = list(dict.fromkeys(float(d) for d in args.doses))
    if any(not math.isfinite(d) or not 0 < d < 100 for d in doses):
        raise ValueError("Invalid trace dose")
    candidates = list(dict.fromkeys(args.candidates))
    models = list(dict.fromkeys(args.models))
    if "raoult" not in models:
        models.insert(0, "raoult")
    props = core.load_reference(args.reference)
    registry = registry_load(args.registry)
    for cas in candidates:
        if cas not in props or cas not in registry:
            raise ValueError(f"Missing reviewed property or registry entry: {cas}")
    products = core.load_product_compositions(args.input, [str(x) for x in args.product_ids])
    engine = core.GammaModels(args.hanna_repo)
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    coverage, status, components, targets = [], [], [], []
    for product, base in products.items():
        for target in candidates:
            exists = any(c == target for c, _ in base)
            coverage.append({"product_id": product, "target_CAS": target,
                             "chemical": registry[target]["chemical_label"],
                             "result": "SKIP_ALREADY_IN_EPA_BASELINE" if exists else "ELIGIBLE_COUNTERFACTUAL",
                             "explanation": "Never overwrite or duplicate an EPA ingredient" if exists
                                            else "Hypothetical addition, not an EPA measured compound"})
            if exists:
                continue
            for dose in doses:
                mixture = make_trace_mixture(base, target, dose)
                relation = cutoff_relation(dose, registry[target])
                cass = [c for c, _ in mixture]
                for model in models:
                    st = {"product_id": product, "target_CAS": target, "target_wt_pct": dose,
                          "model_id": model, "status": "OK", "error_message": "", "elapsed_seconds": None}
                    started = time.perf_counter()
                    try:
                        xs, gamma, pp, ys, vap_wt, terms, shares = compute_one(
                            mixture, props, engine, model, args.total_air_ppm)
                        liquid_rank = pd.Series([-w for _, w in mixture]).rank(method="min").astype(int)
                        screening_rank = pd.Series(-shares).rank(method="min").astype(int)
                        for j, (cas, wt) in enumerate(mixture):
                            row = {"product_id": product, "target_CAS": target,
                                   "target_wt_pct": dose, "model_id": model,
                                   "CAS": cas, "chemical": props[cas]["chemical"],
                                   "composition_source": "COUNTERFACTUAL_INSERTION" if cas == target
                                   else "EPA_REPORTED_RESCALED",
                                   "liquid_wt_pct": wt, "liquid_mol_pct": float(xs[j]*100),
                                   "gamma": float(gamma[j]),
                                   "vapor_pressure_kPa": props[cas]["VP_25C_kPa"],
                                   "partial_pressure_kPa": float(pp[j]),
                                   "vapor_mol_pct": float(ys[j]*100),
                                   "vapor_wt_pct": float(vap_wt[j]*100),
                                   "vapor_wt_enrichment_fold": float(vap_wt[j]*100/wt),
                                   "OEL_NIOSH_ppm": props[cas]["OEL_NIOSH_ppm"],
                                   "hypothetical_component_air_ppm": float(ys[j]*args.total_air_ppm),
                                   "OEL_weighted_term": float(terms[j]),
                                   "OEL_weighted_share_pct": float(shares[j]*100),
                                   "liquid_rank": int(liquid_rank[j]),
                                   "screening_rank": int(screening_rank[j]),
                                   "rank_change": int(liquid_rank[j]-screening_rank[j]),
                                   "interpretation": "COMPARATIVE_ONLY_NOT_EXPOSURE_OR_TOXICITY"}
                            components.append(row)
                            if cas == target:
                                targets.append({"product_id": product, "target_CAS": target,
                                                "chemical": registry[target]["chemical_label"],
                                                "target_wt_pct": dose, "model_id": model,
                                                "gamma": float(gamma[j]),
                                                "vapor_mol_pct": float(ys[j]*100),
                                                "vapor_wt_pct": float(vap_wt[j]*100),
                                                "vapor_wt_enrichment_fold": float(vap_wt[j]*100/wt),
                                                "OEL_NIOSH_ppm": props[cas]["OEL_NIOSH_ppm"],
                                                "OEL_weighted_share_pct": float(shares[j]*100),
                                                "liquid_rank": int(liquid_rank[j]),
                                                "screening_rank": int(screening_rank[j]),
                                                "rank_change": int(liquid_rank[j]-screening_rank[j]),
                                                "screening_cutoff_wt_pct": registry[target]["screening_cutoff_wt_pct"],
                                                "cutoff_relation": relation,
                                                "legal_review_status": registry[target]["legal_review_status"],
                                                "legal_omission_assessment": "NOT_PERFORMED",
                                                "real_product_detection": "NOT_DETERMINED_COUNTERFACTUAL"})
                    except Exception as exc:
                        st["status"] = fail_status(exc)
                        st["error_message"] = f"{type(exc).__name__}: {str(exc)[:400]}"
                    finally:
                        st["elapsed_seconds"] = round(time.perf_counter()-started, 4)
                        status.append(st)
                        print(f"{product} {target} {dose:.3f} wt% {model}: {st['status']}", flush=True)

    out = core.csv_save
    out(pd.DataFrame(coverage), output/"01_scenario_coverage.csv")
    states = pd.DataFrame(status)
    out(states, output/"02_model_status.csv")
    detail = pd.DataFrame(components)
    target = pd.DataFrame(targets)
    if not detail.empty:
        out(detail, output/"03_component_predictions.csv")
    if not target.empty:
        out(target, output/"04_trace_screening.csv")
        group = target.groupby(["product_id","target_CAS","target_wt_pct"], as_index=False).agg(
            models_success=("model_id","nunique"),
            vapor_wt_min=("vapor_wt_pct","min"),vapor_wt_max=("vapor_wt_pct","max"),
            share_min=("OEL_weighted_share_pct","min"),
            share_max=("OEL_weighted_share_pct","max"),
            rank_min=("screening_rank","min"),rank_max=("screening_rank","max"))
        group["vapor_wt_range_pct_points"] = group.vapor_wt_max-group.vapor_wt_min
        group["share_range_pct_points"] = group.share_max-group.share_min
        out(group,output/"05_model_agreement.csv")
        out(target[["product_id","target_CAS","target_wt_pct","model_id","liquid_rank",
                   "screening_rank","rank_change","OEL_weighted_share_pct"]],
            output/"06_rank_discordance.csv")
    meta = {"input":str(args.input),"registry":str(args.registry),
            "temperature_C":25,"temperature_extrapolation":False,
            "model_ids":models,"target_CAS":candidates,"doses_wt_pct":doses,
            "source":"EPA-reported formulas with hypothetical inserted chemicals",
            "total_air_ppm_hypothetical":args.total_air_ppm,
            "oel_basis":"US NIOSH REL, not Korean OEL",
            "worker_exposure_measured":False,"mixed_toxicology_additivity_verified":False,
            "Korean_legal_compliance_evaluated":False,
            "model_validation_against_independent_VLE":False,
            "successful_calculations":int((states.status=="OK").sum()),
            "unsuccessful_calculations":int((states.status!="OK").sum()),
            "interpretation":"screening only; no measured exposure, legal omission, or actual mixture risk"}
    (output/"07_analysis_metadata.json").write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding="utf-8")
    print("OUTPUT",output.resolve(),"OK",meta["successful_calculations"],
          "other",meta["unsuccessful_calculations"])
    if args.strict_models and (states.status!="OK").any():
        raise RuntimeError("At least one model failed; inspect 02_model_status.csv")
    return meta


if __name__=="__main__":
    try:
        run(arguments())
    except Exception as err:
        print(f"[ERROR] {type(err).__name__}: {err}",file=sys.stderr)
        sys.exit(1)
