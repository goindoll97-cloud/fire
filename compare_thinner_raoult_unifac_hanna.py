# -*- coding: utf-8 -*-
"""ChemExpo paint thinner: Raoult vs UNIFAC vs HANNA, 25 deg C pilot.

RESEARCH DESIGN
- Real publicly *reported* product composition from EPA ChemExpo.
- Benzene doses (0, 0.039, 0.040, 0.099, 0.100 wt%) are COUNTERFACTUAL
  additions, NOT measured in those US paint thinner products.
- Raoult, original UNIFAC, optional Dortmund UNIFAC, and optional pre-trained
  multi-component HANNA each predict gamma_i independently.
- y_i = (x_i * gamma_i * P_i^sat) / sum_j(x_j*gamma_j*P_j^sat).
- The relative OEL-weighted contribution is a *screening index*, NOT actual
  occupational exposure, cancer risk, regulatory mixture hazard, or proof of
  lawful MSDS omission. Component toxicological additivity is not established.
- A model that fails is marked FAILED/UNAVAILABLE; no gamma=1 fallback.

Run:
  python compare_thinner_raoult_unifac_hanna.py
  python compare_thinner_raoult_unifac_hanna.py --models raoult unifac
  python compare_thinner_raoult_unifac_hanna.py --input "D:\\fire\\chemexpo_paint_thinner_output\\04_model_input_candidates_REVIEW.csv" --hanna-repo "D:\\fire\\HANNA"

Dependencies: numpy, pandas; optional thermo (UNIFAC); optional official HANNA
GitHub repository, its conda environment, and its trained checkpoints (HANNA).
See accompanying README for setup, source references, and limitations.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = Path(r"D:\fire\chemexpo_paint_thinner_output\04_model_input_candidates_REVIEW.csv")
DEFAULT_REFERENCE = SCRIPT_DIR / "thinner_vle_reference_25C.csv"
DEFAULT_HANNA_REPO = Path(r"D:\fire\HANNA")
RESEARCH_PRODUCT_IDS = ("1208563", "273605", "1183746")
LABELS = {
    "1208563": "EPA-1 Silver paint thinner",
    "273605": "EPA-2 Lacquer thinner",
    "1183746": "EPA-3 Vinyl thinner",
}
BENZENE_CAS = "71-43-2"
DEFAULT_BENZENE_DOSES = (0.0, 0.039, 0.040, 0.099, 0.100)
MODEL_IDS = ("raoult", "unifac", "dortmund", "hanna")
MODEL_NAMES = {
    "raoult": "Raoult_ideal",
    "unifac": "UNIFAC_original",
    "dortmund": "UNIFAC_Dortmund2006",
    "hanna": "HANNA_pretrained_2026",
}


def csv_save(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False, encoding="utf-8-sig")


def valid_float(v, field: str) -> float:
    val = float(v)
    if not math.isfinite(val):
        raise ValueError(f"Invalid {field}: {v!r}")
    return val


def load_reference(path: Path) -> dict[str, dict]:
    if not path.is_file():
        raise FileNotFoundError(
            f"25C reference table absent: {path}\n"
            "Place thinner_vle_reference_25C.csv next to the Python script."
        )
    d = pd.read_csv(path, dtype={"CAS": str}, encoding="utf-8-sig")
    required = ["CAS", "chemical", "MW_g_mol", "VP_25C_kPa", "OEL_NIOSH_ppm", "SMILES"]
    if any(k not in d for k in required):
        raise ValueError(f"Reference file needs columns {required}")
    if d.CAS.duplicated().any():
        raise ValueError("Duplicated reference CAS keys")
    ans = {}
    for _, v in d.iterrows():
        ca = str(v.CAS).strip()
        item = v.to_dict()
        for c in ("MW_g_mol", "VP_25C_kPa", "OEL_NIOSH_ppm"):
            item[c] = valid_float(item[c], f"{ca}:{c}")
            if item[c] <= 0:
                raise ValueError(f"Nonpositive {c} for CAS {ca}")
        item["SMILES"] = str(item["SMILES"]).strip()
        if not item["SMILES"]:
            raise ValueError(f"SMILES missing: {ca}")
        ans[ca] = item
    if BENZENE_CAS not in ans:
        raise ValueError("Benzene (71-43-2) not included in property reference")
    return ans


def load_product_compositions(path: Path, product_ids: list[str]) -> dict[str, list[tuple[str, float]]]:
    if not path.is_file():
        raise FileNotFoundError(f"ChemExpo 04 input not found: {path}")
    d = pd.read_csv(path, dtype={"Data Document ID": str, "CAS": str}, encoding="utf-8-sig")
    mandatory = ["Data Document ID", "CAS", "wt_percent", "Weight Fraction Type"]
    missing = [c for c in mandatory if c not in d.columns]
    if missing:
        raise ValueError(f"Missing ChemExpo columns {missing} in {path}")
    records = {}
    for pid in product_ids:
        s = d.loc[d["Data Document ID"].astype(str).str.strip() == pid].copy()
        if s.empty:
            raise ValueError(f"ChemExpo product document {pid} not found in 04 file")
        s["CAS"] = s["CAS"].astype(str).str.strip()
        s["wt_percent"] = pd.to_numeric(s["wt_percent"], errors="raise")
        if (s["wt_percent"] <= 0).any() or not np.isfinite(s["wt_percent"]).all():
            raise ValueError(f"Invalid liquid fraction for document {pid}")
        if s.CAS.duplicated().any():
            raise ValueError(f"Duplicate CAS entry for document {pid}; manually resolve before predicting")
        if any(str(v).strip().lower() != "reported" for v in s["Weight Fraction Type"]):
            raise ValueError(f"Not all composition rows are 'reported' in document {pid}")
        total = float(s["wt_percent"].sum())
        if abs(total - 100.0) > 0.05:
            raise ValueError(f"Sum of documented fractions is {total:.5g}%, not 100%, for {pid}")
        if (s.CAS == BENZENE_CAS).any():
            raise ValueError(f"EPA product {pid} already has benzene; this dose-insertion policy would duplicate it")
        records[pid] = [(str(r.CAS), float(r.wt_percent)) for _, r in s.iterrows()]
    return records


def make_counterfactual(base: list[tuple[str, float]], benzene_wt_pct: float):
    # Preserve original component proportions while adding benzene.
    result = [(c, w * (100.0 - benzene_wt_pct) / 100.0) for c, w in base]
    if benzene_wt_pct > 0:
        result.append((BENZENE_CAS, float(benzene_wt_pct)))
    return result


def liquid_mole_fractions(components, properties):
    for cas, _ in components:
        if cas not in properties:
            raise ValueError(f"No reference MW/VP/OEL/SMILES for CAS {cas}")
    moles = np.asarray([w / properties[cas]["MW_g_mol"] for cas, w in components], float)
    return moles / moles.sum()


@contextmanager
def cd_to(folder: Path):
    before = Path.cwd()
    os.chdir(folder)
    try:
        yield
    finally:
        os.chdir(before)


class GammaModels:
    """Actual model implementations. No fallback to ideal coefficients on failure."""
    def __init__(self, hanna_repo: Path):
        self.hanna_repo = hanna_repo
        self.predictor = None
        self.unifac_import = None

    def predict(self, model: str, cas_list: list[str], xs: np.ndarray,
                temperature_k: float, properties: dict) -> np.ndarray:
        if model == "raoult":
            return np.ones(len(cas_list), float)
        if model in ("unifac", "dortmund"):
            return self.unifac(model, cas_list, xs, temperature_k)
        if model == "hanna":
            return self.hanna(cas_list, xs, temperature_k, properties)
        raise ValueError(f"Unknown model: {model}")

    @staticmethod
    def unifac(model, cas_list, xs, t):
        try:
            from thermo.unifac import (UNIFAC, UFIP, UFSG, DOUFIP2006, DOUFSG,
                                       UNIFAC_group_assignment_DDBST)
        except ImportError as exc:
            raise RuntimeError("thermo missing: pip install thermo (in the SAME Python environment)") from exc
        kind = "MODIFIED_UNIFAC" if model == "dortmund" else "UNIFAC"
        groups = []
        for cas in cas_list:
            assign = UNIFAC_group_assignment_DDBST(cas, kind)
            if not assign:
                raise RuntimeError(f"DDBST {kind} subgroup mapping unavailable for CAS {cas}")
            groups.append(assign)
        if model == "dortmund":
            obj = UNIFAC.from_subgroups(T=t, xs=xs.tolist(), chemgroups=groups,
                                       version=1, interaction_data=DOUFIP2006, subgroups=DOUFSG)
        else:
            obj = UNIFAC.from_subgroups(T=t, xs=xs.tolist(), chemgroups=groups,
                                       version=0, interaction_data=UFIP, subgroups=UFSG)
        return np.asarray(obj.gammas(), dtype=float).reshape(-1)

    def hanna(self, cas_list, xs, t, properties):
        repo = self.hanna_repo
        if not repo.is_dir() or not (repo / "utils" / "HANNA_predictor.py").is_file():
            raise RuntimeError(f"Official HANNA repository absent in {repo}. Clone it and install environment.yml")
        model_path = repo / "models" / "HANNA" / "ensemble"
        if not model_path.is_dir() or not list(model_path.glob("HANNA_parameters_binary*.pt")):
            raise RuntimeError(f"HANNA trained ensemble checkpoints missing in {model_path}")
        if str(repo) not in sys.path:
            sys.path.insert(0, str(repo))
        with cd_to(repo):
            if self.predictor is None:
                try:
                    from utils.HANNA_predictor import HANNA_Predictor
                except ImportError as exc:
                    raise RuntimeError("HANNA dependencies absent; activate environment.yml's conda environment") from exc
                self.predictor = HANNA_Predictor()
            smiles = [properties[cas]["SMILES"] for cas in cas_list]
            # Official HANNA returns ln(gamma), not gamma.
            ln_gamma, _ = self.predictor.predict(
                smiles_list=smiles,
                molar_fractions=np.asarray([xs], float),
                temperature=float(t),
                verbose=False,
            )
        result = np.asarray(ln_gamma, float)
        if result.size != len(cas_list):
            raise ValueError(f"Unexpected HANNA ln_gamma shape {result.shape}")
        return np.exp(result.reshape(-1))


def evaluate_model(components, properties, gammas, xs, total_air_ppm):
    cas_list = [cas for cas, _ in components]
    gammas = np.asarray(gammas, float).reshape(-1)
    if gammas.shape != (len(cas_list),):
        raise ValueError("Incorrect activity coefficient shape")
    if not np.isfinite(gammas).all() or (gammas <= 0).any():
        raise ValueError(f"Invalid activity coefficients: {gammas}")
    psats = np.asarray([properties[c]["VP_25C_kPa"] for c in cas_list], float)
    partials = xs * gammas * psats
    if not np.isfinite(partials).all() or partials.sum() <= 0:
        raise ValueError("Invalid partial pressures")
    ys = partials / partials.sum()
    oels = np.asarray([properties[c]["OEL_NIOSH_ppm"] for c in cas_list], float)
    rel_terms = ys / oels
    shares = rel_terms / rel_terms.sum()
    return partials, ys, rel_terms, shares


def parse_args():
    p = argparse.ArgumentParser(description="Compare Raoult/UNIFAC/HANNA for ChemExpo paint thinner")
    p.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="04_model_input_candidates_REVIEW.csv path")
    p.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE, help="25C reference table CSV")
    p.add_argument("--hanna-repo", type=Path, default=DEFAULT_HANNA_REPO, help="local clone of official HANNA")
    p.add_argument("--output-dir", type=Path, default=None, help="default: folder beside input CSV")
    p.add_argument("--models", nargs="+", choices=MODEL_IDS, default=["raoult", "unifac", "hanna"])
    p.add_argument("--product-ids", nargs="+", default=list(RESEARCH_PRODUCT_IDS))
    p.add_argument("--benzene-doses", nargs="+", type=float, default=list(DEFAULT_BENZENE_DOSES))
    p.add_argument("--strict-models", action="store_true", help="Fail after export if any selected model was not calculated")
    p.add_argument("--total-air-ppm", type=float, default=100.0,
                   help="Hypothetical total air vapor ppm, NOT measured. Relative share is independent of this value")
    return p.parse_args()


def run(args) -> dict:
    # 25C reference vapor pressures are FIXED. DO NOT vary temperature without updating VP data.
    if not 0 < args.total_air_ppm < 1e6:
        raise ValueError("total-air-ppm must be >0 and <1e6")
    for d in args.benzene_doses:
        if not (0 <= d < 100):
            raise ValueError(f"Invalid benzene dose {d}")
    model_ids = list(dict.fromkeys(args.models))
    if "raoult" not in model_ids:
        model_ids.insert(0, "raoult")  # reference comparison always included
    out = args.output_dir if args.output_dir is not None else (args.input.parent / "thinner_vle_comparison_output")
    out.mkdir(parents=True, exist_ok=True)
    props = load_reference(args.reference)
    products = load_product_compositions(args.input, [str(s) for s in args.product_ids])
    predictor = GammaModels(args.hanna_repo)
    status, detail, benzene = [], [], []
    number = 0
    for pid, base in products.items():
        for dose in args.benzene_doses:
            dose = float(dose)
            components = make_counterfactual(base, dose)
            cas_list = [cas for cas, _ in components]
            xs = liquid_mole_fractions(components, props)
            for model in model_ids:
                number += 1
                started = time.perf_counter()
                entry = dict(product_id=pid, product_name=LABELS.get(pid, pid),
                             benzene_wt_pct_scenario=dose, model=MODEL_NAMES[model],
                             model_id=model, status="OK", error_message="", elapsed_sec=np.nan,
                             n_components=len(components), source_type="counterfactual_not_measured")
                try:
                    gamma = predictor.predict(model, cas_list, xs, 298.15, props)
                    partial, ys, rel_terms, shares = evaluate_model(
                        components, props, gamma, xs, args.total_air_ppm)
                    # Mass balance hard checks, not silent correction.
                    if not np.isclose(sum(w for _, w in components), 100.0, atol=1e-7):
                        raise AssertionError("liquid mass balance failure")
                    if not np.isclose(ys.sum(), 1.0, atol=1e-10):
                        raise AssertionError("vapor composition mass balance failure")
                    for i, (cas, wt) in enumerate(components):
                        detail.append(dict(
                            product_id=pid, product_name=LABELS.get(pid, pid),
                            benzene_wt_pct_scenario=dose, model=MODEL_NAMES[model], model_id=model,
                            chemical=props[cas]["chemical"], CAS=cas,
                            component_basis=("hypothetical_benzene_dose" if cas == BENZENE_CAS
                                             else "EPA_reported_original_scaled"),
                            liquid_wt_pct=wt, liquid_mol_fraction=float(xs[i]),
                            gamma=float(gamma[i]), VP_25C_kPa=props[cas]["VP_25C_kPa"],
                            partial_pressure_kPa=float(partial[i]), vapor_mol_pct=float(ys[i]*100),
                            OEL_NIOSH_ppm=props[cas]["OEL_NIOSH_ppm"],
                            hypothetical_air_ppm=float(args.total_air_ppm*ys[i]),
                            OEL_weighted_relative_index=float(rel_terms[i]),
                            OEL_weighted_share_pct=float(shares[i]*100),
                            risk_interpretation="screening_only_not_actual_worker_exposure",
                        ))
                    if dose > 0:
                        bz = cas_list.index(BENZENE_CAS)
                        benzene.append(dict(
                            product_id=pid, product_name=LABELS.get(pid, pid),
                            benzene_wt_pct_scenario=dose, under_0p1_pct=dose < 0.1,
                            legal_omission_status="NOT_DETERMINED",
                            model=MODEL_NAMES[model], model_id=model,
                            benzene_gamma=float(gamma[bz]),
                            benzene_vapor_mol_pct=float(ys[bz]*100),
                            benzene_OEL_weighted_share_pct=float(shares[bz]*100),
                            benzene_hypothetical_air_ppm=float(args.total_air_ppm*ys[bz]),
                            benzene_hypothetical_air_over_NIOSH_REL=float(args.total_air_ppm*ys[bz]/props[BENZENE_CAS]["OEL_NIOSH_ppm"]),
                            basis="hypothetical_scenario_not_measured",
                        ))
                except Exception as exc:
                    entry["status"] = "UNAVAILABLE" if (
                        isinstance(exc, (ImportError, ModuleNotFoundError)) or
                        any(s in str(exc).lower() for s in ("missing", "absent", "not found", "unavailable"))
                    ) else "FAILED"
                    entry["error_message"] = f"{type(exc).__name__}: {str(exc)[:420]}"
                finally:
                    entry["elapsed_sec"] = round(time.perf_counter() - started, 3)
                    status.append(entry)
                    print(f"[{number:>2}] {pid}, Bz {dose:.3f}%, {model:9s}: {entry['status']}"
                          + (f" ({entry['error_message'][:110]})" if entry["status"] != "OK" else ""), flush=True)

    status_df = pd.DataFrame(status)
    detail_df = pd.DataFrame(detail)
    bz_df = pd.DataFrame(benzene)
    csv_save(status_df, out / "01_model_status.csv")
    if not detail_df.empty:
        csv_save(detail_df, out / "02_component_predictions.csv")
    if not bz_df.empty:
        csv_save(bz_df, out / "03_benzene_comparison_long.csv")
        piv = bz_df.pivot(index=["product_id", "product_name", "benzene_wt_pct_scenario"],
                          columns="model_id", values=["benzene_vapor_mol_pct", "benzene_OEL_weighted_share_pct"])
        piv.columns = ["_".join((metric, model)) for metric, model in piv.columns]
        csv_save(piv.reset_index(), out / "04_benzene_comparison_wide.csv")
        anchor = bz_df.query("model_id == 'raoult'")[
            ["product_id", "benzene_wt_pct_scenario", "benzene_vapor_mol_pct", "benzene_OEL_weighted_share_pct"]
        ].rename(columns={"benzene_vapor_mol_pct": "raoult_benzene_vapor_mol_pct",
                          "benzene_OEL_weighted_share_pct": "raoult_benzene_OEL_share_pct"})
        diffs = bz_df.merge(anchor, on=["product_id", "benzene_wt_pct_scenario"], how="left", validate="many_to_one")
        diffs["delta_vapor_pct_points_vs_raoult"] = diffs["benzene_vapor_mol_pct"] - diffs["raoult_benzene_vapor_mol_pct"]
        diffs["delta_OEL_share_pct_points_vs_raoult"] = diffs["benzene_OEL_weighted_share_pct"] - diffs["raoult_benzene_OEL_share_pct"]
        csv_save(diffs, out / "05_differences_vs_raoult.csv")
        counts = bz_df.groupby(["product_id", "benzene_wt_pct_scenario"], as_index=False).agg(
            models_succeeded=("model_id", "nunique"),
            vapor_min=("benzene_vapor_mol_pct", "min"), vapor_max=("benzene_vapor_mol_pct", "max"),
            share_min=("benzene_OEL_weighted_share_pct", "min"),
            share_max=("benzene_OEL_weighted_share_pct", "max"),
        )
        counts["vapor_range_pct_points"] = counts.vapor_max-counts.vapor_min
        counts["share_range_pct_points"] = counts.share_max-counts.share_min
        csv_save(counts, out / "06_across_model_ranges.csv")
    meta = {
        "input": str(args.input), "reference": str(args.reference), "HANNA_repository": str(args.hanna_repo),
        "temperature_C": 25.0, "model_ids": model_ids,
        "product_ids": list(products), "benzene_wt_pct_scenarios": list(args.benzene_doses),
        "hypothetical_total_air_ppm": args.total_air_ppm,
        "real_EPA_product_benzene_measurements": False,
        "OEL_system": "NIOSH REL (US); Korea OEL not yet evaluated",
        "actual_worker_exposure": False,
        "legal_MSDS_compliance_assessed": False,
        "model_predictions_available": status_df.groupby("model_id")["status"].apply(lambda z: (z == "OK").sum()).to_dict(),
        "reproduction_note": "EPA report-derived original formulas plus hypothetical benzene; do not infer from these outputs legal adequacy or exposed-worker risk",
    }
    (out / "07_run_metadata.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    print("\nOUTPUT:", out.resolve())
    print("Successful calculations by model:", meta["model_predictions_available"])
    if args.strict_models and (status_df["status"] != "OK").any():
        raise RuntimeError("Some model runs were unavailable/failed. Check 01_model_status.csv; files were still saved.")
    return meta


if __name__ == "__main__":
    try:
        run(parse_args())
    except Exception as exc:
        print(f"\n[ERROR] {type(exc).__name__}: {exc}", file=sys.stderr)
        sys.exit(1)
