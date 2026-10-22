# -*- coding: utf-8 -*-
"""Composition-constrained inverse design of *hypothetical* MSDS information gaps.

Finds a local high-Korean-TWA-weighted-vapor-score composition near each
reported ChemExpo thinner formulation, keeping the trace chemical concentration
below the *preliminary hazard-specific* Korean MSDS Annex 6 cutoff.

IMPORTANT: This does NOT determine lawful ingredient omission, workplace
exposure, mixture toxicity or the correctness of an MSDS. No model is silently
substituted after a failure. The original EPA products do NOT contain the
hypothetically inserted trace chemical in these experiments.

Needs: numpy pandas scipy; optional thermo; local official HANNA repository.
  python optimize_msds_gap.py --optimize-model hanna --check-models raoult unifac hanna
  python optimize_msds_gap.py --optimize-model raoult --check-models raoult
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import minimize

import compare_thinner_raoult_unifac_hanna as core
import apply_korean_msds_oel as kr
import analyze_thinner_multichem as multi

ROOT = Path(__file__).resolve().parent
DEFAULT_OUT = ROOT / 'chemexpo_paint_thinner_output' / 'msds_optimized_screening_output'


@dataclass(frozen=True)
class BoundsPolicy:
    """Plausibility constraints relative to the original *nontrace* composition.

    max_reallocation_pp = total variation of percentage weights, not sum of
    absolute changes. 10 pp means up to 10 grams per 100 grams can change hands
    among main solvents. This is an exploration policy, not measured uncertainty.
    """
    max_reallocation_pp: float = 10.0
    relative_change: float = 0.50

    def validate(self) -> None:
        if not (0.0 < self.max_reallocation_pp <= 50.0):
            raise ValueError('max-reallocation-pp must be in (0,50]')
        if not (0.0 < self.relative_change < 1.0):
            raise ValueError('relative-change must be in (0,1)')


def parse_args(argv=None):
    p = argparse.ArgumentParser(description='Inverse design for Korean-MSDS vapor screening (25 C only)')
    p.add_argument('--input', type=Path, default=multi.DEFAULT_INPUT)
    p.add_argument('--reference', type=Path, default=core.DEFAULT_REFERENCE)
    p.add_argument('--registry', type=Path, default=multi.DEFAULT_REGISTRY)
    p.add_argument('--oel-json', type=Path, default=kr.OEL_JSON)
    p.add_argument('--msds-json', type=Path, default=kr.MSDS_JSON)
    p.add_argument('--hanna-repo', type=Path, default=Path(r'D:\fire\HANNA'))
    p.add_argument('--output-dir', type=Path, default=DEFAULT_OUT)
    p.add_argument('--product-ids', nargs='+', default=list(core.RESEARCH_PRODUCT_IDS))
    p.add_argument('--target-cas', nargs='+', default=['71-43-2'])
    p.add_argument('--doses', type=float, nargs='+', default=[0.039, 0.099])
    p.add_argument('--optimize-model', choices=core.MODEL_IDS, default='hanna')
    p.add_argument('--check-models', choices=core.MODEL_IDS, nargs='+', default=['raoult','unifac','hanna'])
    p.add_argument('--max-reallocation-pp', type=float, default=10.0)
    p.add_argument('--relative-change', type=float, default=0.5)
    p.add_argument('--maxiter', type=int, default=45)
    p.add_argument('--starts', type=int, default=2, help='Number of deterministic SLSQP starts per scenario')
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--strict-models', action='store_true', help='Nonzero exit if any cross-model check unavailable')
    p.add_argument('--exploratory-cutoff-wt-pct', type=float, default=None,
                   help='Explicit NONLEGAL what-if cutoff only for unreviewed CAS')
    return p.parse_args(argv)


def trace_cutoff(cas: str, registry: dict, annex6_lines: list[str], exploratory: float | None):
    """A numeric scenario cutoff is *never* proof of lawful omission."""
    item = registry.get(cas)
    if item is None:
        raise ValueError(f'Missing registry entry for CAS {cas}')
    raw = str(item.get('screening_cutoff_wt_pct', '')).strip()
    if raw:
        cutoff = float(raw)
        label = 'PRELIMINARY_HAZARD_SPECIFIC_CUTOFF_NOT_LEGAL_OMISSION'
        if cas == '71-43-2':
            joined = '\n'.join(annex6_lines)
            import re
            if not re.search(r'발암성[^\n]*0\.1\s*%', joined) or not math.isclose(cutoff,0.1,abs_tol=1e-8):
                raise ValueError('Benzene preliminary 0.1% assumption differs from Korean MSDS Annex 6')
        # All other CAS need independently verified GHS classification before
        # being called a regulatory threshold. Nothing above claims approved.
    else:
        if exploratory is None:
            raise ValueError(f'No verified numeric cutoff for CAS {cas}. '
                             'Use --exploratory-cutoff-wt-pct to label it a nonlegal scenario.')
        cutoff = float(exploratory)
        label = 'EXPLORATORY_USER_ASSUMED_CUTOFF_NOT_REGULATORY'
    if not (math.isfinite(cutoff) and 0 < cutoff < 100):
        raise ValueError('Cutoff must be in (0,100)')
    return cutoff, label


def normalized_base(base: list[tuple[str, float]]) -> np.ndarray:
    arr = np.array([float(w) for _, w in base], dtype=float)
    if len(arr) < 2 or not np.isfinite(arr).all() or (arr <= 0).any():
        raise ValueError('Baseline must have two or more finite positive components')
    if not np.isclose(arr.sum(),100.,atol=0.05):
        raise ValueError('Baseline fractions must total 100%')
    return arr / arr.sum()


def feasible(p: np.ndarray, p0: np.ndarray, policy: BoundsPolicy, tol=1e-6) -> bool:
    p, p0 = np.asarray(p,float), np.asarray(p0,float)
    return (p.shape == p0.shape and np.isfinite(p).all() and
            abs(p.sum()-1.) <= tol and
            (p >= p0*(1-policy.relative_change)-tol).all() and
            (p <= p0*(1+policy.relative_change)+tol).all() and
            0.5*np.abs(p-p0).sum() <= policy.max_reallocation_pp/100. + tol)


def mixture(base, target: str, dose: float, p: np.ndarray):
    if len(base) != len(p) or len(set(c for c,_ in base)) != len(base):
        raise ValueError('Invalid base composition')
    if target in {cas for cas,_ in base}:
        raise ValueError('Cannot add a trace chemical already in the original product')
    if not 0 < dose < 100:
        raise ValueError('Invalid dose')
    if not np.isfinite(p).all() or (p <= 0).any() or not np.isclose(sum(p), 1.,atol=1e-6):
        raise ValueError('Candidate composition is not a simplex')
    return [(cas,float(pi)*(100.-dose)) for (cas,_),pi in zip(base,p)] + [(target, float(dose))]


def evaluate(components: list[tuple[str,float]], target: str, model: str,
             predictor, props: dict, twa: dict[str,float]) -> dict[str,Any]:
    cas_list=[x[0] for x in components]
    if target not in cas_list or len(set(cas_list)) != len(cas_list):
        raise ValueError('Target missing or repeated CAS in mixture')
    if any(x not in twa for x in cas_list):
        raise ValueError(f'Korean TWA not mapped: {set(cas_list)-set(twa)}')
    xs=core.liquid_mole_fractions(components, props)
    gammas=np.asarray(predictor.predict(model,cas_list,xs,298.15,props),float).reshape(-1)
    if gammas.size != len(cas_list) or not np.isfinite(gammas).all() or (gammas<=0).any():
        raise ValueError('Invalid model activity coefficients')
    psats=np.array([props[c]['VP_25C_kPa'] for c in cas_list],float)
    partial=xs*gammas*psats
    if not np.isfinite(partial).all() or (partial<=0).any():
        raise ValueError('Invalid vapor-phase partial pressures')
    ys=partial/partial.sum()
    terms=ys/np.array([twa[c] for c in cas_list],float)
    if (terms<=0).any() or not np.isfinite(terms).all():
        raise ValueError('Invalid Korean TWA weighted terms')
    shares=terms/terms.sum()
    molar_mass=np.array([props[c]['MW_g_mol'] for c in cas_list],float)
    vapor_wt=ys*molar_mass/np.dot(ys,molar_mass)
    i=cas_list.index(target)
    target_liquid_wt=dict(components)[target]
    if not (np.isclose(ys.sum(),1.,atol=1e-10) and np.isclose(vapor_wt.sum(),1.,atol=1e-10)
            and np.isclose(shares.sum(),1.,atol=1e-10) and np.isclose(sum(w for _,w in components),100.,atol=1e-6)):
        raise AssertionError('Mass balance failed')
    return {
        'KR_share_pct':float(100*shares[i]),
        'vapor_mol_pct':float(100*ys[i]),
        'vapor_wt_pct':float(100*vapor_wt[i]),
        'vapor_wt_enrichment_fold':float(100*vapor_wt[i]/target_liquid_wt),
        'KR_relative_term':float(ys[i]/twa[target]),
        'KR_TWA_ppm':float(twa[target]),
        'KR_rank':int(1+np.sum(shares>shares[i]+1e-14)),
        'n_components':len(cas_list),
        'temperature_C':25.,
    }


def optimize_one(base, target, dose, model, predictor, props, twa, policy: BoundsPolicy,
                 starts=2, maxiter=45, seed=42):
    """Bounded multistart SLSQP; the best *verified feasible* result is retained.

    This is an inverse-design search, not a proof of global optimality. A
    surrogate AI model can be wrong, so results must be checked against other
    predictors and independent experimental VLE where available.
    """
    policy.validate()
    p0=normalized_base(base)
    if not (isinstance(starts,int) and starts>=1 and maxiter>=1):
        raise ValueError('starts and maxiter must both be positive integers')
    limit=policy.max_reallocation_pp/100.
    lower=p0*(1-policy.relative_change)
    upper=p0*(1+policy.relative_change)
    log=[]; calls=[0]; candidates=[]
    def score(p):
        if not feasible(p,p0,policy,tol=2e-6):
            return None
        result=evaluate(mixture(base,target,dose,p),target,model,predictor,props,twa)
        return result
    initial=score(p0)
    if initial is None:
        raise AssertionError('Baseline infeasible')
    candidates.append((p0.copy(),initial))
    log.append({'trial':'baseline','objective_KR_share_pct':initial['KR_share_pct'],
                'feasible':True,'solver_success':True,'evaluations':1,'solver_message':'baseline reported mixture'})
    rng=np.random.default_rng(seed)
    # Generate true-feasible deterministic starting positions with a small
    # mass transfer. Sampling only within the specified trust region prevents
    # the optimizer from inventing a totally different product.
    starting=[p0]
    for _ in range(starts-1):
        a,b=rng.choice(len(p0),size=2,replace=False)
        q=p0.copy()
        step=min(limit*.5, 0.35*p0[a], 0.35*p0[b],p0[a]-lower[a],upper[b]-p0[b])
        q[a]-=step;q[b]+=step
        starting.append(q)
    constraints=[{'type':'eq','fun':lambda p: float(np.sum(p)-1.)},
                 {'type':'ineq','fun':lambda p: float(limit-0.5*np.sum(np.abs(p-p0)))}]
    for k,start in enumerate(starting,1):
        score0=score(start)
        if score0 is not None:
            candidates.append((start.copy(),score0))
        cache:dict[tuple,float]={}
        def objective(p):
            calls[0]+=1
            # Numerical optimizer may transiently test points outside the
            # admissible region; do not send those to the VLE model.
            if not feasible(p,p0,policy,tol=1e-7):
                return 1000. + 1000.*(abs(p.sum()-1.) + max(0.,0.5*np.abs(p-p0).sum()-limit))
            key=tuple(np.round(p,12))
            if key not in cache:
                cache[key]=-score(p)['KR_share_pct']
            return cache[key]
        t0=time.perf_counter()
        res=minimize(objective,start,method='SLSQP',bounds=list(zip(lower,upper)),
                     constraints=constraints,options={'maxiter':maxiter,'ftol':1e-9,'disp':False})
        valid=feasible(res.x,p0,policy,tol=5e-6)
        best=None
        if valid:
            best=score(res.x)
            candidates.append((res.x.copy(),best))
        log.append({'trial':f'SLSQP_start_{k}','objective_KR_share_pct':best['KR_share_pct'] if best else np.nan,
                    'feasible':bool(valid),'solver_success':bool(res.success),'evaluations':int(res.nfev),
                    'solver_message':str(res.message)[:250],'elapsed_sec':round(time.perf_counter()-t0,3)})
    p_star,result=max(candidates,key=lambda pair:pair[1]['KR_share_pct'])
    if not feasible(p_star,p0,policy):
        raise AssertionError('Best candidate violates compositional trust region')
    return p0,p_star,initial,result,log


def run(args):
    policy=BoundsPolicy(args.max_reallocation_pp,args.relative_change)
    policy.validate()
    if len(args.doses)==0 or len(args.product_ids)==0 or len(args.target_cas)==0:
        raise ValueError('No scenarios requested')
    doses=list(dict.fromkeys(args.doses))
    if any(not math.isfinite(d) or not (0.0<d<100.0) for d in doses):
        raise ValueError('Invalid trace dose')
    props=core.load_reference(args.reference)
    registry=multi.load_registry(args.registry)
    oel_doc,oel_meta=kr.regulatory_document(args.oel_json)
    msds_doc,msds_meta=kr.regulatory_document(args.msds_json)
    annex6=kr.verify_msds_law(msds_doc)
    oel_lines=kr.annex_lines(kr.find_annex(oel_doc,'0001'))
    products=core.load_product_compositions(args.input,[str(p) for p in args.product_ids])
    targets=list(dict.fromkeys(args.target_cas))
    for c in targets:
        if c not in props:raise ValueError(f'CAS {c} has no MW/SMILES/VP reference; unable to evaluate')
    allcas=sorted({c for b in products.values() for c,_ in b}.union(targets))
    twa={c:kr.parse_twa_row(oel_lines,c)['KR_TWA_ppm'] for c in allcas}
    cutoff_data={c:trace_cutoff(c,registry,annex6,args.exploratory_cutoff_wt_pct) for c in targets}
    if args.optimize_model not in args.check_models:
        models=[args.optimize_model]+list(dict.fromkeys(args.check_models))
    else:
        models=list(dict.fromkeys(args.check_models))
    predictor=core.GammaModels(args.hanna_repo)
    statuses=[];comp=[];summary=[];search=[];model_rows=[]
    for pid,base in products.items():
        p0=normalized_base(base)
        for target in targets:
            cutoff,kind=cutoff_data[target]
            if target in dict(base):
                for dose in doses:
                    statuses.append(dict(product_id=pid,target_CAS=target,dose_wt_pct=dose,
                                        status='SKIP_PRESENT_IN_REPORTED_PRODUCT',error_message='',cutoff_type=kind))
                continue
            for dose in doses:
                tag=dict(product_id=pid,product_name=core.LABELS.get(pid,pid),target_CAS=target,
                         dose_wt_pct=dose,cutoff_wt_pct=cutoff,cutoff_type=kind,
                         legal_omission_determination='NOT_PERFORMED')
                if not dose < cutoff - 1e-10:
                    statuses.append({**tag,'status':'SKIP_NOT_STRICTLY_BELOW_CUTOFF','error_message':''})
                    continue
                try:
                    p0,pstar,baseline,best,trials=optimize_one(
                        base,target,dose,args.optimize_model,predictor,props,twa,policy,
                        starts=args.starts,maxiter=args.maxiter,seed=args.seed)
                    statuses.append({**tag,'status':'OK','error_message':''})
                    for j,(cas,orig) in enumerate(base):
                        comp.append({**tag,'CAS':cas,'EPA_original_wt_pct':orig,
                                     'baseline_trace_wt_pct':float(p0[j]*(100.-dose)),
                                     'optimized_trace_wt_pct':float(pstar[j]*(100.-dose)),
                                     'change_wt_percentage_points':float((pstar[j]-p0[j])*(100.-dose))})
                    comp.append({**tag,'CAS':target,'EPA_original_wt_pct':0.,
                                 'baseline_trace_wt_pct':dose,'optimized_trace_wt_pct':dose,
                                 'change_wt_percentage_points':0.})
                    for row in trials:search.append({**tag,'optimize_model':args.optimize_model,**row})
                    summary.append({**tag,'optimize_model':args.optimize_model,
                                    'baseline_KR_share_pct':baseline['KR_share_pct'],
                                    'optimized_KR_share_pct':best['KR_share_pct'],
                                    'gain_percentage_points':best['KR_share_pct']-baseline['KR_share_pct'],
                                    'baseline_vapor_mol_pct':baseline['vapor_mol_pct'],
                                    'optimized_vapor_mol_pct':best['vapor_mol_pct'],
                                    'baseline_vapor_wt_pct':baseline['vapor_wt_pct'],
                                    'optimized_vapor_wt_pct':best['vapor_wt_pct'],
                                    'reallocated_main_components_pp':float(0.5*np.abs(pstar-p0).sum()*100.),
                                    'max_reallocation_pp':policy.max_reallocation_pp,
                                    'relative_change_limit':policy.relative_change,
                                    'global_optimum_proven':False,
                                    'EPA_target_measured':False})
                    for stage,p in [('baseline',p0),('optimized',pstar)]:
                        components=mixture(base,target,dose,p)
                        for model in models:
                            try:
                                r=evaluate(components,target,model,predictor,props,twa)
                                model_rows.append({**tag,'stage':stage,'model_id':model,'status':'OK',
                                                   'error_message':'',**r})
                            except Exception as e:
                                model_rows.append({**tag,'stage':stage,'model_id':model,
                                                   'status':multi.classify_error(e),
                                                   'error_message':f'{type(e).__name__}: {str(e)[:240]}'})
                except Exception as e:
                    statuses.append({**tag,'status':multi.classify_error(e),
                                     'error_message':f'{type(e).__name__}: {str(e)[:400]}'})
    out=Path(args.output_dir)
    out.mkdir(parents=True,exist_ok=True)
    # Always export status even for failed model, avoiding silent fallbacks.
    pd.DataFrame(statuses).to_csv(out/'01_scenario_status.csv',index=False,encoding='utf-8-sig')
    if summary:
        pd.DataFrame(summary).to_csv(out/'02_optimization_summary.csv',index=False,encoding='utf-8-sig')
        pd.DataFrame(comp).to_csv(out/'03_baseline_vs_optimized_compositions.csv',index=False,encoding='utf-8-sig')
        pd.DataFrame(model_rows).to_csv(out/'04_cross_model_checks.csv',index=False,encoding='utf-8-sig')
        pd.DataFrame(search).to_csv(out/'05_search_audit.csv',index=False,encoding='utf-8-sig')
    metadata={'analysis_type':'counterfactual_composition_constrained_inverse_design',
              'temperature_C':25.,'vapor_pressure_validity':'25C reference table only',
              'optimizer':'SLSQP deterministic multistart; best feasible candidate, not global proof',
              'optimize_model':args.optimize_model,'cross_check_models':models,
              'max_reallocation_pp':policy.max_reallocation_pp,
              'relative_component_change_bound':policy.relative_change,
              'optimizer_starts':args.starts,'optimizer_maxiter':args.maxiter,'seed':args.seed,
              'KR_OEL_notice':oel_meta.get('발령번호'),
              'KR_OEL_effective':oel_meta.get('시행일자'),
              'MSDS_notice':msds_meta.get('발령번호'),
              'MSDS_effective':msds_meta.get('시행일자'),
              'legal_omission_assessed':False,'actual_worker_exposure_assessed':False,
              'experimentally_validated':False,'model_validation_required':True,
              'EPA_target_measured':False,'endpoints_toxicologically_additive':False,
              'warning':'Optimizing a predicted weighted share is a SCREENING exercise, not evidence of an illegal SDS or health risk; all trace dosing is counterfactual.',
              'successful_scenarios':sum(z['status']=='OK' for z in statuses),
              'failed_scenarios':sum(z['status'] in ('FAILED','UNAVAILABLE') for z in statuses),
              'skipped_scenarios':sum(z['status'].startswith('SKIP') for z in statuses),
              'cross_model_failed':sum(z['status']!='OK' for z in model_rows)}
    (out/'06_run_metadata.json').write_text(json.dumps(metadata,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps({k:metadata[k] for k in ('successful_scenarios','failed_scenarios','skipped_scenarios','cross_model_failed')},ensure_ascii=False))
    print('Output:',out.resolve())
    if (metadata['failed_scenarios'] or
            (args.strict_models and metadata['cross_model_failed'])):
        raise RuntimeError('Optimization or cross-validation failed; inspect 01 and 04 CSV for details')
    return metadata


if __name__=='__main__':
    try:
        run(parse_args())
    except Exception as exc:
        print(f'[ERROR] {type(exc).__name__}: {exc}',file=sys.stderr)
        sys.exit(1)
