# -*- coding: utf-8 -*-
"""2026 multicomponent HANNA -> thinner-271 optimization adapter.

Install this file at:
    D:\\fire2\\research\\korean_thinners_271\\src\\hanna_gamma_adapter.py
with HANNA cloned at:
    D:\\fire2\\HANNA

The optimizer expects **gamma**, whereas HANNA predicts **ln(gamma)**.
No Raoult fallback: missing checkpoints or invalid prediction raise errors.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Sequence

import numpy as np

_MODEL = None


def _hanna_dir() -> Path:
    """Find sibling HANNA checkout relative to fire repository root."""
    # Adapter lives in fire/research/korean_thinners_271/src/.
    repo_root = Path(__file__).resolve().parents[3]
    return Path(os.environ.get("HANNA_HOME", str(repo_root / "HANNA"))).resolve()


def _get_model():
    global _MODEL
    if _MODEL is not None:
        return _MODEL

    hanna_dir = _hanna_dir()
    expected = hanna_dir / "utils" / "HANNA_predictor.py"
    if not expected.is_file():
        raise FileNotFoundError(
            f"2026 HANNA 소스 파일을 찾지 못했습니다: {expected}. "
            "HANNA_HOME 환경 변수 또는 D:\\fire2\\HANNA 경로를 확인하세요."
        )

    # Upstream predictor loads checkpoints/scalers via relative paths;
    # initialize inside HANNA root and restore the user's working directory.
    hanna_path = str(hanna_dir)
    if hanna_path not in sys.path:
        sys.path.insert(0, hanna_path)

    previous_cwd = Path.cwd()
    try:
        os.chdir(hanna_dir)
        from utils.HANNA_predictor import HANNA_Predictor
        _MODEL = HANNA_Predictor()
    finally:
        os.chdir(previous_cwd)
    return _MODEL


def predict_gammas(smiles: Sequence[str], x: Sequence[float], temperature: float) -> list[float]:
    """Return gamma_i for the existing optimizer's interface.

    Args:
        smiles: SMILES for each verified pure component, ordered as x.
        x: *molar* fractions, not SDS mass percentages; must sum to 1.
        temperature: Kelvin.

    Raises instead of silently falling back to Raoult if HANNA fails.
    """
    smiles = list(smiles)
    fractions = np.asarray(x, dtype=float)
    temperature = float(temperature)
    if len(smiles) < 2 or any(not isinstance(s, str) or not s.strip() for s in smiles):
        raise ValueError("2개 이상의 검증된 SMILES가 필요합니다.")
    if fractions.ndim != 1 or len(fractions) != len(smiles):
        raise ValueError("성분 수와 1차원 몰분율 배열 길이가 일치해야 합니다.")
    if not np.all(np.isfinite(fractions)) or np.any(fractions < 0):
        raise ValueError("몰분율에는 유한한 0 이상의 숫자만 사용할 수 있습니다.")
    if not np.isclose(fractions.sum(), 1.0, atol=1e-8, rtol=0):
        raise ValueError("몰분율 합계가 1이어야 합니다.")
    if not np.isfinite(temperature) or temperature <= 0:
        raise ValueError("온도는 양의 켈빈 값이어야 합니다.")

    ln_gamma, _gE_over_RT = _get_model().predict(
        smiles, [fractions.tolist()], temperature
    )
    ln_gamma = np.asarray(ln_gamma, dtype=float)
    if ln_gamma.shape != (1, len(smiles)) or not np.all(np.isfinite(ln_gamma)):
        raise ValueError("HANNA가 유효한 ln(gamma) 배열을 반환하지 않았습니다.")
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        gammas = np.exp(ln_gamma[0])
    if not np.all(np.isfinite(gammas)) or np.any(gammas <= 0):
        raise ValueError("HANNA의 활동도계수 변환 결과가 유효하지 않습니다.")
    return gammas.tolist()