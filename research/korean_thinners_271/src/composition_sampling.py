# -*- coding: utf-8 -*-
"""Feasible MSDS interval compositions via a bounded-simplex pairwise Gibbs chain.

At each step a random pair (i,j) is uniformly redistributed conditional on the
other coordinates. Uniform distribution on the feasible bounded simplex is an
invariant distribution, but convergence / mixing cannot be assumed or certified.
Use independent seeds, compare against the legacy conditional-rejection sampler,
and never interpret sampled percentiles as confidence intervals of a true formula.
"""
from __future__ import annotations

import math
import random


def check_bounds(comps):
    if len(comps) < 2:
        raise ValueError('At least two components required')
    low, high = [], []
    for r in comps:
        l, u = float(r['함량 하한(%)']), float(r['함량 상한(%)'])
        if not (math.isfinite(l) and math.isfinite(u) and 0 <= l <= u <= 100):
            raise ValueError('Invalid MSDS interval')
        low.append(l)
        high.append(u)
    if sum(low) > 100 + 1e-8 or sum(high) < 100 - 1e-8:
        raise ValueError('MSDS intervals do not admit 100 percent total')
    return low, high


def starting_point(low, high):
    """Construct a deterministic feasible point, including degenerate polytopes."""
    w = low.copy()
    remaining = 100 - sum(w)
    for j in range(len(w)):
        inc = min(remaining, high[j] - w[j])
        w[j] += inc
        remaining -= inc
    if abs(remaining) > 1e-7:
        raise ValueError('Unable to construct feasible bounded simplex')
    return w


def sample_bounded_simplex(comps, count, rng: random.Random, *, burn_in=300, thin=12):
    """Sample by pairwise conditional-uniform updates; preserves sum and bounds.

    Parameters are *Gibbs steps*, not independent samples. Report diagnostics
    from seed/sampler sensitivity; correlated draws have no independent-sample CI.
    """
    low, high = check_bounds(comps)
    if count < 1 or burn_in < 0 or thin < 1:
        raise ValueError('count positive, burn_in nonnegative, thin positive required')
    n = len(low)
    x = starting_point(low, high)
    out = []
    steps = burn_in + count * thin
    for step in range(steps):
        a, b = rng.sample(range(n), 2)
        t = x[a] + x[b]
        left = max(low[a], t - high[b])
        right = min(high[a], t - low[b])
        if right < left - 1e-8:
            raise RuntimeError('Invalid pairwise Gibbs constraint intersection')
        new_a = rng.uniform(left, max(left, right))
        x[a] = new_a
        x[b] = t - new_a
        if step >= burn_in and (step - burn_in) % thin == thin - 1:
            candidate = x.copy()
            if (not math.isclose(sum(candidate), 100, abs_tol=1e-6, rel_tol=0)
                or any(v < l-1e-6 or v > u+1e-6 for v,l,u in zip(candidate, low, high))):
                raise RuntimeError('Gibbs sampler produced infeasible composition')
            out.append(candidate)
    if len(out) != count:
        raise RuntimeError('Sampler draw count mismatch')
    return out