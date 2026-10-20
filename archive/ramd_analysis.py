"""
τRAMD residence-time analysis (Kokh et al., J. Chem. Theory Comput. 2018; tauRAMD.py).

Input: dissociation (stop) times per replica, from any engine. For each replica:
- bootstrap: `rounds` random subsets of `fraction` of the trajectories (drawn
  without replacement, as tauRAMD.py does); the median of each subset is the
  time at which 50% of trajectories have dissociated;
- τ_replica = mean of those medians, SD = their standard deviation;
- KS = max |ECDF(t) − (1 − exp(−t/τ_replica))|, the agreement with a Poisson
  (single-exponential) escape process. tauRAMD.py uses the bootstrap τ as the
  exponential time constant and a binned CDF; here the exact ECDF is used.
  0.1–0.2 is typical.
Over replicas: τ = mean of τ_replica, SD = population SD (ddof=0), as tauRAMD.py.

Trajectories that did not dissociate are kept at their stop time (max_time) as
censored values; they are the longest, so the median is unaffected while fewer
than half are censored. With half or more censored, τ_replica is only a lower
bound and is flagged.
"""

from __future__ import annotations

from typing import Dict, Any, Sequence, Optional

import numpy as np


def ks_exponential(times: Sequence[float], tau: float) -> float:
    t = np.sort(np.asarray(times, dtype=float))
    n = len(t)
    model = 1.0 - np.exp(-t / tau)
    upper = np.arange(1, n + 1) / n - model
    lower = model - np.arange(0, n) / n
    return float(max(upper.max(), lower.max()))


def replica_tau(times: Sequence[float], censored: Optional[Sequence[bool]] = None, rounds: int = 5000,
                fraction: float = 0.8, seed: Optional[int] = 0) -> Dict[str, Any]:
    t = np.asarray(times, dtype=float)
    cens = np.zeros(len(t), bool) if censored is None else np.asarray(censored, bool)
    n = len(t)
    if n < 2:
        raise ValueError("need at least 2 trajectories per replica")
    rng = np.random.default_rng(seed)
    k = max(1, int(fraction * n))
    medians = np.array([np.median(rng.permutation(t)[:k]) for _ in range(rounds)])
    tau, sd = float(medians.mean()), float(medians.std())
    return {
        "n": n,
        "dissociated": int((~cens).sum()),
        "tau_ps": round(tau, 2),
        "sd_ps": round(sd, 2),
        "median_ps": round(float(np.median(t)), 2),
        "ks": round(ks_exponential(t, tau), 3),
        "lower_bound": bool(cens.sum() * 2 >= n),
    }


def tau_ramd(replicas: Dict[Any, Sequence[float]], censored: Optional[Dict[Any, Sequence[bool]]] = None,
             rounds: int = 5000, fraction: float = 0.8, seed: Optional[int] = 0,
             min_trajectories: int = 10) -> Dict[str, Any]:
    """τ per replica and overall, plus warnings about the protocol's own checks."""
    per = {}
    warnings = []
    for r, times in replicas.items():
        if len(times) < 2:
            warnings.append(f"replica {r}: only {len(times)} trajectories, skipped")
            continue
        per[r] = replica_tau(times, (censored or {}).get(r), rounds, fraction, seed)
        if per[r]["n"] < min_trajectories:
            warnings.append(f"replica {r}: {per[r]['n']} trajectories (< {min_trajectories})")
        if per[r]["lower_bound"]:
            warnings.append(f"replica {r}: half or more trajectories did not dissociate; τ is a lower bound "
                            "(longer max_time or a larger force needed)")
    out: Dict[str, Any] = {"replicas": {str(r): v for r, v in per.items()}}
    if per:
        taus = np.array([v["tau_ps"] for v in per.values()])
        sds = np.array([v["sd_ps"] for v in per.values()])
        out.update({
            "tau_ps": round(float(taus.mean()), 2),
            "tau_sd_ps": round(float(taus.std()), 2),  # spread across replicas (ddof=0, as tauRAMD.py)
            "mean_bootstrap_sd_ps": round(float(sds.mean()), 2),
            "n_replicas": len(per),
        })
        if len(per) < 4:
            warnings.append(f"{len(per)} replicas (< 4)")
        if len(per) > 1 and taus.std() > sds.mean():
            warnings.append("replica-to-replica spread exceeds the bootstrap uncertainty: add replicas")
        if taus.min() < 100.0:
            warnings.append(f"shortest replica τ {taus.min():.0f} ps < 100 ps: the force may be too strong")
        ks = [v["ks"] for v in per.values()]
        if max(ks) > 0.3:
            warnings.append(f"KS up to {max(ks):.2f} (> 0.3): stop times are not single-exponential")
    out["warnings"] = warnings
    return out


def format_report(res: Dict[str, Any]) -> str:
    lines = ["replica  n  dissoc  tau_ps   sd_ps  median_ps  KS"]
    for r, v in res["replicas"].items():
        lines.append(f"{r:>7} {v['n']:>2} {v['dissociated']:>7} {v['tau_ps']:>7.1f} {v['sd_ps']:>7.1f} "
                     f"{v['median_ps']:>10.1f} {v['ks']:>5.2f}" + ("  (lower bound)" if v["lower_bound"] else ""))
    if "tau_ps" in res:
        lines.append(f"τ = {res['tau_ps']:.1f} ± {res['tau_sd_ps']:.1f} ps over {res['n_replicas']} replicas "
                     f"({res['tau_ps'] / 1000:.3f} ns); mean bootstrap SD {res['mean_bootstrap_sd_ps']:.1f} ps")
    lines += [f"[WARNING] {w}" for w in res["warnings"]]
    return "\n".join(lines)
