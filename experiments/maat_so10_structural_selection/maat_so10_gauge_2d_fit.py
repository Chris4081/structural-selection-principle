#!/usr/bin/env python3
"""
MAAT-SO(10) 2D parameter fit:
fit alpha_GUT and M_GUT simultaneously.

Proof-of-concept only:
- one-loop SM running
- no threshold corrections
- no two-loop RG
- no SUSY/non-SUSY splitting
"""

import numpy as np
import json
from pathlib import Path
from scipy.optimize import differential_evolution, minimize

OUTPUT_DIR = Path("outputs")
OUTPUT_DIR.mkdir(exist_ok=True)

# ----------------------------
# Experimental low-energy input
# ----------------------------

MZ = 91.1876

alpha_exp = {
    "alpha1": 0.01695,
    "alpha2": 0.03382,
    "alpha3": 0.11840,
}

sigma = {
    "alpha1": 0.00020,
    "alpha2": 0.00030,
    "alpha3": 0.00100,
}

# One-loop SM beta coefficients, GUT-normalized
b = {
    "alpha1": 41 / 10,
    "alpha2": -19 / 6,
    "alpha3": -7,
}


# ----------------------------
# RG running
# ----------------------------

def run_alpha_down(alpha_gut, log10_mgut, b_i):
    mgut = 10.0 ** log10_mgut

    if alpha_gut <= 0 or mgut <= MZ:
        return np.nan

    inv_alpha_low = (1.0 / alpha_gut) + (b_i / (2 * np.pi)) * np.log(mgut / MZ)

    if inv_alpha_low <= 0:
        return np.nan

    return 1.0 / inv_alpha_low


# ----------------------------
# SM chi-square
# ----------------------------

def sm_fit_cost(alpha_gut, log10_mgut):
    cost = 0.0

    for key in alpha_exp:
        pred = run_alpha_down(alpha_gut, log10_mgut, b[key])

        if not np.isfinite(pred) or pred <= 0:
            return 1e12

        cost += ((pred - alpha_exp[key]) / sigma[key]) ** 2

    return cost


# ----------------------------
# MAAT supports
# ----------------------------

def maat_supports(alpha_gut, log10_mgut):
    preds = np.array([
        run_alpha_down(alpha_gut, log10_mgut, b["alpha1"]),
        run_alpha_down(alpha_gut, log10_mgut, b["alpha2"]),
        run_alpha_down(alpha_gut, log10_mgut, b["alpha3"]),
    ])

    obs = np.array([
        alpha_exp["alpha1"],
        alpha_exp["alpha2"],
        alpha_exp["alpha3"],
    ])

    if np.any(~np.isfinite(preds)) or np.any(preds <= 0):
        return {
            "H": 1e-12,
            "B": 1e-12,
            "S": 1e-12,
            "V": 1e-12,
            "R_resp": 1e-12,
            "R_rob": 1e-12,
            "d_H": 1e12,
            "d_B": 1e12,
            "d_S": 1e12,
            "d_V": 1e12,
        }

    # H: closeness to observed low-energy couplings
    d_H = np.mean(((preds - obs) / obs) ** 2)

    # B: balance of predicted low-energy sectors
    ratio = preds / np.mean(preds)
    d_B = np.var(ratio)

    # S: controlled running activity
    flow_strength = np.mean(np.abs(preds - alpha_gut) / alpha_gut)
    target_activity = 1.0
    d_S = (flow_strength - target_activity) ** 2

    # V: connectedness / sector convergence
    d_V = (np.max(preds) - np.min(preds)) / np.mean(preds)

    H = 1.0 / (1.0 + d_H)
    B = 1.0 / (1.0 + d_B)
    S = 1.0 / (1.0 + d_S)
    V = 1.0 / (1.0 + d_V)

    R_resp = (H * B * V) ** (1.0 / 3.0)
    R_rob = min(R_resp, (H * B * S * V) ** 0.25)

    return {
        "H": H,
        "B": B,
        "S": S,
        "V": V,
        "R_resp": R_resp,
        "R_rob": R_rob,
        "d_H": d_H,
        "d_B": d_B,
        "d_S": d_S,
        "d_V": d_V,
    }


def maat_cost(alpha_gut, log10_mgut):
    eps = 1e-12

    lam = {
        "H": 1.0,
        "B": 1.0,
        "S": 1.0,
        "V": 1.0,
        "cl": 2.0,
    }

    s = maat_supports(alpha_gut, log10_mgut)

    F = 0.0
    for key in ["H", "B", "S", "V"]:
        F += -lam[key] * np.log(eps + s[key])

    F += -lam["cl"] * np.log(eps + s["R_rob"])

    return F


# ----------------------------
# Total objective
# ----------------------------

def total_cost_vec(x):
    alpha_gut, log10_mgut = x

    chi2 = sm_fit_cost(alpha_gut, log10_mgut)
    fmaat = maat_cost(alpha_gut, log10_mgut)

    # Weak prior: plausible SO(10)/GUT scale window
    # Keeps runaway low-scale fake unification suppressed.
    mgut_prior = ((log10_mgut - 16.0) / 2.0) ** 2

    return chi2 + fmaat + mgut_prior


# ----------------------------
# Fit
# ----------------------------

bounds = [
    (0.005, 0.08),   # alpha_GUT
    (12.0, 18.5),    # log10(M_GUT / GeV)
]

global_result = differential_evolution(
    total_cost_vec,
    bounds=bounds,
    tol=1e-10,
    polish=False,
    seed=42,
)

local_result = minimize(
    total_cost_vec,
    global_result.x,
    method="Nelder-Mead",
    options={
        "xatol": 1e-13,
        "fatol": 1e-9,
        "maxiter": 20000,
    },
)

alpha_best, log10_mgut_best = local_result.x
mgut_best = 10.0 ** log10_mgut_best
g_best = np.sqrt(4.0 * np.pi * alpha_best)

print("\n=== MAAT-SO(10) 2D Parameter Fit ===")
print(f"success          : {local_result.success}")
print(f"alpha_GUT        : {alpha_best:.12f}")
print(f"g_GUT            : {g_best:.12f}")
print(f"log10(M_GUT/GeV) : {log10_mgut_best:.12f}")
print(f"M_GUT [GeV]      : {mgut_best:.6e}")
print(f"total cost       : {total_cost_vec(local_result.x):.6f}")
print(f"SM chi2          : {sm_fit_cost(alpha_best, log10_mgut_best):.6f}")
print(f"MAAT cost        : {maat_cost(alpha_best, log10_mgut_best):.6f}")

print("\n--- Predicted low-energy couplings ---")
for key in alpha_exp:
    pred = run_alpha_down(alpha_best, log10_mgut_best, b[key])
    pull = (pred - alpha_exp[key]) / sigma[key]
    print(
        f"{key}: pred={pred:.8f} | obs={alpha_exp[key]:.8f} | pull={pull:+.3f}"
    )

print("\n--- MAAT supports ---")
s = maat_supports(alpha_best, log10_mgut_best)
for key, value in s.items():
    print(f"{key:8s}: {value:.8f}")

summary = {
    "success": bool(local_result.success),
    "alpha_GUT": float(alpha_best),
    "g_GUT": float(g_best),
    "log10_M_GUT_GeV": float(log10_mgut_best),
    "M_GUT_GeV": float(mgut_best),
    "total_cost": float(total_cost_vec(local_result.x)),
    "SM_chi2": float(sm_fit_cost(alpha_best, log10_mgut_best)),
    "MAAT_cost": float(maat_cost(alpha_best, log10_mgut_best)),
    "predicted_couplings": {
        key: float(run_alpha_down(alpha_best, log10_mgut_best, b[key]))
        for key in alpha_exp
    },
    "observed_couplings": {key: float(value) for key, value in alpha_exp.items()},
    "supports": {key: float(value) for key, value in s.items()},
    "status": "one-loop diagnostic; no threshold corrections; no precision unification claim",
}

with (OUTPUT_DIR / "gauge_2d_summary.json").open("w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)

print(f"\nWrote {OUTPUT_DIR / 'gauge_2d_summary.json'}")
