#!/usr/bin/env python3
"""
MAAT-SO(10) Yukawa Fit v0.4
Fuller coupled 1-loop Yukawa RG prototype with QCD-enhanced bottom running.

Boundary:
    y_b(MGUT) = y_tau(MGUT) = y_btau

Fit parameters:
    y_t(MGUT)
    y_btau(MGUT)
    log10(MGUT/GeV)
    Delta_b

Low-energy threshold:
    y_b_eff(MZ) = y_b_RG(MZ) * (1 + Delta_b)

Status:
    Still a prototype:
    - one-loop gauge background
    - one-loop third-generation Yukawa RG
    - no SUSY thresholds
    - no two-loop running
"""

import json
import os
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution, minimize

OUTPUT_DIR = Path("outputs")
OUTPUT_DIR.mkdir(exist_ok=True)

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

v = 246.22
MZ = 91.1876

m_obs = {
    "top": 172.76,
    "bottom": 4.18,
    "tau": 1.77686,
}

m_sigma = {
    "top": 1.0,
    "bottom": 0.15,
    "tau": 0.002,
}

y_obs = {k: np.sqrt(2.0) * m_obs[k] / v for k in m_obs}
y_sigma = {k: np.sqrt(2.0) * m_sigma[k] / v for k in m_sigma}

alpha_gut_ref = 0.022676468338
log10_mgut_ref = 15.822115615160

# SM one-loop beta coefficients, GUT-normalized
b_gauge = {
    "alpha1": 41 / 10,
    "alpha2": -19 / 6,
    "alpha3": -7,
}

sigma_delta_b = 0.35


# ------------------------------------------------------------
# Gauge running
# ------------------------------------------------------------

def alpha_at_scale(alpha_gut, log10_mgut, b_i, mu):
    """
    One-loop:
        1/alpha(mu) = 1/alpha(MGUT) + b/(2pi) ln(MGUT/mu)
    """
    mgut = 10.0 ** log10_mgut

    if alpha_gut <= 0 or mgut <= 0 or mu <= 0:
        return np.nan

    inv_alpha = (1.0 / alpha_gut) + (b_i / (2.0 * np.pi)) * np.log(mgut / mu)

    if inv_alpha <= 0:
        return np.nan

    return 1.0 / inv_alpha


def gauge_couplings_at_t(t, log10_mgut):
    """
    t = ln(mu)
    Returns g1,g2,g3 at scale mu = exp(t).
    """
    mu = np.exp(t)

    a1 = alpha_at_scale(alpha_gut_ref, log10_mgut, b_gauge["alpha1"], mu)
    a2 = alpha_at_scale(alpha_gut_ref, log10_mgut, b_gauge["alpha2"], mu)
    a3 = alpha_at_scale(alpha_gut_ref, log10_mgut, b_gauge["alpha3"], mu)

    if not np.all(np.isfinite([a1, a2, a3])):
        return None
    if a1 <= 0 or a2 <= 0 or a3 <= 0:
        return None

    return (
        np.sqrt(4.0 * np.pi * a1),
        np.sqrt(4.0 * np.pi * a2),
        np.sqrt(4.0 * np.pi * a3),
    )


# ------------------------------------------------------------
# Coupled one-loop third-generation Yukawa RG
# ------------------------------------------------------------

def yukawa_beta(t, y, log10_mgut):
    """
    y = [yt, yb, ytau]
    Running downward is handled by solve_ivp from t_GUT to t_Z.
    """

    yt, yb, ytau = y

    if yt <= 0 or yb <= 0 or ytau <= 0:
        return [0.0, 0.0, 0.0]

    gs = gauge_couplings_at_t(t, log10_mgut)
    if gs is None:
        return [0.0, 0.0, 0.0]

    g1, g2, g3 = gs
    loop = 16.0 * np.pi**2

    beta_t = yt * (
        (9.0 / 2.0) * yt**2
        + (3.0 / 2.0) * yb**2
        - ((17.0 / 20.0) * g1**2 + (9.0 / 4.0) * g2**2 + 8.0 * g3**2)
    ) / loop

    beta_b = yb * (
        (3.0 / 2.0) * yt**2
        + (9.0 / 2.0) * yb**2
        + ytau**2
        - ((1.0 / 4.0) * g1**2 + (9.0 / 4.0) * g2**2 + 8.0 * g3**2)
    ) / loop

    beta_tau = ytau * (
        3.0 * yb**2
        + (5.0 / 2.0) * ytau**2
        - ((9.0 / 4.0) * g1**2 + (9.0 / 4.0) * g2**2)
    ) / loop

    return [beta_t, beta_b, beta_tau]


def run_yukawas_down(params):
    """
    Runs from MGUT to MZ.

    Because t decreases from high to low scale, solve_ivp integrates
    from ln(MGUT) to ln(MZ). This automatically applies the sign correctly.
    """
    yt_gut, ybtau_gut, log10_mgut, delta_b = params

    if yt_gut <= 0 or ybtau_gut <= 0:
        return None

    mgut = 10.0 ** log10_mgut
    if mgut <= MZ:
        return None

    t_gut = np.log(mgut)
    t_z = np.log(MZ)

    y0 = [yt_gut, ybtau_gut, ybtau_gut]

    sol = solve_ivp(
        fun=lambda t, y: yukawa_beta(t, y, log10_mgut),
        t_span=(t_gut, t_z),
        y0=y0,
        method="RK45",
        rtol=1e-7,
        atol=1e-10,
        max_step=0.25,
    )

    if not sol.success:
        return None

    yt_z, yb_z_rg, ytau_z = sol.y[:, -1]

    if not np.all(np.isfinite([yt_z, yb_z_rg, ytau_z])):
        return None
    if yt_z <= 0 or yb_z_rg <= 0 or ytau_z <= 0:
        return None

    yb_z_eff = yb_z_rg * (1.0 + delta_b)

    if yb_z_eff <= 0:
        return None

    return {
        "top": yt_z,
        "bottom_rg": yb_z_rg,
        "bottom": yb_z_eff,
        "tau": ytau_z,
    }


# ------------------------------------------------------------
# MAAT supports
# ------------------------------------------------------------

def bad_supports():
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


def maat_supports(params):
    pred = run_yukawas_down(params)

    if pred is None:
        return bad_supports()

    yt = pred["top"]
    yb = pred["bottom"]
    ytau = pred["tau"]

    p = np.array([yt, yb, ytau])
    obs = np.array([y_obs["top"], y_obs["bottom"], y_obs["tau"]])

    if np.any(p <= 0) or np.any(~np.isfinite(p)):
        return bad_supports()

    # H: low-energy fit consistency
    d_H = np.mean(((p - obs) / obs) ** 2)

    # B: effective low-energy b/tau ratio consistency
    ratio_pred = yb / ytau
    ratio_obs = y_obs["bottom"] / y_obs["tau"]
    d_B = (np.log(ratio_pred / ratio_obs)) ** 2

    # S: hierarchy consistency
    hierarchy = np.log(yt / np.sqrt(yb * ytau))
    hierarchy_target = np.log(y_obs["top"] / np.sqrt(y_obs["bottom"] * y_obs["tau"]))
    d_S = (hierarchy - hierarchy_target) ** 2

    # V: connectedness before threshold.
    # This measures how far RG-only b/tau relation is from observed b/tau.
    yb_rg = pred["bottom_rg"]
    d_V = (np.log((yb_rg / ytau) / ratio_obs)) ** 2

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


# ------------------------------------------------------------
# Cost functions
# ------------------------------------------------------------

def chi2_yukawa(params):
    pred = run_yukawas_down(params)

    if pred is None:
        return 1e12

    cost = 0.0
    for k in ["top", "bottom", "tau"]:
        if pred[k] <= 0 or not np.isfinite(pred[k]):
            return 1e12
        cost += ((pred[k] - y_obs[k]) / y_sigma[k]) ** 2

    return cost


def maat_cost(params):
    eps = 1e-12

    lam = {
        "H": 1.0,
        "B": 1.0,
        "S": 1.0,
        "V": 1.0,
        "cl": 2.0,
    }

    s = maat_supports(params)

    F = 0.0
    for key in ["H", "B", "S", "V"]:
        F += -lam[key] * np.log(eps + s[key])

    F += -lam["cl"] * np.log(eps + s["R_rob"])

    return F


def threshold_penalty(params):
    _, _, _, delta_b = params
    return (delta_b / sigma_delta_b) ** 2


def complexity_penalty(params):
    yt_gut, ybtau_gut, log10_mgut, delta_b = params

    scale_penalty = ((log10_mgut - log10_mgut_ref) / 2.0) ** 2
    yuk_penalty = 0.01 * (np.log(yt_gut) ** 2 + np.log(ybtau_gut) ** 2)

    return scale_penalty + yuk_penalty


def total_cost(params):
    yt_gut, ybtau_gut, log10_mgut, delta_b = params

    if not (0.2 <= yt_gut <= 3.0):
        return 1e12
    if not (0.002 <= ybtau_gut <= 0.2):
        return 1e12
    if not (12.0 <= log10_mgut <= 18.5):
        return 1e12
    if not (-0.8 <= delta_b <= 3.0):
        return 1e12
    if 1.0 + delta_b <= 0:
        return 1e12

    return (
        chi2_yukawa(params)
        + maat_cost(params)
        + threshold_penalty(params)
        + complexity_penalty(params)
    )


# ------------------------------------------------------------
# Fit
# ------------------------------------------------------------

bounds = [
    (0.2, 3.0),       # y_t(MGUT)
    (0.002, 0.2),     # y_btau(MGUT)
    (12.0, 18.5),     # log10(MGUT/GeV)
    (-0.8, 3.0),      # Delta_b
]

global_result = differential_evolution(
    total_cost,
    bounds=bounds,
    seed=42,
    tol=1e-7,
    polish=False,
    workers=1,
    updating="immediate",
    maxiter=300,
    popsize=18,
)

local_result = minimize(
    total_cost,
    global_result.x,
    method="L-BFGS-B",
    bounds=bounds,
    options={
        "maxiter": 50000,
        "ftol": 1e-12,
        "gtol": 1e-10,
    },
)

params_best = local_result.x
yt_gut, ybtau_gut, log10_mgut, delta_b = params_best
mgut = 10.0 ** log10_mgut

pred = run_yukawas_down(params_best)
s = maat_supports(params_best)

print("\n=== MAAT-SO(10) Yukawa Fit v0.4 Coupled 1-loop RG ===")
print(f"success global       : {global_result.success}")
print(f"success local        : {local_result.success}")
print(f"message local        : {local_result.message}")
print(f"total cost           : {total_cost(params_best):.8f}")
print(f"chi2_yukawa          : {chi2_yukawa(params_best):.8f}")
print(f"MAAT cost            : {maat_cost(params_best):.8f}")
print(f"threshold penalty    : {threshold_penalty(params_best):.8f}")
print(f"complexity penalty   : {complexity_penalty(params_best):.8f}")
print(f"y_t(MGUT)            : {yt_gut:.10f}")
print(f"y_btau(MGUT)         : {ybtau_gut:.10f}")
print(f"log10(MGUT/GeV)      : {log10_mgut:.10f}")
print(f"MGUT [GeV]           : {mgut:.6e}")
print(f"Delta_b              : {delta_b:.10f}")
print(f"bottom boost 1+Delta : {1.0 + delta_b:.10f}")

print("\n--- Predicted low-energy Yukawas and masses ---")
for k in ["top", "bottom", "tau"]:
    y = pred[k]
    mass = y * v / np.sqrt(2.0)
    pull = (y - y_obs[k]) / y_sigma[k]
    print(
        f"{k:7s}: y_pred={y:.10f} | y_obs={y_obs[k]:.10f} | "
        f"pull={pull:+.4f} | m_pred={mass:.6f} GeV"
    )

print("\n--- Bottom threshold diagnostic ---")
print(f"bottom RG before threshold   : {pred['bottom_rg']:.10f}")
print(f"bottom after threshold       : {pred['bottom']:.10f}")
print(f"required multiplicative boost: {(pred['bottom'] / pred['bottom_rg']):.6f}")

print("\n--- MAAT supports ---")
for key, value in s.items():
    print(f"{key:8s}: {value:.8f}")

print("\n--- SO(10) unification diagnostics ---")
print(f"b/tau at MGUT        : {1.0:.6f}  (imposed)")
print(f"b/tau at MZ RG-only  : {pred['bottom_rg'] / pred['tau']:.6f}")
print(f"b/tau at MZ effective: {pred['bottom'] / pred['tau']:.6f}")
print(f"b/tau at MZ observed : {y_obs['bottom'] / y_obs['tau']:.6f}")
print(f"top hierarchy pred   : {pred['top'] / np.sqrt(pred['bottom'] * pred['tau']):.6f}")

# ============================================================
# Outputs and figure
# ============================================================

summary = {
    "success_global": bool(global_result.success),
    "success_local": bool(local_result.success),
    "message_local": str(local_result.message),
    "total_cost": float(total_cost(params_best)),
    "chi2_yukawa": float(chi2_yukawa(params_best)),
    "MAAT_cost": float(maat_cost(params_best)),
    "threshold_penalty": float(threshold_penalty(params_best)),
    "complexity_penalty": float(complexity_penalty(params_best)),
    "y_t_MGUT": float(yt_gut),
    "y_btau_MGUT": float(ybtau_gut),
    "log10_MGUT_GeV": float(log10_mgut),
    "MGUT_GeV": float(mgut),
    "Delta_b": float(delta_b),
    "bottom_boost_1_plus_Delta": float(1.0 + delta_b),
    "predicted_yukawas": {key: float(pred[key]) for key in ["top", "bottom", "tau"]},
    "predicted_masses_GeV": {
        key: float(pred[key] * v / np.sqrt(2.0)) for key in ["top", "bottom", "tau"]
    },
    "observed_yukawas": {key: float(value) for key, value in y_obs.items()},
    "supports": {key: float(value) for key, value in s.items()},
    "b_tau_MZ_RG_only": float(pred["bottom_rg"] / pred["tau"]),
    "b_tau_MZ_effective": float(pred["bottom"] / pred["tau"]),
    "b_tau_MZ_observed": float(y_obs["bottom"] / y_obs["tau"]),
    "status": "SO(10)-motivated b-tau boundary with one-loop SM running and fitted bottom-threshold correction",
}

with (OUTPUT_DIR / "yukawa_v04_summary.json").open("w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)

print(f"\nWrote {OUTPUT_DIR / 'yukawa_v04_summary.json'}")

os.environ.setdefault("MPLCONFIGDIR", "/tmp/codex-mpl-cache")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

log10_mgut_vals = np.linspace(14, 18, 80)
ybtau_vals = np.linspace(0.002, 0.02, 80)

Z = np.zeros((len(ybtau_vals), len(log10_mgut_vals)))

# best-fit Werte aus deinem Run:
yt_gut_fixed = yt_gut
delta_b_fixed = delta_b

for i, ybtau in enumerate(ybtau_vals):
    for j, log10_mgut in enumerate(log10_mgut_vals):
        params = [yt_gut_fixed, ybtau, log10_mgut, delta_b_fixed]
        val = total_cost(params)
        Z[i, j] = np.log10(val + 1e-6)

plt.figure(figsize=(7, 5))

im = plt.imshow(
    Z,
    origin="lower",
    aspect="auto",
    extent=[
        log10_mgut_vals.min(),
        log10_mgut_vals.max(),
        ybtau_vals.min(),
        ybtau_vals.max(),
    ],
)

plt.colorbar(im, label=r"$\log_{10}(\mathcal{L}_{\mathrm{total}})$")

best_log10_mgut = 16.269
best_ybtau = 0.00774

plt.scatter(
    best_log10_mgut,
    best_ybtau,
    color="red",
    s=70,
    edgecolor="black",
    label="Best fit"
)

plt.xlabel(r"$\log_{10}(M_{\mathrm{GUT}} / \mathrm{GeV})$")
plt.ylabel(r"$y_{b\tau}(M_{\mathrm{GUT}})$")

plt.title("MAAT–SO(10) Selection Landscape")

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "maat_selection_landscape.png", dpi=300)
print(f"Wrote {OUTPUT_DIR / 'maat_selection_landscape.png'}")
