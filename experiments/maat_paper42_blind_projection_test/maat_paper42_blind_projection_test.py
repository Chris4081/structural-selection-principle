#!/usr/bin/env python3
"""
MAAT Paper 42 — Blind Projection Test
=====================================
Derives C_proj(z) from response weights lambda_a without free projection
parameters.

Pipeline:
    defects -> covariance C -> lambda
    -> gamma_lambda, Bstar_lambda, alpha_lambda
    -> C_proj_lambda(z)
    -> residual / CCI / null tests

Run:
    python3 maat_paper42_blind_projection_test.py
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
import numpy as np

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "maat_paper42_matplotlib"),
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent
OUT = BASE_DIR / "outputs_paper42"
OUT.mkdir(exist_ok=True)

EPS = 1e-12
RNG_SEED = 42
N_PERM = 20000

# ============================================================
# Reference cosmology
# ============================================================

H0 = 67.4
Omega_m0 = 0.315
Omega_L0 = 0.685
sigma8_0 = 0.811
gamma_growth = 0.55


def E_lcdm(z):
    return np.sqrt(Omega_m0 * (1 + z) ** 3 + Omega_L0)


def Omega_m_z(z):
    Ez = E_lcdm(z)
    return Omega_m0 * (1 + z) ** 3 / (Ez ** 2)


def growth_factor_approx(z):
    """Carroll-Press-Turner style growth approximation, normalized to D(0)=1."""
    Om = Omega_m_z(z)
    Ol = Omega_L0 / E_lcdm(z) ** 2

    g = (5 * Om / 2) / (
        Om ** (4 / 7)
        - Ol
        + (1 + Om / 2) * (1 + Ol / 70)
    )

    g0 = (5 * Omega_m0 / 2) / (
        Omega_m0 ** (4 / 7)
        - Omega_L0
        + (1 + Omega_m0 / 2) * (1 + Omega_L0 / 70)
    )

    return g / (g0 * (1 + z))


def fsigma8_lcdm(z):
    D = growth_factor_approx(z)
    f = Omega_m_z(z) ** gamma_growth
    return f * D * sigma8_0


# ============================================================
# Compact fσ8 dataset
# ============================================================

growth_data = np.array([
    [0.020, 0.428, 0.046],
    [0.067, 0.423, 0.055],
    [0.100, 0.370, 0.130],
    [0.170, 0.510, 0.060],
    [0.220, 0.420, 0.070],
    [0.250, 0.351, 0.058],
    [0.370, 0.460, 0.038],
    [0.410, 0.450, 0.040],
    [0.570, 0.444, 0.038],
    [0.600, 0.430, 0.040],
    [0.780, 0.380, 0.040],
    [0.800, 0.470, 0.080],
    [1.400, 0.482, 0.116],
])

z_data = growth_data[:, 0]
fs8_obs = growth_data[:, 1]
fs8_err = growth_data[:, 2]
fs8_ref = fsigma8_lcdm(z_data)

residual = fs8_obs - fs8_ref
residual_sigma = residual / fs8_err
abs_residual_sigma = np.abs(residual_sigma)


# ============================================================
# Basic statistics
# ============================================================

def rankdata(x):
    """Simple average-rank implementation."""
    x = np.asarray(x)
    order = np.argsort(x)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(len(x), dtype=float)

    # average ties
    vals = x[order]
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and vals[j + 1] == vals[i]:
            j += 1
        if j > i:
            avg = 0.5 * (i + j)
            ranks[order[i:j+1]] = avg
        i = j + 1
    return ranks + 1.0


def pearsonr(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x - np.mean(x)
    y = y - np.mean(y)
    den = np.sqrt(np.sum(x*x) * np.sum(y*y))
    if den < EPS:
        return 0.0
    return float(np.sum(x*y) / den)


def spearmanr(x, y):
    return pearsonr(rankdata(x), rankdata(y))


def permutation_pvalue(x, y, stat_fn=spearmanr, n_perm=N_PERM, seed=RNG_SEED):
    rng = np.random.default_rng(seed)
    obs = stat_fn(x, y)
    vals = np.empty(n_perm)
    y = np.asarray(y)
    for i in range(n_perm):
        vals[i] = stat_fn(x, rng.permutation(y))
    p = (np.sum(np.abs(vals) >= abs(obs)) + 1) / (n_perm + 1)
    return float(obs), float(p), vals


# ============================================================
# Build blind defects from theory-side quantities
# ============================================================

def latent_depth(z, zc=1.1, sharpness=3.5, floor=0.20):
    return floor + (1 - floor) / (1 + np.exp(sharpness * (z - zc)))


def normalized(x):
    x = np.asarray(x, dtype=float)
    lo, hi = np.min(x), np.max(x)
    if abs(hi - lo) < EPS:
        return np.zeros_like(x)
    return (x - lo) / (hi - lo + EPS)


def build_blind_defects(z):
    """
    Blind defect construction:
    No observed residuals are used here.

    H: smoothness / curvature defect of growth
    B: expansion-balance defect from dE/dz
    S: structural activity from breadth growth
    V: coherence-depth imbalance
    """
    z = np.asarray(z)
    E = E_lcdm(z)
    D = growth_factor_approx(z)
    L = latent_depth(z)

    B_exp = E * (1 + z)
    D_acc = D * L

    # Numerical derivatives on grid
    dD_dz = np.gradient(D, z)
    d2D_dz2 = np.gradient(dD_dz, z)
    dE_dz = np.gradient(E, z)
    dBexp_dz = np.gradient(B_exp, z)

    # Defects: non-negative and dimensionless after normalization
    d_H = normalized(np.abs(d2D_dz2))
    d_B = normalized(np.abs(dE_dz / (E + EPS)))
    d_S = normalized(np.abs(dBexp_dz / (B_exp + EPS)))
    imbalance = np.abs(B_exp - D_acc) / (B_exp + D_acc + EPS)
    d_V = normalized(imbalance)

    return {
        "z": z,
        "E": E,
        "D": D,
        "L": L,
        "B_exp": B_exp,
        "D_acc": D_acc,
        "d_H": d_H,
        "d_B": d_B,
        "d_S": d_S,
        "d_V": d_V,
    }


def response_lambda_from_defects(defects, eta=0.05):
    Dmat = np.column_stack([
        defects["d_H"],
        defects["d_B"],
        defects["d_S"],
        defects["d_V"],
    ])

    mean_d = np.mean(Dmat, axis=0)

    # Target: lower-than-mean defect regime, not zero to avoid overreaction
    d_star = np.percentile(Dmat, 20, axis=0)

    C = np.cov(Dmat, rowvar=False)
    ridge = eta * (np.trace(C) / C.shape[0] + EPS) * np.eye(C.shape[0])

    lam = np.linalg.solve(C + ridge, mean_d - d_star)

    return lam, C, mean_d, d_star


def lambda_projection_params(lam):
    shares = np.abs(lam) / (np.sum(np.abs(lam)) + EPS)
    pi_H, pi_B, pi_S, pi_V = shares

    gamma_lam = (pi_V + pi_B) / (pi_S + EPS)
    Bstar_lam = 1.0 / (pi_H + pi_B + EPS)
    alpha_lam = 1.0 + pi_S / (pi_V + EPS)

    return {
        "pi_H": float(pi_H),
        "pi_B": float(pi_B),
        "pi_S": float(pi_S),
        "pi_V": float(pi_V),
        "gamma_lambda": float(gamma_lam),
        "Bstar_lambda": float(Bstar_lam),
        "alpha_lambda": float(alpha_lam),
    }


def cproj_from_lambda(defects, params):
    B_exp = defects["B_exp"]
    D_acc = defects["D_acc"]

    gamma_lam = params["gamma_lambda"]
    Bstar_lam = params["Bstar_lambda"]
    alpha_lam = params["alpha_lambda"]

    S_proj = B_exp * np.tanh((B_exp / (Bstar_lam + EPS)) ** alpha_lam)
    C_proj = S_proj / (1.0 + gamma_lam * D_acc + EPS)

    C_proj_norm = C_proj / (C_proj[0] + EPS)

    R_proj = np.abs(S_proj - gamma_lam * D_acc) / (
        S_proj + gamma_lam * D_acc + EPS
    )

    return S_proj, C_proj_norm, R_proj


def estimate_transition(z, R_proj):
    d1 = np.gradient(R_proj, z)
    d2 = np.gradient(d1, z)
    n = len(z)
    lo = int(0.08 * n)
    hi = int(0.92 * n)
    idx = lo + np.argmax(np.abs(d2[lo:hi]))
    return float(z[idx]), float(abs(d2[idx]))


# ============================================================
# v1.2.1 closure
# ============================================================

def support_from_defect(d):
    return 1.0 / (1.0 + np.maximum(d, 0.0))


def build_supports_from_blind_projection(z_grid, defects_grid, C_proj_grid, R_proj_grid):
    # Interpolate blind defects to data points
    d_H = np.interp(z_data, z_grid, defects_grid["d_H"])
    d_B = np.interp(z_data, z_grid, defects_grid["d_B"])
    d_S = np.interp(z_data, z_grid, defects_grid["d_S"])
    d_V = np.interp(z_data, z_grid, defects_grid["d_V"])

    H = support_from_defect(d_H)
    B = support_from_defect(d_B)
    S = support_from_defect(d_S)
    V = support_from_defect(d_V)

    R_resp = (H * B * V) ** (1 / 3)
    full_bound = (H * B * S * V) ** (1 / 4)
    R_rob = np.minimum(R_resp, full_bound)

    CCI_min = S / (H + B + V + EPS)
    CCI_diag = S / (H + B + V + R_rob + EPS)

    C_proj_data = np.interp(z_data, z_grid, C_proj_grid)
    R_proj_data = np.interp(z_data, z_grid, R_proj_grid)

    return {
        "H": H,
        "B": B,
        "S": S,
        "V": V,
        "R_resp": R_resp,
        "R_rob": R_rob,
        "CCI_min": CCI_min,
        "CCI_diag": CCI_diag,
        "C_proj_lambda": C_proj_data,
        "R_proj_lambda": R_proj_data,
    }


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 72)
    print("MAAT Paper 42 — Blind Projection Test")
    print("=" * 72)

    z_grid = np.linspace(0.0, 2.0, 900)

    defects_grid = build_blind_defects(z_grid)

    lam, C, mean_d, d_star = response_lambda_from_defects(defects_grid, eta=0.05)
    params = lambda_projection_params(lam)

    S_proj, C_proj_lambda, R_proj_lambda = cproj_from_lambda(defects_grid, params)
    z_tr, curvature_peak = estimate_transition(z_grid, R_proj_lambda)

    supports = build_supports_from_blind_projection(
        z_grid, defects_grid, C_proj_lambda, R_proj_lambda
    )

    diagnostics = {
        "C_proj_lambda": supports["C_proj_lambda"],
        "R_proj_lambda": supports["R_proj_lambda"],
        "H": supports["H"],
        "B": supports["B"],
        "S": supports["S"],
        "V": supports["V"],
        "R_resp": supports["R_resp"],
        "R_rob": supports["R_rob"],
        "CCI_min": supports["CCI_min"],
        "CCI_diag": supports["CCI_diag"],
    }

    # Correlation table
    rows = []
    for name, x in diagnostics.items():
        pear_s = pearsonr(x, residual_sigma)
        spear_s, p_s, _ = permutation_pvalue(x, residual_sigma, spearmanr, seed=RNG_SEED)

        pear_abs = pearsonr(x, abs_residual_sigma)
        spear_abs, p_abs, _ = permutation_pvalue(
            x, abs_residual_sigma, spearmanr, seed=RNG_SEED + 1
        )

        rows.append([name, pear_s, spear_s, p_s, pear_abs, spear_abs, p_abs])

    # Main null tests for CCI_diag
    cci_obs, cci_p, cci_null = permutation_pvalue(
        diagnostics["CCI_diag"], abs_residual_sigma, spearmanr, seed=RNG_SEED + 2
    )

    cproj_obs, cproj_p, cproj_null = permutation_pvalue(
        diagnostics["C_proj_lambda"], abs_residual_sigma, spearmanr, seed=RNG_SEED + 3
    )

    # Save CSV tables
    table = np.array(rows, dtype=object)
    with open(OUT / "paper42_correlation_table.csv", "w", encoding="utf-8") as f:
        f.write("diagnostic,pearson_signed,spearman_signed,p_signed,pearson_abs,spearman_abs,p_abs\n")
        for row in rows:
            f.write(",".join([row[0]] + [f"{v:.8f}" for v in row[1:]]) + "\n")

    data_out = np.column_stack([
        z_data,
        fs8_obs,
        fs8_err,
        fs8_ref,
        residual_sigma,
        abs_residual_sigma,
        diagnostics["C_proj_lambda"],
        diagnostics["R_proj_lambda"],
        diagnostics["H"],
        diagnostics["B"],
        diagnostics["S"],
        diagnostics["V"],
        diagnostics["R_resp"],
        diagnostics["R_rob"],
        diagnostics["CCI_min"],
        diagnostics["CCI_diag"],
    ])
    np.savetxt(
        OUT / "paper42_blind_signature_table.csv",
        data_out,
        delimiter=",",
        header=(
            "z,fsigma8_obs,sigma,fsigma8_LCDM,residual_sigma,abs_residual_sigma,"
            "C_proj_lambda,R_proj_lambda,H,B,S,V,R_resp,R_rob,CCI_min,CCI_diag"
        ),
        comments=""
    )

    projection_out = np.column_stack([
        z_grid,
        defects_grid["B_exp"],
        defects_grid["D_acc"],
        S_proj,
        C_proj_lambda,
        R_proj_lambda,
        defects_grid["d_H"],
        defects_grid["d_B"],
        defects_grid["d_S"],
        defects_grid["d_V"],
    ])
    np.savetxt(
        OUT / "paper42_lambda_projection_curve.csv",
        projection_out,
        delimiter=",",
        header="z,B_exp,D_acc,S_proj_lambda,C_proj_lambda,R_proj_lambda,d_H,d_B,d_S,d_V",
        comments=""
    )

    summary = {
        "model": "MAAT Paper 42 Blind Projection Test",
        "n_growth_points": int(len(z_data)),
        "lambda": {
            "lambda_H": float(lam[0]),
            "lambda_B": float(lam[1]),
            "lambda_S": float(lam[2]),
            "lambda_V": float(lam[3]),
        },
        "lambda_shares": {
            "pi_H": params["pi_H"],
            "pi_B": params["pi_B"],
            "pi_S": params["pi_S"],
            "pi_V": params["pi_V"],
        },
        "derived_projection_parameters": {
            "gamma_lambda": params["gamma_lambda"],
            "Bstar_lambda": params["Bstar_lambda"],
            "alpha_lambda": params["alpha_lambda"],
        },
        "transition": {
            "z_transition": z_tr,
            "curvature_peak": curvature_peak,
        },
        "mean_d": mean_d.tolist(),
        "d_star": d_star.tolist(),
        "correlations": {
            row[0]: {
                "pearson_signed": float(row[1]),
                "spearman_signed": float(row[2]),
                "p_signed": float(row[3]),
                "pearson_abs": float(row[4]),
                "spearman_abs": float(row[5]),
                "p_abs": float(row[6]),
            }
            for row in rows
        },
        "main_null_tests": {
            "CCI_diag_vs_abs_residual_sigma": {
                "spearman": float(cci_obs),
                "p_perm": float(cci_p),
                "null_mean": float(np.mean(cci_null)),
                "null_std": float(np.std(cci_null)),
            },
            "C_proj_lambda_vs_abs_residual_sigma": {
                "spearman": float(cproj_obs),
                "p_perm": float(cproj_p),
                "null_mean": float(np.mean(cproj_null)),
                "null_std": float(np.std(cproj_null)),
            },
        },
        "note": (
            "Blind projection test: projection parameters are derived from "
            "lambda response weights. No epsilon fit, no gamma scan, no template tuning."
        ),
    }

    with open(OUT / "paper42_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    # ========================================================
    # Figures
    # ========================================================

    plt.figure(figsize=(7.2, 4.6))
    plt.plot(z_grid, C_proj_lambda, label=r"$C_{\rm proj,\lambda}$")
    plt.plot(z_grid, R_proj_lambda, label=r"$R_{\rm proj,\lambda}$")
    plt.axvline(z_tr, linestyle="--", label=fr"$z_{{tr}}\approx {z_tr:.3f}$")
    plt.xlabel("z")
    plt.ylabel("blind projection diagnostics")
    plt.title("MAAT Paper 42 — Response-Derived Projection")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / "fig1_blind_projection_curve.png", dpi=220)
    plt.close()

    plt.figure(figsize=(7.2, 4.6))
    plt.bar(["H", "B", "S", "V"], [params["pi_H"], params["pi_B"], params["pi_S"], params["pi_V"]])
    plt.ylabel(r"response share $\pi_a$")
    plt.title("Response Shares Derived from Defect Covariance")
    plt.tight_layout()
    plt.savefig(OUT / "fig2_lambda_response_shares.png", dpi=220)
    plt.close()

    plt.figure(figsize=(7.2, 4.6))
    plt.scatter(diagnostics["CCI_diag"], abs_residual_sigma)
    plt.xlabel(r"$CCI_{\rm diag,\lambda}$")
    plt.ylabel(r"$|r_\sigma|$")
    plt.title(fr"Blind CCI vs Residual Magnitude: $\rho_S={cci_obs:.3f}$, p={cci_p:.4f}")
    plt.tight_layout()
    plt.savefig(OUT / "fig3_blind_cci_vs_abs_residual.png", dpi=220)
    plt.close()

    plt.figure(figsize=(7.2, 4.6))
    plt.scatter(diagnostics["C_proj_lambda"], abs_residual_sigma)
    plt.xlabel(r"$C_{\rm proj,\lambda}$")
    plt.ylabel(r"$|r_\sigma|$")
    plt.title(fr"Blind Projection vs Residual Magnitude: $\rho_S={cproj_obs:.3f}$, p={cproj_p:.4f}")
    plt.tight_layout()
    plt.savefig(OUT / "fig4_blind_cproj_vs_abs_residual.png", dpi=220)
    plt.close()

    plt.figure(figsize=(7.2, 4.6))
    plt.hist(cci_null, bins=40, alpha=0.75, label="permutation null")
    plt.axvline(cci_obs, linestyle="--", linewidth=2, label=fr"observed $\rho_S={cci_obs:.3f}$")
    plt.axvline(-cci_obs, linestyle=":", linewidth=2)
    plt.xlabel(r"Spearman $\rho_S$")
    plt.ylabel("count")
    plt.title("Null Test: Blind CCI vs Residual Magnitude")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / "fig5_null_test_blind_cci.png", dpi=220)
    plt.close()

    # Console output
    print("\n--- Response-derived lambda ---")
    print(f"lambda_H: {lam[0]: .6f}")
    print(f"lambda_B: {lam[1]: .6f}")
    print(f"lambda_S: {lam[2]: .6f}")
    print(f"lambda_V: {lam[3]: .6f}")

    print("\n--- Response shares ---")
    print(f"pi_H: {params['pi_H']:.6f}")
    print(f"pi_B: {params['pi_B']:.6f}")
    print(f"pi_S: {params['pi_S']:.6f}")
    print(f"pi_V: {params['pi_V']:.6f}")

    print("\n--- Derived projection parameters ---")
    print(f"gamma_lambda : {params['gamma_lambda']:.6f}")
    print(f"Bstar_lambda : {params['Bstar_lambda']:.6f}")
    print(f"alpha_lambda : {params['alpha_lambda']:.6f}")
    print(f"z_transition : {z_tr:.6f}")

    print("\n--- Main blind signature tests ---")
    print(f"CCI_diag vs |residual_sigma|     : Spearman={cci_obs:.4f}, p={cci_p:.4f}")
    print(f"C_proj_lambda vs |residual_sigma|: Spearman={cproj_obs:.4f}, p={cproj_p:.4f}")

    print("\n--- Full correlation table ---")
    print(f"{'diagnostic':>14} {'Sp signed':>10} {'p':>9} {'Sp abs':>10} {'p':>9}")
    for row in rows:
        print(f"{row[0]:>14} {row[2]:10.4f} {row[3]:9.4f} {row[5]:10.4f} {row[6]:9.4f}")

    print(f"\nSaved outputs to: {OUT.resolve()}")


if __name__ == "__main__":
    main()
