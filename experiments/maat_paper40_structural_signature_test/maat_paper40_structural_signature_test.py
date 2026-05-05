#!/usr/bin/env python3
"""
MAAT Paper 40 — Structural Signature Test of Projection CCI
==========================================================

Tests whether MAAT v1.2.1 structural diagnostics correlate with
growth-data residuals.

Core question:
    Does CCI_diag(z_i) align with residual structure in fσ8 data?

Residual definitions:
    residual_raw   = fσ8_obs - fσ8_LCDM
    residual_sigma = (fσ8_obs - fσ8_LCDM) / sigma_i

Tested diagnostics:
    CCI_min
    CCI_diag
    R_resp
    R_rob
    H, B, S, V
    C_proj

Scientific status:
    Diagnostic signature test only.
    Not a full cosmological likelihood.
    Not a Boltzmann-code perturbation calculation.

Run:
    python3 maat_paper40_cci_signature_test.py
"""

import json
import os
from pathlib import Path

import numpy as np
os.environ.setdefault("MPLCONFIGDIR", "/tmp/codex-mpl-cache")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("outputs_paper40")
OUT.mkdir(exist_ok=True)

eps = 1e-12

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
    return Omega_m0 * (1 + z) ** 3 / np.maximum(Ez ** 2, eps)


def growth_factor_approx(z):
    """
    Carroll-Press-Turner style approximate growth factor,
    normalized to D(0)=1.
    """
    Om = Omega_m_z(z)
    Ol = Omega_L0 / np.maximum(E_lcdm(z) ** 2, eps)

    g = (5 * Om / 2) / np.maximum(
        Om ** (4 / 7)
        - Ol
        + (1 + Om / 2) * (1 + Ol / 70),
        eps,
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
# Projection CCI model
# ============================================================

def latent_depth(z, zc=1.1, sharpness=3.5, floor=0.20):
    return floor + (1 - floor) / (1 + np.exp(sharpness * (z - zc)))


def projection_template(z, alpha=1.0, gamma_proj=1.0):
    """
    Minimal C_proj model from Papers 33/34/39.
    """
    Ez = E_lcdm(z)
    D = growth_factor_approx(z)

    B_exp = Ez * (1 + z)
    Lz = latent_depth(z)
    D_acc = D * Lz

    B_norm = B_exp / np.maximum(B_exp[0], eps)
    Q = np.tanh((B_norm ** alpha) / np.maximum(np.max(B_norm ** alpha), eps))

    S_proj = B_exp * Q
    C_proj = S_proj / np.maximum(1 + gamma_proj * D_acc, eps)
    C_norm = C_proj / np.maximum(C_proj[0], eps)

    R_proj = np.abs(S_proj - gamma_proj * D_acc) / np.maximum(
        S_proj + gamma_proj * D_acc,
        eps,
    )

    return C_norm, R_proj


def estimate_transition(z, R_proj):
    d1 = np.gradient(R_proj, z)
    d2 = np.gradient(d1, z)

    n = len(z)
    lo = int(0.08 * n)
    hi = int(0.92 * n)

    idx_local = int(np.argmax(np.abs(d2[lo:hi])))
    idx = lo + idx_local

    return float(z[idx]), float(np.abs(d2[idx]))


# ============================================================
# MAAT v1.2.1 closure
# ============================================================

def support_from_defect(d):
    return 1.0 / (1.0 + np.maximum(d, 0.0))


def r_resp(H, B, V):
    return np.power(np.maximum(H * B * V, eps), 1.0 / 3.0)


def r_rob(H, B, S, V):
    Rresp = r_resp(H, B, V)
    full_bound = np.power(np.maximum(H * B * S * V, eps), 1.0 / 4.0)
    return np.minimum(Rresp, full_bound)


def cci_v121(H, B, S, V, use_robustness=True):
    if use_robustness:
        Rrob = r_rob(H, B, S, V)
        return S / np.maximum(H + B + V + Rrob + eps, eps)
    return S / np.maximum(H + B + V + eps, eps)


def build_v121_projection_fields(z, C_proj, R_proj, fs8_lcdm, fs8_obs=None):
    """
    Builds proxy supports.

    H: projection smoothness/coherence support
    B: data-balance support from fσ8 residuals when data are available
    S: projection activity support
    V: connectedness support from bounded projection residual R_proj
    """
    Cn = C_proj / np.maximum(np.max(C_proj), eps)

    dC = np.gradient(Cn, z)
    dH = np.abs(dC) / np.maximum(np.max(np.abs(dC)), eps)

    S = support_from_defect(1.0 - Cn)

    dV = R_proj / np.maximum(np.max(R_proj), eps)

    if fs8_obs is not None:
        frac_res = np.abs(fs8_obs - fs8_lcdm) / np.maximum(np.abs(fs8_lcdm), eps)
        dB = frac_res / np.maximum(np.max(frac_res), eps)
    else:
        dB = np.abs(Cn - np.mean(Cn)) / np.maximum(np.std(Cn), eps)

    H = support_from_defect(dH)
    B = support_from_defect(dB)
    V = support_from_defect(dV)

    Rresp = r_resp(H, B, V)
    Rrob = r_rob(H, B, S, V)

    return {
        "H": H,
        "B": B,
        "S": S,
        "V": V,
        "R_resp": Rresp,
        "R_rob": Rrob,
        "CCI_min": cci_v121(H, B, S, V, use_robustness=False),
        "CCI_diag": cci_v121(H, B, S, V, use_robustness=True),
    }


# ============================================================
# Correlation tools
# ============================================================

def pearson_corr(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    if len(x) < 3:
        return np.nan

    x0 = x - np.mean(x)
    y0 = y - np.mean(y)

    denom = np.sqrt(np.sum(x0 ** 2) * np.sum(y0 ** 2))
    if denom <= eps:
        return np.nan

    return float(np.sum(x0 * y0) / denom)


def rankdata_average(x):
    """
    Average-rank implementation for Spearman correlation without scipy.
    """
    x = np.asarray(x)
    order = np.argsort(x)
    ranks = np.empty(len(x), dtype=float)

    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and x[order[j + 1]] == x[order[i]]:
            j += 1

        avg_rank = 0.5 * (i + j) + 1.0
        ranks[order[i:j + 1]] = avg_rank
        i = j + 1

    return ranks


def spearman_corr(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    if len(x) < 3:
        return np.nan

    rx = rankdata_average(x)
    ry = rankdata_average(y)

    return pearson_corr(rx, ry)


def permutation_pvalue(x, y, corr_fn, n_perm=20000, seed=42):
    """
    Two-sided permutation p-value for correlation.

    Uses deterministic RNG seed for reproducibility.
    """
    rng = np.random.default_rng(seed)

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    observed = corr_fn(x, y)

    if not np.isfinite(observed):
        return np.nan

    count = 0
    abs_obs = abs(observed)

    for _ in range(n_perm):
        yp = rng.permutation(y)
        r = corr_fn(x, yp)
        if np.isfinite(r) and abs(r) >= abs_obs:
            count += 1

    return float((count + 1) / (n_perm + 1))


def correlation_report(name, x, residuals_raw, residuals_sigma, n_perm=20000):
    pear_raw = pearson_corr(x, residuals_raw)
    spear_raw = spearman_corr(x, residuals_raw)

    pear_sig = pearson_corr(x, residuals_sigma)
    spear_sig = spearman_corr(x, residuals_sigma)

    p_pear_raw = permutation_pvalue(x, residuals_raw, pearson_corr, n_perm=n_perm, seed=101)
    p_spear_raw = permutation_pvalue(x, residuals_raw, spearman_corr, n_perm=n_perm, seed=102)
    p_pear_sig = permutation_pvalue(x, residuals_sigma, pearson_corr, n_perm=n_perm, seed=103)
    p_spear_sig = permutation_pvalue(x, residuals_sigma, spearman_corr, n_perm=n_perm, seed=104)

    return {
        "name": name,
        "pearson_raw": pear_raw,
        "spearman_raw": spear_raw,
        "pearson_sigma": pear_sig,
        "spearman_sigma": spear_sig,
        "p_pearson_raw_perm": p_pear_raw,
        "p_spearman_raw_perm": p_spear_raw,
        "p_pearson_sigma_perm": p_pear_sig,
        "p_spearman_sigma_perm": p_spear_sig,
    }




# ============================================================
# Robustness / null-model tests (Maatis v1.2.1)
# ============================================================

def random_field_null_test(diagnostic, target, corr_fn=None, n_perm=50000, seed=777):
    if corr_fn is None:
        corr_fn = spearman_corr
    rng = np.random.default_rng(seed)
    diagnostic = np.asarray(diagnostic, dtype=float)
    target = np.asarray(target, dtype=float)
    observed = corr_fn(diagnostic, target)
    null_corrs = np.array([corr_fn(rng.permutation(diagnostic), target) for _ in range(n_perm)])
    p_value = (np.sum(np.abs(null_corrs) >= abs(observed)) + 1) / (n_perm + 1)
    return {"observed_corr": float(observed), "null_mean": float(np.mean(null_corrs)),
            "null_std": float(np.std(null_corrs)), "p_value": float(p_value), "n_perm": int(n_perm)}


def redshift_shuffle_test(z, diagnostic, target, corr_fn=None, n_perm=50000, seed=888):
    if corr_fn is None:
        corr_fn = spearman_corr
    rng = np.random.default_rng(seed)
    diagnostic = np.asarray(diagnostic, dtype=float)
    target = np.asarray(target, dtype=float)
    observed = corr_fn(diagnostic, target)
    null_corrs = np.array([corr_fn(diagnostic[rng.permutation(len(z))], target) for _ in range(n_perm)])
    p_value = (np.sum(np.abs(null_corrs) >= abs(observed)) + 1) / (n_perm + 1)
    return {"observed_corr": float(observed), "shuffle_mean": float(np.mean(null_corrs)),
            "shuffle_std": float(np.std(null_corrs)), "p_value": float(p_value), "n_perm": int(n_perm)}

# ============================================================
# Data
# ============================================================

growth_data = np.array([
    [0.02, 0.428, 0.046],
    [0.067, 0.423, 0.055],
    [0.10, 0.370, 0.130],
    [0.17, 0.510, 0.060],
    [0.22, 0.420, 0.070],
    [0.25, 0.351, 0.058],
    [0.37, 0.460, 0.038],
    [0.41, 0.450, 0.040],
    [0.57, 0.444, 0.038],
    [0.60, 0.430, 0.040],
    [0.78, 0.380, 0.040],
    [0.80, 0.470, 0.080],
    [1.40, 0.482, 0.116],
])

z_data = growth_data[:, 0]
fs8_obs = growth_data[:, 1]
fs8_err = growth_data[:, 2]


# ============================================================
# Main pipeline
# ============================================================

def main():
    z_grid = np.linspace(0.0, 2.0, 800)

    C_proj_grid, R_proj_grid = projection_template(z_grid)
    z_tr, curvature_peak = estimate_transition(z_grid, R_proj_grid)

    C_data = np.interp(z_data, z_grid, C_proj_grid)
    R_data = np.interp(z_data, z_grid, R_proj_grid)
    fs8_lcdm_data = fsigma8_lcdm(z_data)

    fields = build_v121_projection_fields(
        z_data,
        C_data,
        R_data,
        fs8_lcdm_data,
        fs8_obs=fs8_obs,
    )

    residual_raw = fs8_obs - fs8_lcdm_data
    residual_sigma = residual_raw / np.maximum(fs8_err, eps)
    abs_residual_sigma = np.abs(residual_sigma)

    diagnostics = {
        "C_proj": C_data,
        "R_proj": R_data,
        "H": fields["H"],
        "B": fields["B"],
        "S": fields["S"],
        "V": fields["V"],
        "R_resp": fields["R_resp"],
        "R_rob": fields["R_rob"],
        "CCI_min": fields["CCI_min"],
        "CCI_diag": fields["CCI_diag"],
    }

    reports = []
    for name, values in diagnostics.items():
        reports.append(
            correlation_report(
                name,
                values,
                residual_raw,
                residual_sigma,
                n_perm=20000,
            )
        )

    # Additional test: correlation with residual magnitude
    magnitude_reports = []
    for name, values in diagnostics.items():
        magnitude_reports.append({
            "name": name,
            "pearson_abs_sigma": pearson_corr(values, abs_residual_sigma),
            "spearman_abs_sigma": spearman_corr(values, abs_residual_sigma),
            "p_pearson_abs_sigma_perm": permutation_pvalue(
                values,
                abs_residual_sigma,
                pearson_corr,
                n_perm=20000,
                seed=201,
            ),
            "p_spearman_abs_sigma_perm": permutation_pvalue(
                values,
                abs_residual_sigma,
                spearman_corr,
                n_perm=20000,
                seed=202,
            ),
        })

    # Save data table
    table = np.column_stack([
        z_data,
        fs8_obs,
        fs8_err,
        fs8_lcdm_data,
        residual_raw,
        residual_sigma,
        abs_residual_sigma,
        C_data,
        R_data,
        fields["H"],
        fields["B"],
        fields["S"],
        fields["V"],
        fields["R_resp"],
        fields["R_rob"],
        fields["CCI_min"],
        fields["CCI_diag"],
    ])

    np.savetxt(
        OUT / "paper40_signature_table.csv",
        table,
        delimiter=",",
        header=(
            "z,fsigma8_obs,sigma,fsigma8_lcdm,"
            "residual_raw,residual_sigma,abs_residual_sigma,"
            "C_proj,R_proj,H,B,S,V,R_resp,R_rob,CCI_min,CCI_diag"
        ),
        comments="",
    )


    # Breakthrough / robustness tests (Maatis v1.2.1)
    breakthrough_tests = {}
    primary_signature = fields["CCI_diag"]
    target_signature = abs_residual_sigma
    breakthrough_tests["CCI_diag_vs_abs_residual_sigma"] = {
        "signature": "corr(CCI_diag, |residual_sigma|)",
        "original_spearman": float(spearman_corr(primary_signature, target_signature)),
        "original_pearson": float(pearson_corr(primary_signature, target_signature)),
        "random_field_null_spearman": random_field_null_test(primary_signature, target_signature, n_perm=50000, seed=777),
        "redshift_shuffle_spearman": redshift_shuffle_test(z_data, primary_signature, target_signature, n_perm=50000, seed=888),
    }
    breakthrough_tests["CCI_min_vs_abs_residual_sigma"] = {
        "signature": "corr(CCI_min, |residual_sigma|)",
        "original_spearman": float(spearman_corr(fields["CCI_min"], target_signature)),
        "original_pearson": float(pearson_corr(fields["CCI_min"], target_signature)),
        "random_field_null_spearman": random_field_null_test(fields["CCI_min"], target_signature, n_perm=50000, seed=779),
        "redshift_shuffle_spearman": redshift_shuffle_test(z_data, fields["CCI_min"], target_signature, n_perm=50000, seed=889),
    }

    # Save reports
    summary = {
        "model": "MAAT Paper 40 Structural Signature Test",
        "n_growth_points": int(len(z_data)),
        "z_transition_estimate": z_tr,
        "transition_curvature_peak": curvature_peak,
        "tested_hypothesis": (
            "Structural diagnostics, especially CCI_diag, should correlate "
            "with the residual structure of f_sigma8 observations relative to LCDM."
        ),
        "residual_definitions": {
            "residual_raw": "fsigma8_obs - fsigma8_LCDM",
            "residual_sigma": "(fsigma8_obs - fsigma8_LCDM)/sigma",
            "abs_residual_sigma": "abs(residual_sigma)",
        },
        "correlations_signed_residuals": reports,
        "correlations_residual_magnitude": magnitude_reports,
        "breakthrough_tests": breakthrough_tests,
        "note": (
            "Diagnostic signature test only. Small sample size N=13. "
            "Permutation p-values are exploratory and not a full cosmological likelihood."
        ),
    }

    with open(OUT / "paper40_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    # Console output
    print("=== MAAT Paper 40 Structural Signature Test ===")
    print(f"Growth points         : {len(z_data)}")
    print(f"z_transition estimate : {z_tr:.6f}")
    print("\n--- Signed residual correlations ---")
    print(f"{'diagnostic':>10}  {'Pearson σ':>10}  {'p':>9}  {'Spearman σ':>11}  {'p':>9}")

    for r in reports:
        print(
            f"{r['name']:>10}  "
            f"{r['pearson_sigma']:10.4f}  {r['p_pearson_sigma_perm']:9.4f}  "
            f"{r['spearman_sigma']:11.4f}  {r['p_spearman_sigma_perm']:9.4f}"
        )

    print("\n--- Residual magnitude correlations |residual_sigma| ---")
    print(f"{'diagnostic':>10}  {'Pearson |σ|':>12}  {'p':>9}  {'Spearman |σ|':>13}  {'p':>9}")

    for r in magnitude_reports:
        print(
            f"{r['name']:>10}  "
            f"{r['pearson_abs_sigma']:12.4f}  {r['p_pearson_abs_sigma_perm']:9.4f}  "
            f"{r['spearman_abs_sigma']:13.4f}  {r['p_spearman_abs_sigma_perm']:9.4f}"
        )

    # Figures
    plt.figure(figsize=(7, 4.8))
    plt.axhline(0, color="black", linewidth=0.8)
    plt.errorbar(
        z_data,
        residual_sigma,
        yerr=np.ones_like(residual_sigma),
        fmt="o",
        label=r"$r_\sigma=(f\sigma_{8,obs}-f\sigma_{8,\Lambda CDM})/\sigma$",
    )
    plt.plot(z_data, fields["CCI_diag"], "s-", label=r"$CCI_{\rm diag}$")
    plt.xlabel("z")
    plt.ylabel("scaled value")
    plt.title("CCI_diag vs signed residual structure")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / "fig1_cci_diag_vs_signed_residual.png", dpi=220)
    plt.close()

    plt.figure(figsize=(7, 4.8))
    plt.scatter(fields["CCI_diag"], residual_sigma, s=65)
    for zi, x, y in zip(z_data, fields["CCI_diag"], residual_sigma):
        plt.annotate(f"{zi:.2f}", (x, y), fontsize=8, alpha=0.75)
    plt.axhline(0, color="black", linewidth=0.8)
    plt.xlabel(r"$CCI_{\rm diag}$")
    plt.ylabel(r"signed residual $r_\sigma$")
    plt.title("Signature scatter: CCI_diag vs residual")
    plt.tight_layout()
    plt.savefig(OUT / "fig2_scatter_cci_diag_signed_residual.png", dpi=220)
    plt.close()

    plt.figure(figsize=(7, 4.8))
    plt.scatter(fields["CCI_diag"], abs_residual_sigma, s=65)
    for zi, x, y in zip(z_data, fields["CCI_diag"], abs_residual_sigma):
        plt.annotate(f"{zi:.2f}", (x, y), fontsize=8, alpha=0.75)
    plt.xlabel(r"$CCI_{\rm diag}$")
    plt.ylabel(r"$|r_\sigma|$")
    plt.title("Signature scatter: CCI_diag vs residual magnitude")
    plt.tight_layout()
    plt.savefig(OUT / "fig3_scatter_cci_diag_abs_residual.png", dpi=220)
    plt.close()

    plt.figure(figsize=(9, 4.8))
    names = [r["name"] for r in reports]
    vals = [r["spearman_sigma"] for r in reports]
    plt.bar(names, vals)
    plt.axhline(0, color="black", linewidth=0.8)
    plt.xticks(rotation=35, ha="right")
    plt.ylabel("Spearman correlation with signed residual")
    plt.title("Structural diagnostics vs residual structure")
    plt.tight_layout()
    plt.savefig(OUT / "fig4_spearman_signed_residuals.png", dpi=220)
    plt.close()

    plt.figure(figsize=(9, 4.8))
    names = [r["name"] for r in magnitude_reports]
    vals = [r["spearman_abs_sigma"] for r in magnitude_reports]
    plt.bar(names, vals)
    plt.axhline(0, color="black", linewidth=0.8)
    plt.xticks(rotation=35, ha="right")
    plt.ylabel("Spearman correlation with |residual_sigma|")
    plt.title("Structural diagnostics vs residual magnitude")
    plt.tight_layout()
    plt.savefig(OUT / "fig5_spearman_abs_residuals.png", dpi=220)
    plt.close()


    # Null distribution figure
    rng_plot = np.random.default_rng(999)
    null_corrs_plot = np.array([
        spearman_corr(rng_plot.permutation(fields["CCI_diag"]), abs_residual_sigma)
        for _ in range(10000)
    ])
    obs_val = spearman_corr(fields["CCI_diag"], abs_residual_sigma)
    plt.figure(figsize=(7, 4.8))
    plt.hist(null_corrs_plot, bins=40, alpha=0.75, color="steelblue", label="randomized null")
    plt.axvline(obs_val, color="red", linewidth=2, label=f"observed r={obs_val:.3f}")
    plt.axvline(-obs_val, color="red", linestyle="--", linewidth=1)
    plt.xlabel("Spearman r")
    plt.ylabel("count")
    plt.title("Null test: CCI_diag vs |r_sigma|")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / "fig6_breakthrough_null_test_cci_diag.png", dpi=220)
    plt.close()

    print("\n--- Breakthrough robustness tests ---")
    for name, test in breakthrough_tests.items():
        print(f"\n{name}")
        print(f"  Original Spearman : {test['original_spearman']:.4f}")
        rn = test['random_field_null_spearman']
        rs = test['redshift_shuffle_spearman']
        print(f"  Random-field null : mean={rn['null_mean']:.4f}, std={rn['null_std']:.4f}, p={rn['p_value']:.4f}")
        print(f"  Redshift shuffle  : mean={rs['shuffle_mean']:.4f}, std={rs['shuffle_std']:.4f}, p={rs['p_value']:.4f}")

    print(f"\nSaved outputs to: {OUT.resolve()}")


if __name__ == "__main__":
    main()
