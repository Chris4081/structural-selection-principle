#!/usr/bin/env python3
"""
MAAT Paper 39 — Observable Signature Proxy with v1.2.1 Robustness Closure
=========================================================================

fσ8(z) modulated by C_proj(z), extended with MAAT v1.2.1:

    R_resp = (H B V)^(1/3)
    R_rob  = min(R_resp, (H B S V)^(1/4))
    Stability = R_rob

CCI:
    CCI_min  = S / (H + B + V + eps)
    CCI_diag = S / (H + B + V + R_rob + eps)

R_resp is not separately included in the CCI denominator to avoid double counting.

Run:
    python3 maat_paper39_observable_signature_v121.py
"""

import numpy as np
import json
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/codex-mpl-cache")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path("outputs_paper39")
OUT.mkdir(exist_ok=True)

# -----------------------------
# Reference cosmology
# -----------------------------
H0 = 67.4
Omega_m0 = 0.315
Omega_L0 = 0.685
sigma8_0 = 0.811
gamma_growth = 0.55

eps = 1e-12


def E_lcdm(z):
    return np.sqrt(Omega_m0 * (1 + z) ** 3 + Omega_L0)


def Omega_m_z(z):
    Ez = E_lcdm(z)
    return Omega_m0 * (1 + z) ** 3 / (Ez ** 2)


def growth_factor_approx(z):
    """
    Carroll-Press-Turner style approximate growth factor,
    normalized to D(0)=1.
    """
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


# -----------------------------
# Projection CCI model
# -----------------------------
def latent_depth(z, zc=1.1, sharpness=3.5, floor=0.20):
    return floor + (1 - floor) / (1 + np.exp(sharpness * (z - zc)))


def projection_template(z, alpha=1.0, gamma_proj=1.0):
    """
    Minimal C_proj model from Papers 33/34.
    """
    Ez = E_lcdm(z)
    D = growth_factor_approx(z)

    B_exp = Ez * (1 + z)
    Lz = latent_depth(z)
    D_acc = D * Lz

    B_norm = B_exp / B_exp[0]
    Q = np.tanh((B_norm ** alpha) / np.max(B_norm ** alpha))

    S_proj = B_exp * Q
    C_proj = S_proj / (1 + gamma_proj * D_acc + eps)
    C_norm = C_proj / C_proj[0]

    R_proj = np.abs(S_proj - gamma_proj * D_acc) / (
        S_proj + gamma_proj * D_acc + eps
    )

    return C_norm, R_proj


def estimate_transition(z, R_proj):
    """
    z_tr = argmax |d²R_proj/dz²| excluding edges.
    """
    d1 = np.gradient(R_proj, z)
    d2 = np.gradient(d1, z)

    n = len(z)
    lo = int(0.08 * n)
    hi = int(0.92 * n)

    idx_local = np.argmax(np.abs(d2[lo:hi]))
    idx = lo + idx_local

    return float(z[idx]), float(np.abs(d2[idx]))


# -----------------------------
# MAAT v1.2.1 Respect / Robustness closure
# -----------------------------
def support_from_defect(d):
    """Γ = 1/(1+d), support score in (0,1]."""
    return 1.0 / (1.0 + np.maximum(d, 0.0))


def r_resp(H, B, V):
    """Implicit respect: R_resp = (H B V)^(1/3)."""
    return np.power(np.maximum(H * B * V, eps), 1.0 / 3.0)


def r_rob(H, B, S, V):
    """Emergent robustness: R_rob = min(R_resp, (H B S V)^(1/4))."""
    Rresp = r_resp(H, B, V)
    full_bound = np.power(np.maximum(H * B * S * V, eps), 1.0 / 4.0)
    return np.minimum(Rresp, full_bound)


def cci_v121(H, B, S, V, use_robustness=True):
    """
    v1.2.1 CCI.

    Minimal:
        CCI = S / (H+B+V+eps)

    Diagnostic:
        CCI = S / (H+B+V+R_rob+eps)

    R_resp is not included separately to avoid double counting.
    """
    if use_robustness:
        Rrob = r_rob(H, B, S, V)
        return S / (H + B + V + Rrob + eps)

    return S / (H + B + V + eps)


def build_v121_projection_fields(z, C_proj, R_proj, fs8_lcdm, fs8_obs=None):
    """
    Builds H,B,S,V proxy supports for Paper 39.

    H: projection smoothness/coherence support
    B: data-balance support from fσ8 residuals when data are available
    S: projection activity support
    V: connectedness support from bounded projection residual R_proj

    These are diagnostic proxy supports, not fundamental field measurements.
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

    CCI_min = cci_v121(H, B, S, V, use_robustness=False)
    CCI_diag = cci_v121(H, B, S, V, use_robustness=True)

    return {
        "H": H,
        "B": B,
        "S": S,
        "V": V,
        "R_resp": Rresp,
        "R_rob": Rrob,
        "CCI_min": CCI_min,
        "CCI_diag": CCI_diag,
        "dH": dH,
        "dB": dB,
        "dV": dV,
    }


# -----------------------------
# Compact fσ8 dataset
# -----------------------------
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


# -----------------------------
# Build model grid
# -----------------------------
z_grid = np.linspace(0.0, 2.0, 800)

C_proj_grid, R_proj_grid = projection_template(
    z_grid,
    alpha=1.0,
    gamma_proj=1.0,
)

z_tr, curvature_peak = estimate_transition(z_grid, R_proj_grid)

v121_grid = build_v121_projection_fields(
    z_grid,
    C_proj_grid,
    R_proj_grid,
    fsigma8_lcdm(z_grid),
    fs8_obs=None,
)

C_data = np.interp(z_data, z_grid, C_proj_grid)
R_data = np.interp(z_data, z_grid, R_proj_grid)

C_log = np.log1p(C_data)
C_hat = (C_log - np.mean(C_log)) / (np.std(C_log) + eps)
C_hat = np.tanh(C_hat)

fs8_lcdm_data = fsigma8_lcdm(z_data)

v121_data = build_v121_projection_fields(
    z_data,
    C_data,
    R_data,
    fs8_lcdm_data,
    fs8_obs=fs8_obs,
)

R_resp_data = v121_data["R_resp"]
R_rob_data = v121_data["R_rob"]
CCI_min_data = v121_data["CCI_min"]
CCI_diag_data = v121_data["CCI_diag"]


def fsigma8_maat_proxy(epsilon):
    return fs8_lcdm_data * (1 + epsilon * C_hat)


def chi2(epsilon):
    model = fsigma8_maat_proxy(epsilon)
    return float(np.sum(((model - fs8_obs) / fs8_err) ** 2))


# -----------------------------
# epsilon scan
# -----------------------------
eps_values = np.linspace(-0.01, 0.01, 1001)
chi2_values = np.array([chi2(e) for e in eps_values])

best_idx = np.argmin(chi2_values)
eps_best = float(eps_values[best_idx])
chi2_best = float(chi2_values[best_idx])
chi2_lcdm = chi2(0.0)
delta_chi2 = chi2_best - chi2_lcdm

fs8_best = fsigma8_maat_proxy(eps_best)
frac_dev = (fs8_best - fs8_lcdm_data) / np.maximum(fs8_lcdm_data, eps)

max_abs_frac_dev = float(np.max(np.abs(frac_dev)))
mean_abs_frac_dev = float(np.mean(np.abs(frac_dev)))


# -----------------------------
# Save CSV
# -----------------------------
csv_data = np.column_stack([
    z_data,
    fs8_obs,
    fs8_err,
    fs8_lcdm_data,
    fs8_best,
    C_data,
    C_hat,
    frac_dev,
    v121_data["H"],
    v121_data["B"],
    v121_data["S"],
    v121_data["V"],
    R_resp_data,
    R_rob_data,
    CCI_min_data,
    CCI_diag_data,
])

np.savetxt(
    OUT / "paper39_growth_signature_results.csv",
    csv_data,
    delimiter=",",
    header=(
        "z,fsigma8_obs,sigma,fsigma8_lcdm,fsigma8_maat,"
        "C_proj,C_hat,frac_dev,H,B,S,V,R_resp,R_rob,CCI_min,CCI_diag"
    ),
    comments=""
)

scan_data = np.column_stack([eps_values, chi2_values])
np.savetxt(
    OUT / "paper39_epsilon_chi2_scan.csv",
    scan_data,
    delimiter=",",
    header="epsilon,chi2",
    comments=""
)

projection_data = np.column_stack([
    z_grid,
    C_proj_grid,
    R_proj_grid,
    v121_grid["H"],
    v121_grid["B"],
    v121_grid["S"],
    v121_grid["V"],
    v121_grid["R_resp"],
    v121_grid["R_rob"],
    v121_grid["CCI_min"],
    v121_grid["CCI_diag"],
])

np.savetxt(
    OUT / "paper39_projection_template_v121.csv",
    projection_data,
    delimiter=",",
    header="z,C_proj_norm,R_proj,H,B,S,V,R_resp,R_rob,CCI_min,CCI_diag",
    comments=""
)


# -----------------------------
# Save JSON summary
# -----------------------------
summary = {
    "model": "MAAT Paper 39 Observable Signature Proxy v1.2.1",
    "definition": "O_MAAT(z)=fsigma8_LCDM(z)*(1+epsilon*C_hat_proj(z))",
    "n_growth_points": int(len(z_data)),
    "epsilon_scan_range": [-0.01, 0.01],
    "epsilon_best": eps_best,
    "chi2_LCDM": chi2_lcdm,
    "chi2_MAAT_best": chi2_best,
    "delta_chi2": delta_chi2,
    "max_abs_frac_dev": max_abs_frac_dev,
    "mean_abs_frac_dev": mean_abs_frac_dev,
    "z_transition_estimate": z_tr,
    "transition_curvature_peak": curvature_peak,
    "maat_v121_closure": {
        "R_resp": "(H*B*V)^(1/3)",
        "R_rob": "min(R_resp, (H*B*S*V)^(1/4))",
        "CCI_min": "S/(H+B+V+eps)",
        "CCI_diag": "S/(H+B+V+R_rob+eps)",
        "R_resp_mean": float(np.mean(R_resp_data)),
        "R_rob_mean": float(np.mean(R_rob_data)),
        "CCI_min_mean": float(np.mean(CCI_min_data)),
        "CCI_diag_mean": float(np.mean(CCI_diag_data)),
        "note": "R_resp is not included separately in the CCI denominator to avoid double counting."
    },
    "note": (
        "Diagnostic proxy only. Not a full perturbation calculation, "
        "not a Boltzmann-code result, and not a precision cosmological fit."
    ),
}

with open(OUT / "paper39_summary.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)


# -----------------------------
# Figures
# -----------------------------
plt.figure(figsize=(7, 4.5))
plt.plot(z_grid, C_proj_grid, label=r"$C_{\rm proj}(z)$")
plt.axvline(z_tr, linestyle="--", label=fr"$z_{{tr}}\approx {z_tr:.3f}$")
plt.xlabel("z")
plt.ylabel(r"normalised $C_{\rm proj}$")
plt.title("MAAT Projection Template")
plt.legend()
plt.tight_layout()
plt.savefig(OUT / "fig1_projection_template.png", dpi=220)
plt.close()

plt.figure(figsize=(7, 4.5))
plt.errorbar(z_data, fs8_obs, yerr=fs8_err, fmt="o", label="observed")
plt.plot(z_data, fs8_lcdm_data, "s-", label=r"$\Lambda$CDM")
plt.plot(z_data, fs8_best, "^-", label=fr"MAAT proxy $\epsilon={eps_best:.4f}$")
plt.xlabel("z")
plt.ylabel(r"$f\sigma_8(z)$")
plt.title("Growth Observable with MAAT Projection Modulation")
plt.legend()
plt.tight_layout()
plt.savefig(OUT / "fig2_growth_comparison.png", dpi=220)
plt.close()

plt.figure(figsize=(7, 4.5))
res_lcdm = fs8_obs - fs8_lcdm_data
res_maat = fs8_obs - fs8_best
plt.axhline(0, color="black", linewidth=0.8)
plt.errorbar(z_data, res_lcdm, yerr=fs8_err, fmt="o", label=r"$\Lambda$CDM residual")
plt.errorbar(z_data, res_maat, yerr=fs8_err, fmt="s", label="MAAT proxy residual")
plt.xlabel("z")
plt.ylabel(r"residual in $f\sigma_8$")
plt.title("Residual Structure")
plt.legend()
plt.tight_layout()
plt.savefig(OUT / "fig3_residuals.png", dpi=220)
plt.close()

plt.figure(figsize=(7, 4.5))
plt.plot(eps_values, chi2_values)
plt.axvline(eps_best, linestyle="--", label=fr"best $\epsilon={eps_best:.4f}$")
plt.axhline(chi2_lcdm, linestyle=":", label=fr"$\Lambda$CDM $\chi^2={chi2_lcdm:.3f}$")
plt.xlabel(r"$\epsilon$")
plt.ylabel(r"$\chi^2$")
plt.title(r"$\chi^2(\epsilon)$ Scan")
plt.legend()
plt.tight_layout()
plt.savefig(OUT / "fig4_chi2_scan.png", dpi=220)
plt.close()

plt.figure(figsize=(7, 4.5))
plt.plot(z_grid, v121_grid["R_resp"], label=r"$R_{\rm resp}$")
plt.plot(z_grid, v121_grid["R_rob"], "--", label=r"$R_{\rm rob}$")
plt.xlabel("z")
plt.ylabel("closure score")
plt.title("MAAT v1.2.1 Respect / Robustness Closure")
plt.legend()
plt.tight_layout()
plt.savefig(OUT / "fig5_v121_robustness_closure.png", dpi=220)
plt.close()

plt.figure(figsize=(7, 4.5))
plt.plot(z_grid, v121_grid["CCI_min"], label=r"$CCI_{\rm min}$")
plt.plot(z_grid, v121_grid["CCI_diag"], "--", label=r"$CCI_{\rm diag}$")
plt.xlabel("z")
plt.ylabel("CCI")
plt.title("MAAT v1.2.1 CCI Diagnostics")
plt.legend()
plt.tight_layout()
plt.savefig(OUT / "fig6_v121_cci.png", dpi=220)
plt.close()

plt.figure(figsize=(7, 4.5))
plt.plot(z_grid, v121_grid["H"], label="H")
plt.plot(z_grid, v121_grid["B"], label="B")
plt.plot(z_grid, v121_grid["S"], label="S")
plt.plot(z_grid, v121_grid["V"], label="V")
plt.xlabel("z")
plt.ylabel("support")
plt.title("MAAT v1.2.1 Projection Support Fields")
plt.legend()
plt.tight_layout()
plt.savefig(OUT / "fig7_v121_support_fields.png", dpi=220)
plt.close()


# -----------------------------
# Console output
# -----------------------------
print("=== MAAT Paper 39 Observable Signature Proxy v1.2.1 ===")
print(f"Growth points              : {len(z_data)}")
print(f"z_transition estimate      : {z_tr:.6f}")
print(f"epsilon_best               : {eps_best:.6f}")
print(f"chi2_LCDM                  : {chi2_lcdm:.6f}")
print(f"chi2_MAAT_best             : {chi2_best:.6f}")
print(f"delta_chi2                 : {delta_chi2:.6f}")
print(f"max |Δfσ8/fσ8|             : {100*max_abs_frac_dev:.4f}%")
print(f"mean |Δfσ8/fσ8|            : {100*mean_abs_frac_dev:.4f}%")
print(f"mean R_resp                : {np.mean(R_resp_data):.6f}")
print(f"mean R_rob                 : {np.mean(R_rob_data):.6f}")
print(f"mean CCI_min               : {np.mean(CCI_min_data):.6f}")
print(f"mean CCI_diag              : {np.mean(CCI_diag_data):.6f}")
print(f"Saved outputs to           : {OUT.resolve()}")
