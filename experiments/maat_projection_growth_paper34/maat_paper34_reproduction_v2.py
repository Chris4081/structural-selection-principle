#!/usr/bin/env python3
"""
MAAT Paper 34 — Full Reproduction Script v2
Christof Krieg, 2026, maat-research.com
"""

from __future__ import annotations
import json
from pathlib import Path
from dataclasses import dataclass, asdict
import numpy as np
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 11, "figure.dpi": 150, "savefig.dpi": 150})

OUTDIR = Path("paper34_outputs")
OUTDIR.mkdir(exist_ok=True)
EPS = 1e-10

@dataclass(frozen=True)
class LCDMParams:
    H0: float = 67.4
    Omm: float = 0.315
    OmL: float = 0.685
    sigma8_0: float = 0.811
    gamma_growth: float = 0.55

@dataclass(frozen=True)
class ProjectionParams:
    gamma: float = 2.0
    alpha: float = 1.0
    zc: float = 1.1
    s: float = 3.5
    L_floor: float = 0.20
    edge_frac: float = 0.08

COSMO = LCDMParams()
BASE = ProjectionParams()

def E(z, p=COSMO):
    z = np.asarray(z, dtype=float)
    return np.sqrt(p.Omm*(1+z)**3 + p.OmL)

def omega_m_z(z, p=COSMO):
    z = np.asarray(z, dtype=float)
    return p.Omm*(1+z)**3 / E(z,p)**2

def omega_lambda_z(z, p=COSMO):
    z = np.asarray(z, dtype=float)
    return p.OmL / E(z,p)**2

def g_cpt(z, p=COSMO):
    om = omega_m_z(z,p); ol = omega_lambda_z(z,p)
    return (5*om/2) / (om**(4/7) - ol + (1+om/2)*(1+ol/70))

def D_growth(z, p=COSMO):
    z = np.asarray(z, dtype=float)
    return (g_cpt(z,p)/(1+z)) / g_cpt(0.0,p)

def f_growth(z, p=COSMO):
    return omega_m_z(z,p)**p.gamma_growth

def fsigma8_lcdm(z, p=COSMO):
    return f_growth(z,p) * p.sigma8_0 * D_growth(z,p)

def L_latent(z, pp=BASE):
    z = np.asarray(z, dtype=float)
    return pp.L_floor + (1-pp.L_floor)/(1+np.exp(pp.s*(z-pp.zc)))

def B_exp(z, p=COSMO):
    z = np.asarray(z, dtype=float)
    return E(z,p)*(1+z)

def D_acc(z, p=COSMO, pp=BASE):
    return D_growth(z,p) * L_latent(z,pp)

def S_proj(z_arr, p=COSMO, pp=BASE):
    z_arr = np.asarray(z_arr, dtype=float)
    be = B_exp(z_arr,p)
    b0 = float(B_exp(0.0,p))
    bnorm = (be/b0)**pp.alpha
    q_alpha = np.tanh(bnorm/np.max(bnorm))
    return be * q_alpha

def C_proj(z_arr, p=COSMO, pp=BASE):
    z_arr = np.asarray(z_arr, dtype=float)
    return S_proj(z_arr,p,pp) / (1 + pp.gamma*D_acc(z_arr,p,pp))

def R_proj(z_arr, p=COSMO, pp=BASE):
    """v0.11-consistent: |S_proj - gamma*D_acc| / (S_proj + gamma*D_acc + eps)"""
    sp = S_proj(z_arr,p,pp)
    da = pp.gamma * D_acc(z_arr,p,pp)
    return np.abs(sp - da) / (sp + da + EPS)

def curvature_abs(y, x):
    return np.abs(np.gradient(np.gradient(y,x),x))

def z_transition(z_arr, p=COSMO, pp=BASE):
    rp = R_proj(z_arr,p,pp)
    curv = curvature_abs(rp, z_arr)
    i1 = int(pp.edge_frac*len(z_arr))
    i2 = int((1-pp.edge_frac)*len(z_arr))
    local_idx = int(np.argmax(curv[i1:i2]))
    idx = i1 + local_idx
    return float(z_arr[idx]), idx

def chi2(y_obs, y_model, y_err):
    return float(np.sum(((y_obs-y_model)/y_err)**2))

def weighted_best_scale(y_obs, template, y_err):
    w = 1.0/np.maximum(y_err,EPS)**2
    return float(np.sum(w*y_obs*template)/np.sum(w*template**2))

FSIG8_DATA = np.array([
    [0.067,0.423,0.055],[0.17,0.510,0.060],[0.22,0.420,0.070],
    [0.25,0.3512,0.0583],[0.37,0.4602,0.0378],[0.41,0.450,0.040],
    [0.44,0.413,0.080],[0.57,0.441,0.043],[0.60,0.390,0.063],
    [0.60,0.550,0.120],[0.78,0.380,0.040],[0.80,0.470,0.080],
    [1.40,0.482,0.116],
], dtype=float)

HZ_DATA = np.array([
    [0.07,69.0,19.6],[0.09,69.0,12.0],[0.12,68.6,26.2],
    [0.17,83.0,8.0],[0.179,75.0,4.0],[0.199,75.0,5.0],
    [0.20,72.9,29.6],[0.27,77.0,14.0],[0.28,88.8,36.6],
    [0.352,83.0,14.0],[0.38,83.0,13.5],[0.40,95.0,17.0],
    [0.40,77.0,10.2],[0.42,87.1,11.2],[0.45,92.8,12.9],
    [0.47,89.0,50.0],[0.48,80.9,9.0],[0.48,97.0,62.0],
    [0.59,104.0,13.0],[0.68,92.0,8.0],[0.78,105.0,12.0],
    [0.88,125.0,17.0],[0.88,90.0,40.0],[0.90,117.0,23.0],
    [1.04,154.0,20.0],[1.30,168.0,17.0],[1.36,160.0,33.6],
    [1.43,177.0,18.0],[1.53,140.0,14.0],[1.75,202.0,40.0],
    [1.97,186.5,50.4],
], dtype=float)

def main():
    print("="*72)
    print("MAAT Paper 34 — Full Reproduction v2")
    print("="*72)

    z_arr = np.linspace(0.01, 3.0, 2000)
    c = C_proj(z_arr)
    c_norm = c / c[0]
    rp = R_proj(z_arr)
    curv = curvature_abs(rp, z_arr)
    ztr, idx_tr = z_transition(z_arr)

    print(f"\n--- Baseline (gamma={BASE.gamma}, alpha={BASE.alpha}) ---")
    print(f"z_tr              = {ztr:.5f}")
    print(f"C_norm(z_tr)      = {c_norm[idx_tr]:.5f}")
    print(f"R_proj(z_tr)      = {rp[idx_tr]:.5e}")
    print(f"C_norm(z=3)       = {c_norm[-1]:.5f}")
    print(f"Note: z_tr differs from Paper 33 (~0.85) because R_proj")
    print(f"now correctly uses gamma*D_acc (v0.11-consistent residual)")

    # fσ8 diagnostic
    z_obs=FSIG8_DATA[:,0]; y_obs=FSIG8_DATA[:,1]; y_err=FSIG8_DATA[:,2]
    y_lcdm = fsigma8_lcdm(z_obs)
    chi2_lcdm = chi2(y_obs, y_lcdm, y_err)
    dof_lcdm = len(z_obs)

    interp_c = interp1d(z_arr, c_norm, kind="linear", fill_value="extrapolate")
    template = interp_c(z_obs)
    scale = weighted_best_scale(y_obs, template, y_err)
    y_cproj = scale * template
    chi2_cproj = chi2(y_obs, y_cproj, y_err)
    dof_cproj = len(z_obs) - 1

    print(f"\n--- Diagnostic: C_proj vs f*sigma8 ({len(z_obs)} points) ---")
    print(f"LCDM:   chi2/dof = {chi2_lcdm/dof_lcdm:.5f}")
    print(f"C_proj: chi2/dof = {chi2_cproj/dof_cproj:.5f}")
    print(f"Scale (WLS):      = {scale:.6f}")
    print(f"Interpretation: large chi2 expected — C_proj != f*sigma8")

    # H(z) reference
    z_h=HZ_DATA[:,0]; H_obs=HZ_DATA[:,1]; H_err=HZ_DATA[:,2]
    H_th = COSMO.H0 * E(z_h)
    chi2_Hz = chi2(H_obs, H_th, H_err)
    dof_Hz = len(z_h)
    print(f"\n--- H(z) reference ({len(z_h)} points) ---")
    print(f"LCDM: chi2/dof = {chi2_Hz/dof_Hz:.5f}")

    # Parameter scan 25x25
    print(f"\n--- Parameter scan: gamma x alpha (25x25 = 625 points) ---")
    gammas = np.linspace(0.5, 5.0, 25)
    alphas = np.linspace(0.5, 2.0, 25)
    ztr_grid = np.zeros((len(gammas), len(alphas)))
    for i, gam in enumerate(gammas):
        for j, alp in enumerate(alphas):
            pp = ProjectionParams(gamma=float(gam), alpha=float(alp),
                                  zc=BASE.zc, s=BASE.s,
                                  L_floor=BASE.L_floor, edge_frac=BASE.edge_frac)
            ztr_grid[i,j], _ = z_transition(z_arr, COSMO, pp)

    print(f"z_tr min    = {np.min(ztr_grid):.5f}")
    print(f"z_tr max    = {np.max(ztr_grid):.5f}")
    print(f"z_tr median = {np.median(ztr_grid):.5f}")
    print(f"z_tr std    = {np.std(ztr_grid):.5f}")

    # Export CSVs
    core = np.column_stack([z_arr, E(z_arr), D_growth(z_arr), L_latent(z_arr),
                            B_exp(z_arr), D_acc(z_arr), S_proj(z_arr),
                            C_proj(z_arr), c_norm, rp, curv])
    np.savetxt(OUTDIR/"paper34_core_curves.csv", core, delimiter=",",
               header="z,E,D_growth,L_latent,B_exp,D_acc,S_proj,C_proj,C_proj_norm,R_proj,abs_d2R_dz2",
               comments="")

    comp = np.column_stack([z_obs, y_obs, y_err, y_lcdm, y_cproj])
    np.savetxt(OUTDIR/"paper34_fsigma8_comparison.csv", comp, delimiter=",",
               header="z,fsigma8_obs,sigma,fsigma8_LCDM,Cproj_rescaled", comments="")

    rows = [(g,a,ztr_grid[i,j]) for i,g in enumerate(gammas)
                                  for j,a in enumerate(alphas)]
    np.savetxt(OUTDIR/"paper34_scan_gamma_alpha.csv", np.array(rows), delimiter=",",
               header="gamma,alpha,z_tr", comments="")

    # Figures
    fig, axes = plt.subplots(1, 2, figsize=(13,5))
    ax = axes[0]
    ax.errorbar(z_obs, y_obs, yerr=y_err, fmt="o", color="steelblue",
                capsize=3, label="fσ₈ data", zorder=3)
    zp = np.linspace(0.05, 1.5, 400)
    ax.plot(zp, fsigma8_lcdm(zp), "k-", lw=1.8, label="ΛCDM fσ₈")
    ax.plot(zp, interp_c(zp)*scale, "r--", lw=2.0, label="C_proj (WLS rescaled)")
    ax.axvline(ztr, color="darkorange", ls=":", lw=2.0, label=f"z_tr={ztr:.3f}")
    ax.set_xlabel("z"); ax.set_ylabel("fσ₈(z)")
    ax.set_title("Projection diagnostic vs growth data")
    ax.legend(fontsize=9); ax.grid(alpha=0.25)

    ax2 = axes[1]
    ax2.plot(z_arr, rp, "b-", lw=2.0, label="R_proj(z)")
    ax2r = ax2.twinx()
    ax2r.plot(z_arr, curv, "r-", lw=1.5, alpha=0.72, label="|d²R/dz²|")
    ax2.axvline(ztr, color="darkorange", ls=":", lw=2.0, label=f"z_tr={ztr:.3f}")
    ax2.set_xlabel("z"); ax2.set_ylabel("R_proj(z)", color="b")
    ax2r.set_ylabel("|d²R_proj/dz²|", color="r")
    ax2.set_title("Transition marker (v0.11 residual)")
    lines1,labs1 = ax2.get_legend_handles_labels()
    lines2,labs2 = ax2r.get_legend_handles_labels()
    ax2.legend(lines1+lines2, labs1+labs2, fontsize=9, loc="upper left")
    ax2.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUTDIR/"fig_paper34_main.png", bbox_inches="tight")
    plt.close(fig)

    # Components figure
    figc, axc = plt.subplots(figsize=(8,5))
    axc.plot(z_arr, B_exp(z_arr)/B_exp(z_arr)[0], lw=2, label="B_exp / B_exp(0)")
    axc.plot(z_arr, D_acc(z_arr)/D_acc(z_arr)[0], lw=2, label="D_acc / D_acc(0)")
    axc.plot(z_arr, c_norm, lw=2, label="C_proj / C_proj(0)")
    axc.axvline(ztr, color="darkorange", ls=":", lw=2.0, label=f"z_tr={ztr:.3f}")
    axc.set_xlabel("z"); axc.set_ylabel("normalised value")
    axc.set_title("Projection components")
    axc.legend(fontsize=9); axc.grid(alpha=0.25)
    figc.tight_layout()
    figc.savefig(OUTDIR/"fig_paper34_components.png", bbox_inches="tight")
    plt.close(figc)

    # Scan heatmap
    fig2, ax = plt.subplots(figsize=(7,5))
    cf = ax.contourf(alphas, gammas, ztr_grid, levels=22, cmap="RdYlBu_r")
    plt.colorbar(cf, ax=ax, label="z_tr")
    ax.axhline(BASE.gamma, color="white", ls="--", lw=1.8, label="baseline γ=2")
    ax.axvline(BASE.alpha, color="white", ls=":", lw=1.8, label="baseline α=1")
    ax.set_xlabel("α"); ax.set_ylabel("γ")
    ax.set_title("Transition redshift z_tr(γ, α) — v0.11")
    ax.legend(fontsize=9)
    fig2.tight_layout()
    fig2.savefig(OUTDIR/"fig_paper34_scan.png", bbox_inches="tight")
    plt.close(fig2)

    # JSON summary
    summary = {
        "cosmology": asdict(COSMO),
        "projection_baseline": asdict(BASE),
        "baseline_results": {
            "z_tr": ztr, "C_norm_z_tr": float(c_norm[idx_tr]),
            "R_proj_z_tr": float(rp[idx_tr]),
            "curvature_z_tr": float(curv[idx_tr]),
            "C_norm_z3": float(c_norm[-1])
        },
        "fsigma8_diagnostic": {
            "n_points": int(len(z_obs)),
            "chi2_lcdm": float(chi2_lcdm),
            "dof_lcdm": int(dof_lcdm),
            "chi2_dof_lcdm": float(chi2_lcdm/dof_lcdm),
            "chi2_cproj": float(chi2_cproj),
            "dof_cproj": int(dof_cproj),
            "chi2_dof_cproj": float(chi2_cproj/dof_cproj),
            "scale_wls": float(scale)
        },
        "Hz_reference": {
            "n_points": int(len(z_h)),
            "chi2_lcdm": float(chi2_Hz),
            "dof_lcdm": int(dof_Hz),
            "chi2_dof_lcdm": float(chi2_Hz/dof_Hz)
        },
        "scan_gamma_alpha": {
            "n_points": int(ztr_grid.size),
            "z_tr_min": float(np.min(ztr_grid)),
            "z_tr_max": float(np.max(ztr_grid)),
            "z_tr_median": float(np.median(ztr_grid)),
            "z_tr_mean": float(np.mean(ztr_grid)),
            "z_tr_std": float(np.std(ztr_grid))
        },
        "interpretation_note": (
            "C_proj is a projection-stress observable, not a physical f_sigma8 "
            "growth model. z_tr=1.044 (vs Paper 33 ~0.85) reflects the v0.11-"
            "consistent R_proj formulation using gamma*D_acc in the residual."
        )
    }
    with open(OUTDIR/"paper34_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*72}")
    print("PAPER 34 KEY RESULTS")
    print(f"{'='*72}")
    print(f"Baseline z_tr            = {ztr:.5f}")
    print(f"chi2/dof ΛCDM fσ₈        = {chi2_lcdm/dof_lcdm:.5f}")
    print(f"chi2/dof C_proj fσ₈      = {chi2_cproj/dof_cproj:.5f}")
    print(f"chi2/dof ΛCDM H(z)       = {chi2_Hz/dof_Hz:.5f}")
    print(f"z_tr scan median         = {np.median(ztr_grid):.5f}")
    print(f"z_tr scan range          = [{np.min(ztr_grid):.3f}, {np.max(ztr_grid):.3f}]")
    print(f"Outputs: {OUTDIR.resolve()}")
    print("Run: python3 maat_paper34_reproduction_v2.py")

if __name__ == "__main__":
    main()
