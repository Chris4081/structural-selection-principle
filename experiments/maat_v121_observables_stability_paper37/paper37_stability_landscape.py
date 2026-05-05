#!/usr/bin/env python3
# MAAT Paper 37 — Stability & Robustness Landscape of MAAT v1.2.1
# ----------------------------------------------------------------
# Purpose:
#   Systematically scan projection strength eta_proj and damping gamma_proj
#   to test whether MAAT v1.2.1 is robust or fine-tuned.
#
# Outputs:
#   - CSV with scan results
#   - stability heatmap
#   - robustness heatmap
#   - max ΔH/H heatmap
#   - max Δfσ8 heatmap
#
# Scientific status:
#   Toy/proxy cosmological scan. Not a Boltzmann-code calculation.

from pathlib import Path
import os
import json
import numpy as np
import pandas as pd

if "MPLCONFIGDIR" not in os.environ:
    mpl_cache = Path(".matplotlib-cache")
    mpl_cache.mkdir(exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(mpl_cache.resolve())

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR = Path("stability_landscape_outputs")
OUTDIR.mkdir(exist_ok=True)

# -----------------------------
# Baseline cosmology
# -----------------------------
rho_m0 = 0.90
rho_r0 = 3.0e-4
rho_L = 2.10

# MAAT field parameters
epsilon = 1.0e-6
mu = 27.8256
m_lambda = 0.55
Z_lambda = 1.0

lambda_H = 1.0
lambda_B = 1.0
lambda_S = 1.0
lambda_V = 1.0

kappa_struct = 0.35

# Integration controls
a0 = 0.30
a1 = 1.00
phi_dot0 = 2.37815
lambda0 = 0.05
lambda_dot0 = 0.0
N = 4000
dt = 0.002

# Scan grid
eta_vals = np.linspace(0.00, 0.16, 33)
gamma_vals = np.linspace(0.05, 1.20, 33)
OMEGA_TRUST_THRESHOLD = 0.02
W_MAX_TRUST = 5.0


def safe_clip01(x):
    return np.clip(x, 0.0, 1.0)


def U_of_X(X):
    Gamma = 1.0 / (1.0 + X)
    return -np.log((epsilon + Gamma) / (epsilon + 1.0))


def U_X(X):
    Gamma = 1.0 / (1.0 + X)
    dGamma = -1.0 / (1.0 + X) ** 2
    return -(1.0 / (epsilon + Gamma)) * dGamma


def U_XX(X):
    h = 1e-5 * max(1.0, abs(X))
    return (U_X(X + h) - U_X(max(X - h, 0.0))) / (2.0 * h)


def preliminary_projection(a):
    z = 1.0 / a - 1.0
    return 1.0 / (1.0 + 0.25 * z)


def structural_fields(a, phi_dot, lam, lam_dot):
    X = 0.5 * phi_dot**2
    U = U_of_X(X)

    activity = abs(phi_dot) + abs(lam_dot) + abs(lam)
    normalized_activity = activity / (1.0 + activity)

    d_H = U / (1.0 + U)
    H = 1.0 / (1.0 + d_H)

    balance_load = X + 0.5 * lam_dot**2 + 0.5 * m_lambda**2 * lam**2
    d_B = balance_load / (1.0 + balance_load)
    B = 1.0 / (1.0 + d_B)

    S = 4.0 * normalized_activity * (1.0 - normalized_activity)
    S = safe_clip01(S)

    coupling = abs(lam * phi_dot)
    coupling_coherence = coupling / (1.0 + coupling)
    smoothness = (H * B) ** 0.5
    accessibility = preliminary_projection(a)

    V = safe_clip01(0.45 * smoothness + 0.35 * coupling_coherence + 0.20 * accessibility)

    respect_implicit = (H * B * V) ** (1.0 / 3.0)
    base_robustness = (H * B * S * V) ** 0.25
    R_emergent = min(respect_implicit, base_robustness)

    # Linear diagnostic proxy only. The canonical MAAT functional in Paper 36
    # uses log-supports; this proxy is kept for lightweight monitoring.
    F_struct_proxy = (
        lambda_H * (1.0 - H)
        + lambda_B * (1.0 - B)
        + lambda_S * (1.0 - S)
        + lambda_V * (1.0 - V)
    )

    return {
        "X": X,
        "U": U,
        "H_struct": H,
        "B_struct": B,
        "S_struct": S,
        "V_struct": V,
        "Respect_implicit": respect_implicit,
        "R_emergent": R_emergent,
        "F_struct_proxy": F_struct_proxy,
    }


def projection_layer(a, s, eta_proj, gamma_proj):
    z = 1.0 / a - 1.0
    E_proxy = np.sqrt((rho_m0 / a**3 + rho_r0 / a**4 + rho_L) / 3.0)
    B_exp = E_proxy * (1.0 + z)
    D_acc = ((s["H_struct"] * s["B_struct"] * s["V_struct"] * s["R_emergent"]) ** 0.25) * (1.0 / (1.0 + z))
    C_proj = np.tanh(B_exp) / (1.0 + gamma_proj * D_acc)
    C_proj = safe_clip01(C_proj)
    F_proj = eta_proj * (1.0 - C_proj) ** 2
    return C_proj, F_proj, B_exp, D_acc


def densities(a, phi_dot, lam, lam_dot, eta_proj, gamma_proj):
    s = structural_fields(a, phi_dot, lam, lam_dot)
    X = s["X"]
    U = s["U"]
    UX = U_X(X)

    rho_m = rho_m0 / a**3
    rho_r = rho_r0 / a**4

    rho_raw = (
        0.5 * Z_lambda * lam_dot**2
        + 0.5 * m_lambda**2 * lam**2
        + mu * lam * (2.0 * X * UX - U)
    )

    p_raw = (
        0.5 * Z_lambda * lam_dot**2
        - 0.5 * m_lambda**2 * lam**2
        + mu * lam * U
    )

    C_proj, F_proj, B_exp, D_acc = projection_layer(a, s, eta_proj, gamma_proj)

    selection_factor = 1.0 + kappa_struct * (1.0 - s["R_emergent"]) + F_proj

    rho_maat = rho_raw * selection_factor
    p_maat = p_raw * selection_factor - F_proj * abs(rho_raw) / (abs(rho_raw) + 1.0)

    rho_total = rho_m + rho_r + rho_L + rho_maat
    Hubble = np.sqrt(max(rho_total, 1e-12) / 3.0)

    extra = {
        **s,
        "C_proj": C_proj,
        "F_proj": F_proj,
        "B_exp": B_exp,
        "D_acc": D_acc,
        "rho_raw_MAAT": rho_raw,
        "p_raw_MAAT": p_raw,
        "selection_factor": selection_factor,
    }

    return X, rho_m, rho_r, rho_maat, p_maat, Hubble, extra


def rhs(y, eta_proj, gamma_proj):
    a, phi, phi_dot, lam, lam_dot = y
    X, rho_m, rho_r, rho_maat, p_maat, Hubble, extra = densities(a, phi_dot, lam, lam_dot, eta_proj, gamma_proj)

    UX = U_X(X)
    UXX = U_XX(X)

    PX = 1.0 + mu * lam * UX
    denom = PX + 2.0 * X * mu * lam * UXX

    if abs(denom) < 1e-10:
        denom = np.sign(denom) * 1e-10 if denom != 0 else 1e-10

    a_dot = a * Hubble
    phi_ddot = -(mu * UX * lam_dot * phi_dot + 3.0 * Hubble * PX * phi_dot) / denom
    lam_ddot = (mu * U_of_X(X) - 3.0 * Hubble * Z_lambda * lam_dot - m_lambda**2 * lam) / Z_lambda

    return np.array([a_dot, phi_dot, phi_ddot, lam_dot, lam_ddot], dtype=float)


def rk4_step(y, dt, eta_proj, gamma_proj):
    k1 = rhs(y, eta_proj, gamma_proj)
    k2 = rhs(y + 0.5 * dt * k1, eta_proj, gamma_proj)
    k3 = rhs(y + 0.5 * dt * k2, eta_proj, gamma_proj)
    k4 = rhs(y + dt * k3, eta_proj, gamma_proj)
    return y + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0


def simulate_one(eta_proj, gamma_proj):
    y = np.array([a0, 0.0, phi_dot0, lambda0, lambda_dot0], dtype=float)
    rows = []
    t = 0.0

    stable = True

    for _ in range(N):
        a, phi, phi_dot, lam, lam_dot = y
        if a >= a1:
            break

        X, rho_m, rho_r, rho_maat, p_maat, Hubble, extra = densities(a, phi_dot, lam, lam_dot, eta_proj, gamma_proj)

        H_lcdm = np.sqrt((rho_m + rho_r + rho_L) / 3.0)
        Omega_maat = rho_maat / (rho_m + rho_r + rho_L + rho_maat)
        w_raw = p_maat / rho_maat if abs(rho_maat) > 1e-12 else np.nan

        UX = U_X(X)
        UXX = U_XX(X)
        PX = 1.0 + mu * lam * UX
        denom = PX + 2.0 * X * mu * lam * UXX
        cs2 = PX / denom if abs(denom) > 1e-12 else np.nan

        z = 1.0 / a - 1.0

        rows.append({
            "z": z,
            "a": a,
            "H_MAAT": Hubble,
            "H_LCDM": H_lcdm,
            "delta_H": (Hubble - H_lcdm) / H_lcdm,
            "rho_MAAT": rho_maat,
            "Omega_MAAT": Omega_maat,
            "w_raw": w_raw,
            "P_X": PX,
            "cs2": cs2,
            "H_struct": extra["H_struct"],
            "B_struct": extra["B_struct"],
            "S_struct": extra["S_struct"],
            "V_struct": extra["V_struct"],
            "Respect_implicit": extra["Respect_implicit"],
            "R_emergent": extra["R_emergent"],
            "C_proj": extra["C_proj"],
            "F_proj": extra["F_proj"],
        })

        y = rk4_step(y, dt, eta_proj, gamma_proj)
        t += dt

        if not np.all(np.isfinite(y)):
            stable = False
            break

    df = pd.DataFrame(rows)

    if len(df) < 50:
        stable = False
        return None, stable

    df = df.sort_values("z")

    rho_m = rho_m0 * (1.0 + df["z"]) ** 3
    sigma8_0 = 0.811
    gamma_growth = 0.55
    Om_maat = rho_m / (3.0 * df["H_MAAT"] ** 2)
    Om_lcdm = rho_m / (3.0 * df["H_LCDM"] ** 2)
    f_maat = Om_maat.clip(lower=0) ** gamma_growth
    f_lcdm = Om_lcdm.clip(lower=0) ** gamma_growth
    D = 1.0 / (1.0 + df["z"])

    df["fsigma8_MAAT_proxy"] = f_maat * sigma8_0 * D
    df["fsigma8_LCDM_proxy"] = f_lcdm * sigma8_0 * D
    df["delta_fsigma8"] = (df["fsigma8_MAAT_proxy"] - df["fsigma8_LCDM_proxy"]) / df["fsigma8_LCDM_proxy"]

    return df, stable


def evaluate_model(df, stable, eta_proj, gamma_proj):
    if df is None or not stable:
        return {
            "eta_proj": eta_proj,
            "gamma_proj": gamma_proj,
            "stable": False,
        }

    trusted_w_mask = (
        np.abs(df["Omega_MAAT"]) >= OMEGA_TRUST_THRESHOLD
    ) & np.isfinite(df["w_raw"]) & (np.abs(df["w_raw"]) <= W_MAX_TRUST)

    trusted_w = df.loc[trusted_w_mask, "w_raw"]

    max_abs_delta_H = float(np.nanmax(np.abs(df["delta_H"])))
    max_abs_delta_fsigma8 = float(np.nanmax(np.abs(df["delta_fsigma8"])))
    max_Omega_MAAT = float(np.nanmax(df["Omega_MAAT"]))
    min_cs2 = float(np.nanmin(df["cs2"]))
    max_cs2 = float(np.nanmax(df["cs2"]))
    min_P_X = float(np.nanmin(df["P_X"]))

    mean_H = float(np.nanmean(df["H_struct"]))
    mean_B = float(np.nanmean(df["B_struct"]))
    mean_S = float(np.nanmean(df["S_struct"]))
    mean_V = float(np.nanmean(df["V_struct"]))
    mean_respect = float(np.nanmean(df["Respect_implicit"]))
    mean_R = float(np.nanmean(df["R_emergent"]))
    mean_C_proj = float(np.nanmean(df["C_proj"]))

    physical_stable = (
        min_P_X > 0
        and min_cs2 > 0
        and np.isfinite(max_abs_delta_H)
        and max_abs_delta_H < 0.08
        and max_abs_delta_fsigma8 < 0.08
        and max_Omega_MAAT < 0.10
    )

    # Diagnostic scan score only. This is not a new v1.2.1 selection
    # functional; it ranks stable proxy runs by small observable deviations
    # and high derived structural supports.
    obs_score = max(0.0, 1.0 - max_abs_delta_H / 0.08) * max(0.0, 1.0 - max_abs_delta_fsigma8 / 0.08)
    structure_score = (mean_H * mean_B * mean_S * mean_V * mean_respect * mean_R) ** (1.0 / 6.0)
    scan_diagnostic_score = obs_score * structure_score

    return {
        "eta_proj": eta_proj,
        "gamma_proj": gamma_proj,
        "stable": bool(physical_stable),
        "max_abs_delta_H": max_abs_delta_H,
        "max_abs_delta_fsigma8": max_abs_delta_fsigma8,
        "max_Omega_MAAT": max_Omega_MAAT,
        "min_cs2": min_cs2,
        "max_cs2": max_cs2,
        "min_P_X": min_P_X,
        "mean_H_struct": mean_H,
        "mean_B_struct": mean_B,
        "mean_S_struct": mean_S,
        "mean_V_struct": mean_V,
        "mean_Respect_implicit": mean_respect,
        "mean_R_emergent": mean_R,
        "mean_C_proj": mean_C_proj,
        "trusted_w_fraction": float(trusted_w_mask.mean()),
        "trusted_w_min": float(np.nanmin(trusted_w)) if len(trusted_w) else np.nan,
        "trusted_w_max": float(np.nanmax(trusted_w)) if len(trusted_w) else np.nan,
        "scan_diagnostic_score": float(scan_diagnostic_score),
    }


def run_scan():
    results = []
    total = len(eta_vals) * len(gamma_vals)
    count = 0

    for eta_proj in eta_vals:
        for gamma_proj in gamma_vals:
            count += 1
            df, stable = simulate_one(eta_proj, gamma_proj)
            result = evaluate_model(df, stable, eta_proj, gamma_proj)
            results.append(result)
            if count % 50 == 0:
                print(f"Progress: {count}/{total}")

    scan = pd.DataFrame(results)
    scan.to_csv(OUTDIR / "paper37_maat_v121_landscape_scan.csv", index=False)
    return scan


def heatmap(scan, value_col, title, filename, cmap="viridis"):
    pivot = scan.pivot(index="gamma_proj", columns="eta_proj", values=value_col)
    plt.figure(figsize=(9, 6))
    plt.imshow(
        pivot.values,
        origin="lower",
        aspect="auto",
        extent=[eta_vals.min(), eta_vals.max(), gamma_vals.min(), gamma_vals.max()],
        cmap=cmap,
    )
    plt.colorbar(label=value_col)
    plt.xlabel("eta_proj")
    plt.ylabel("gamma_proj")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(OUTDIR / filename, dpi=240)
    plt.close()


def save_plots(scan):
    scan_plot = scan.copy()
    scan_plot["stable_numeric"] = scan_plot["stable"].astype(float)

    heatmap(
        scan_plot,
        "stable_numeric",
        "MAAT v1.2.1 Stability Map",
        "fig1_stability_map.png",
        cmap="Greens",
    )

    heatmap(
        scan_plot,
        "scan_diagnostic_score",
        "MAAT v1.2.1 Diagnostic Scan Score",
        "fig2_robustness_heatmap.png",
        cmap="viridis",
    )

    heatmap(
        scan_plot,
        "max_abs_delta_H",
        "Maximum Relative Expansion Deviation",
        "fig3_max_delta_H_heatmap.png",
        cmap="magma",
    )

    heatmap(
        scan_plot,
        "max_abs_delta_fsigma8",
        "Maximum Growth Proxy Deviation",
        "fig4_max_delta_fsigma8_heatmap.png",
        cmap="magma",
    )

    heatmap(
        scan_plot,
        "mean_R_emergent",
        "Mean Emergent Robustness",
        "fig5_mean_R_emergent_heatmap.png",
        cmap="plasma",
    )

    heatmap(
        scan_plot,
        "trusted_w_fraction",
        "Trusted w Diagnostic Fraction",
        "fig6_trusted_w_fraction_heatmap.png",
        cmap="cividis",
    )


def save_summary(scan):
    stable_scan = scan[scan["stable"] == True].copy()

    best = None
    if len(stable_scan):
        best_row = stable_scan.sort_values("scan_diagnostic_score", ascending=False).iloc[0]
        best = best_row.to_dict()

    summary = {
        "model": "MAAT Paper 37 Stability & Robustness Landscape v1.2.1",
        "n_models": int(len(scan)),
        "n_stable": int(scan["stable"].sum()) if "stable" in scan else 0,
        "stable_fraction": float(scan["stable"].mean()) if "stable" in scan else 0.0,
        "eta_range": [float(eta_vals.min()), float(eta_vals.max())],
        "gamma_range": [float(gamma_vals.min()), float(gamma_vals.max())],
        "best_model": best,
        "note": "Toy/proxy parameter scan. Growth is not Boltzmann-code based. The scan diagnostic score is a ranking aid, not a proof."
    }

    with open(OUTDIR / "paper37_landscape_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    return summary


def main():
    print("=== Paper 37: MAAT v1.2.1 Stability & Robustness Landscape ===")
    scan = run_scan()
    save_plots(scan)
    summary = save_summary(scan)

    print("\nSummary:")
    for k, v in summary.items():
        print(f"{k}: {v}")

    print(f"\nSaved outputs to: {OUTDIR.resolve()}")
    print("Figures:")
    print(" - fig1_stability_map.png")
    print(" - fig2_robustness_heatmap.png")
    print(" - fig3_max_delta_H_heatmap.png")
    print(" - fig4_max_delta_fsigma8_heatmap.png")
    print(" - fig5_mean_R_emergent_heatmap.png")
    print(" - fig6_trusted_w_fraction_heatmap.png")


if __name__ == "__main__":
    main()
