#!/usr/bin/env python3
# MAAT v1.2.1 — Observable Predictions Layer with Emergent Robustness
# --------------------------------------------------------------------
# Safe rewritten version of the original MAAT v0.10 observable script.
#
# Conceptual update:
#   - H, B, S, V are primary structural fields.
#   - Respect is implicit, inferred from H, B, V.
#   - Robustness is emergent, inferred from H, B, S, V and implicit respect.
#   - MAAT remains a subleading structural correction close to ΛCDM.
#   - w_MAAT is treated as a diagnostic, not as a classical fluid equation of state.
#
# This script DOES NOT overwrite the v0.10 output directory.

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

OUTDIR = Path("observable_outputs")
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

# Structural sector weights
lambda_H = 1.0
lambda_B = 1.0
lambda_S = 1.0
lambda_V = 1.0

# Structural correction strength
kappa_struct = 0.35

# Projection-layer parameters
eta_proj = 0.08
gamma_proj = 0.35

# Diagnostic thresholds
OMEGA_TRUST_THRESHOLD = 0.02
W_MAX_TRUST = 5.0
W_CLIP = 2.0

# Initial state
# y = [a, phi, phi_dot, lambda, lambda_dot]
a0 = 0.30
a1 = 1.00
phi_dot0 = 2.37815
lambda0 = 0.05
lambda_dot0 = 0.0

# Integration controls
N = 4000
dt = 0.002


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
    """Simple redshift accessibility prior."""
    z = 1.0 / a - 1.0
    return 1.0 / (1.0 + 0.25 * z)


def structural_fields(a, phi_dot, lam, lam_dot):
    """
    Primary MAAT fields:
      H = Harmony / consistency proxy
      B = Balance / energetic control proxy
      S = Creative activity / productive dynamics
      V = Connectedness / coupling + smoothness + accessibility

    Respect is implicit:
      Respect_implicit = (H * B * V)^(1/3)

    Robustness is emergent:
      R_emergent = min(Respect_implicit, (H * B * S * V)^(1/4))
    """
    X = 0.5 * phi_dot**2
    U = U_of_X(X)

    activity = abs(phi_dot) + abs(lam_dot) + abs(lam)
    normalized_activity = activity / (1.0 + activity)

    # H: bounded cost consistency
    d_H = U / (1.0 + U)
    H = 1.0 / (1.0 + d_H)

    # B: penalizes excessive kinetic/lambda load
    balance_load = X + 0.5 * lam_dot**2 + 0.5 * m_lambda**2 * lam**2
    d_B = balance_load / (1.0 + balance_load)
    B = 1.0 / (1.0 + d_B)

    # S: activity should be present but not runaway
    S = 4.0 * normalized_activity * (1.0 - normalized_activity)
    S = safe_clip01(S)

    # V: improved connectedness sector
    coupling = abs(lam * phi_dot)
    coupling_component = coupling / (1.0 + coupling)
    smoothness_component = (H * B) ** 0.5
    accessibility_component = preliminary_projection(a)

    V = safe_clip01(
        0.45 * smoothness_component
        + 0.35 * coupling_component
        + 0.20 * accessibility_component
    )

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

    # CCI compatible diagnostic: activity over coherence/stabilizers
    CCI_v121 = S / (H + B + V + respect_implicit + R_emergent + epsilon)

    return {
        "X": X,
        "U": U,
        "H_struct": H,
        "B_struct": B,
        "S_struct": S,
        "V_struct": V,
        "V_coupling": coupling_component,
        "V_smoothness": smoothness_component,
        "V_accessibility": accessibility_component,
        "Respect_implicit": respect_implicit,
        "R_emergent": R_emergent,
        "F_struct_proxy": F_struct_proxy,
        "CCI_v121": CCI_v121,
    }


def projection_layer(a, s):
    z = 1.0 / a - 1.0
    E_proxy = np.sqrt((rho_m0 / a**3 + rho_r0 / a**4 + rho_L) / 3.0)
    B_exp = E_proxy * (1.0 + z)

    D_acc = ((s["H_struct"] * s["B_struct"] * s["V_struct"] * s["R_emergent"]) ** 0.25) * (1.0 / (1.0 + z))
    C_proj = np.tanh(B_exp) / (1.0 + gamma_proj * D_acc)
    C_proj = safe_clip01(C_proj)

    F_proj = eta_proj * (1.0 - C_proj) ** 2
    return C_proj, F_proj, B_exp, D_acc


def densities(a, phi_dot, lam, lam_dot):
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

    C_proj, F_proj, B_exp, D_acc = projection_layer(a, s)

    selection_factor = 1.0 + kappa_struct * (1.0 - s["R_emergent"]) + F_proj

    rho_maat = rho_raw * selection_factor

    # Stabilized pressure proxy: avoids projection pressure dominating near rho_raw≈0.
    p_maat = p_raw * selection_factor - F_proj * abs(rho_raw) / (abs(rho_raw) + 1.0)

    rho_total = rho_m + rho_r + rho_L + rho_maat
    H = np.sqrt(max(rho_total, 1e-12) / 3.0)

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

    return X, rho_m, rho_r, rho_maat, p_maat, H, extra


def stable_w_diagnostic(p_maat, rho_maat, Omega_maat):
    """
    w_raw is kept for diagnostics.
    w_stable is trusted only if MAAT fraction is relevant and raw w is not pathological.
    MAAT is interpreted as a structural correction, not as a classical fluid.
    """
    w_raw = p_maat / rho_maat if abs(rho_maat) > 1e-12 else np.nan

    trusted = (
        abs(Omega_maat) >= OMEGA_TRUST_THRESHOLD
        and np.isfinite(w_raw)
        and abs(w_raw) <= W_MAX_TRUST
    )

    w_stable = w_raw if trusted else np.nan
    w_plot = np.clip(w_stable, -W_CLIP, W_CLIP) if np.isfinite(w_stable) else np.nan

    return w_raw, w_stable, w_plot, trusted


def rhs(y):
    a, phi, phi_dot, lam, lam_dot = y

    X, rho_m, rho_r, rho_maat, p_maat, H, extra = densities(a, phi_dot, lam, lam_dot)

    UX = U_X(X)
    UXX = U_XX(X)

    PX = 1.0 + mu * lam * UX
    denom = PX + 2.0 * X * mu * lam * UXX

    if abs(denom) < 1e-10:
        denom = np.sign(denom) * 1e-10 if denom != 0 else 1e-10

    a_dot = a * H
    phi_ddot = -(mu * UX * lam_dot * phi_dot + 3.0 * H * PX * phi_dot) / denom
    lam_ddot = (mu * U_of_X(X) - 3.0 * H * Z_lambda * lam_dot - m_lambda**2 * lam) / Z_lambda

    return np.array([a_dot, phi_dot, phi_ddot, lam_dot, lam_ddot], dtype=float)


def rk4_step(y, dt):
    k1 = rhs(y)
    k2 = rhs(y + 0.5 * dt * k1)
    k3 = rhs(y + 0.5 * dt * k2)
    k4 = rhs(y + dt * k3)
    return y + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0


def simulate():
    y = np.array([a0, 0.0, phi_dot0, lambda0, lambda_dot0], dtype=float)

    rows = []
    t = 0.0

    for _ in range(N):
        a, phi, phi_dot, lam, lam_dot = y

        if a >= a1:
            break

        X, rho_m, rho_r, rho_maat, p_maat, H, extra = densities(a, phi_dot, lam, lam_dot)

        H_lcdm = np.sqrt((rho_m + rho_r + rho_L) / 3.0)
        Omega_maat = rho_maat / (rho_m + rho_r + rho_L + rho_maat)

        w_raw, w_stable, w_plot, w_trusted = stable_w_diagnostic(p_maat, rho_maat, Omega_maat)

        UX = U_X(X)
        UXX = U_XX(X)
        PX = 1.0 + mu * lam * UX
        denom = PX + 2.0 * X * mu * lam * UXX
        cs2 = PX / denom if abs(denom) > 1e-12 else np.nan

        z = 1.0 / a - 1.0

        rows.append({
            "t": t,
            "a": a,
            "z": z,
            "phi_dot": phi_dot,
            "lambda": lam,
            "lambda_dot": lam_dot,
            "X": X,
            "H_MAAT": H,
            "H_LCDM": H_lcdm,
            "delta_H": (H - H_lcdm) / H_lcdm,
            "rho_MAAT": rho_maat,
            "rho_raw_MAAT": extra["rho_raw_MAAT"],
            "p_MAAT": p_maat,
            "p_raw_MAAT": extra["p_raw_MAAT"],
            "Omega_MAAT": Omega_maat,
            "w_MAAT_raw": w_raw,
            "w_MAAT_stable": w_stable,
            "w_MAAT_plot": w_plot,
            "w_MAAT_trusted": w_trusted,
            "P_X": PX,
            "cs2": cs2,
            "H_struct": extra["H_struct"],
            "B_struct": extra["B_struct"],
            "S_struct": extra["S_struct"],
            "V_struct": extra["V_struct"],
            "V_coupling": extra["V_coupling"],
            "V_smoothness": extra["V_smoothness"],
            "V_accessibility": extra["V_accessibility"],
            "Respect_implicit": extra["Respect_implicit"],
            "R_emergent": extra["R_emergent"],
            "F_struct_proxy": extra["F_struct_proxy"],
            "CCI_v121": extra["CCI_v121"],
            "C_proj": extra["C_proj"],
            "F_proj": extra["F_proj"],
            "B_exp": extra["B_exp"],
            "D_acc": extra["D_acc"],
            "selection_factor": extra["selection_factor"],
        })

        y = rk4_step(y, dt)
        t += dt

        if not np.all(np.isfinite(y)):
            break

    df = pd.DataFrame(rows)
    df = df.sort_values("z")
    return df


def growth_proxy(df):
    """
    Simple GR-like growth proxy using MAAT-modified expansion.
    This is not a Boltzmann-code calculation.
    """
    sigma8_0 = 0.811
    gamma = 0.55

    rho_m = rho_m0 * (1.0 + df["z"])**3

    Om_maat = rho_m / (3.0 * df["H_MAAT"]**2)
    Om_lcdm = rho_m / (3.0 * df["H_LCDM"]**2)

    f_maat = Om_maat.clip(lower=0)**gamma
    f_lcdm = Om_lcdm.clip(lower=0)**gamma

    D = 1.0 / (1.0 + df["z"])

    df["fsigma8_MAAT_proxy"] = f_maat * sigma8_0 * D
    df["fsigma8_LCDM_proxy"] = f_lcdm * sigma8_0 * D
    df["delta_fsigma8"] = (df["fsigma8_MAAT_proxy"] - df["fsigma8_LCDM_proxy"]) / df["fsigma8_LCDM_proxy"]

    return df


def save_single_plot(x, y, xlabel, ylabel, title, filename, label=None, y2=None, label2=None):
    plt.figure(figsize=(8, 5))
    plt.plot(x, y, label=label)
    if y2 is not None:
        plt.plot(x, y2, "--", label=label2)
    if label or label2:
        plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / filename, dpi=240)
    plt.close()


def save_plots(df):
    save_single_plot(
        df["z"], df["H_MAAT"],
        "z", "H(z) [dimensionless]",
        "MAAT v1.2.1 — Expansion History",
        "fig1_Hz_MAAT_vs_LCDM.png",
        label="MAAT v1.2.1", y2=df["H_LCDM"], label2="ΛCDM reference"
    )

    save_single_plot(
        df["z"], df["delta_H"],
        "z", "(H_MAAT - H_LCDM) / H_LCDM",
        "Relative Expansion Deviation",
        "fig2_delta_H.png"
    )

    save_single_plot(
        df["z"], df["Omega_MAAT"],
        "z", "Ω_MAAT(z)",
        "MAAT Density Fraction",
        "fig3_Omega_MAAT.png"
    )

    save_single_plot(
        df["z"], df["w_MAAT_plot"],
        "z", "w_MAAT(z), trusted and clipped to [-2,2]",
        "Effective Equation of State — Trusted Diagnostic",
        "fig4_w_MAAT_stable.png"
    )

    save_single_plot(
        df["z"], df["w_MAAT_raw"].clip(-10, 10),
        "z", "w_MAAT raw, clipped to [-10,10]",
        "Raw Equation of State Diagnostic",
        "fig4b_w_MAAT_raw.png"
    )

    save_single_plot(
        df["z"], df["fsigma8_MAAT_proxy"],
        "z", "fσ8(z) proxy",
        "Growth Observable Proxy",
        "fig5_fsigma8_proxy.png",
        label="MAAT proxy", y2=df["fsigma8_LCDM_proxy"], label2="ΛCDM proxy"
    )

    save_single_plot(
        df["z"], df["delta_fsigma8"],
        "z", "Δ fσ8 / fσ8",
        "Relative Growth Proxy Deviation",
        "fig6_delta_fsigma8.png"
    )

    plt.figure(figsize=(8, 5))
    plt.plot(df["z"], df["H_struct"], label="H: Harmony")
    plt.plot(df["z"], df["B_struct"], label="B: Balance")
    plt.plot(df["z"], df["S_struct"], label="S: Creative activity")
    plt.plot(df["z"], df["V_struct"], label="V: Connectedness")
    plt.xlabel("z")
    plt.ylabel("Structural value [0,1]")
    plt.title("Primary MAAT Structural Fields")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig7_structural_fields_HBSV.png", dpi=240)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(df["z"], df["V_struct"], label="V total")
    plt.plot(df["z"], df["V_coupling"], label="coupling component")
    plt.plot(df["z"], df["V_smoothness"], label="smoothness component")
    plt.plot(df["z"], df["V_accessibility"], label="accessibility component")
    plt.xlabel("z")
    plt.ylabel("Value [0,1]")
    plt.title("Connectedness Sector Decomposition")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig7b_connectedness_decomposition.png", dpi=240)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(df["z"], df["Respect_implicit"], label="Implicit Respect")
    plt.plot(df["z"], df["R_emergent"], label="Emergent Robustness")
    plt.xlabel("z")
    plt.ylabel("Value [0,1]")
    plt.title("Implicit Respect and Emergent Robustness")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig8_respect_robustness.png", dpi=240)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(df["z"], df["C_proj"], label="C_proj")
    plt.plot(df["z"], df["F_proj"], label="F_proj")
    plt.xlabel("z")
    plt.ylabel("Projection value")
    plt.title("MAAT Projection Layer")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig9_projection_layer.png", dpi=240)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(df["z"], df["CCI_v121"])
    plt.xlabel("z")
    plt.ylabel("CCI_v1.2.1")
    plt.title("Critical Coherence Index Diagnostic")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig10_CCI_v121.png", dpi=240)
    plt.close()


def save_summary(df):
    trusted_w = df.loc[df["w_MAAT_trusted"], "w_MAAT_stable"]

    summary = {
        "model": "MAAT v1.2.1 Observable Predictions with Emergent Robustness",
        "n_points": int(len(df)),
        "z_min": float(df["z"].min()),
        "z_max": float(df["z"].max()),
        "max_abs_delta_H": float(np.nanmax(np.abs(df["delta_H"]))),
        "max_Omega_MAAT": float(np.nanmax(df["Omega_MAAT"])),
        "max_abs_delta_fsigma8": float(np.nanmax(np.abs(df["delta_fsigma8"]))),
        "min_P_X": float(np.nanmin(df["P_X"])),
        "min_cs2": float(np.nanmin(df["cs2"])),
        "max_cs2": float(np.nanmax(df["cs2"])),
        "mean_H_struct": float(np.nanmean(df["H_struct"])),
        "mean_B_struct": float(np.nanmean(df["B_struct"])),
        "mean_S_struct": float(np.nanmean(df["S_struct"])),
        "mean_V_struct": float(np.nanmean(df["V_struct"])),
        "mean_Respect_implicit": float(np.nanmean(df["Respect_implicit"])),
        "mean_R_emergent": float(np.nanmean(df["R_emergent"])),
        "mean_CCI_v121": float(np.nanmean(df["CCI_v121"])),
        "max_CCI_v121": float(np.nanmax(df["CCI_v121"])),
        "mean_C_proj": float(np.nanmean(df["C_proj"])),
        "trusted_w_points": int(df["w_MAAT_trusted"].sum()),
        "trusted_w_min": float(np.nanmin(trusted_w)) if len(trusted_w) else None,
        "trusted_w_max": float(np.nanmax(trusted_w)) if len(trusted_w) else None,
        "trusted_w_final": float(trusted_w.iloc[-1]) if len(trusted_w) else None,
        "note": (
            "Growth is a simple proxy, not a Boltzmann-code calculation. "
            "Respect is implicit; robustness is emergent. "
            "w_MAAT_stable is trusted only where |Omega_MAAT| exceeds threshold."
        ),
    }

    with open(OUTDIR / "maat_v121_observable_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    return summary


def main():
    df = simulate()
    df = growth_proxy(df)

    df.to_csv(OUTDIR / "maat_v121_observable_predictions.csv", index=False)
    save_plots(df)
    summary = save_summary(df)

    print("\n=== MAAT v1.2.1 Observable Predictions: Emergent Robustness ===")
    for k, v in summary.items():
        print(f"{k}: {v}")

    print(f"\nSaved outputs to: {OUTDIR.resolve()}")
    print("Main figures:")
    print(" - fig1_Hz_MAAT_vs_LCDM.png")
    print(" - fig2_delta_H.png")
    print(" - fig3_Omega_MAAT.png")
    print(" - fig4_w_MAAT_stable.png")
    print(" - fig4b_w_MAAT_raw.png")
    print(" - fig5_fsigma8_proxy.png")
    print(" - fig6_delta_fsigma8.png")
    print(" - fig7_structural_fields_HBSV.png")
    print(" - fig7b_connectedness_decomposition.png")
    print(" - fig8_respect_robustness.png")
    print(" - fig9_projection_layer.png")
    print(" - fig10_CCI_v121.png")


if __name__ == "__main__":
    main()
