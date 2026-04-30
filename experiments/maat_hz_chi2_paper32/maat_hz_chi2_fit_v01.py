#!/usr/bin/env python3
# MAAT Paper 32 — χ² parameter scan against Cosmic Chronometer H(z)
# Scans μ and phidot0 for the v0.10 MAAT scalar trajectory and compares H(z) to data.

from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUTDIR = Path("maat_hz_chi2_fit_v01")
OUTDIR.mkdir(exist_ok=True)

# -----------------------------
# Background reference
# -----------------------------
H0 = 67.4
OMEGA_M = 0.315
OMEGA_L = 0.685

rho_m0 = 0.90
rho_r0 = 3.0e-4
rho_L = 2.10

epsilon = 1.0e-6
m_lambda = 0.55
Z_lambda = 1.0

a0 = 0.30
a1 = 1.00
lambda0 = 0.05
lambda_dot0 = 0.0

N_STEPS = 2500
DT = 0.002

# Compact Cosmic Chronometer H(z) table: z, H_obs, sigma_H
CHRONOMETER_DATA = np.array([
    (0.070, 69.0, 19.6),
    (0.090, 69.0, 12.0),
    (0.120, 68.6, 26.2),
    (0.170, 83.0, 8.0),
    (0.179, 75.0, 4.0),
    (0.199, 75.0, 5.0),
    (0.200, 72.9, 29.6),
    (0.270, 77.0, 14.0),
    (0.280, 88.8, 36.6),
    (0.352, 83.0, 14.0),
    (0.3802, 83.0, 13.5),
    (0.400, 95.0, 17.0),
    (0.4004, 77.0, 10.2),
    (0.4247, 87.1, 11.2),
    (0.4497, 92.8, 12.9),
    (0.470, 89.0, 49.6),
    (0.4783, 80.9, 9.0),
    (0.480, 97.0, 62.0),
    (0.593, 104.0, 13.0),
    (0.680, 92.0, 8.0),
    (0.781, 105.0, 12.0),
    (0.875, 125.0, 17.0),
    (0.880, 90.0, 40.0),
    (0.900, 117.0, 23.0),
    (1.037, 154.0, 20.0),
    (1.300, 168.0, 17.0),
    (1.363, 160.0, 33.6),
    (1.430, 177.0, 18.0),
    (1.530, 140.0, 14.0),
    (1.750, 202.0, 40.0),
    (1.965, 186.5, 50.4),
])

z_data = CHRONOMETER_DATA[:, 0]
H_obs = CHRONOMETER_DATA[:, 1]
sigma_H = CHRONOMETER_DATA[:, 2]


def lcdm_H(z):
    return H0 * np.sqrt(OMEGA_M * (1.0 + z) ** 3 + OMEGA_L)


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


def densities(a, phi_dot, lam, lam_dot, mu):
    X = 0.5 * phi_dot**2
    U = U_of_X(X)
    UX = U_X(X)

    rho_m = rho_m0 / a**3
    rho_r = rho_r0 / a**4

    rho_maat = (
        0.5 * Z_lambda * lam_dot**2
        + 0.5 * m_lambda**2 * lam**2
        + mu * lam * (2.0 * X * UX - U)
    )

    p_maat = (
        0.5 * Z_lambda * lam_dot**2
        - 0.5 * m_lambda**2 * lam**2
        + mu * lam * U
    )

    rho_total = rho_m + rho_r + rho_L + rho_maat
    if rho_total <= 0 or not np.isfinite(rho_total):
        return None

    H = np.sqrt(rho_total / 3.0)
    return X, rho_m, rho_r, rho_maat, p_maat, H


def rhs(y, mu):
    a, phi, phi_dot, lam, lam_dot = y
    den = densities(a, phi_dot, lam, lam_dot, mu)

    if den is None:
        return None

    X, rho_m, rho_r, rho_maat, p_maat, H = den

    UX = U_X(X)
    UXX = U_XX(X)

    PX = 1.0 + mu * lam * UX
    denom = PX + 2.0 * X * mu * lam * UXX

    if not np.isfinite(PX) or not np.isfinite(denom) or abs(denom) < 1e-12:
        return None

    a_dot = a * H
    phi_ddot = -(mu * UX * lam_dot * phi_dot + 3.0 * H * PX * phi_dot) / denom
    lam_ddot = (mu * U_of_X(X) - 3.0 * H * Z_lambda * lam_dot - m_lambda**2 * lam) / Z_lambda

    return np.array([a_dot, phi_dot, phi_ddot, lam_dot, lam_ddot], dtype=float)


def rk4_step(y, dt, mu):
    k1 = rhs(y, mu)
    if k1 is None:
        return None
    k2 = rhs(y + 0.5 * dt * k1, mu)
    if k2 is None:
        return None
    k3 = rhs(y + 0.5 * dt * k2, mu)
    if k3 is None:
        return None
    k4 = rhs(y + dt * k3, mu)
    if k4 is None:
        return None

    y_new = y + dt * (k1 + 2*k2 + 2*k3 + k4) / 6.0
    if not np.all(np.isfinite(y_new)):
        return None

    return y_new


def simulate(mu, phidot0):
    y = np.array([a0, 0.0, phidot0, lambda0, lambda_dot0], dtype=float)
    rows = []

    for _ in range(N_STEPS):
        a, phi, phi_dot, lam, lam_dot = y

        if a >= a1:
            break

        den = densities(a, phi_dot, lam, lam_dot, mu)
        if den is None:
            return None

        X, rho_m, rho_r, rho_maat, p_maat, H = den

        UX = U_X(X)
        UXX = U_XX(X)
        PX = 1.0 + mu * lam * UX
        denom = PX + 2.0 * X * mu * lam * UXX
        cs2 = PX / denom if abs(denom) > 1e-12 else np.nan

        if PX <= 0 or denom <= 0 or cs2 <= 0 or rho_maat < -1e-8:
            return None

        z = 1.0 / a - 1.0

        rows.append({
            "a": a,
            "z": z,
            "H_dimless": H,
            "rho_MAAT": rho_maat,
            "Omega_MAAT": rho_maat / (rho_m + rho_r + rho_L + rho_maat),
            "w_MAAT": p_maat / rho_maat if abs(rho_maat) > 1e-12 else np.nan,
            "P_X": PX,
            "cs2": cs2,
        })

        y = rk4_step(y, DT, mu)
        if y is None:
            return None

    df = pd.DataFrame(rows)
    if len(df) < 20:
        return None

    return df.sort_values("z")


def chi2_for_model(df):
    z_grid = df["z"].to_numpy()
    H_dimless = df["H_dimless"].to_numpy()

    # Ensure data range is covered
    if z_grid.min() > z_data.min() or z_grid.max() < z_data.max():
        return np.inf, None

    H0_model_dimless = np.interp(0.0, z_grid, H_dimless)
    H_model = H0 * H_dimless / H0_model_dimless
    H_interp = np.interp(z_data, z_grid, H_model)

    chi2 = float(np.sum(((H_interp - H_obs) / sigma_H) ** 2))
    return chi2, H_interp


def main():
    # Conservative scan: enough to test without exploding runtime.
    mu_values = np.logspace(-2, 2.2, 22)       # 0.01 to ~158
    phidot_values = np.linspace(0.15, 3.0, 28)

    rows = []

    lcdm_chi2 = float(np.sum(((lcdm_H(z_data) - H_obs) / sigma_H) ** 2))
    n_data = len(z_data)
    n_scan_parameters = 2
    dof_scan = n_data - n_scan_parameters

    best = None
    best_df = None
    best_H_interp = None

    total = len(mu_values) * len(phidot_values)
    count = 0

    for mu in mu_values:
        for phidot0 in phidot_values:
            count += 1
            df = simulate(mu, phidot0)

            if df is None:
                rows.append({
                    "mu": mu,
                    "phidot0": phidot0,
                    "stable": False,
                    "chi2": np.nan,
                    "chi2_per_point": np.nan,
                    "max_Omega_MAAT": np.nan,
                    "max_abs_delta_H": np.nan,
                    "min_P_X": np.nan,
                    "min_cs2": np.nan,
                    "max_cs2": np.nan,
                })
                continue

            chi2, H_interp = chi2_for_model(df)
            if not np.isfinite(chi2):
                stable = False
            else:
                stable = True

            z_grid = df["z"].to_numpy()
            H_dimless = df["H_dimless"].to_numpy()
            H0_model_dimless = np.interp(0.0, z_grid, H_dimless)
            H_model = H0 * H_dimless / H0_model_dimless
            H_lcdm_grid = lcdm_H(z_grid)
            max_abs_delta_H = float(np.nanmax(np.abs((H_model - H_lcdm_grid) / H_lcdm_grid)))

            row = {
                "mu": mu,
                "phidot0": phidot0,
                "stable": stable,
                "chi2": chi2 if stable else np.nan,
                "chi2_per_point": chi2 / n_data if stable else np.nan,
                "chi2_reduced_k2": chi2 / dof_scan if stable else np.nan,
                "delta_chi2_vs_LCDM": chi2 - lcdm_chi2 if stable else np.nan,
                "max_Omega_MAAT": float(np.nanmax(df["Omega_MAAT"])),
                "max_abs_delta_H": max_abs_delta_H,
                "min_P_X": float(np.nanmin(df["P_X"])),
                "min_cs2": float(np.nanmin(df["cs2"])),
                "max_cs2": float(np.nanmax(df["cs2"])),
            }
            rows.append(row)

            if stable and (best is None or chi2 < best["chi2"]):
                best = row
                best_df = df.copy()
                best_H_interp = H_interp.copy()

    res = pd.DataFrame(rows)
    res.to_csv(OUTDIR / "maat_hz_chi2_scan_results.csv", index=False)

    summary = {
        "lcdm_chi2": lcdm_chi2,
        "lcdm_chi2_per_point": lcdm_chi2 / n_data,
        "lcdm_reduced_chi2_fixed_reference": lcdm_chi2 / n_data,
        "n_data_points": int(n_data),
        "n_scan_parameters": int(n_scan_parameters),
        "degrees_of_freedom_for_scan": int(dof_scan),
        "n_scan_total": int(total),
        "n_stable": int(res["stable"].sum()),
        "best": best,
        "note": "First diagnostic parameter scan; not a full cosmological inference."
    }

    with open(OUTDIR / "maat_hz_chi2_fit_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\n=== MAAT H(z) χ² Scan v0.1 ===")
    print(f"ΛCDM chi2: {lcdm_chi2:.6f}")
    print(f"ΛCDM chi2/N: {lcdm_chi2 / n_data:.6f}")
    print(f"MAAT scan degrees of freedom (N-k, k=2): {dof_scan}")
    print(f"Stable models: {int(res['stable'].sum())}/{total}")

    if best is None:
        print("No stable model covered the H(z) data range.")
        return

    print("\n--- Best MAAT model ---")
    for k, v in best.items():
        print(f"{k}: {v}")

    # Save best trajectory
    best_df.to_csv(OUTDIR / "best_maat_trajectory.csv", index=False)

    # Smooth best model curve
    z_smooth = np.linspace(z_data.min(), z_data.max(), 400)

    z_grid = best_df["z"].to_numpy()
    H_dimless = best_df["H_dimless"].to_numpy()
    H0_model_dimless = np.interp(0.0, z_grid, H_dimless)
    H_model_grid = H0 * H_dimless / H0_model_dimless

    H_best_smooth = np.interp(z_smooth, z_grid, H_model_grid)
    H_lcdm_smooth = lcdm_H(z_smooth)

    # Plot 1: H(z) data comparison
    plt.figure(figsize=(8.5, 5.5))
    plt.errorbar(z_data, H_obs, yerr=sigma_H, fmt="o", capsize=3, label="Cosmic Chronometers")
    plt.plot(z_smooth, H_lcdm_smooth, "--", label="ΛCDM reference")
    plt.plot(z_smooth, H_best_smooth, label="Best MAAT scan model")
    plt.xlabel("z")
    plt.ylabel("H(z) [km s$^{-1}$ Mpc$^{-1}$]")
    plt.title("Best MAAT H(z) Fit vs Cosmic Chronometer Data")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig1_best_MAAT_vs_Hz_data.png", dpi=260)
    plt.close()

    # Plot 2: chi2 heatmap
    pivot = res.pivot(index="phidot0", columns="mu", values="chi2")
    plt.figure(figsize=(9, 5.8))
    plt.imshow(
        pivot.values,
        origin="lower",
        aspect="auto",
        extent=[
            np.log10(mu_values.min()),
            np.log10(mu_values.max()),
            phidot_values.min(),
            phidot_values.max(),
        ],
    )
    plt.colorbar(label=r"$\chi^2$")
    plt.xlabel(r"$\log_{10}\mu$")
    plt.ylabel(r"$\dot{\phi}_0$")
    plt.title(r"MAAT H(z) $\chi^2$ Scan")
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig2_chi2_heatmap.png", dpi=260)
    plt.close()

    # Plot 3: relative deviation best model
    delta_best = (H_best_smooth - H_lcdm_smooth) / H_lcdm_smooth
    plt.figure(figsize=(8.5, 5.5))
    plt.axhline(0, linestyle="--")
    plt.plot(z_smooth, delta_best)
    plt.xlabel("z")
    plt.ylabel("(H_MAAT - H_LCDM) / H_LCDM")
    plt.title("Best-Fit MAAT Relative Expansion Deviation")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig3_best_relative_H_deviation.png", dpi=260)
    plt.close()

    print(f"\nSaved outputs to: {OUTDIR.resolve()}")
    print("Main files:")
    print(" - maat_hz_chi2_scan_results.csv")
    print(" - maat_hz_chi2_fit_summary.json")
    print(" - best_maat_trajectory.csv")
    print(" - fig1_best_MAAT_vs_Hz_data.png")
    print(" - fig2_chi2_heatmap.png")
    print(" - fig3_best_relative_H_deviation.png")


if __name__ == "__main__":
    main()
