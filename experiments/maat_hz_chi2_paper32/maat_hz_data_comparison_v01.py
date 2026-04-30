#!/usr/bin/env python3
# MAAT Paper 32 — First H(z) Data Comparison
# Compares MAAT v0.10 H(z) prediction with Cosmic Chronometer H(z) data.

from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUTDIR = Path("maat_hz_data_comparison_v01")
OUTDIR.mkdir(exist_ok=True)

V10_CSV = (
    Path(__file__).resolve().parents[1]
    / "maat_observable_predictions_v10"
    / "outputs"
    / "maat_v10_observables.csv"
)

H0 = 67.4
OMEGA_M = 0.315
OMEGA_L = 0.685


# Compact Cosmic Chronometer H(z) table.
# Columns: z, H_obs, sigma_H
CHRONOMETER_DATA = [
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
]


def lcdm_H(z):
    return H0 * np.sqrt(OMEGA_M * (1 + z) ** 3 + OMEGA_L)


def load_maat_prediction():
    if not V10_CSV.exists():
        raise FileNotFoundError(
            f"Could not find {V10_CSV}. "
            "Run maat_observable_predictions_v10.py first."
        )

    df = pd.read_csv(V10_CSV)

    required = {"z"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in v0.10 CSV: {missing}")

    # Support both the early TOE-column names and the repository v0.10 names.
    if "H_MAAT" not in df.columns and "E_maat" in df.columns:
        df["H_MAAT"] = df["E_maat"]
    if "H_LCDM" not in df.columns and "E_lcdm" in df.columns:
        df["H_LCDM"] = df["E_lcdm"]

    required_model = {"H_MAAT", "H_LCDM"}
    missing_model = required_model - set(df.columns)
    if missing_model:
        raise ValueError(f"Missing model columns in v0.10 CSV: {missing_model}")

    df = df.sort_values("z").dropna(subset=["z", "H_MAAT"])
    return df


def interpolate_model(df_model, z_data):
    z_grid = df_model["z"].to_numpy()
    H_maat_dimless = df_model["H_MAAT"].to_numpy()

    # Convert dimensionless model H to km/s/Mpc by matching H0 normalization.
    H_maat = H0 * H_maat_dimless / np.interp(0.0, z_grid, H_maat_dimless)

    H_interp = np.interp(z_data, z_grid, H_maat)
    return H_interp


def chi_square(H_model, H_obs, sigma):
    return float(np.sum(((H_model - H_obs) / sigma) ** 2))


def main():
    df_model = load_maat_prediction()

    data = pd.DataFrame(
        CHRONOMETER_DATA,
        columns=["z", "H_obs", "sigma_H"]
    )

    data["H_LCDM"] = lcdm_H(data["z"])
    data["H_MAAT"] = interpolate_model(df_model, data["z"])

    data["pull_LCDM"] = (data["H_LCDM"] - data["H_obs"]) / data["sigma_H"]
    data["pull_MAAT"] = (data["H_MAAT"] - data["H_obs"]) / data["sigma_H"]

    chi2_lcdm = chi_square(data["H_LCDM"], data["H_obs"], data["sigma_H"])
    chi2_maat = chi_square(data["H_MAAT"], data["H_obs"], data["sigma_H"])

    n = len(data)
    # No parameters are fitted in this fixed-branch diagnostic.
    k_fixed = 0
    dof_fixed = n - k_fixed

    summary = {
        "n_data_points": n,
        "n_fitted_parameters_fixed_branch": k_fixed,
        "degrees_of_freedom_fixed_branch": dof_fixed,
        "chi2_LCDM": chi2_lcdm,
        "chi2_MAAT": chi2_maat,
        "chi2_per_point_LCDM": chi2_lcdm / n,
        "chi2_per_point_MAAT": chi2_maat / n,
        "reduced_chi2_LCDM_fixed": chi2_lcdm / dof_fixed,
        "reduced_chi2_MAAT_fixed": chi2_maat / dof_fixed,
        "delta_chi2_MAAT_minus_LCDM": chi2_maat - chi2_lcdm,
        "max_abs_pull_LCDM": float(np.max(np.abs(data["pull_LCDM"]))),
        "max_abs_pull_MAAT": float(np.max(np.abs(data["pull_MAAT"]))),
        "note": "This is a first diagnostic comparison, not a parameter fit."
    }

    data.to_csv(OUTDIR / "maat_vs_chronometer_Hz_table.csv", index=False)

    with open(OUTDIR / "maat_vs_chronometer_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Smooth model curves
    z_smooth = np.linspace(data["z"].min(), data["z"].max(), 400)
    H_lcdm_smooth = lcdm_H(z_smooth)

    z_grid = df_model["z"].to_numpy()
    H_dimless = df_model["H_MAAT"].to_numpy()
    H_maat_grid = H0 * H_dimless / np.interp(0.0, z_grid, H_dimless)
    H_maat_smooth = np.interp(z_smooth, z_grid, H_maat_grid)

    # Figure 1: H(z)
    plt.figure(figsize=(8.5, 5.5))
    plt.errorbar(
        data["z"], data["H_obs"], yerr=data["sigma_H"],
        fmt="o", capsize=3, label="Cosmic Chronometers"
    )
    plt.plot(z_smooth, H_lcdm_smooth, "--", label="ΛCDM reference")
    plt.plot(z_smooth, H_maat_smooth, label="MAAT v0.10")
    plt.xlabel("z")
    plt.ylabel("H(z) [km s$^{-1}$ Mpc$^{-1}$]")
    plt.title("MAAT v0.10 vs Cosmic Chronometer H(z) Data")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig1_MAAT_vs_Hz_data.png", dpi=260)
    plt.close()

    # Figure 2: pulls
    plt.figure(figsize=(8.5, 5.5))
    plt.axhline(0, linestyle="--")
    plt.scatter(data["z"], data["pull_LCDM"], label="ΛCDM pull")
    plt.scatter(data["z"], data["pull_MAAT"], marker="x", label="MAAT pull")
    plt.xlabel("z")
    plt.ylabel("(H_model - H_obs) / σ_H")
    plt.title("H(z) Pulls Relative to Cosmic Chronometer Data")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig2_Hz_pulls.png", dpi=260)
    plt.close()

    # Figure 3: MAAT - LCDM model deviation
    delta_model = (H_maat_smooth - H_lcdm_smooth) / H_lcdm_smooth

    plt.figure(figsize=(8.5, 5.5))
    plt.axhline(0, linestyle="--")
    plt.plot(z_smooth, delta_model)
    plt.xlabel("z")
    plt.ylabel("(H_MAAT - H_LCDM) / H_LCDM")
    plt.title("Relative MAAT Expansion Deviation")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig3_MAAT_relative_H_deviation.png", dpi=260)
    plt.close()

    print("\n=== MAAT v0.10 First H(z) Data Comparison ===")
    for k, v in summary.items():
        print(f"{k}: {v}")

    print(f"\nSaved outputs to: {OUTDIR.resolve()}")
    print("Main files:")
    print(" - maat_vs_chronometer_Hz_table.csv")
    print(" - maat_vs_chronometer_summary.json")
    print(" - fig1_MAAT_vs_Hz_data.png")
    print(" - fig2_Hz_pulls.png")
    print(" - fig3_MAAT_relative_H_deviation.png")


if __name__ == "__main__":
    main()
