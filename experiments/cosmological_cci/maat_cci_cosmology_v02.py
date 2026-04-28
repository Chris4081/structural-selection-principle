#!/usr/bin/env python3
"""Reproduce the Cosmological Critical Coherence Index v0.2 artifacts.

The script generates:

- maat_cci_cosmology_v02.csv
- maat_cci_cosmology_v02_data_comparison.csv
- maat_cci_cosmology_v02_plot.png
- maat_cci_cosmology_v02_Hz_comparison.png
- maat_cci_cosmology_v02_data_comparison.png
- maat_cci_cosmology_v02_residuals.png

It is a compact diagnostic pipeline, not a precision cosmology likelihood.
The chronometer table is a small literature-derived compilation used only
for the observational projection shown in the accompanying note.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


H0 = 67.4
OMEGA_M = 0.315
OMEGA_L = 0.685
SIGMA8_0 = 0.811
EPS = 1.0e-9


def e_lcdm(z: np.ndarray | float) -> np.ndarray | float:
    return np.sqrt(OMEGA_M * (1.0 + z) ** 3 + OMEGA_L)


def omega_m_z(z: np.ndarray | float) -> np.ndarray | float:
    ez2 = e_lcdm(z) ** 2
    return OMEGA_M * (1.0 + z) ** 3 / ez2


def omega_l_z(z: np.ndarray | float) -> np.ndarray | float:
    ez2 = e_lcdm(z) ** 2
    return OMEGA_L / ez2


def growth_suppression(z: np.ndarray | float) -> np.ndarray | float:
    om = omega_m_z(z)
    ol = omega_l_z(z)
    denominator = (
        om ** (4.0 / 7.0)
        - ol
        + (1.0 + om / 2.0) * (1.0 + ol / 70.0)
    )
    return (5.0 * om / 2.0) / denominator


def growth_factor(z: np.ndarray | float) -> np.ndarray | float:
    return growth_suppression(z) / (growth_suppression(0.0) * (1.0 + z))


def cci_from_e_and_growth(
    z: np.ndarray | float,
    e_value: np.ndarray | float,
    d_growth: np.ndarray | float,
) -> np.ndarray | float:
    gamma_inst = e_value
    gamma_prod = 1.0 + z
    gamma_coh = d_growth
    u_struct = np.abs(gamma_inst - gamma_coh) / (gamma_inst + gamma_coh + EPS)
    return gamma_inst * gamma_prod * (1.0 + u_struct) / (1.0 + gamma_coh + EPS)


def write_model_grid(out_path: Path) -> dict[str, np.ndarray]:
    z = np.linspace(0.0, 10.0, 900)
    e = e_lcdm(z)
    omz = omega_m_z(z)
    olz = omega_l_z(z)
    d = growth_factor(z)
    sigma8_z = SIGMA8_0 * d
    gamma_inst = e
    gamma_prod = 1.0 + z
    gamma_coh = d
    u_struct = np.abs(e - d) / (e + d + EPS)
    cci = cci_from_e_and_growth(z, e, d)
    cci_norm = cci / cci[0]

    columns = [
        "z",
        "E_H_over_H0",
        "Omega_m_z",
        "Omega_L_z",
        "D_growth_norm",
        "sigma8_z",
        "Gamma_inst",
        "Gamma_prod",
        "Gamma_coh",
        "U_struct",
        "CCI_cos_v02",
        "CCI_cos_v02_normalized",
    ]
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for row in zip(
            z,
            e,
            omz,
            olz,
            d,
            sigma8_z,
            gamma_inst,
            gamma_prod,
            gamma_coh,
            u_struct,
            cci,
            cci_norm,
        ):
            writer.writerow([f"{value:.18e}" for value in row])

    return {
        "z": z,
        "e": e,
        "d": d,
        "cci": cci,
        "cci_norm": cci_norm,
    }


def read_chronometers(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    z_values: list[float] = []
    h_values: list[float] = []
    h_errors: list[float] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            z_values.append(float(row["z"]))
            h_values.append(float(row["H_obs_km_s_Mpc"]))
            h_errors.append(float(row["H_err"]))
    return np.array(z_values), np.array(h_values), np.array(h_errors)


def write_data_comparison(raw_path: Path, out_path: Path) -> dict[str, np.ndarray]:
    z, h_obs, h_err = read_chronometers(raw_path)
    h_lcdm = H0 * e_lcdm(z)
    h_residual = h_obs - h_lcdm
    h_pull = h_residual / h_err

    d = growth_factor(z)
    sigma8_z = SIGMA8_0 * d
    gamma_coh = d

    e_obs = h_obs / H0
    e_lcdm_values = h_lcdm / H0
    u_obs = np.abs(e_obs - gamma_coh) / (e_obs + gamma_coh + EPS)
    cci_obs = cci_from_e_and_growth(z, e_obs, d)
    cci_lcdm = cci_from_e_and_growth(z, e_lcdm_values, d)
    cci0 = cci_from_e_and_growth(0.0, 1.0, 1.0)
    cci_obs_norm = cci_obs / cci0
    cci_lcdm_norm = cci_lcdm / cci0
    cci_residual = cci_obs_norm - cci_lcdm_norm

    cci_hi = cci_from_e_and_growth(z, (h_obs + h_err) / H0, d) / cci0
    cci_lo = cci_from_e_and_growth(z, np.maximum(h_obs - h_err, EPS) / H0, d) / cci0
    cci_err = 0.5 * np.abs(cci_hi - cci_lo)

    columns = [
        "z",
        "H_obs_km_s_Mpc",
        "H_err",
        "H_LCDM_Planck",
        "H_residual",
        "H_pull_sigma",
        "D_growth_LCDM",
        "sigma8_z_proxy",
        "Gamma_coh",
        "U_struct_obs",
        "CCI_obs_norm",
        "CCI_LCDM_norm",
        "CCI_residual",
        "CCI_err_norm_approx",
    ]
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for row in zip(
            z,
            h_obs,
            h_err,
            h_lcdm,
            h_residual,
            h_pull,
            d,
            sigma8_z,
            gamma_coh,
            u_obs,
            cci_obs_norm,
            cci_lcdm_norm,
            cci_residual,
            cci_err,
        ):
            writer.writerow([f"{value:.16g}" for value in row])

    return {
        "z": z,
        "h_obs": h_obs,
        "h_err": h_err,
        "h_lcdm": h_lcdm,
        "cci_obs_norm": cci_obs_norm,
        "cci_lcdm_norm": cci_lcdm_norm,
        "cci_residual": cci_residual,
        "cci_err": cci_err,
    }


def make_model_plot(model: dict[str, np.ndarray], out_path: Path) -> None:
    z = model["z"]
    plt.figure(figsize=(12, 7))
    plt.semilogy(z, model["cci_norm"], lw=3, label="MAAT CCI_cos v0.2 (normalised)")
    plt.semilogy(z, model["e"], "--", lw=2.5, label="Expansion stress E(z)")
    plt.semilogy(z, model["d"], ":", lw=3, label="Linear growth coherence D(z)")
    plt.axvspan(0.0, 0.5, color="#dbeaf2", alpha=0.45, label="Late structured era")
    plt.axvspan(0.5, 3.0, color="#dbeaf2", alpha=0.30, label="Structure-formation era")
    plt.axvspan(3.0, 10.0, color="#dbeaf2", alpha=0.20, label="Early high-stress regime")
    plt.xlabel("Redshift z")
    plt.ylabel("Normalised value")
    plt.title("MAAT CCI-Cosmology v0.2\nExpansion / activity stress vs. linear growth coherence")
    plt.grid(alpha=0.25, which="both")
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.savefig(out_path, dpi=220)
    plt.close()


def make_hz_plot(data: dict[str, np.ndarray], out_path: Path) -> None:
    plt.figure(figsize=(10, 6))
    z_grid = np.linspace(0.0, 2.1, 500)
    plt.plot(z_grid, H0 * e_lcdm(z_grid), lw=2.5, label="Planck LCDM")
    plt.errorbar(
        data["z"],
        data["h_obs"],
        yerr=data["h_err"],
        fmt="o",
        ms=4,
        capsize=2,
        label="Cosmic Chronometers",
    )
    plt.xlabel("Redshift z")
    plt.ylabel(r"$H(z)$ [km s$^{-1}$ Mpc$^{-1}$]")
    plt.title("Cosmic Chronometer H(z) data vs. Planck LCDM")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=220)
    plt.close()


def make_cci_data_plot(data: dict[str, np.ndarray], out_path: Path) -> None:
    plt.figure(figsize=(10, 6))
    plt.plot(data["z"], data["cci_lcdm_norm"], lw=2.5, label="Planck LCDM projection")
    plt.errorbar(
        data["z"],
        data["cci_obs_norm"],
        yerr=data["cci_err"],
        fmt="o",
        ms=4,
        capsize=2,
        label="Chronometer projection",
    )
    plt.xlabel("Redshift z")
    plt.ylabel(r"Normalised $\mathrm{CCI}_{\mathrm{cos}}$")
    plt.title("Cosmological CCI projection from H(z) data")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=220)
    plt.close()


def make_residual_plot(data: dict[str, np.ndarray], out_path: Path) -> None:
    plt.figure(figsize=(10, 5))
    plt.axhline(0.0, color="black", lw=1.0)
    plt.errorbar(
        data["z"],
        data["cci_residual"],
        yerr=data["cci_err"],
        fmt="o",
        ms=4,
        capsize=2,
    )
    plt.xlabel("Redshift z")
    plt.ylabel(r"$\Delta \widehat{\mathrm{CCI}}_{\mathrm{cos}}$")
    plt.title("CCI residuals relative to Planck LCDM")
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(out_path, dpi=220)
    plt.close()


def main() -> None:
    base = Path(__file__).resolve().parent
    plot_dir = base / "plots"
    plot_dir.mkdir(exist_ok=True)
    raw_path = base / "maat_cci_cosmology_v02_chronometers.csv"

    model = write_model_grid(base / "maat_cci_cosmology_v02.csv")
    data = write_data_comparison(raw_path, base / "maat_cci_cosmology_v02_data_comparison.csv")

    make_model_plot(model, plot_dir / "maat_cci_cosmology_v02_plot.png")
    make_hz_plot(data, plot_dir / "maat_cci_cosmology_v02_Hz_comparison.png")
    make_cci_data_plot(data, plot_dir / "maat_cci_cosmology_v02_data_comparison.png")
    make_residual_plot(data, plot_dir / "maat_cci_cosmology_v02_residuals.png")

    print("Wrote maat_cci_cosmology_v02.csv")
    print("Wrote maat_cci_cosmology_v02_data_comparison.csv")
    print("Wrote four PNG figures to plots/.")


if __name__ == "__main__":
    main()
