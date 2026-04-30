#!/usr/bin/env python3
"""MAAT v0.10 observable predictions from the v0.9 FLRW trajectory.

This script converts the representative v0.9 stable FLRW trajectory into a
minimal observable projection:

1. H(z) / H0 compared with an internally normalised LCDM reference.
2. Relative Hubble deviation Delta H / H.
3. MAAT equation of state w_MAAT(z).
4. MAAT density fraction Omega_MAAT(z).
5. A simple growth proxy f sigma_8(z).

The calculation is intentionally a diagnostic benchmark, not a fit to
cosmological data and not a full perturbation solver.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_INPUT = (
    Path(__file__).resolve().parents[1]
    / "maat_dynamic_fields_v05_v09"
    / "v09_flrw_stability_scan"
    / "outputs"
    / "representative_flrw_trajectory.csv"
)


def cumulative_trapezoid(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Small scipy-free cumulative trapezoid helper."""

    out = np.zeros_like(y, dtype=float)
    if len(y) < 2:
        return out
    dx = np.diff(x)
    area = 0.5 * (y[1:] + y[:-1]) * dx
    out[1:] = np.cumsum(area)
    return out


def prepare_observables(input_csv: Path, sigma8_0: float = 0.811) -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(input_csv)

    # Use the a ~= 1 slice as the effective present epoch. This avoids using
    # the far-future tail of the v0.9 stability trajectory as z=0.
    ref_idx = int((df["a"] - 1.0).abs().idxmin())
    ref = df.loc[ref_idx]
    a0 = float(ref["a"])
    h0 = float(ref["H"])

    branch = df[df["a"] <= a0].copy()
    branch["z"] = a0 / branch["a"] - 1.0
    branch = branch.sort_values("z").reset_index(drop=True)

    z = branch["z"].to_numpy(dtype=float)
    e_maat = branch["H"].to_numpy(dtype=float) / h0

    rho_m0 = float(ref["rho_m"])
    rho_r0 = float(ref["rho_r"])
    rho_l0 = float(ref["rho_lambda_cosmo"])
    rho_lcdm0 = rho_m0 + rho_r0 + rho_l0

    omega_m0 = rho_m0 / rho_lcdm0
    omega_r0 = rho_r0 / rho_lcdm0
    omega_l0 = rho_l0 / rho_lcdm0

    e_lcdm = np.sqrt(
        omega_m0 * (1.0 + z) ** 3
        + omega_r0 * (1.0 + z) ** 4
        + omega_l0
    )
    delta_h = (e_maat - e_lcdm) / e_lcdm

    rho_maat = branch["rho_maat"].to_numpy(dtype=float)
    p_maat = branch["p_maat"].to_numpy(dtype=float)
    rho_total = branch["rho_total"].to_numpy(dtype=float)
    w_maat = np.divide(
        p_maat,
        rho_maat,
        out=np.full_like(p_maat, np.nan, dtype=float),
        where=np.abs(rho_maat) > 1e-14,
    )
    omega_maat = rho_maat / rho_total

    omega_m_lcdm = omega_m0 * (1.0 + z) ** 3 / (e_lcdm**2)
    omega_m_maat = branch["rho_m"].to_numpy(dtype=float) / rho_total
    omega_m_lcdm = np.clip(omega_m_lcdm, 0.0, 1.0)
    omega_m_maat = np.clip(omega_m_maat, 0.0, 1.0)

    gamma_gr = 0.55
    f_lcdm = omega_m_lcdm**gamma_gr
    f_maat = omega_m_maat**gamma_gr

    # Approximate D(z) from f=d ln D/d ln a, normalised to D(0)=1.
    integrand_lcdm = f_lcdm / (1.0 + z)
    integrand_maat = f_maat / (1.0 + z)
    d_lcdm = np.exp(-cumulative_trapezoid(integrand_lcdm, z))
    d_maat = np.exp(-cumulative_trapezoid(integrand_maat, z))

    fs8_lcdm = sigma8_0 * f_lcdm * d_lcdm
    fs8_maat = sigma8_0 * f_maat * d_maat
    delta_fs8 = (fs8_maat - fs8_lcdm) / fs8_lcdm

    obs = pd.DataFrame(
        {
            "z": z,
            "a": branch["a"].to_numpy(dtype=float),
            "E_maat": e_maat,
            "E_lcdm": e_lcdm,
            "delta_H_over_H": delta_h,
            "w_maat": w_maat,
            "Omega_maat": omega_maat,
            "Omega_m_maat": omega_m_maat,
            "Omega_m_lcdm": omega_m_lcdm,
            "f_maat": f_maat,
            "f_lcdm": f_lcdm,
            "D_maat": d_maat,
            "D_lcdm": d_lcdm,
            "fsigma8_maat": fs8_maat,
            "fsigma8_lcdm": fs8_lcdm,
            "delta_fsigma8": delta_fs8,
        }
    )

    max_delta_idx = int(np.nanargmax(np.abs(delta_h)))
    max_fs8_idx = int(np.nanargmax(np.abs(delta_fs8)))

    summary = {
        "input_csv": str(input_csv),
        "reference_index": ref_idx,
        "a_reference": a0,
        "H_reference": h0,
        "z_min": float(np.nanmin(z)),
        "z_max": float(np.nanmax(z)),
        "Omega_m0_lcdm_reference": omega_m0,
        "Omega_r0_lcdm_reference": omega_r0,
        "Omega_lambda0_lcdm_reference": omega_l0,
        "Omega_maat_at_reference": float(ref["omega_maat"]),
        "w_maat_at_reference": float(ref["w_maat"]),
        "max_abs_delta_H_over_H": float(np.nanmax(np.abs(delta_h))),
        "z_at_max_abs_delta_H_over_H": float(z[max_delta_idx]),
        "mean_abs_delta_H_over_H": float(np.nanmean(np.abs(delta_h))),
        "max_Omega_maat": float(np.nanmax(omega_maat)),
        "z_at_max_Omega_maat": float(z[int(np.nanargmax(omega_maat))]),
        "w_maat_min": float(np.nanmin(w_maat)),
        "w_maat_max": float(np.nanmax(w_maat)),
        "max_abs_delta_fsigma8": float(np.nanmax(np.abs(delta_fs8))),
        "z_at_max_abs_delta_fsigma8": float(z[max_fs8_idx]),
        "sigma8_0": sigma8_0,
        "growth_index_gamma": gamma_gr,
        "status": "diagnostic observable projection, not data fit",
    }

    return obs, summary


def plot_observables(obs: pd.DataFrame, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    z = obs["z"].to_numpy()

    plt.figure(figsize=(7.2, 4.6))
    plt.plot(z, obs["E_lcdm"], label=r"$\Lambda$CDM reference", lw=2.4)
    plt.plot(z, obs["E_maat"], label="MAAT v0.10 projection", lw=2.2)
    plt.xlabel("redshift z")
    plt.ylabel(r"$E(z)=H(z)/H_0$")
    plt.title("Expansion history")
    plt.legend()
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(outdir / "hubble_history_vs_lcdm.png", dpi=220)
    plt.close()

    plt.figure(figsize=(7.2, 4.2))
    plt.axhline(0.0, color="black", lw=1.0)
    plt.plot(z, obs["delta_H_over_H"], color="#9a3412", lw=2.3)
    plt.xlabel("redshift z")
    plt.ylabel(r"$(H_{\rm MAAT}-H_{\Lambda{\rm CDM}})/H_{\Lambda{\rm CDM}}$")
    plt.title("Relative Hubble deviation")
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(outdir / "relative_hubble_deviation.png", dpi=220)
    plt.close()

    plt.figure(figsize=(7.2, 4.2))
    plt.axhline(-1.0, color="black", lw=1.0, ls="--", label=r"$w=-1$")
    plt.plot(z, obs["w_maat"], color="#1d4ed8", lw=2.3, label=r"$w_{\rm MAAT}$")
    plt.xlabel("redshift z")
    plt.ylabel(r"$w_{\rm MAAT}(z)$")
    plt.title("MAAT equation of state")
    plt.legend()
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(outdir / "maat_equation_of_state.png", dpi=220)
    plt.close()

    plt.figure(figsize=(7.2, 4.2))
    plt.plot(z, obs["Omega_maat"], color="#047857", lw=2.3)
    plt.xlabel("redshift z")
    plt.ylabel(r"$\Omega_{\rm MAAT}(z)$")
    plt.title("MAAT density fraction")
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(outdir / "maat_density_fraction.png", dpi=220)
    plt.close()

    fig, ax = plt.subplots(2, 1, figsize=(7.2, 6.4), sharex=True)
    ax[0].plot(z, obs["fsigma8_lcdm"], label=r"$\Lambda$CDM proxy", lw=2.2)
    ax[0].plot(z, obs["fsigma8_maat"], label="MAAT proxy", lw=2.2)
    ax[0].set_ylabel(r"$f\sigma_8(z)$")
    ax[0].legend()
    ax[0].grid(alpha=0.25)
    ax[1].axhline(0.0, color="black", lw=1.0)
    ax[1].plot(z, obs["delta_fsigma8"], color="#7c3aed", lw=2.2)
    ax[1].set_xlabel("redshift z")
    ax[1].set_ylabel(r"$\Delta f\sigma_8/f\sigma_8$")
    ax[1].grid(alpha=0.25)
    fig.suptitle("Growth proxy")
    fig.tight_layout()
    fig.savefig(outdir / "growth_proxy_fsigma8.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(2, 2, figsize=(10.5, 7.8))
    ax = ax.ravel()
    ax[0].plot(z, obs["E_lcdm"], label=r"$\Lambda$CDM", lw=2.1)
    ax[0].plot(z, obs["E_maat"], label="MAAT", lw=2.1)
    ax[0].set_title("A: expansion")
    ax[0].set_ylabel(r"$E(z)$")
    ax[0].legend(fontsize=8)
    ax[1].axhline(0.0, color="black", lw=0.9)
    ax[1].plot(z, obs["delta_H_over_H"], color="#9a3412", lw=2.1)
    ax[1].set_title("B: Hubble deviation")
    ax[1].set_ylabel(r"$\Delta H/H$")
    ax[2].axhline(-1.0, color="black", lw=0.9, ls="--")
    ax[2].plot(z, obs["w_maat"], color="#1d4ed8", lw=2.1)
    ax[2].set_title("C: equation of state")
    ax[2].set_xlabel("redshift z")
    ax[2].set_ylabel(r"$w_{\rm MAAT}$")
    ax[3].plot(z, obs["Omega_maat"], color="#047857", lw=2.1)
    ax[3].set_title("D: density fraction")
    ax[3].set_xlabel("redshift z")
    ax[3].set_ylabel(r"$\Omega_{\rm MAAT}$")
    for axis in ax:
        axis.grid(alpha=0.22)
    fig.suptitle("MAAT v0.10 observable prediction layer")
    fig.tight_layout()
    fig.savefig(outdir / "observable_predictions_summary.png", dpi=240)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--sigma8", type=float, default=0.811)
    args = parser.parse_args()

    obs, summary = prepare_observables(args.input, sigma8_0=args.sigma8)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    obs.to_csv(args.output_dir / "maat_v10_observables.csv", index=False)
    with (args.output_dir / "maat_v10_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    plot_observables(obs, args.output_dir)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

