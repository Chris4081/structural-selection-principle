#!/usr/bin/env python3
"""
Extra Phenomenological Paper - MAAT v1.6 dark matter interpretation.

This script builds a minimal galaxy-rotation toy model:

* an exponential baryonic disk proxy;
* a smooth observed rotation curve with a flat outer part;
* an inferred effective dark residual from v_obs^2 - v_baryon^2;
* MAAT v1.6 structural supports along radius.

It does not fit real galaxy data and does not claim to replace particle dark
matter.  It operationalises the sentence:

    Dark matter is the gravitationally inferred structural residual required
    to keep the observed dynamical trajectory coherent.
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

_cache_dir = Path(tempfile.gettempdir()) / "maat_dark_matter_matplotlib_cache"
_cache_dir.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_cache_dir))
os.environ.setdefault("XDG_CACHE_HOME", str(_cache_dir))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


OUTDIR = Path("outputs_phenomenological")
EPS = 1e-12


def smooth_rotation_curve(r: np.ndarray, v_flat: float = 205.0, r_turn: float = 2.5) -> np.ndarray:
    """Observed toy rotation curve in km/s."""
    return v_flat * (1.0 - np.exp(-r / r_turn))


def baryonic_rotation_curve(
    r: np.ndarray, v_peak: float = 130.0, r_scale: float = 3.0, gas_floor: float = 35.0
) -> np.ndarray:
    """Simple exponential-disk-like baryonic contribution in km/s."""
    disk = v_peak * (r / r_scale) * np.exp(0.5 * (1.0 - r / r_scale))
    gas = gas_floor * (1.0 - np.exp(-r / 8.0))
    return np.sqrt(np.maximum(disk * disk + gas * gas, 0.0))


def gradient_log(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    return np.gradient(np.log(np.maximum(y, EPS)), np.log(np.maximum(x, EPS)))


def compute_supports(r: np.ndarray, v_obs: np.ndarray, v_bar: np.ndarray) -> pd.DataFrame:
    residual_v2 = np.maximum(v_obs * v_obs - v_bar * v_bar, 0.0)
    v_dm = np.sqrt(residual_v2)
    total_v2 = np.maximum(v_obs * v_obs, EPS)
    residual_fraction = residual_v2 / total_v2

    # Spherical-equivalent inferred mass in arbitrary units: M ~ r v^2.
    m_obs = r * v_obs * v_obs
    m_bar = r * v_bar * v_bar
    m_res = np.maximum(m_obs - m_bar, 0.0)
    dm_dr = np.gradient(m_res, r)
    rho_residual = np.maximum(dm_dr / (4.0 * np.pi * r * r + EPS), 0.0)
    rho_norm = rho_residual / (np.nanmax(rho_residual) + EPS)

    # H: Poisson/rotation consistency after adding residual.
    recon_error = np.abs(v_obs * v_obs - (v_bar * v_bar + v_dm * v_dm)) / total_v2
    H = 1.0 / (1.0 + 100.0 * recon_error)

    # B: radial force-balance smoothness.  Coherent outer disks prefer low
    # d ln v / d ln r magnitude, not necessarily zero everywhere.
    slope = gradient_log(v_obs, r)
    B = 1.0 / (1.0 + np.abs(slope))

    # S_eff: controlled activity in the transition region where residual support
    # emerges.  Too little residual is baryon-dominated; too much abruptness is
    # overactive.
    transition_activity = np.abs(np.gradient(residual_fraction, r))
    act_norm = transition_activity / (np.nanmedian(transition_activity) + EPS)
    S_eff = np.exp(-0.5 * (np.log(np.maximum(act_norm, EPS)) / 1.25) ** 2)

    # V: baryon-residual coupling.  We avoid the tautological statement that
    # v_dark^2/v_obs^2 is connected to 1 - v_baryon^2/v_obs^2 by comparing a
    # baryonic density proxy to the residual-density proxy instead.
    baryon_proxy = np.exp(-r / 3.0) + 0.22 * np.exp(-r / 8.0)
    baryon_proxy = baryon_proxy / (np.nanmax(baryon_proxy) + EPS)
    res_slope = np.abs(gradient_log(rho_norm + 0.03, r))
    bar_slope = np.abs(gradient_log(baryon_proxy + 0.03, r))
    shape_mismatch = np.abs(res_slope - bar_slope) / (res_slope + bar_slope + EPS)
    radial_overlap = np.sqrt(np.maximum(rho_norm, 0.0) * np.maximum(baryon_proxy, 0.0))
    V = np.clip(0.65 / (1.0 + shape_mismatch) + 0.35 * radial_overlap, 0.0, 1.0)

    R_resp = np.power(np.maximum(H * B * V, 0.0), 1.0 / 3.0)
    R_rob = np.minimum(R_resp, np.power(np.maximum(H * B * S_eff * V, 0.0), 0.25))

    # Structural dark score: high when residual is dynamically important and
    # robustly connected, not merely large.
    structural_dark_support = residual_fraction * R_rob * V
    warning = residual_fraction * (1.0 - R_rob) * (1.0 + np.abs(slope))

    return pd.DataFrame(
        {
            "r_kpc": r,
            "v_obs_km_s": v_obs,
            "v_baryon_km_s": v_bar,
            "v_dark_residual_km_s": v_dm,
            "residual_fraction": residual_fraction,
            "rho_residual_norm": rho_norm,
            "slope_dlnv_dlnr": slope,
            "H": H,
            "B": B,
            "S_eff": S_eff,
            "V": V,
            "R_resp": R_resp,
            "R_rob": R_rob,
            "structural_dark_support": structural_dark_support,
            "structural_warning": warning,
        }
    )


def plot_curves(df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    ax.plot(df["r_kpc"], df["v_obs_km_s"], lw=2.4, label="observed toy curve", color="#1b4965")
    ax.plot(df["r_kpc"], df["v_baryon_km_s"], lw=2.2, label="baryonic contribution", color="#ca6702")
    ax.plot(df["r_kpc"], df["v_dark_residual_km_s"], lw=2.2, label="effective dark residual", color="#5e548e")
    ax.set_xlabel("radius [kpc]")
    ax.set_ylabel("rotation velocity [km/s]")
    ax.set_title("Toy rotation curve and inferred dark residual")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig1_rotation_curve_residual.png", dpi=180)
    plt.close(fig)


def plot_supports(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for key, color in [
        ("H", "#0077b6"),
        ("B", "#2a9d8f"),
        ("S_eff", "#e9c46a"),
        ("V", "#9b5de5"),
        ("R_rob", "#ef476f"),
    ]:
        axes[0].plot(df["r_kpc"], df[key], lw=2.0, label=key, color=color)
    axes[0].set_ylabel("support")
    axes[0].set_title("MAAT v1.6 supports along radius")
    axes[0].grid(alpha=0.25)
    axes[0].legend(ncol=5, fontsize=8)
    axes[1].plot(df["r_kpc"], df["residual_fraction"], lw=2.2, label="residual fraction", color="#5e548e")
    axes[1].plot(
        df["r_kpc"],
        df["structural_dark_support"],
        lw=2.2,
        label="structural dark support",
        color="#006d77",
    )
    axes[1].plot(df["r_kpc"], df["rho_residual_norm"], lw=1.8, label="normalised residual density", color="#bb3e03")
    axes[1].set_xlabel("radius [kpc]")
    axes[1].set_ylabel("diagnostic")
    axes[1].grid(alpha=0.25)
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig2_structural_supports.png", dpi=180)
    plt.close(fig)


def plot_phase(df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    sc = ax.scatter(
        df["residual_fraction"],
        df["R_rob"],
        c=df["r_kpc"],
        s=34 + 80 * df["structural_dark_support"],
        cmap="viridis",
        alpha=0.85,
    )
    ax.set_xlabel("dark residual fraction")
    ax.set_ylabel("MAAT robustness")
    ax.set_title("Dark residual as structural support, not only missing mass")
    ax.grid(alpha=0.25)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("radius [kpc]")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig3_residual_robustness_phase.png", dpi=180)
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    r = np.linspace(0.4, 28.0, 240)
    v_obs = smooth_rotation_curve(r)
    v_bar = baryonic_rotation_curve(r)
    df = compute_supports(r, v_obs, v_bar)

    df.to_csv(OUTDIR / "dark_matter_structural_coherence_profile.csv", index=False)
    plot_curves(df)
    plot_supports(df)
    plot_phase(df)

    outer = df[df["r_kpc"] >= 12.0]
    transition = df[(df["r_kpc"] >= 5.0) & (df["r_kpc"] <= 12.0)]
    summary = {
        "paper": "extra_phenomenological",
        "title": "Dark Matter as Structural Coherence Residual in MAAT v1.6",
        "status": "toy phenomenological diagnostic, not a replacement for particle dark matter",
        "max_observed_velocity_km_s": float(df["v_obs_km_s"].max()),
        "outer_mean_residual_fraction": float(outer["residual_fraction"].mean()),
        "outer_mean_R_rob": float(outer["R_rob"].mean()),
        "transition_peak_structural_dark_support": float(transition["structural_dark_support"].max()),
        "radius_peak_structural_dark_support_kpc": float(
            df.loc[df["structural_dark_support"].idxmax(), "r_kpc"]
        ),
        "mean_H": float(df["H"].mean()),
        "mean_B": float(df["B"].mean()),
        "mean_S_eff": float(df["S_eff"].mean()),
        "mean_V": float(df["V"].mean()),
        "mean_R_rob": float(df["R_rob"].mean()),
    }
    with open(OUTDIR / "dark_matter_structural_coherence_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
