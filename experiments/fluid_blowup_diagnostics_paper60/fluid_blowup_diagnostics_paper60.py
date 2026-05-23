#!/usr/bin/env python3
"""
Paper 60 - MAAT fluid blow-up diagnostics.

This script is deliberately conservative.  It does not attempt to solve the
Navier-Stokes regularity problem.  It tests whether MAAT-style structural
supports can act as early-warning diagnostics in two controlled fluid settings:

1. viscous Burgers shock formation as a finite-gradient warning case;
2. unforced 2D incompressible Navier-Stokes as a regularity/control case.

Outputs are written to outputs_paper60/.
"""

from __future__ import annotations

import json
import math
import os
import tempfile
from pathlib import Path

_cache_dir = Path(tempfile.gettempdir()) / "maat_paper60_matplotlib_cache"
_cache_dir.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_cache_dir))
os.environ.setdefault("XDG_CACHE_HOME", str(_cache_dir))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


EPS = 1e-12
OUTDIR = Path("outputs_paper60")


def rankdata(x: np.ndarray) -> np.ndarray:
    """Small tie-aware rank implementation to avoid scipy dependency."""
    x = np.asarray(x)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and x[order[j + 1]] == x[order[i]]:
            j += 1
        avg = 0.5 * (i + j) + 1.0
        ranks[order[i : j + 1]] = avg
        i = j + 1
    return ranks


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx = rankdata(np.asarray(x))
    ry = rankdata(np.asarray(y))
    rx -= rx.mean()
    ry -= ry.mean()
    denom = np.sqrt(np.sum(rx * rx) * np.sum(ry * ry))
    return float(np.sum(rx * ry) / denom) if denom > 0 else float("nan")


def supports_from_defects(
    d_h: float,
    d_b: float,
    activity_pressure: float,
    spectral_coherence: float,
    activity_center: float = 1.0,
    activity_width: float = 0.85,
) -> dict[str, float]:
    """Convert operational fluid defects into MAAT support scores."""
    h = 1.0 / (1.0 + max(d_h, 0.0))
    b = 1.0 / (1.0 + max(d_b, 0.0))
    log_a = math.log(max(activity_pressure, EPS))
    log_c = math.log(max(activity_center, EPS))
    s_eff = math.exp(-0.5 * ((log_a - log_c) / activity_width) ** 2)
    v = float(np.clip(spectral_coherence, 0.0, 1.0))
    r_resp = (h * b * v) ** (1.0 / 3.0)
    r_rob = min(r_resp, (h * b * s_eff * v) ** 0.25)
    # Warning grows when activity/tail pressure rises while robustness degrades.
    warning = (1.0 - r_rob) * math.sqrt(max(activity_pressure, 0.0))
    return {
        "H": h,
        "B": b,
        "S_eff": s_eff,
        "V": v,
        "R_resp": r_resp,
        "R_rob": r_rob,
        "MAAT_warning": warning,
    }


def spectral_derivative_1d(u: np.ndarray, k: np.ndarray, order: int = 1) -> np.ndarray:
    return np.fft.ifft((1j * k) ** order * np.fft.fft(u)).real


def low_mode_fraction_1d(uhat: np.ndarray, k: np.ndarray, cutoff: int = 8) -> float:
    power = np.abs(uhat) ** 2
    total = float(np.sum(power)) + EPS
    return float(np.sum(power[np.abs(k) <= cutoff]) / total)


def run_burgers() -> tuple[pd.DataFrame, dict[str, float]]:
    n = 256
    viscosity = 5e-4
    dt = 5e-4
    t_end = 1.20
    sample_every = 10
    x = 2.0 * np.pi * np.arange(n) / n
    k = np.fft.fftfreq(n, d=1.0 / n)
    dealias = np.abs(k) <= n / 3

    u = np.sin(x) + 0.25 * np.sin(2.0 * x)
    ux0 = np.cos(x) + 0.5 * np.cos(2.0 * x)
    shock_time_inviscid = -1.0 / float(np.min(ux0))
    enstrophy_ref = float(np.mean(ux0**2))

    rows: list[dict[str, float]] = []
    prev_u = u.copy()
    prev_energy = 0.5 * float(np.mean(u * u))
    prev_t = 0.0

    steps = int(round(t_end / dt))
    for step in range(1, steps + 1):
        uhat = np.fft.fft(u)
        ux = spectral_derivative_1d(u, k, 1)
        nonlinear = -u * ux
        nhat = np.fft.fft(nonlinear)
        nhat[~dealias] = 0.0
        uhat_new = (uhat + dt * nhat) / (1.0 + dt * viscosity * k**2)
        uhat_new[~dealias] = 0.0
        u = np.fft.ifft(uhat_new).real

        if step % sample_every != 0 and step != steps:
            continue

        t = step * dt
        sample_dt = t - prev_t
        ux = spectral_derivative_1d(u, k, 1)
        uxx = spectral_derivative_1d(u, k, 2)
        rhs = -u * ux + viscosity * uxx
        ut = (u - prev_u) / sample_dt
        residual = ut - rhs
        rhs_rms = float(np.sqrt(np.mean(rhs * rhs))) + EPS
        d_h = float(np.sqrt(np.mean(residual * residual)) / rhs_rms)

        energy = 0.5 * float(np.mean(u * u))
        d_energy = (energy - prev_energy) / sample_dt
        diss = viscosity * float(np.mean(ux * ux))
        d_b = abs(d_energy + diss) / (abs(diss) + EPS)

        enstrophy = float(np.mean(ux * ux))
        max_grad = float(np.max(np.abs(ux)))
        activity_pressure = enstrophy / (enstrophy_ref + EPS)
        uhat = np.fft.fft(u)
        low_frac = low_mode_fraction_1d(uhat, k, cutoff=8)
        high_frac = 1.0 - low_frac
        supports = supports_from_defects(
            d_h=d_h,
            d_b=d_b,
            activity_pressure=activity_pressure,
            spectral_coherence=low_frac,
            activity_center=1.0,
            activity_width=0.80,
        )
        # High-k pressure captures shock-like spectral cascade.
        supports["MAAT_warning"] *= 1.0 + 2.0 * high_frac

        rows.append(
            {
                "case": "burgers_shock_warning",
                "t": t,
                "shock_time_inviscid": shock_time_inviscid,
                "energy": energy,
                "dissipation": diss,
                "equation_residual_defect": d_h,
                "energy_balance_defect": d_b,
                "enstrophy": enstrophy,
                "max_gradient": max_grad,
                "activity_pressure": activity_pressure,
                "low_mode_fraction": low_frac,
                "high_mode_fraction": high_frac,
                **supports,
            }
        )
        prev_u = u.copy()
        prev_energy = energy
        prev_t = t

    df = pd.DataFrame(rows)
    pre = df[df["t"] <= shock_time_inviscid]
    max_warning = float(df["MAAT_warning"].max())
    threshold = 0.7 * max_warning
    first_alert = df.loc[df["MAAT_warning"] >= threshold, "t"].min()
    lead_time = float(shock_time_inviscid - first_alert) if np.isfinite(first_alert) else float("nan")
    summary = {
        "shock_time_inviscid": shock_time_inviscid,
        "max_gradient_final": float(df["max_gradient"].iloc[-1]),
        "max_warning": max_warning,
        "first_70pct_warning_time": float(first_alert),
        "warning_lead_time_before_inviscid_shock": lead_time,
        "spearman_warning_vs_max_gradient": spearman(df["MAAT_warning"].values, df["max_gradient"].values),
        "spearman_pre_shock_warning_vs_time_to_shock_inverse": spearman(
            pre["MAAT_warning"].values,
            1.0 / np.maximum(pre["shock_time_inviscid"].values - pre["t"].values, 1e-3),
        ),
        "min_R_rob": float(df["R_rob"].min()),
    }
    return df, summary


def make_initial_vorticity(n: int, seed: int = 60) -> np.ndarray:
    rng = np.random.default_rng(seed)
    noise = rng.normal(size=(n, n))
    kh = np.fft.fftfreq(n, d=1.0 / n)
    kx, ky = np.meshgrid(kh, kh, indexing="ij")
    k2 = kx * kx + ky * ky
    filt = np.exp(-0.5 * k2 / 9.0)
    what = np.fft.fft2(noise) * filt
    what[0, 0] = 0.0
    omega = np.fft.ifft2(what).real
    omega /= np.std(omega) + EPS
    return omega


def low_mode_fraction_2d(fhat: np.ndarray, kx: np.ndarray, ky: np.ndarray, cutoff: int = 8) -> float:
    power = np.abs(fhat) ** 2
    total = float(np.sum(power)) + EPS
    return float(np.sum(power[(kx * kx + ky * ky) <= cutoff**2]) / total)


def run_2d_navier_stokes() -> tuple[pd.DataFrame, dict[str, float]]:
    n = 64
    viscosity = 2e-3
    dt = 0.004
    t_end = 4.0
    sample_every = 10
    k = np.fft.fftfreq(n, d=1.0 / n)
    kx, ky = np.meshgrid(k, k, indexing="ij")
    k2 = kx * kx + ky * ky
    k2_safe = k2.copy()
    k2_safe[0, 0] = 1.0
    dealias = (np.abs(kx) <= n / 3) & (np.abs(ky) <= n / 3)

    omega = make_initial_vorticity(n)
    omega_hat = np.fft.fft2(omega)
    omega_hat[~dealias] = 0.0
    grad_omega_x = np.fft.ifft2(1j * kx * omega_hat).real
    grad_omega_y = np.fft.ifft2(1j * ky * omega_hat).real
    pal_ref = 0.5 * float(np.mean(grad_omega_x**2 + grad_omega_y**2))

    rows: list[dict[str, float]] = []
    prev_omega = omega.copy()
    prev_z = 0.5 * float(np.mean(omega * omega))
    prev_t = 0.0

    steps = int(round(t_end / dt))
    for step in range(1, steps + 1):
        omega_hat = np.fft.fft2(omega)
        psi_hat = -omega_hat / k2_safe
        psi_hat[0, 0] = 0.0
        u = np.fft.ifft2(1j * ky * psi_hat).real
        v = np.fft.ifft2(-1j * kx * psi_hat).real
        ox = np.fft.ifft2(1j * kx * omega_hat).real
        oy = np.fft.ifft2(1j * ky * omega_hat).real
        jac = u * ox + v * oy
        jac_hat = np.fft.fft2(jac)
        jac_hat[~dealias] = 0.0
        omega_hat_new = (omega_hat - dt * jac_hat) / (1.0 + dt * viscosity * k2)
        omega_hat_new[~dealias] = 0.0
        omega_hat_new[0, 0] = 0.0
        omega = np.fft.ifft2(omega_hat_new).real

        if step % sample_every != 0 and step != steps:
            continue

        t = step * dt
        sample_dt = t - prev_t
        omega_hat = np.fft.fft2(omega)
        psi_hat = -omega_hat / k2_safe
        psi_hat[0, 0] = 0.0
        u = np.fft.ifft2(1j * ky * psi_hat).real
        v = np.fft.ifft2(-1j * kx * psi_hat).real
        ox = np.fft.ifft2(1j * kx * omega_hat).real
        oy = np.fft.ifft2(1j * ky * omega_hat).real
        lap_omega = np.fft.ifft2(-k2 * omega_hat).real
        rhs = -(u * ox + v * oy) + viscosity * lap_omega
        ot = (omega - prev_omega) / sample_dt
        residual = ot - rhs
        rhs_rms = float(np.sqrt(np.mean(rhs * rhs))) + EPS
        d_h = float(np.sqrt(np.mean(residual * residual)) / rhs_rms)

        energy = 0.5 * float(np.mean(u * u + v * v))
        enstrophy = 0.5 * float(np.mean(omega * omega))
        palinstrophy = 0.5 * float(np.mean(ox * ox + oy * oy))
        d_z = (enstrophy - prev_z) / sample_dt
        ens_diss = viscosity * float(np.mean(ox * ox + oy * oy))
        d_b = abs(d_z + ens_diss) / (abs(ens_diss) + EPS)

        activity_pressure = palinstrophy / (pal_ref + EPS)
        low_frac = low_mode_fraction_2d(omega_hat, kx, ky, cutoff=8)
        high_frac = 1.0 - low_frac
        supports = supports_from_defects(
            d_h=d_h,
            d_b=d_b,
            activity_pressure=activity_pressure,
            spectral_coherence=low_frac,
            activity_center=1.0,
            activity_width=0.90,
        )
        supports["MAAT_warning"] *= 1.0 + high_frac

        rows.append(
            {
                "case": "2d_navier_stokes_regular_control",
                "t": t,
                "energy": energy,
                "enstrophy": enstrophy,
                "palinstrophy": palinstrophy,
                "equation_residual_defect": d_h,
                "enstrophy_balance_defect": d_b,
                "activity_pressure": activity_pressure,
                "low_mode_fraction": low_frac,
                "high_mode_fraction": high_frac,
                **supports,
            }
        )
        prev_omega = omega.copy()
        prev_z = enstrophy
        prev_t = t

    df = pd.DataFrame(rows)
    summary = {
        "initial_enstrophy": float(df["enstrophy"].iloc[0]),
        "final_enstrophy": float(df["enstrophy"].iloc[-1]),
        "max_palinstrophy": float(df["palinstrophy"].max()),
        "max_warning": float(df["MAAT_warning"].max()),
        "min_R_rob": float(df["R_rob"].min()),
        "spearman_warning_vs_palinstrophy": spearman(df["MAAT_warning"].values, df["palinstrophy"].values),
        "false_high_warning_fraction_gt_1": float(np.mean(df["MAAT_warning"].values > 1.0)),
    }
    return df, summary


def plot_burgers(df: pd.DataFrame) -> None:
    fig, ax1 = plt.subplots(figsize=(10, 5.6))
    ax1.plot(df["t"], df["max_gradient"], color="#b23a48", lw=2.0, label="max |u_x|")
    ax1.set_xlabel("time")
    ax1.set_ylabel("max gradient", color="#b23a48")
    ax1.tick_params(axis="y", labelcolor="#b23a48")
    ax1.axvline(df["shock_time_inviscid"].iloc[0], color="black", ls="--", lw=1.3, label="inviscid shock time")
    ax2 = ax1.twinx()
    ax2.plot(df["t"], df["MAAT_warning"], color="#1d7874", lw=2.0, label="MAAT warning")
    ax2.plot(df["t"], 1.0 - df["R_rob"], color="#f4a261", lw=1.8, alpha=0.9, label="1 - R_rob")
    ax2.set_ylabel("warning / robustness loss", color="#1d7874")
    ax2.tick_params(axis="y", labelcolor="#1d7874")
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2, loc="upper left")
    ax1.grid(alpha=0.25)
    fig.suptitle("Burgers shock warning: gradient growth and MAAT diagnostics")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig1_burgers_warning_timeseries.png", dpi=180)
    plt.close(fig)


def plot_ns2d(df: pd.DataFrame) -> None:
    fig, ax1 = plt.subplots(figsize=(10, 5.6))
    ax1.plot(df["t"], df["enstrophy"], color="#31572c", lw=2.0, label="enstrophy")
    ax1.plot(df["t"], df["palinstrophy"], color="#90a955", lw=2.0, label="palinstrophy")
    ax1.set_xlabel("time")
    ax1.set_ylabel("fluid activity")
    ax1.grid(alpha=0.25)
    ax2 = ax1.twinx()
    ax2.plot(df["t"], df["MAAT_warning"], color="#5e548e", lw=2.0, label="MAAT warning")
    ax2.plot(df["t"], df["R_rob"], color="#00a6a6", lw=1.8, label="R_rob")
    ax2.set_ylabel("warning / robustness")
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2, loc="upper right")
    fig.suptitle("2D Navier-Stokes regular control: decaying activity and bounded warning")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig2_ns2d_regular_control.png", dpi=180)
    plt.close(fig)


def plot_scatter(burgers: pd.DataFrame, ns2d: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    axes[0].scatter(burgers["max_gradient"], burgers["MAAT_warning"], s=18, color="#b23a48", alpha=0.8)
    axes[0].set_xlabel("max |u_x|")
    axes[0].set_ylabel("MAAT warning")
    axes[0].set_title("Burgers: warning tracks steepening")
    axes[0].grid(alpha=0.25)
    axes[1].scatter(ns2d["palinstrophy"], ns2d["MAAT_warning"], s=18, color="#5e548e", alpha=0.8)
    axes[1].set_xlabel("palinstrophy")
    axes[1].set_ylabel("MAAT warning")
    axes[1].set_title("2D NS: bounded control response")
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig3_warning_scatter.png", dpi=180)
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    burgers, burgers_summary = run_burgers()
    ns2d, ns2d_summary = run_2d_navier_stokes()

    burgers.to_csv(OUTDIR / "paper60_burgers_diagnostics.csv", index=False)
    ns2d.to_csv(OUTDIR / "paper60_ns2d_diagnostics.csv", index=False)

    plot_burgers(burgers)
    plot_ns2d(ns2d)
    plot_scatter(burgers, ns2d)

    summary = {
        "paper": 60,
        "title": "MAAT Fluid Coherence and Blow-up Diagnostics",
        "status": "toy diagnostic benchmark, not a Navier-Stokes regularity proof",
        "burgers": burgers_summary,
        "navier_stokes_2d": ns2d_summary,
    }
    with open(OUTDIR / "paper60_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
