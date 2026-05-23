#!/usr/bin/env python3
"""
Extra Phenomenological Paper - MAAT v1.6 Navier-Stokes regularity attempt.

This is not a proof of the Navier-Stokes regularity problem.  It is a bounded
numerical programme:

* solve 3D incompressible Navier-Stokes on a periodic box;
* monitor Beale-Kato-Majda (BKM) risk via integral ||omega||_inf dt;
* compute MAAT v1.6 structural supports H, B, S, V, R along the trajectory;
* evaluate the symbolic ToE_MAAT-inspired structural action

  ToE_MAAT = integral [ (H+B+S+V+R) * Z * cascade_resistance ]
                      / [ DeltaE + DeltaQ + D0 ] dt.

The cascade-resistance factor is operationalised as a spectral-tail and
vortex-stretching resistance term.  Higher action means more coherent bounded
evolution; higher warning means stronger structural risk.
"""

from __future__ import annotations

import json
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

_cache_dir = Path(tempfile.gettempdir()) / "maat_navier_stokes_phenomenological_cache"
_cache_dir.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_cache_dir))
os.environ.setdefault("XDG_CACHE_HOME", str(_cache_dir))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


EPS = 1e-12
OUTDIR = Path("outputs_phenomenological")
FIXED_BASELINE_PENALTY = 11.0


@dataclass(frozen=True)
class Scenario:
    name: str
    viscosity: float
    dt: float
    t_end: float
    amplitude: float


def rankdata(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and x[order[j + 1]] == x[order[i]]:
            j += 1
        ranks[order[i : j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return ranks


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx = rankdata(np.asarray(x, dtype=float))
    ry = rankdata(np.asarray(y, dtype=float))
    rx -= rx.mean()
    ry -= ry.mean()
    denom = np.sqrt(np.sum(rx * rx) * np.sum(ry * ry))
    return float(np.sum(rx * ry) / denom) if denom > 0 else float("nan")


def make_grid(n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x = 2.0 * np.pi * np.arange(n) / n
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
    k = np.fft.fftfreq(n, d=1.0 / n)
    return xx, yy, zz, k


def taylor_green_initial(n: int, amplitude: float) -> np.ndarray:
    x, y, z, _ = make_grid(n)
    u = np.empty((3, n, n, n), dtype=float)
    u[0] = amplitude * np.sin(x) * np.cos(y) * np.cos(z)
    u[1] = -amplitude * np.cos(x) * np.sin(y) * np.cos(z)
    u[2] = 0.0
    return u


def spectral_setup(n: int) -> dict[str, np.ndarray]:
    k = np.fft.fftfreq(n, d=1.0 / n)
    kx, ky, kz = np.meshgrid(k, k, k, indexing="ij")
    k2 = kx * kx + ky * ky + kz * kz
    k2_safe = k2.copy()
    k2_safe[0, 0, 0] = 1.0
    dealias = (np.abs(kx) <= n / 3) & (np.abs(ky) <= n / 3) & (np.abs(kz) <= n / 3)
    kmag = np.sqrt(k2)
    return {
        "kx": kx,
        "ky": ky,
        "kz": kz,
        "k2": k2,
        "k2_safe": k2_safe,
        "dealias": dealias,
        "kmag": kmag,
    }


def fft_vec(u: np.ndarray) -> np.ndarray:
    return np.stack([np.fft.fftn(u[i]) for i in range(3)], axis=0)


def ifft_vec(uh: np.ndarray) -> np.ndarray:
    return np.stack([np.fft.ifftn(uh[i]).real for i in range(3)], axis=0)


def project_incompressible(uh: np.ndarray, spec: dict[str, np.ndarray]) -> np.ndarray:
    kx, ky, kz, k2 = spec["kx"], spec["ky"], spec["kz"], spec["k2_safe"]
    divh = kx * uh[0] + ky * uh[1] + kz * uh[2]
    out = uh.copy()
    out[0] -= kx * divh / k2
    out[1] -= ky * divh / k2
    out[2] -= kz * divh / k2
    out[:, 0, 0, 0] = 0.0
    return out


def deriv_scalar(fh: np.ndarray, axis: int, spec: dict[str, np.ndarray]) -> np.ndarray:
    kk = [spec["kx"], spec["ky"], spec["kz"]][axis]
    return np.fft.ifftn(1j * kk * fh).real


def gradient_vector(u: np.ndarray, spec: dict[str, np.ndarray]) -> np.ndarray:
    uh = fft_vec(u)
    grad = np.empty((3, 3) + u.shape[1:], dtype=float)
    for comp in range(3):
        for axis in range(3):
            grad[comp, axis] = deriv_scalar(uh[comp], axis, spec)
    return grad


def curl(u: np.ndarray, spec: dict[str, np.ndarray]) -> np.ndarray:
    uh = fft_vec(u)
    ux_y = deriv_scalar(uh[0], 1, spec)
    ux_z = deriv_scalar(uh[0], 2, spec)
    uy_x = deriv_scalar(uh[1], 0, spec)
    uy_z = deriv_scalar(uh[1], 2, spec)
    uz_x = deriv_scalar(uh[2], 0, spec)
    uz_y = deriv_scalar(uh[2], 1, spec)
    return np.stack([uz_y - uy_z, ux_z - uz_x, uy_x - ux_y], axis=0)


def divergence(u: np.ndarray, spec: dict[str, np.ndarray]) -> np.ndarray:
    uh = fft_vec(u)
    divh = 1j * (spec["kx"] * uh[0] + spec["ky"] * uh[1] + spec["kz"] * uh[2])
    return np.fft.ifftn(divh).real


def rhs_navier_stokes(u: np.ndarray, viscosity: float, spec: dict[str, np.ndarray]) -> np.ndarray:
    grad = gradient_vector(u, spec)
    adv = np.empty_like(u)
    for comp in range(3):
        adv[comp] = u[0] * grad[comp, 0] + u[1] * grad[comp, 1] + u[2] * grad[comp, 2]
    nh = fft_vec(-adv)
    nh[:, ~spec["dealias"]] = 0.0
    nh = project_incompressible(nh, spec)
    uh = fft_vec(u)
    visc = ifft_vec(-viscosity * spec["k2"][None, ...] * uh)
    return ifft_vec(nh) + visc


def step(u: np.ndarray, viscosity: float, dt: float, spec: dict[str, np.ndarray]) -> np.ndarray:
    grad = gradient_vector(u, spec)
    adv = np.empty_like(u)
    for comp in range(3):
        adv[comp] = u[0] * grad[comp, 0] + u[1] * grad[comp, 1] + u[2] * grad[comp, 2]
    nh = fft_vec(-adv)
    nh[:, ~spec["dealias"]] = 0.0
    nh = project_incompressible(nh, spec)
    uh = fft_vec(u)
    uh_new = (uh + dt * nh) / (1.0 + dt * viscosity * spec["k2"][None, ...])
    uh_new[:, ~spec["dealias"]] = 0.0
    uh_new = project_incompressible(uh_new, spec)
    return ifft_vec(uh_new)


def low_mode_fraction(uh: np.ndarray, spec: dict[str, np.ndarray], cutoff: int = 5) -> float:
    power = np.sum(np.abs(uh) ** 2, axis=0)
    total = float(np.sum(power)) + EPS
    return float(np.sum(power[spec["kmag"] <= cutoff]) / total)


def high_tail_fraction(uh: np.ndarray, spec: dict[str, np.ndarray], cutoff_ratio: float = 0.45) -> float:
    power = np.sum(np.abs(uh) ** 2, axis=0)
    total = float(np.sum(power)) + EPS
    kmax = float(np.max(spec["kmag"][spec["dealias"]]))
    return float(np.sum(power[spec["kmag"] >= cutoff_ratio * kmax]) / total)


def compute_diagnostics(
    u: np.ndarray,
    u_prev: np.ndarray,
    dt_sample: float,
    viscosity: float,
    spec: dict[str, np.ndarray],
    prev_energy: float,
    pal_ref: float,
) -> dict[str, float]:
    rhs = rhs_navier_stokes(u, viscosity, spec)
    ut = (u - u_prev) / dt_sample
    residual = ut - rhs
    d_eq = float(np.sqrt(np.mean(residual * residual)) / (np.sqrt(np.mean(rhs * rhs)) + EPS))
    div = divergence(u, spec)
    grad = gradient_vector(u, spec)
    grad_norm = float(np.sqrt(np.mean(grad * grad))) + EPS
    d_div = float(np.sqrt(np.mean(div * div)) / grad_norm)

    energy = 0.5 * float(np.mean(np.sum(u * u, axis=0)))
    dissipation = viscosity * float(np.mean(np.sum(grad * grad, axis=(0, 1))))
    d_energy = (energy - prev_energy) / dt_sample
    delta_e = abs(d_energy + dissipation) / (abs(dissipation) + EPS)

    omega = curl(u, spec)
    omega_mag = np.sqrt(np.sum(omega * omega, axis=0))
    omega_inf = float(np.max(omega_mag))
    enstrophy = 0.5 * float(np.mean(omega_mag * omega_mag))

    # Vortex stretching: omega dot grad(u).  This is the true 3D risk channel.
    stretch_vec = np.empty_like(u)
    for comp in range(3):
        stretch_vec[comp] = (
            omega[0] * grad[comp, 0] + omega[1] * grad[comp, 1] + omega[2] * grad[comp, 2]
        )
    stretch_mag = np.sqrt(np.sum(stretch_vec * stretch_vec, axis=0))
    stretch_mean = float(np.mean(stretch_mag))
    stretch_inf = float(np.max(stretch_mag))

    # Palinstrophy-like gradient of vorticity.
    wh = fft_vec(omega)
    wx = ifft_vec(1j * spec["kx"][None, ...] * wh)
    wy = ifft_vec(1j * spec["ky"][None, ...] * wh)
    wz = ifft_vec(1j * spec["kz"][None, ...] * wh)
    pal = 0.5 * float(np.mean(np.sum(wx * wx + wy * wy + wz * wz, axis=0)))
    activity_pressure = pal / (pal_ref + EPS)

    uh = fft_vec(u)
    low_frac = low_mode_fraction(uh, spec, cutoff=5)
    tail_frac = high_tail_fraction(uh, spec)

    h = 1.0 / (1.0 + d_eq + d_div)
    b = 1.0 / (1.0 + delta_e)
    # Bounded activity: healthy around A=1, stressed when activity is too small or too large.
    s_eff = math.exp(-0.5 * (math.log(max(activity_pressure, EPS)) / 1.05) ** 2)
    v = float(np.clip(low_frac, 0.0, 1.0))
    r_resp = (h * b * v) ** (1.0 / 3.0)
    r_rob = min(r_resp, (h * b * s_eff * v) ** 0.25)

    stretch_pressure = stretch_mean / (enstrophy + EPS)
    delta_q = tail_frac + 0.1 * d_div
    z_partition = math.exp(-(1.0 - r_rob)) * (0.5 + 0.5 * v)
    cascade_resistance = 1.0 / (1.0 + 2.0 * tail_frac + 0.15 * stretch_pressure)
    toe_integrand = (
        (h + b + s_eff + v + r_rob) * z_partition * cascade_resistance
    ) / (delta_e + delta_q + FIXED_BASELINE_PENALTY + EPS)
    maat_warning = (1.0 - r_rob) * math.sqrt(max(activity_pressure, 0.0)) * (1.0 + tail_frac) * (
        1.0 + 0.15 * stretch_pressure
    )

    return {
        "energy": energy,
        "dissipation": dissipation,
        "equation_residual_defect": d_eq,
        "divergence_defect": d_div,
        "energy_balance_defect": delta_e,
        "omega_inf": omega_inf,
        "enstrophy": enstrophy,
        "palinstrophy": pal,
        "vortex_stretch_mean": stretch_mean,
        "vortex_stretch_inf": stretch_inf,
        "stretch_pressure": stretch_pressure,
        "activity_pressure": activity_pressure,
        "low_mode_fraction": low_frac,
        "spectral_tail_fraction": tail_frac,
        "DeltaE": delta_e,
        "DeltaQ": delta_q,
        "H": h,
        "B": b,
        "S_eff": s_eff,
        "V": v,
        "R_resp": r_resp,
        "R_rob": r_rob,
        "Z_partition": z_partition,
        "cascade_resistance": cascade_resistance,
        "ToE_MAAT_integrand": toe_integrand,
        "MAAT_warning": maat_warning,
    }


def simulate_scenario(scenario: Scenario, n: int = 24, sample_every: int = 5) -> tuple[pd.DataFrame, dict[str, float]]:
    spec = spectral_setup(n)
    u = taylor_green_initial(n, scenario.amplitude)
    uh = project_incompressible(fft_vec(u), spec)
    uh[:, ~spec["dealias"]] = 0.0
    u = ifft_vec(uh)

    omega0 = curl(u, spec)
    wh0 = fft_vec(omega0)
    wx0 = ifft_vec(1j * spec["kx"][None, ...] * wh0)
    wy0 = ifft_vec(1j * spec["ky"][None, ...] * wh0)
    wz0 = ifft_vec(1j * spec["kz"][None, ...] * wh0)
    pal_ref = 0.5 * float(np.mean(np.sum(wx0 * wx0 + wy0 * wy0 + wz0 * wz0, axis=0))) + EPS

    prev_u = u.copy()
    prev_energy = 0.5 * float(np.mean(np.sum(u * u, axis=0)))
    prev_t = 0.0
    bkm_integral = 0.0
    toe_action = 0.0
    warning_action = 0.0
    rows: list[dict[str, float]] = []

    steps = int(round(scenario.t_end / scenario.dt))
    for step_idx in range(1, steps + 1):
        u = step(u, scenario.viscosity, scenario.dt, spec)
        if step_idx % sample_every != 0 and step_idx != steps:
            continue
        t = step_idx * scenario.dt
        dt_sample = t - prev_t
        diag = compute_diagnostics(u, prev_u, dt_sample, scenario.viscosity, spec, prev_energy, pal_ref)
        bkm_integral += diag["omega_inf"] * dt_sample
        toe_action += diag["ToE_MAAT_integrand"] * dt_sample
        warning_action += diag["MAAT_warning"] * dt_sample
        rows.append(
            {
                "scenario": scenario.name,
                "t": t,
                "viscosity": scenario.viscosity,
                "amplitude": scenario.amplitude,
                "BKM_integral": bkm_integral,
                "ToE_MAAT_action": toe_action,
                "warning_action": warning_action,
                **diag,
            }
        )
        prev_u = u.copy()
        prev_energy = diag["energy"]
        prev_t = t

    df = pd.DataFrame(rows)
    summary = {
        "scenario": scenario.name,
        "viscosity": scenario.viscosity,
        "amplitude": scenario.amplitude,
        "final_energy": float(df["energy"].iloc[-1]),
        "max_omega_inf": float(df["omega_inf"].max()),
        "final_BKM_integral": float(df["BKM_integral"].iloc[-1]),
        "min_R_rob": float(df["R_rob"].min()),
        "max_warning": float(df["MAAT_warning"].max()),
        "final_ToE_MAAT_action": float(df["ToE_MAAT_action"].iloc[-1]),
        "final_warning_action": float(df["warning_action"].iloc[-1]),
        "max_spectral_tail_fraction": float(df["spectral_tail_fraction"].max()),
        "spearman_warning_vs_omega_inf": spearman(df["MAAT_warning"].values, df["omega_inf"].values),
        "spearman_warning_vs_stretch": spearman(df["MAAT_warning"].values, df["vortex_stretch_mean"].values),
        "spearman_toe_integrand_vs_Rrob": spearman(df["ToE_MAAT_integrand"].values, df["R_rob"].values),
    }
    return df, summary


def plot_timeseries(all_df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for name, sub in all_df.groupby("scenario"):
        axes[0, 0].plot(sub["t"], sub["omega_inf"], lw=2, label=name)
        axes[0, 1].plot(sub["t"], sub["BKM_integral"], lw=2, label=name)
        axes[1, 0].plot(sub["t"], sub["MAAT_warning"], lw=2, label=name)
        axes[1, 1].plot(sub["t"], sub["R_rob"], lw=2, label=name)
    axes[0, 0].set_title("Vorticity infinity norm")
    axes[0, 1].set_title("BKM integral proxy")
    axes[1, 0].set_title("MAAT warning")
    axes[1, 1].set_title("Robustness closure")
    for ax in axes.flat:
        ax.set_xlabel("time")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig1_bkm_warning_timeseries.png", dpi=180)
    plt.close(fig)


def plot_action(all_df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for name, sub in all_df.groupby("scenario"):
        axes[0].plot(sub["t"], sub["ToE_MAAT_action"], lw=2, label=name)
        axes[1].plot(sub["t"], sub["ToE_MAAT_integrand"], lw=2, label=name)
    axes[0].set_title("Integrated ToE_MAAT structural action")
    axes[1].set_title("Instantaneous structural-action integrand")
    for ax in axes:
        ax.set_xlabel("time")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig2_toe_maat_action.png", dpi=180)
    plt.close(fig)


def plot_phase(all_df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7, 5.5))
    for name, sub in all_df.groupby("scenario"):
        ax.scatter(sub["omega_inf"], sub["MAAT_warning"], s=18, alpha=0.75, label=name)
    ax.set_xlabel("||omega||_inf")
    ax.set_ylabel("MAAT warning")
    ax.set_title("Warning vs BKM risk channel")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig3_warning_vs_vorticity.png", dpi=180)
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    scenarios = [
        Scenario("taylor_green_moderate_viscosity", viscosity=0.02, dt=0.006, t_end=2.4, amplitude=1.0),
        Scenario("taylor_green_low_viscosity_high_stress", viscosity=0.006, dt=0.003, t_end=1.8, amplitude=1.25),
    ]

    frames: list[pd.DataFrame] = []
    summaries: list[dict[str, float]] = []
    for scenario in scenarios:
        df, summary = simulate_scenario(scenario)
        frames.append(df)
        summaries.append(summary)

    all_df = pd.concat(frames, ignore_index=True)
    all_df.to_csv(OUTDIR / "navier_stokes_maat_diagnostics.csv", index=False)
    pd.DataFrame(summaries).to_csv(OUTDIR / "summary_by_scenario.csv", index=False)

    plot_timeseries(all_df)
    plot_action(all_df)
    plot_phase(all_df)

    summary = {
        "paper": "extra_phenomenological",
        "title": "BKM-Aware MAAT Diagnostics for Navier-Stokes",
        "status": "3D spectral diagnostic programme, not a proof",
        "symbolic_formula": "ToE_MAAT = integral [ (H+B+S+V+R)*Z*cascade_resistance ] / (DeltaE+DeltaQ+D0) dt",
        "D0_fixed_baseline_penalty": FIXED_BASELINE_PENALTY,
        "scenarios": summaries,
    }
    with open(OUTDIR / "navier_stokes_maat_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
