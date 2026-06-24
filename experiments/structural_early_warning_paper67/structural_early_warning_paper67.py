#!/usr/bin/env python3
"""Paper 67: structural early warning as a refinement trigger.

This benchmark asks a utility question rather than a diagnostic-only question:
at the same false-alarm rate, does a MAAT-style structural warning functional
trigger adaptive refinement earlier than standard vorticity monitors?

The implemented stage is a local forced 2D Navier--Stokes vorticity simulation.
JHTDB is treated as a predeclared external replication route; no JHTDB data are
downloaded or redistributed by this script.
"""

from __future__ import annotations

import json
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

_CACHE = Path(tempfile.gettempdir()) / "maat_paper67_matplotlib_cache"
_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


EPS = 1.0e-12
SEED = 67067
OUTDIR = Path("outputs_paper67")
BOOTSTRAP_SAMPLES = 10000


def rankdata(x: np.ndarray) -> np.ndarray:
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
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(x) < 4 or len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
        return float("nan")
    rx = rankdata(x)
    ry = rankdata(y)
    rx -= rx.mean()
    ry -= ry.mean()
    denom = float(np.sqrt(np.sum(rx * rx) * np.sum(ry * ry)))
    return float(np.sum(rx * ry) / denom) if denom > 0 else float("nan")


def normalize01(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    lo = float(np.nanmin(x))
    hi = float(np.nanmax(x))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo < EPS:
        return np.zeros_like(x)
    return (x - lo) / (hi - lo)


def to_jsonable(obj):
    """Convert numpy scalars and NaN/Inf values to strict JSON-compatible data."""
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, tuple):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        val = float(obj)
        return val if np.isfinite(val) else None
    return obj


def rolling_future_max(values: np.ndarray, steps: int) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.empty_like(values)
    for i in range(len(values)):
        j = min(len(values), i + steps + 1)
        out[i] = float(np.max(values[i:j]))
    return out


def event_onsets(mask: np.ndarray) -> list[int]:
    mask = np.asarray(mask, dtype=bool)
    onsets: list[int] = []
    prev = False
    for i, val in enumerate(mask):
        if val and not prev:
            onsets.append(i)
        prev = bool(val)
    return onsets


def first_alert_leads(
    times: np.ndarray,
    alert: np.ndarray,
    onsets: list[int],
    lookback_steps: int,
) -> list[float]:
    leads: list[float] = []
    for idx in onsets:
        lo = max(0, idx - lookback_steps)
        alert_idx = np.where(alert[lo : idx + 1])[0]
        if len(alert_idx) == 0:
            leads.append(float("nan"))
            continue
        first = lo + int(alert_idx[0])
        leads.append(float(times[idx] - times[first]))
    return leads


def threshold_at_false_alarm(
    values: np.ndarray,
    non_event_mask: np.ndarray,
    false_alarm_rate: float,
) -> float:
    pool = np.asarray(values, dtype=float)[non_event_mask]
    if len(pool) == 0:
        return float("inf")
    return float(np.quantile(pool, 1.0 - false_alarm_rate))


@dataclass
class Grid:
    n: int
    kx: np.ndarray
    ky: np.ndarray
    k2: np.ndarray
    k2_safe: np.ndarray
    dealias: np.ndarray
    x: np.ndarray
    y: np.ndarray


def make_grid(n: int) -> Grid:
    k = np.fft.fftfreq(n, d=1.0 / n)
    kx, ky = np.meshgrid(k, k, indexing="ij")
    k2 = kx * kx + ky * ky
    k2_safe = k2.copy()
    k2_safe[0, 0] = 1.0
    dealias = (np.abs(kx) <= n / 3) & (np.abs(ky) <= n / 3)
    x = 2.0 * np.pi * np.arange(n) / n
    y = 2.0 * np.pi * np.arange(n) / n
    return Grid(n=n, kx=kx, ky=ky, k2=k2, k2_safe=k2_safe, dealias=dealias, x=x, y=y)


def make_initial_vorticity(grid: Grid, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    noise = rng.normal(size=(grid.n, grid.n))
    filt = np.exp(-0.5 * grid.k2 / 12.0)
    what = np.fft.fft2(noise) * filt
    what[0, 0] = 0.0
    omega = np.fft.ifft2(what).real
    omega /= np.std(omega) + EPS
    return 0.25 * omega


def forcing_field(grid: Grid, t: float) -> np.ndarray:
    """Low-mode deterministic forcing with intermittent amplitude pulses."""
    x, y = np.meshgrid(grid.x, grid.y, indexing="ij")
    base = (
        np.sin(4 * x + 2 * y + 0.7 * t)
        + 0.75 * np.cos(3 * x - 4 * y - 0.4 * t)
        + 0.55 * np.sin(5 * x - y + 0.2 * t)
    )
    pulses = (
        1.0
        + 0.9 * np.exp(-0.5 * ((t - 8.0) / 0.90) ** 2)
        + 0.8 * np.exp(-0.5 * ((t - 16.0) / 1.00) ** 2)
        + 0.9 * np.exp(-0.5 * ((t - 24.0) / 0.95) ** 2)
    )
    return 0.050 * pulses * base


def low_mode_fraction(fhat: np.ndarray, grid: Grid, cutoff: int = 8) -> float:
    power = np.abs(fhat) ** 2
    total = float(np.sum(power)) + EPS
    return float(np.sum(power[grid.k2 <= cutoff**2]) / total)


def high_mode_fraction(fhat: np.ndarray, grid: Grid, cutoff: int = 16) -> float:
    power = np.abs(fhat) ** 2
    total = float(np.sum(power)) + EPS
    return float(np.sum(power[grid.k2 >= cutoff**2]) / total)


def supports_from_snapshot(
    d_h: float,
    d_b: float,
    activity_pressure: float,
    low_fraction: float,
    high_fraction: float,
    tail_growth: float,
) -> dict[str, float]:
    h = 1.0 / (1.0 + max(d_h, 0.0))
    b = 1.0 / (1.0 + max(d_b, 0.0))
    log_a = math.log(max(activity_pressure, EPS))
    s_eff = math.exp(-0.5 * (log_a / 0.95) ** 2)
    v = float(np.clip(low_fraction, 0.0, 1.0))
    r_resp = (h * b * v) ** (1.0 / 3.0)
    r_rob = min(r_resp, (h * b * s_eff * v) ** 0.25)
    high_tail_pressure = max(high_fraction, 0.0) / 1.0e-4
    # W_MAAT is intentionally current-state only: it uses the present support
    # state plus a backward finite-difference tail-growth pressure.
    w_maat = (1.0 - r_rob) * math.sqrt(max(activity_pressure, 0.0))
    w_maat *= 1.0 + 0.7 * high_tail_pressure
    w_maat *= 1.0 + 2.0 * max(tail_growth, 0.0)
    return {
        "H": h,
        "B": b,
        "S_eff": s_eff,
        "V": v,
        "R_resp": r_resp,
        "R_rob": r_rob,
        "high_tail_pressure": high_tail_pressure,
        "W_MAAT": w_maat,
    }


def add_warning_ablation_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add causal warning variants used for the component ablation.

    The ablations keep the same calibration protocol as the main monitor.  They
    are intentionally simple: the question is which part of the structural
    warning is carrying the utility in this pilot.
    """
    work = df.copy()
    robustness_loss = np.clip(1.0 - work["R_rob"].values, 0.0, None)
    activity = np.sqrt(np.clip(work["activity_pressure"].values, 0.0, None))
    tail = np.clip(work["high_tail_pressure"].values, 0.0, None)
    growth = np.clip(work["tail_growth_pressure"].values, 0.0, None)
    work["W_robustness_only"] = robustness_loss
    work["W_activity_only"] = activity
    work["W_tail_only"] = tail
    work["W_growth_only"] = growth
    work["W_no_tail"] = robustness_loss * activity * (1.0 + 2.0 * growth)
    work["W_no_growth"] = robustness_loss * activity * (1.0 + 0.7 * tail)
    work["W_no_robustness"] = activity * (1.0 + 0.7 * tail) * (1.0 + 2.0 * growth)
    return work


def simulate_forced_2d(seed: int = SEED, run_id: int = 0) -> pd.DataFrame:
    grid = make_grid(64)
    viscosity = 2.0e-3
    drag = 4.0e-2
    dt = 0.005
    t_end = 30.0
    sample_every = 10
    omega = make_initial_vorticity(grid, seed)
    prev_omega = omega.copy()
    prev_enstrophy = 0.5 * float(np.mean(omega * omega))
    prev_high = 0.0
    prev_t = 0.0

    rows: list[dict[str, float]] = []
    steps = int(round(t_end / dt))
    for step in range(1, steps + 1):
        t = step * dt
        omega_hat = np.fft.fft2(omega)
        psi_hat = omega_hat / grid.k2_safe
        psi_hat[0, 0] = 0.0
        u = np.fft.ifft2(1j * grid.ky * psi_hat).real
        v = -np.fft.ifft2(1j * grid.kx * psi_hat).real
        ox = np.fft.ifft2(1j * grid.kx * omega_hat).real
        oy = np.fft.ifft2(1j * grid.ky * omega_hat).real
        adv_hat = np.fft.fft2(u * ox + v * oy)
        force = forcing_field(grid, t)
        force_hat = np.fft.fft2(force)
        rhs_hat = -adv_hat + force_hat
        rhs_hat[~grid.dealias] = 0.0
        next_hat = (omega_hat + dt * rhs_hat) / (1.0 + dt * (viscosity * grid.k2 + drag))
        next_hat[~grid.dealias] = 0.0
        next_hat[0, 0] = 0.0
        omega = np.fft.ifft2(next_hat).real

        if step % sample_every != 0 and step != steps:
            continue

        sample_dt = t - prev_t
        omega_hat = np.fft.fft2(omega)
        psi_hat = omega_hat / grid.k2_safe
        psi_hat[0, 0] = 0.0
        u = np.fft.ifft2(1j * grid.ky * psi_hat).real
        v = -np.fft.ifft2(1j * grid.kx * psi_hat).real
        ox = np.fft.ifft2(1j * grid.kx * omega_hat).real
        oy = np.fft.ifft2(1j * grid.ky * omega_hat).real
        lap = np.fft.ifft2(-grid.k2 * omega_hat).real
        force = forcing_field(grid, t)
        rhs = -(u * ox + v * oy) + viscosity * lap - drag * omega + force
        ot = (omega - prev_omega) / sample_dt
        residual = ot - rhs
        rhs_rms = float(np.sqrt(np.mean(rhs * rhs))) + EPS
        d_h = float(np.sqrt(np.mean(residual * residual)) / rhs_rms)

        enstrophy = 0.5 * float(np.mean(omega * omega))
        palinstrophy = 0.5 * float(np.mean(ox * ox + oy * oy))
        injection = float(np.mean(omega * force))
        dZ = (enstrophy - prev_enstrophy) / sample_dt
        balance_resid = dZ - injection + 2.0 * viscosity * palinstrophy + 2.0 * drag * enstrophy
        d_b = abs(balance_resid) / (abs(injection) + 2.0 * viscosity * palinstrophy + 2.0 * drag * enstrophy + EPS)

        low = low_mode_fraction(omega_hat, grid, cutoff=8)
        high = high_mode_fraction(omega_hat, grid, cutoff=16)
        activity_pressure = palinstrophy / 35.0
        tail_growth = max(0.0, (high - prev_high) / max(sample_dt, EPS)) / 0.05
        supports = supports_from_snapshot(d_h, d_b, activity_pressure, low, high, tail_growth)
        high_tail_pressure = max(high, 0.0) / 1.0e-4
        refinement_oracle = math.sqrt(max(palinstrophy, 0.0)) * (1.0 + 3.0 * high_tail_pressure)
        rows.append(
            {
                "run_id": run_id,
                "t": t,
                "tau": t + run_id * (t_end + 3.0),
                "enstrophy": enstrophy,
                "palinstrophy": palinstrophy,
                "max_abs_vorticity": float(np.max(np.abs(omega))),
                "rms_vorticity": float(np.sqrt(np.mean(omega * omega))),
                "low_mode_fraction": low,
                "high_mode_fraction": high,
                "tail_growth_pressure": tail_growth,
                "activity_pressure": activity_pressure,
                "equation_residual_defect": d_h,
                "enstrophy_balance_defect": d_b,
                "refinement_oracle": refinement_oracle,
                **supports,
            }
        )
        prev_omega = omega.copy()
        prev_enstrophy = enstrophy
        prev_high = high
        prev_t = t

    df = pd.DataFrame(rows)
    return df


def evaluate_monitor_rows(
    work: pd.DataFrame,
    monitors: dict[str, str],
    event_mask: np.ndarray,
    onset_idx: list[int],
    non_event: np.ndarray,
    false_alarm_rate: float,
    lookback_steps: int,
    label: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    lead_rows: list[dict[str, object]] = []
    for name, col in monitors.items():
        vals = work[col].values
        threshold = threshold_at_false_alarm(vals, non_event, false_alarm_rate)
        alert = vals >= threshold
        false_alarm_observed = float(np.mean(alert[non_event])) if np.any(non_event) else float("nan")
        leads = first_alert_leads(work["tau"].values, alert, onset_idx, lookback_steps)
        finite = np.asarray([x for x in leads if np.isfinite(x)], dtype=float)
        detection_rate = float(len(finite) / max(len(onset_idx), 1))
        median_lead = float(np.median(finite)) if len(finite) else float("nan")
        mean_lead = float(np.mean(finite)) if len(finite) else float("nan")
        utility = detection_rate * (median_lead if np.isfinite(median_lead) else 0.0)
        rows.append(
            {
                "analysis": label,
                "monitor": name,
                "source_column": col,
                "threshold": threshold,
                "false_alarm_rate_target": false_alarm_rate,
                "false_alarm_rate_observed": false_alarm_observed,
                "event_count": int(len(onset_idx)),
                "detection_rate": detection_rate,
                "median_lead_time": median_lead,
                "mean_lead_time": mean_lead,
                "lead_coverage_utility": utility,
                "spearman_monitor_vs_oracle": spearman(vals, work["refinement_oracle"].values),
            }
        )
        for event_id, idx in enumerate(onset_idx):
            lead = float(leads[event_id]) if event_id < len(leads) else float("nan")
            lead_rows.append(
                {
                    "analysis": label,
                    "monitor": name,
                    "event_id": event_id,
                    "run_id": int(work.iloc[idx]["run_id"]),
                    "event_time": float(work.iloc[idx]["t"]),
                    "event_tau": float(work.iloc[idx]["tau"]),
                    "lead_time": lead,
                    "detected": int(np.isfinite(lead)),
                }
            )
        work[f"alert_{name}"] = alert.astype(int)

    result = pd.DataFrame(rows).sort_values(["lead_coverage_utility", "detection_rate"], ascending=False)
    lead_df = pd.DataFrame(lead_rows)
    return result, lead_df


def utility_from_leads(leads: np.ndarray) -> float:
    leads = np.asarray(leads, dtype=float)
    finite = leads[np.isfinite(leads)]
    if len(leads) == 0 or len(finite) == 0:
        return 0.0
    return float(len(finite) / len(leads) * np.median(finite))


def bootstrap_utility_deltas(
    lead_df: pd.DataFrame,
    result: pd.DataFrame,
    reference_monitor: str = "W_MAAT",
    samples: int = BOOTSTRAP_SAMPLES,
    seed: int = SEED + 17,
) -> pd.DataFrame:
    primary = lead_df[lead_df["analysis"] == "primary"].copy()
    pivot = primary.pivot(index="event_id", columns="monitor", values="lead_time")
    rng = np.random.default_rng(seed)
    event_ids = pivot.index.to_numpy()
    rows: list[dict[str, object]] = []
    observed = result.set_index("monitor")["lead_coverage_utility"].to_dict()
    if reference_monitor not in pivot.columns:
        return pd.DataFrame(rows)
    for monitor in pivot.columns:
        if monitor == reference_monitor:
            continue
        deltas = np.empty(samples, dtype=float)
        for i in range(samples):
            sample_ids = rng.choice(event_ids, size=len(event_ids), replace=True)
            ref_u = utility_from_leads(pivot.loc[sample_ids, reference_monitor].to_numpy())
            mon_u = utility_from_leads(pivot.loc[sample_ids, monitor].to_numpy())
            deltas[i] = ref_u - mon_u
        rows.append(
            {
                "reference_monitor": reference_monitor,
                "comparison_monitor": monitor,
                "observed_delta_utility": float(observed[reference_monitor] - observed[monitor]),
                "bootstrap_mean_delta": float(np.mean(deltas)),
                "ci95_low": float(np.quantile(deltas, 0.025)),
                "ci95_high": float(np.quantile(deltas, 0.975)),
                "p_delta_gt_0": float(np.mean(deltas > 0.0)),
                "bootstrap_samples": int(samples),
                "bootstrap_unit": "refinement event onsets",
            }
        )
    return pd.DataFrame(rows).sort_values("observed_delta_utility", ascending=False)


def evaluate_triggers(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    work = add_warning_ablation_columns(df)
    dt_sample = float(np.median(np.diff(np.sort(work["t"].unique()))))
    horizon_steps = int(round(1.0 / dt_sample))
    lookback_steps = int(round(2.0 / dt_sample))
    false_alarm_rate = 0.05
    spinup_time = 4.0
    eval_mask = work["t"].values >= spinup_time
    event_threshold = float(np.quantile(work.loc[eval_mask, "refinement_oracle"].values, 0.90))
    event_mask = work["refinement_oracle"].values >= event_threshold
    event_mask &= eval_mask
    onset_idx = event_onsets(event_mask)
    non_event = (~event_mask) & eval_mask

    monitors = {
        "W_MAAT": "W_MAAT",
        "max_abs_vorticity": "max_abs_vorticity",
        "rms_vorticity": "rms_vorticity",
        "palinstrophy": "palinstrophy",
        "high_mode_fraction": "high_mode_fraction",
    }
    result, primary_leads = evaluate_monitor_rows(
        work, monitors, event_mask, onset_idx, non_event, false_alarm_rate, lookback_steps, "primary"
    )
    ablation_monitors = {
        "W_MAAT": "W_MAAT",
        "robustness_loss_only": "W_robustness_only",
        "activity_only": "W_activity_only",
        "tail_only": "W_tail_only",
        "growth_only": "W_growth_only",
        "W_no_tail": "W_no_tail",
        "W_no_growth": "W_no_growth",
        "W_no_robustness": "W_no_robustness",
    }
    ablation_result, ablation_leads = evaluate_monitor_rows(
        work, ablation_monitors, event_mask, onset_idx, non_event, false_alarm_rate, lookback_steps, "ablation"
    )
    bootstrap = bootstrap_utility_deltas(primary_leads, result)
    summary = {
        "event_threshold_quantile": 0.90,
        "event_threshold": event_threshold,
        "event_count": int(len(onset_idx)),
        "false_alarm_rate_target": false_alarm_rate,
        "spinup_time_excluded": spinup_time,
        "horizon_time": float(horizon_steps * dt_sample),
        "lookback_time": float(lookback_steps * dt_sample),
        "best_by_utility": result.iloc[0].to_dict(),
        "best_by_median_lead": result.sort_values(["median_lead_time", "detection_rate"], ascending=False).iloc[0].to_dict(),
        "monitor_results": result.to_dict(orient="records"),
        "ablation_results": ablation_result.to_dict(orient="records"),
        "bootstrap_delta_results": bootstrap.to_dict(orient="records"),
    }
    lead_df = pd.concat([primary_leads, ablation_leads], ignore_index=True)
    return result, summary, ablation_result, bootstrap, lead_df


def plot_timeseries(df: pd.DataFrame, result: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(9.0, 7.4), sharex=True)
    axes[0].plot(df["t"], normalize01(df["refinement_oracle"].values), label="refinement oracle", color="#111827")
    axes[0].plot(df["t"], normalize01(df["W_MAAT"].values), label="W_MAAT", color="#0f766e")
    axes[0].plot(df["t"], normalize01(df["max_abs_vorticity"].values), label="max |omega|", color="#dc2626", alpha=0.8)
    axes[0].legend(frameon=False, ncol=3, fontsize=8)
    axes[0].set_ylabel("normalised")
    axes[0].set_title("Forced 2D turbulence: warning signals")

    axes[1].plot(df["t"], df["R_rob"], label="R_rob", color="#2563eb")
    axes[1].plot(df["t"], df["V"], label="V low-mode coherence", color="#7c3aed")
    axes[1].legend(frameon=False, ncol=2, fontsize=8)
    axes[1].set_ylabel("support")

    axes[2].plot(df["t"], df["high_mode_fraction"], label="high-k fraction", color="#f97316")
    axes[2].plot(df["t"], df["tail_growth_pressure"], label="tail growth pressure", color="#b45309", alpha=0.8)
    axes[2].legend(frameon=False, ncol=2, fontsize=8)
    axes[2].set_ylabel("tail")
    axes[2].set_xlabel("time")
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_trigger_comparison(result: pd.DataFrame, out: Path) -> None:
    order = ["W_MAAT", "max_abs_vorticity", "rms_vorticity", "palinstrophy", "high_mode_fraction"]
    work = result.set_index("monitor").loc[order].reset_index()
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    colors = ["#0f766e" if m == "W_MAAT" else "#64748b" for m in work["monitor"]]
    ax.bar(work["monitor"], work["lead_coverage_utility"], color=colors)
    ax.set_ylabel("detection rate x median lead time")
    ax.set_title("Adaptive-refinement utility at matched false alarm")
    ax.grid(axis="y", alpha=0.25)
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_false_alarm(result: pd.DataFrame, out: Path) -> None:
    order = ["W_MAAT", "max_abs_vorticity", "rms_vorticity", "palinstrophy", "high_mode_fraction"]
    work = result.set_index("monitor").loc[order].reset_index()
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    ax.bar(work["monitor"], work["false_alarm_rate_observed"], color="#94a3b8")
    ax.axhline(float(work["false_alarm_rate_target"].iloc[0]), color="#dc2626", linestyle="--", label="target")
    ax.set_ylabel("false alarm rate on non-event windows")
    ax.set_title("False alarm calibration")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_ablation_comparison(ablation: pd.DataFrame, out: Path) -> None:
    order = [
        "W_MAAT",
        "W_no_tail",
        "W_no_growth",
        "W_no_robustness",
        "tail_only",
        "growth_only",
        "robustness_loss_only",
        "activity_only",
    ]
    work = ablation.set_index("monitor").loc[order].reset_index()
    fig, ax = plt.subplots(figsize=(9.4, 5.2))
    colors = ["#0f766e" if m == "W_MAAT" else "#64748b" for m in work["monitor"]]
    ax.bar(work["monitor"], work["lead_coverage_utility"], color=colors)
    ax.set_ylabel("detection rate x median lead time")
    ax.set_title("W_MAAT component ablation at matched false alarm")
    ax.grid(axis="y", alpha=0.25)
    ax.tick_params(axis="x", rotation=28)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_bootstrap_ci(bootstrap: pd.DataFrame, out: Path) -> None:
    if bootstrap.empty:
        return
    order = [
        "high_mode_fraction",
        "palinstrophy",
        "max_abs_vorticity",
        "rms_vorticity",
    ]
    work = bootstrap.set_index("comparison_monitor").loc[order].reset_index()
    y = np.arange(len(work))
    x = work["observed_delta_utility"].values
    xerr = np.vstack(
        [
            x - work["ci95_low"].values,
            work["ci95_high"].values - x,
        ]
    )
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    ax.errorbar(x, y, xerr=xerr, fmt="o", color="#0f766e", ecolor="#94a3b8", capsize=4)
    ax.axvline(0.0, color="#dc2626", linestyle="--", linewidth=1.2)
    ax.set_yticks(y)
    ax.set_yticklabels(work["comparison_monitor"])
    ax.set_xlabel("utility difference: W_MAAT - monitor")
    ax.set_title("Event-bootstrap utility differences")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def write_jhtdb_protocol(out: Path) -> None:
    protocol = {
        "status": "predeclared external replication route, not executed in Paper 67",
        "source": "Johns Hopkins Turbulence Database (JHTDB)",
        "url": "https://turbulence.idies.jhu.edu/",
        "candidate_datasets": [
            "forced isotropic turbulence / isotropic1024",
            "JHTDB v2.0 Zarr-backed isotropic turbulence datasets where available",
        ],
        "required_fields": [
            "velocity or vorticity sub-cubes",
            "spatial gradients or enough local stencil data to compute them",
            "time-resolved samples for lead-time evaluation",
        ],
        "predeclared_test": [
            "compute W_MAAT and vorticity/palinstrophy monitors on equal sub-cubes",
            "define refinement events by future high-k stress or gradient-growth oracle",
            "calibrate thresholds at the same false-alarm rate on non-event windows",
            "compare median lead time and detection rate using bootstrap confidence intervals",
        ],
        "license_note": "Do not redistribute raw JHTDB cutouts unless permitted by the current JHTDB terms. Cite JHTDB and the relevant dataset publications.",
    }
    out.write_text(json.dumps(protocol, indent=2), encoding="utf-8")


def write_readme(outdir: Path, summary: dict[str, object]) -> None:
    best = summary["best_by_utility"]
    highk_delta = next(
        (
            row
            for row in summary.get("bootstrap_delta_results", [])
            if row.get("comparison_monitor") == "high_mode_fraction"
        ),
        None,
    )
    lines = [
        "# Paper 67 -- Structural Early Warning as a Refinement Trigger",
        "",
        "Utility benchmark for adaptive-refinement triggers in forced 2D turbulence.",
        "",
        "## Headline",
        "",
        f"- Best lead-coverage utility monitor: `{best['monitor']}`.",
        f"- Utility: `{best['lead_coverage_utility']:.4f}` at false-alarm target `{summary['false_alarm_rate_target']:.3f}`.",
        f"- Event count: `{summary['event_count']}`.",
    ]
    if highk_delta:
        lines.extend(
            [
                f"- Observed utility margin over high-k: `{highk_delta['observed_delta_utility']:.4f}`.",
                f"- Event-bootstrap 95% CI for that margin: `[{highk_delta['ci95_low']:.4f}, {highk_delta['ci95_high']:.4f}]`.",
            ]
        )
    lines.extend(
        [
            "",
            "## Important status",
            "",
            "This first implementation uses a local forced-2D turbulence simulation.",
            "JHTDB is specified as an external replication route but no JHTDB data are downloaded or redistributed.",
            "",
            "## Outputs",
            "",
            "- `paper67_forced2d_timeseries.csv`",
            "- `paper67_trigger_results.csv`",
            "- `paper67_ablation_results.csv`",
            "- `paper67_bootstrap_ci.csv`",
            "- `paper67_event_leads.csv`",
            "- `paper67_summary.json`",
            "- `paper67_jhtdb_protocol.json`",
            "- `fig1_forced2d_warning_timeseries.png`",
            "- `fig2_trigger_leadtime_comparison.png`",
            "- `fig3_false_alarm_calibration.png`",
            "- `fig4_ablation_utility.png`",
            "- `fig5_bootstrap_delta_ci.png`",
        ]
    )
    outdir.joinpath("README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    dfs = [simulate_forced_2d(seed=SEED + 101 * i, run_id=i) for i in range(6)]
    df = pd.concat(dfs, ignore_index=True)
    result, summary, ablation, bootstrap, leads = evaluate_triggers(df)
    df.to_csv(OUTDIR / "paper67_forced2d_timeseries.csv", index=False)
    result.to_csv(OUTDIR / "paper67_trigger_results.csv", index=False)
    ablation.to_csv(OUTDIR / "paper67_ablation_results.csv", index=False)
    bootstrap.to_csv(OUTDIR / "paper67_bootstrap_ci.csv", index=False)
    leads.to_csv(OUTDIR / "paper67_event_leads.csv", index=False)
    json_summary = to_jsonable(summary)
    (OUTDIR / "paper67_summary.json").write_text(
        json.dumps(json_summary, indent=2, allow_nan=False), encoding="utf-8"
    )
    write_jhtdb_protocol(OUTDIR / "paper67_jhtdb_protocol.json")
    write_readme(OUTDIR, summary)
    plot_timeseries(df, result, OUTDIR / "fig1_forced2d_warning_timeseries.png")
    plot_trigger_comparison(result, OUTDIR / "fig2_trigger_leadtime_comparison.png")
    plot_false_alarm(result, OUTDIR / "fig3_false_alarm_calibration.png")
    plot_ablation_comparison(ablation, OUTDIR / "fig4_ablation_utility.png")
    plot_bootstrap_ci(bootstrap, OUTDIR / "fig5_bootstrap_delta_ci.png")
    print(json.dumps(json_summary, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
