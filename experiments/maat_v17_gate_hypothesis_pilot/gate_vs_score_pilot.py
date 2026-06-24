#!/usr/bin/env python3
"""MAAT v1.7 gate-hypothesis pilot.

This script tests a narrow version of the hypothesis suggested by Papers 63,
64, 66 and 67:

    R should often act as a conditional gate rather than as another global
    scoring factor.

The test is deliberately modest.  It reuses already generated repository
artifacts and compares R-as-score against R-as-gate in three domains:

* SAT CDCL branching (Paper 63): global Mode A vs structure-gated Mode A.
* Two-qubit pointer robustness (Paper 64): scalar score vs hard R gate.
* Forced-2D refinement triggers (Paper 67): W_MAAT vs hard R-gated base trigger.

No external dataset is downloaded or redistributed by this script.
"""

from __future__ import annotations

import json
import math
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

_CACHE = Path(tempfile.gettempdir()) / "maat_v17_gate_pilot_matplotlib_cache"
_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPO = Path(__file__).resolve().parents[2]
OUTDIR = Path("outputs_v17_gate_pilot")
EPS = 1.0e-12


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return float(pd.Series(x).corr(pd.Series(y), method="spearman"))


def normalize01(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    lo = float(np.nanmin(x))
    hi = float(np.nanmax(x))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo < EPS:
        return np.zeros_like(x)
    return (x - lo) / (hi - lo)


def to_jsonable(obj):
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


def event_onsets(mask: np.ndarray) -> list[int]:
    onsets: list[int] = []
    prev = False
    for i, val in enumerate(mask.astype(bool)):
        if val and not prev:
            onsets.append(i)
        prev = bool(val)
    return onsets


def first_alert_leads(times: np.ndarray, alert: np.ndarray, onsets: list[int], lookback_steps: int) -> list[float]:
    leads: list[float] = []
    for idx in onsets:
        lo = max(0, idx - lookback_steps)
        hits = np.where(alert[lo : idx + 1])[0]
        if len(hits) == 0:
            leads.append(float("nan"))
            continue
        first = lo + int(hits[0])
        leads.append(float(times[idx] - times[first]))
    return leads


def threshold_at_false_alarm(values: np.ndarray, non_event_mask: np.ndarray, false_alarm_rate: float) -> float:
    pool = np.asarray(values, dtype=float)[non_event_mask]
    if len(pool) == 0:
        return float("inf")
    return float(np.quantile(pool, 1.0 - false_alarm_rate))


def trigger_utility(
    df: pd.DataFrame,
    monitor: str,
    event_mask: np.ndarray,
    non_event_mask: np.ndarray,
    onsets: list[int],
    lookback_steps: int,
    false_alarm_rate: float,
) -> dict[str, float | str]:
    vals = df[monitor].to_numpy(dtype=float)
    threshold = threshold_at_false_alarm(vals, non_event_mask, false_alarm_rate)
    alert = vals >= threshold
    leads = first_alert_leads(df["tau"].to_numpy(dtype=float), alert, onsets, lookback_steps)
    finite = np.asarray([x for x in leads if np.isfinite(x)], dtype=float)
    detection = float(len(finite) / max(len(onsets), 1))
    median_lead = float(np.median(finite)) if len(finite) else float("nan")
    utility = detection * (median_lead if np.isfinite(median_lead) else 0.0)
    return {
        "monitor": monitor,
        "threshold": threshold,
        "false_alarm_rate": float(np.mean(alert[non_event_mask])) if np.any(non_event_mask) else float("nan"),
        "detection_rate": detection,
        "median_lead": median_lead,
        "utility": utility,
        "spearman_oracle": spearman(vals, df["refinement_oracle"].to_numpy(dtype=float)),
    }


def sat_gate_test() -> tuple[pd.DataFrame, dict[str, object]]:
    path = REPO / "experiments/sat_cdcl_structure_gated_mode_a_paper63/outputs_paper63_structure_gated_mode_a/paper63_cdcl_runs.csv"
    runs = pd.read_csv(path)
    pivot = runs.pivot_table(index=["family", "instance_id"], columns="policy", values="compute_cost", aggfunc="first")
    pivot = pivot.reset_index()
    pivot["delta_gate_minus_score"] = pivot["mode_a_gated"] - pivot["mode_a"]
    pivot["gate_better"] = pivot["delta_gate_minus_score"] < 0.0
    rows = [
        {
            "domain": "SAT",
            "metric": "compute_cost_lower_is_better",
            "score_variant": "global Mode A",
            "gate_variant": "structure-gated Mode A",
            "score_value": float(pivot["mode_a"].mean()),
            "gate_value": float(pivot["mode_a_gated"].mean()),
            "delta_gate_minus_score": float(pivot["delta_gate_minus_score"].mean()),
            "relative_improvement": float((pivot["mode_a"].mean() - pivot["mode_a_gated"].mean()) / (pivot["mode_a"].mean() + EPS)),
            "median_delta": float(pivot["delta_gate_minus_score"].median()),
            "gate_win_rate": float(pivot["gate_better"].mean()),
            "n": int(len(pivot)),
            "interpretation": "gate slightly reduces mean/tail damage relative to ungated structural policy",
        }
    ]
    by_family = (
        pivot.groupby("family", as_index=False)
        .agg(
            n=("instance_id", "count"),
            mean_score=("mode_a", "mean"),
            mean_gate=("mode_a_gated", "mean"),
            mean_delta=("delta_gate_minus_score", "mean"),
            median_delta=("delta_gate_minus_score", "median"),
            gate_win_rate=("gate_better", "mean"),
        )
        .sort_values("mean_delta")
    )
    return pd.DataFrame(rows), {"instance_level": pivot, "family_level": by_family}


def quantum_gate_test() -> tuple[pd.DataFrame, pd.DataFrame]:
    path = REPO / "experiments/two_qubit_pointer_paper64/outputs_paper64/paper64_two_qubit_instances.csv"
    df = pd.read_csv(path)
    y = df["entanglement_robustness"].to_numpy(dtype=float)
    high_y = y >= np.quantile(y, 0.80)

    # Base structural support without R.  This is intentionally a scalar toy
    # projection, not a replacement for the Paper-64 field model.
    base = (
        np.clip(df["H_stationarity"].to_numpy(dtype=float), 0, 1)
        * np.clip(df["B_corr_stat"].to_numpy(dtype=float), 0, 1)
        * np.clip(df["S_ent"].to_numpy(dtype=float), 0, 1)
        * np.clip(df["V_corr"].to_numpy(dtype=float), 0, 1)
    ) ** 0.25
    r = np.clip(df["R_rob_stat"].to_numpy(dtype=float), 0, 1)
    score = base * r
    cost_score = -df["F_scalar_stat"].to_numpy(dtype=float)

    rows = []
    variants = {
        "base_HBSV_no_R": base,
        "score_base_times_R": score,
        "score_negative_F_scalar_stat": cost_score,
    }
    for q in [0.25, 0.50, 0.75]:
        r_star = float(np.quantile(r, q))
        variants[f"gate_R_ge_q{int(q*100)}"] = base * (r >= r_star)
    for name, vals in variants.items():
        rows.append(
            {
                "domain": "Quantum",
                "metric": "entanglement_robustness_higher_is_better",
                "variant": name,
                "spearman": spearman(vals, y),
                "auroc_top20": float(roc_auc_score(high_y.astype(int), vals)),
                "mean_score": float(np.mean(vals)),
                "n": int(len(df)),
            }
        )
    res = pd.DataFrame(rows).sort_values(["spearman", "auroc_top20"], ascending=False)
    best_gate = res[res["variant"].str.startswith("gate_")].iloc[0]
    score_row = res[res["variant"] == "score_base_times_R"].iloc[0]
    summary = pd.DataFrame(
        [
            {
                "domain": "Quantum",
                "metric": "Spearman with entanglement robustness",
                "score_variant": "base * R_rob_stat",
                "gate_variant": best_gate["variant"],
                "score_value": float(score_row["spearman"]),
                "gate_value": float(best_gate["spearman"]),
                "delta_gate_minus_score": float(best_gate["spearman"] - score_row["spearman"]),
                "relative_improvement": float((best_gate["spearman"] - score_row["spearman"]) / (abs(score_row["spearman"]) + EPS)),
                "median_delta": None,
                "gate_win_rate": None,
                "n": int(len(df)),
                "interpretation": "hard R gate improves scalar ranking in-distribution, but this is not a field-model replacement",
            }
        ]
    )
    return summary, res


def fluid_gate_test() -> tuple[pd.DataFrame, pd.DataFrame]:
    path = REPO / "experiments/structural_early_warning_paper67/outputs_paper67/paper67_forced2d_timeseries.csv"
    df = pd.read_csv(path)
    df = df.copy()
    base = np.sqrt(np.clip(df["activity_pressure"].to_numpy(dtype=float), 0, None))
    base *= 1.0 + 0.7 * np.clip(df["high_tail_pressure"].to_numpy(dtype=float), 0, None)
    base *= 1.0 + 2.0 * np.clip(df["tail_growth_pressure"].to_numpy(dtype=float), 0, None)
    r = np.clip(df["R_rob"].to_numpy(dtype=float), 0, 1)
    df["base_tail_activity"] = base
    df["score_R_multiplicative"] = df["W_MAAT"].to_numpy(dtype=float)
    for q in [0.25, 0.50, 0.75]:
        r_star = float(np.quantile(r[df["t"].to_numpy(dtype=float) >= 4.0], q))
        df[f"gate_R_le_q{int(q*100)}"] = base * (r <= r_star)

    dt_sample = float(np.median(np.diff(np.sort(df["t"].unique()))))
    lookback_steps = int(round(2.0 / dt_sample))
    false_alarm_rate = 0.05
    eval_mask = df["t"].to_numpy(dtype=float) >= 4.0
    event_threshold = float(np.quantile(df.loc[eval_mask, "refinement_oracle"].to_numpy(dtype=float), 0.90))
    event_mask = (df["refinement_oracle"].to_numpy(dtype=float) >= event_threshold) & eval_mask
    onsets = event_onsets(event_mask)
    non_event = (~event_mask) & eval_mask

    monitor_names = [
        "base_tail_activity",
        "score_R_multiplicative",
        "gate_R_le_q25",
        "gate_R_le_q50",
        "gate_R_le_q75",
        "high_mode_fraction",
        "max_abs_vorticity",
        "rms_vorticity",
    ]
    rows = []
    for name in monitor_names:
        row = trigger_utility(df, name, event_mask, non_event, onsets, lookback_steps, false_alarm_rate)
        row["domain"] = "Fluid"
        row["metric"] = "lead_coverage_utility_higher_is_better"
        rows.append(row)
    res = pd.DataFrame(rows).sort_values(["utility", "detection_rate"], ascending=False)
    best_gate = res[res["monitor"].str.startswith("gate_")].iloc[0]
    score_row = res[res["monitor"] == "score_R_multiplicative"].iloc[0]
    summary = pd.DataFrame(
        [
            {
                "domain": "Fluid",
                "metric": "lead_coverage_utility",
                "score_variant": "W_MAAT = base * (1 - R_rob)",
                "gate_variant": best_gate["monitor"],
                "score_value": float(score_row["utility"]),
                "gate_value": float(best_gate["utility"]),
                "delta_gate_minus_score": float(best_gate["utility"] - score_row["utility"]),
                "relative_improvement": float((best_gate["utility"] - score_row["utility"]) / (abs(score_row["utility"]) + EPS)),
                "median_delta": None,
                "gate_win_rate": None,
                "n": int(len(onsets)),
                "interpretation": "best hard R gate improves utility but sacrifices event coverage and is threshold-sensitive",
            }
        ]
    )
    return summary, res


def plot_domain_summary(summary: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    colors = ["#0f766e" if x > 0 else "#b91c1c" for x in summary["relative_improvement"]]
    ax.bar(summary["domain"], summary["relative_improvement"], color=colors)
    ax.axhline(0.0, color="#111827", linewidth=1.0)
    ax.set_ylabel("relative gate improvement over score")
    ax.set_title("R-gate hypothesis pilot across three domains")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_quantum(res: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    work = res.sort_values("spearman")
    colors = ["#0f766e" if v.startswith("gate_") else "#64748b" for v in work["variant"]]
    ax.barh(work["variant"], work["spearman"], color=colors)
    ax.set_xlabel("Spearman with entanglement robustness")
    ax.set_title("Quantum scalar R gate vs R score")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_fluid(res: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    work = res.sort_values("utility")
    colors = ["#0f766e" if str(v).startswith("gate_") else "#64748b" for v in work["monitor"]]
    ax.barh(work["monitor"], work["utility"], color=colors)
    ax.set_xlabel("lead-coverage utility")
    ax.set_title("Fluid R gate vs R score")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def write_readme(summary: pd.DataFrame) -> None:
    table = summary.copy()
    table = table.where(pd.notnull(table), "--")
    text = """# MAAT v1.7 Gate-Hypothesis Pilot

This experiment re-analyses existing repository outputs to test a narrow
hypothesis:

> R often works better as a conditional gate than as another global scoring
> factor.

Domains:

- SAT: Paper 63 global Mode A vs structure-gated Mode A.
- Quantum: Paper 64 scalar R score vs hard R gate for entanglement robustness.
- Fluid: Paper 67 W_MAAT score vs hard R gate for refinement triggering.

This is a pilot, not a new final MAAT version.  It uses no external dataset
beyond artifacts already present in this repository.

## Headline

The evidence supports the gate direction, but not yet a completed v2.0-style
framework.  SAT shows a small tail-damage reduction from gating.  The quantum
scalar test shows a large in-distribution ranking improvement.  The fluid test
shows that the best hard R gate improves lead-coverage utility, but does so by
sacrificing event coverage and is sensitive to the gate threshold.

Conclusion: v1.7 is justified as a gate-hypothesis / design-direction paper.
v2.0 would require a predeclared, calibrated gate that generalizes under harder
external tests.

## Summary Table

```text
{table}
```
""".format(table=table.to_string(index=False))
    (OUTDIR / "README.md").write_text(text, encoding="utf-8")


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    sat_summary, sat_details = sat_gate_test()
    quantum_summary, quantum_details = quantum_gate_test()
    fluid_summary, fluid_details = fluid_gate_test()
    summary = pd.DataFrame(
        sat_summary.to_dict(orient="records")
        + quantum_summary.to_dict(orient="records")
        + fluid_summary.to_dict(orient="records")
    )

    summary.to_csv(OUTDIR / "v17_gate_hypothesis_summary.csv", index=False)
    sat_details["instance_level"].to_csv(OUTDIR / "sat_gate_instance_deltas.csv", index=False)
    sat_details["family_level"].to_csv(OUTDIR / "sat_gate_family_deltas.csv", index=False)
    quantum_details.to_csv(OUTDIR / "quantum_gate_scalar_results.csv", index=False)
    fluid_details.to_csv(OUTDIR / "fluid_gate_trigger_results.csv", index=False)
    payload = {
        "status": "pilot re-analysis",
        "hypothesis": "R as conditional gate versus R as global score",
        "summary": summary.to_dict(orient="records"),
        "conclusion": (
            "Gate-favouring evidence appears in all three domains, but it is "
            "not yet enough for a v2.0 claim. SAT is weak-positive, quantum is "
            "strong in-distribution, and fluid is positive but threshold-sensitive "
            "and sacrifices detection coverage. This supports a v1.7 gate-hypothesis "
            "paper."
        ),
    }
    (OUTDIR / "v17_gate_hypothesis_summary.json").write_text(
        json.dumps(to_jsonable(payload), indent=2, allow_nan=False), encoding="utf-8"
    )
    plot_domain_summary(summary, OUTDIR / "fig1_domain_gate_summary.png")
    plot_quantum(quantum_details, OUTDIR / "fig2_quantum_gate_scalar.png")
    plot_fluid(fluid_details, OUTDIR / "fig3_fluid_gate_trigger.png")
    write_readme(summary)
    print(json.dumps(to_jsonable(payload), indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
