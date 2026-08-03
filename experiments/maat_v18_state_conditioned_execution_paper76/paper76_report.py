#!/usr/bin/env python3
"""Generate frozen Paper-76 tables, figures, and secondary cost decomposition."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from paper76_execution import PRIMARY, family_bootstrap


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs_paper76_state_conditioned_execution"
_CACHE = Path(os.environ.get("TMPDIR", "/tmp")) / "maat_paper76_report"
(_CACHE / "matplotlib").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE / "matplotlib"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


COLORS = {
    "moms": "#1f4e5f",
    "jeroslow_wang": "#49796b",
    "vsids": "#8b6f47",
    "score_with_R_Paper69": "#a65f46",
    "static_gate_v17_Paper69": "#d08c60",
    PRIMARY: "#c33c32",
    "state_v18_gate_signal_shuffled": "#64748b",
    "random_activation_matched_q": "#94a3b8",
}


def main() -> None:
    records = pd.read_csv(OUT / "paper76_solve_records.csv")
    comparisons = pd.read_csv(OUT / "paper76_primary_comparisons.csv")
    test = records[records.split == "test"].copy()

    selected = list(COLORS)
    means = test[test.policy.isin(selected)].groupby("policy").total_cost.mean().reindex(selected)
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    ax.bar(range(len(means)), means.values, color=[COLORS[key] for key in means.index])
    ax.set_xticks(range(len(means)))
    ax.set_xticklabels([
        "MOMS", "JW", "VSIDS", "score+R", "static v1.7", "state v1.8", "shuffled", "random"
    ], rotation=24, ha="right")
    ax.set_ylabel("Mean total cost (descriptive, instance-weighted)")
    ax.set_title("Paper 76: frozen policy execution")
    ax.grid(axis="y", alpha=0.22)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_policy_total_cost.png", dpi=200)
    plt.close(fig)

    comp = comparisons.copy()
    labels = (
        comp.comparison.str.replace(f"{PRIMARY}_vs_", "", regex=False)
        .str.replace("_Paper69", "", regex=False)
        .str.replace("state_v18_", "", regex=False)
        .str.replace("_", " ", regex=False)
    )
    y = np.arange(len(comp))
    fig, ax = plt.subplots(figsize=(8.7, 5.1))
    xerr = np.vstack([comp["mean"] - comp.ci95_low, comp.ci95_high - comp["mean"]])
    ax.errorbar(comp["mean"], y, xerr=xerr, fmt="o", color="#c33c32", ecolor="#495057", capsize=4)
    ax.axvline(0, color="black", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Family-balanced ΔU (positive favors state v1.8)")
    ax.set_title("Primary paired comparisons with 95% bootstrap intervals")
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_primary_comparisons.png", dpi=200)
    plt.close(fig)

    pair = test[test.policy.isin([PRIMARY, "moms"])].pivot(
        index=["dataset_id", "family"], columns="policy",
        values=["search_cost", "total_cost", "structural_overhead_s", "activation_rate"],
    ).dropna()
    pair.columns = ["_".join(column) for column in pair.columns]
    pair["delta_search"] = pair.search_cost_moms - pair[f"search_cost_{PRIMARY}"]
    pair["delta_total"] = pair.total_cost_moms - pair[f"total_cost_{PRIMARY}"]
    pair["overhead_penalty"] = pair[f"total_cost_{PRIMARY}"] - pair[f"search_cost_{PRIMARY}"]
    pair = pair.reset_index()
    family = pair.groupby("family", as_index=False).agg(
        delta_search=("delta_search", "mean"), delta_total=("delta_total", "mean"),
        overhead_penalty=("overhead_penalty", "mean"),
        activation_rate=(f"activation_rate_{PRIMARY}", "mean"),
    )
    family.to_csv(OUT / "paper76_cost_decomposition_by_family.csv", index=False)
    search_boot = family_bootstrap(pair.delta_search, pair.family)
    total_boot = family_bootstrap(pair.delta_total, pair.family)
    diagnostics = {
        "status": "secondary_descriptive_analysis_no_reclassification",
        "search_only_delta_vs_MOMS": search_boot,
        "total_cost_delta_vs_MOMS": total_boot,
        "mean_primary_structural_overhead_seconds": float(
            test.loc[test.policy == PRIMARY, "structural_overhead_s"].mean()
        ),
        "max_incremental_audit_error": float(records.audit_max_abs_error.max()),
        "note": "Positive delta favors state v1.8. Paper-75 status labels remain unchanged.",
    }
    (OUT / "paper76_cost_decomposition.json").write_text(json.dumps(diagnostics, indent=2), encoding="utf-8")

    positions = np.arange(len(family))
    width = 0.35
    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    ax.bar(positions - width / 2, family.delta_search, width, label="search cost only", color="#437f78")
    ax.bar(positions + width / 2, family.delta_total, width, label="total cost", color="#c33c32")
    ax.axhline(0, color="black", lw=1)
    ax.set_xticks(positions)
    ax.set_xticklabels([
        name.replace("satlib_", "").replace("dimacs_", "").replace("random_3sat_", "")
        for name in family.family
    ], rotation=22, ha="right")
    ax.set_ylabel("ΔU vs MOMS (positive favors state v1.8)")
    ax.set_title("Search degradation and charged overhead by family")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(OUT / "fig3_family_cost_decomposition.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 4.6))
    ax.bar(range(len(family)), family.activation_rate, color="#2f6690")
    ax.axhline(0.25, color="#c33c32", linestyle="--", label="calibration budget q=0.25")
    ax.axhline(0.125, color="#111827", linestyle=":", label="minimum q_test,min=0.125")
    ax.set_xticks(range(len(family)))
    ax.set_xticklabels([
        name.replace("satlib_", "").replace("dimacs_", "").replace("random_3sat_", "")
        for name in family.family
    ], rotation=22, ha="right")
    ax.set_ylabel("Mean decision activation rate")
    ax.set_title("State-conditioned activation under the frozen threshold")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(OUT / "fig4_activation_by_family.png", dpi=200)
    plt.close(fig)

    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()

