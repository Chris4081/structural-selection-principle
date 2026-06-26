#!/usr/bin/env python3
"""Post-hoc diagnostics for Paper 69 SATLIB smoke outputs.

This script does not modify gate parameters and should not be read as a tuning
step. It analyses why the frozen v1.7 gate can lose to MOMS/Jeroslow-Wang on
external SATLIB subsets by comparing instance-level compute cost, family-level
regret, short-clause pressure proxies, and gate activation levels.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "outputs_paper69_sat_gate_challenge"
ANALYSIS_JSON = OUTDIR / "paper69_moms_vs_gate_diagnostics.json"
ANALYSIS_CSV = OUTDIR / "paper69_moms_vs_gate_instance_diagnostics.csv"
FAMILY_CSV = OUTDIR / "paper69_moms_vs_gate_family_diagnostics.csv"


def load_outputs() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    records = pd.read_csv(OUTDIR / "paper69_solve_records.csv")
    features = pd.read_csv(OUTDIR / "paper69_gate_features.csv")
    summary = json.loads((OUTDIR / "paper69_summary.json").read_text(encoding="utf-8"))
    return records, features, summary


def pivot_costs(records: pd.DataFrame) -> pd.DataFrame:
    pivot = records.pivot_table(
        index=["dataset_id", "family", "split"],
        columns="policy",
        values="compute_cost",
        aggfunc="first",
    ).reset_index()
    return pivot


def build_instance_diagnostics(records: pd.DataFrame, features: pd.DataFrame) -> pd.DataFrame:
    pivot = pivot_costs(records)
    keep_features = [
        "dataset_id",
        "alpha",
        "H",
        "B",
        "S",
        "V",
        "R_rob",
        "G_gate",
        "degree_cv",
        "literal_imbalance",
        "clustering_mean",
        "largest_component_frac",
    ]
    df = pivot.merge(features[keep_features], on="dataset_id", how="left")
    for baseline in ["score_with_R", "moms", "jeroslow_wang", "vsids", "progress_only", "gate_shuffled_R"]:
        if baseline in df.columns and "gate_v17" in df.columns:
            df[f"delta_gate_vs_{baseline}"] = df[baseline] - df["gate_v17"]
            df[f"regret_gate_to_{baseline}"] = df["gate_v17"] - df[baseline]
    if "moms" in df.columns and "score_with_R" in df.columns:
        df["delta_moms_vs_score_with_R"] = df["score_with_R"] - df["moms"]
    return df


def family_diagnostics(inst: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "gate_v17",
        "score_with_R",
        "moms",
        "jeroslow_wang",
        "vsids",
        "G_gate",
        "R_rob",
        "S",
        "V",
        "alpha",
        "delta_gate_vs_score_with_R",
        "delta_gate_vs_moms",
        "regret_gate_to_moms",
    ]
    existing = [m for m in metrics if m in inst.columns]
    fam = inst.groupby("family")[existing].agg(["mean", "median", "count"])
    fam.columns = ["_".join(col).strip("_") for col in fam.columns]
    fam = fam.reset_index()
    return fam


def correlations(inst: pd.DataFrame) -> dict[str, dict[str, float]]:
    targets = ["regret_gate_to_moms", "delta_gate_vs_score_with_R", "gate_v17", "moms"]
    predictors = ["G_gate", "R_rob", "S", "V", "H", "B", "alpha", "degree_cv", "literal_imbalance", "clustering_mean"]
    out: dict[str, dict[str, float]] = {}
    for target in targets:
        if target not in inst.columns:
            continue
        vals: dict[str, float] = {}
        for pred in predictors:
            if pred not in inst.columns:
                continue
            sub = inst[[target, pred]].replace([np.inf, -np.inf], np.nan).dropna()
            if len(sub) < 4 or sub[target].nunique() < 2 or sub[pred].nunique() < 2:
                vals[pred] = float("nan")
            else:
                vals[pred] = float(sub[target].corr(sub[pred], method="spearman"))
        out[target] = vals
    return out


def plot_family_regret(fam: pd.DataFrame) -> None:
    if "regret_gate_to_moms_mean" not in fam.columns:
        return
    data = fam.sort_values("regret_gate_to_moms_mean", ascending=False)
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    ax.axhline(0, color="#444", lw=1)
    ax.bar(data["family"], data["regret_gate_to_moms_mean"], color="#b23a48", edgecolor="#1f2933")
    ax.set_ylabel("mean regret: gate cost - MOMS cost")
    ax.set_title("Why MOMS wins: family-level gate regret")
    ax.tick_params(axis="x", rotation=25)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig3_moms_regret_by_family.png", dpi=180)
    plt.close(fig)


def plot_gate_vs_moms(inst: pd.DataFrame) -> None:
    if not {"G_gate", "regret_gate_to_moms", "family"}.issubset(inst.columns):
        return
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    families = sorted(inst["family"].unique())
    cmap = plt.get_cmap("tab10")
    for i, family in enumerate(families):
        sub = inst[inst["family"] == family]
        ax.scatter(sub["G_gate"], sub["regret_gate_to_moms"], s=54, color=cmap(i % 10), label=family, alpha=0.82)
    ax.axhline(0, color="#444", lw=1)
    ax.set_xlabel("frozen v1.7 gate activation G")
    ax.set_ylabel("regret to MOMS (gate cost - MOMS cost)")
    ax.set_title("Gate activation does not capture MOMS short-clause advantage")
    ax.legend(fontsize=7, frameon=False, loc="best")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig4_gate_activation_vs_moms_regret.png", dpi=180)
    plt.close(fig)


def interpretation(summary: dict, inst: pd.DataFrame, fam: pd.DataFrame) -> list[str]:
    lines = [
        "This is a post-hoc diagnostic, not a retuning step.",
        "MOMS is a local short-clause heuristic; the frozen gate is a root-level robustness activation.",
        "A MOMS advantage therefore indicates that short-clause pressure is carrying decision-local information not captured by the current root gate.",
    ]
    if "delta_gate_vs_score_with_R" in inst.columns:
        mean_delta = float(inst["delta_gate_vs_score_with_R"].mean())
        lines.append(
            f"Gate versus score-with-R mean delta in this run is {mean_delta:.4g}; this is not support unless the preregistered CI lower bound is positive."
        )
    if "regret_gate_to_moms_mean" in fam.columns:
        worst = fam.sort_values("regret_gate_to_moms_mean", ascending=False).head(3)
        lines.append(
            "Largest mean gate-to-MOMS regret families: "
            + ", ".join(f"{r.family} ({r.regret_gate_to_moms_mean:.3g})" for r in worst.itertuples())
            + "."
        )
    lines.append(
        "The appropriate next scientific question is why MOMS's short-clause signal dominates, not how to tune the gate until it wins."
    )
    return lines


def main() -> None:
    records, features, summary = load_outputs()
    inst = build_instance_diagnostics(records, features)
    fam = family_diagnostics(inst)
    corr = correlations(inst)
    inst.to_csv(ANALYSIS_CSV, index=False)
    fam.to_csv(FAMILY_CSV, index=False)
    plot_family_regret(fam)
    plot_gate_vs_moms(inst)

    analysis = {
        "paper": 69,
        "analysis_type": "post_hoc_moms_vs_gate_diagnostic_not_retuning",
        "source_summary_status": summary.get("status"),
        "source_execution_scope": summary.get("execution_scope"),
        "n_instances": int(len(inst)),
        "families": sorted(inst["family"].unique().tolist()),
        "core_question": "Which structural information is used by MOMS that the frozen v1.7 gate does not use?",
        "correlations_spearman": corr,
        "interpretation": interpretation(summary, inst, fam),
        "outputs": {
            "instance_diagnostics_csv": str(ANALYSIS_CSV.relative_to(ROOT)),
            "family_diagnostics_csv": str(FAMILY_CSV.relative_to(ROOT)),
            "family_regret_figure": "outputs_paper69_sat_gate_challenge/fig3_moms_regret_by_family.png",
            "gate_activation_figure": "outputs_paper69_sat_gate_challenge/fig4_gate_activation_vs_moms_regret.png",
        },
    }
    ANALYSIS_JSON.write_text(json.dumps(analysis, indent=2), encoding="utf-8")
    print(json.dumps(analysis, indent=2))


if __name__ == "__main__":
    main()

