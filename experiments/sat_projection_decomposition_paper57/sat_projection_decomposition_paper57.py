#!/usr/bin/env python3
"""
Paper 57: SAT defect-field projection decomposition.

This script is intentionally downstream of Paper 55. It does not generate new
SAT instances and does not rerun CDCL solvers. Instead, it takes the canonical
SAT-family CDCL table from Paper 55 and asks where the predictive signal lives:
mean, tail, cluster, or scale projections of the defect-field geometry.

Run from this directory:
    python3 sat_projection_decomposition_paper57.py

Or from the repository root:
    python3 experiments/sat_projection_decomposition_paper57/sat_projection_decomposition_paper57.py
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold, LeaveOneGroupOut, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

SEED = 57057
OUTDIR = Path("outputs_paper57")
OUTDIR.mkdir(parents=True, exist_ok=True)
(OUTDIR / ".mplconfig").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUTDIR / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def locate_paper55_csv() -> Path:
    env = os.environ.get("PAPER55_CSV")
    candidates = []
    if env:
        candidates.append(Path(env))
    here = Path(__file__).resolve().parent
    candidates.extend(
        [
            here.parent / "sat_cdcl_external_paper55" / "outputs_paper55" / "paper55_sat_cdcl_instances.csv",
            Path.cwd().parent / "sat_cdcl_external_paper55" / "outputs_paper55" / "paper55_sat_cdcl_instances.csv",
            Path.cwd()
            / "experiments"
            / "sat_cdcl_external_paper55"
            / "outputs_paper55"
            / "paper55_sat_cdcl_instances.csv",
            Path("/Volumes/MAATSSD/structural-selection-principle/experiments/sat_cdcl_external_paper55/outputs_paper55/paper55_sat_cdcl_instances.csv"),
        ]
    )
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(
        "Could not find paper55_sat_cdcl_instances.csv. Run Paper 55 first or set PAPER55_CSV=/path/to/file."
    )


FEATURE_GROUPS: dict[str, list[str]] = {
    "density_only": ["n_vars", "n_clauses", "alpha", "density_peak"],
    "macro_scalar": ["H", "B", "S", "V", "R_resp", "R_rob", "F_struct_scalar"],
    "mean_projection": ["frust_mean", "F_mean"],
    "tail_projection": [
        "frust_std",
        "frust_q75",
        "frust_q90",
        "frust_q95",
        "frust_q99",
        "frust_max",
        "frust_tail95_mean",
        "frust_tail_mass",
        "frust_gini",
        "R_tail",
        "F_tail",
    ],
    "cluster_projection": [
        "hotspot_fraction",
        "hotspot_cluster_frac",
        "hotspot_edge_density",
        "R_cluster",
        "F_cluster",
    ],
    "scale_projection": [
        "covariance_condition_log",
        "expansion_defect",
        "propagation_depth_mean",
        "contradiction_rate_mean",
        "backdoor_density",
        "F_scale",
    ],
    "standard_graph": [
        "n_vars",
        "n_clauses",
        "alpha",
        "degree_cv",
        "max_degree_norm",
        "literal_imbalance",
        "mean_var_graph_degree",
        "clustering_mean",
        "component_frac",
        "largest_component_frac",
        "expansion_defect",
        "backdoor_density",
    ],
}

FEATURE_GROUPS["mean_tail"] = FEATURE_GROUPS["mean_projection"] + FEATURE_GROUPS["tail_projection"]
FEATURE_GROUPS["tail_cluster"] = FEATURE_GROUPS["tail_projection"] + FEATURE_GROUPS["cluster_projection"]
FEATURE_GROUPS["tail_scale"] = FEATURE_GROUPS["tail_projection"] + FEATURE_GROUPS["scale_projection"]
FEATURE_GROUPS["cluster_scale"] = FEATURE_GROUPS["cluster_projection"] + FEATURE_GROUPS["scale_projection"]
FEATURE_GROUPS["all_defect_projections"] = (
    FEATURE_GROUPS["mean_projection"]
    + FEATURE_GROUPS["tail_projection"]
    + FEATURE_GROUPS["cluster_projection"]
    + FEATURE_GROUPS["scale_projection"]
    + ["R_sat", "F_struct_multi"]
)

INCREMENTAL_GROUPS = {
    "mean_only": FEATURE_GROUPS["mean_projection"],
    "mean_plus_tail": FEATURE_GROUPS["mean_projection"] + FEATURE_GROUPS["tail_projection"],
    "mean_tail_cluster": FEATURE_GROUPS["mean_projection"] + FEATURE_GROUPS["tail_projection"] + FEATURE_GROUPS["cluster_projection"],
    "mean_tail_cluster_scale": FEATURE_GROUPS["all_defect_projections"],
}


def score_predictions(y: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    rho = spearmanr(y, pred).correlation
    return {
        "r2": float(r2_score(y, pred)),
        "rmse": float(math.sqrt(mean_squared_error(y, pred))),
        "spearman": float(0.0 if np.isnan(rho) else rho),
    }


def evaluate_cv(df: pd.DataFrame, feature_sets: dict[str, list[str]], rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame]:
    y = df["cdcl_hardness"].to_numpy()
    kfold = KFold(n_splits=5, shuffle=True, random_state=SEED)
    rows = []
    pred_rows = df[["instance_id", "family", "cdcl_hardness"]].copy()
    for name, cols in feature_sets.items():
        X = df[cols].to_numpy()
        models = {
            "ridge": make_pipeline(StandardScaler(), Ridge(alpha=2.0)),
            "rf": RandomForestRegressor(n_estimators=280, min_samples_leaf=4, random_state=SEED, n_jobs=1),
        }
        for model_name, model in models.items():
            pred = cross_val_predict(model, X, y, cv=kfold)
            metrics = score_predictions(y, pred)
            rows.append({"feature_set": name, "model": model_name, "n_features": len(cols), **metrics})
            if model_name == "rf":
                pred_rows[f"pred_{name}"] = pred

    # Shuffled null for the all-projection field. This tests whether feature
    # coherence matters beyond marginal distributions.
    cols = FEATURE_GROUPS["all_defect_projections"]
    shuffle_rows = []
    for rep in range(20):
        X = df[cols].to_numpy().copy()
        for j in range(X.shape[1]):
            rng.shuffle(X[:, j])
        model = RandomForestRegressor(n_estimators=220, min_samples_leaf=4, random_state=SEED + rep, n_jobs=1)
        pred = cross_val_predict(model, X, y, cv=kfold)
        shuffle_rows.append(score_predictions(y, pred))
    rows.append(
        {
            "feature_set": "shuffled_all_projection_null",
            "model": "rf",
            "n_features": len(cols),
            "r2": float(np.mean([r["r2"] for r in shuffle_rows])),
            "rmse": float(np.mean([r["rmse"] for r in shuffle_rows])),
            "spearman": float(np.mean([r["spearman"] for r in shuffle_rows])),
        }
    )
    return pd.DataFrame(rows).sort_values(["model", "r2"], ascending=[True, False]), pred_rows


def evaluate_lfo(df: pd.DataFrame, feature_sets: dict[str, list[str]]) -> pd.DataFrame:
    y = df["cdcl_hardness"].to_numpy()
    groups = df["family"].to_numpy()
    logo = LeaveOneGroupOut()
    rows = []
    for name, cols in feature_sets.items():
        X = df[cols].to_numpy()
        pred = np.full_like(y, fill_value=np.nan, dtype=float)
        for train, test in logo.split(X, y, groups):
            model = RandomForestRegressor(n_estimators=280, min_samples_leaf=4, random_state=SEED, n_jobs=1)
            model.fit(X[train], y[train])
            pred[test] = model.predict(X[test])
            rows.append({"feature_set": name, "held_out_family": groups[test][0], "n_test": len(test), **score_predictions(y[test], pred[test])})
        rows.append({"feature_set": name, "held_out_family": "ALL_LFO", "n_test": len(y), **score_predictions(y, pred)})
    return pd.DataFrame(rows)


def single_feature_correlations(df: pd.DataFrame) -> pd.DataFrame:
    cols = []
    for group in ["mean_projection", "tail_projection", "cluster_projection", "scale_projection", "macro_scalar", "standard_graph"]:
        for col in FEATURE_GROUPS[group]:
            if col not in cols:
                cols.append(col)
    y = df["cdcl_hardness"].to_numpy()
    rows = []
    for col in cols:
        rho = spearmanr(df[col].to_numpy(), y).correlation
        rows.append(
            {
                "feature": col,
                "projection": next(k for k, v in FEATURE_GROUPS.items() if col in v and k in ["mean_projection", "tail_projection", "cluster_projection", "scale_projection", "macro_scalar", "standard_graph"]),
                "spearman": float(0.0 if np.isnan(rho) else rho),
                "abs_spearman": float(abs(0.0 if np.isnan(rho) else rho)),
            }
        )
    return pd.DataFrame(rows).sort_values("abs_spearman", ascending=False)


def group_permutation_importance(df: pd.DataFrame, groups: dict[str, list[str]], rng: np.random.Generator) -> pd.DataFrame:
    y = df["cdcl_hardness"].to_numpy()
    cols = FEATURE_GROUPS["all_defect_projections"]
    X = df[cols].to_numpy()
    model = RandomForestRegressor(n_estimators=420, min_samples_leaf=4, random_state=SEED, n_jobs=1)
    model.fit(X, y)
    base = r2_score(y, model.predict(X))
    rows = []
    col_index = {c: i for i, c in enumerate(cols)}
    for group, group_cols in groups.items():
        idx = [col_index[c] for c in group_cols if c in col_index]
        drops = []
        for _ in range(30):
            Xp = X.copy()
            perm = rng.permutation(len(Xp))
            Xp[:, idx] = Xp[perm][:, idx]
            drops.append(base - r2_score(y, model.predict(Xp)))
        rows.append({"projection": group, "train_r2_base": float(base), "importance_drop_mean": float(np.mean(drops)), "importance_drop_std": float(np.std(drops))})
    return pd.DataFrame(rows).sort_values("importance_drop_mean", ascending=False)


def make_plots(
    cv: pd.DataFrame,
    lfo: pd.DataFrame,
    single: pd.DataFrame,
    group_imp: pd.DataFrame,
    pred: pd.DataFrame,
) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")

    rf = cv[cv["model"] == "rf"].copy()
    order = rf.sort_values("r2", ascending=False)["feature_set"].tolist()
    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    ax.bar(rf.set_index("feature_set").loc[order].index, rf.set_index("feature_set").loc[order, "r2"], color="#2f6f9f")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("5-fold CV R2")
    ax.set_title("Paper 57: Where does SAT/CDCL hardness signal live?")
    ax.tick_params(axis="x", rotation=35)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig1_projection_cv_r2.png", dpi=180)
    plt.close(fig)

    lfo_all = lfo[lfo["held_out_family"] == "ALL_LFO"].copy()
    order = lfo_all.sort_values("spearman", ascending=False)["feature_set"].tolist()
    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    ax.bar(lfo_all.set_index("feature_set").loc[order].index, lfo_all.set_index("feature_set").loc[order, "spearman"], color="#ef6f6c")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Leave-family-out Spearman")
    ax.set_title("Family-transfer stress test by projection")
    ax.tick_params(axis="x", rotation=35)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig2_projection_lfo_spearman.png", dpi=180)
    plt.close(fig)

    top = single.head(14).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.6, 5.8))
    colors = {
        "mean_projection": "#4c78a8",
        "tail_projection": "#f58518",
        "cluster_projection": "#54a24b",
        "scale_projection": "#b279a2",
        "macro_scalar": "#e45756",
        "standard_graph": "#72b7b2",
    }
    ax.barh(top["feature"], top["spearman"], color=[colors.get(x, "#777777") for x in top["projection"]])
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Spearman correlation with CDCL hardness")
    ax.set_title("Top single-feature correlations")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig3_top_single_feature_correlations.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 4.7))
    ax.bar(group_imp["projection"], group_imp["importance_drop_mean"], yerr=group_imp["importance_drop_std"], color="#7a5195")
    ax.set_ylabel("Grouped permutation R2 drop")
    ax.set_title("Projection-level importance inside all-defect model")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig4_group_permutation_importance.png", dpi=180)
    plt.close(fig)

    best = rf.sort_values("r2", ascending=False).iloc[0]["feature_set"]
    pred_col = f"pred_{best}"
    fig, ax = plt.subplots(figsize=(6, 5.5))
    for fam, sub in pred.groupby("family"):
        ax.scatter(sub["cdcl_hardness"], sub[pred_col], s=18, alpha=0.75, label=fam)
    lo = min(pred["cdcl_hardness"].min(), pred[pred_col].min())
    hi = max(pred["cdcl_hardness"].max(), pred[pred_col].max())
    ax.plot([lo, hi], [lo, hi], color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Observed CDCL hardness")
    ax.set_ylabel(f"Predicted hardness ({best})")
    ax.set_title("Best projection model prediction scatter")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig5_best_projection_scatter.png", dpi=180)
    plt.close(fig)


def main() -> None:
    rng = np.random.default_rng(SEED)
    input_csv = locate_paper55_csv()
    df = pd.read_csv(input_csv)

    feature_sets = dict(FEATURE_GROUPS)
    feature_sets.update(INCREMENTAL_GROUPS)
    cv, pred = evaluate_cv(df, feature_sets, rng)
    lfo = evaluate_lfo(df, feature_sets)
    single = single_feature_correlations(df)
    group_imp = group_permutation_importance(
        df,
        {
            "mean": FEATURE_GROUPS["mean_projection"],
            "tail": FEATURE_GROUPS["tail_projection"],
            "cluster": FEATURE_GROUPS["cluster_projection"],
            "scale": FEATURE_GROUPS["scale_projection"],
        },
        rng,
    )

    cv.to_csv(OUTDIR / "paper57_projection_cv_results.csv", index=False)
    lfo.to_csv(OUTDIR / "paper57_projection_leave_family_out.csv", index=False)
    single.to_csv(OUTDIR / "paper57_single_feature_correlations.csv", index=False)
    group_imp.to_csv(OUTDIR / "paper57_group_permutation_importance.csv", index=False)
    pred.to_csv(OUTDIR / "paper57_projection_predictions.csv", index=False)
    make_plots(cv, lfo, single, group_imp, pred)

    rf = cv[cv["model"] == "rf"].set_index("feature_set")
    lfo_all = lfo[lfo["held_out_family"] == "ALL_LFO"].set_index("feature_set")
    summary = {
        "seed": SEED,
        "input_csv": str(input_csv),
        "n_instances": int(len(df)),
        "families": sorted(df["family"].unique().tolist()),
        "target": "Paper 55 cdcl_hardness = log1p(mean_conflicts + 0.10*mean_decisions + 0.01*mean_propagations)",
        "best_rf_feature_set": str(rf["r2"].idxmax()),
        "best_rf_r2": float(rf["r2"].max()),
        "best_rf_spearman": float(rf.loc[rf["r2"].idxmax(), "spearman"]),
        "projection_rf_results": rf[["r2", "rmse", "spearman"]].to_dict(orient="index"),
        "lfo_overall": lfo_all[["r2", "rmse", "spearman"]].to_dict(orient="index"),
        "top_single_features": single.head(10).to_dict(orient="records"),
        "group_importance": group_imp.to_dict(orient="records"),
    }
    with (OUTDIR / "paper57_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
