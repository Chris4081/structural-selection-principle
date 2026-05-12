#!/usr/bin/env python3
"""Paper 53: Robustness of frozen structural-mode transfer.

This script extends the Paper 52 protocol in three ways:

1. It tests several source-only fitting rules instead of only NNLS.
2. It evaluates the frozen architectures on more than two domains.
3. It records whether the Paper 52 V/connectivity dominance is stable under
   fit-rule changes, source-domain changes, and shuffled-defect nulls.

All target-domain score directions are predeclared as "higher score = higher
structural quality".  No sign flips or target-domain retuning are allowed.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import pearsonr, spearmanr
from sklearn.feature_selection import mutual_info_regression
from sklearn.linear_model import Lasso


EPS = 1.0e-9
SEED = 4081
SECTORS = ["H", "B", "S", "V", "R_rob"]
FIT_RULES = ["equal", "nnls", "ridge_pos", "lasso_pos", "response_top20", "v_only"]


@dataclass
class DomainData:
    name: str
    df: pd.DataFrame
    supports: np.ndarray
    defects: np.ndarray
    quality: np.ndarray
    quality_label: str
    notes: str


def normalize01(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    lo = np.nanmin(x)
    hi = np.nanmax(x)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo < EPS:
        return np.zeros_like(x, dtype=float)
    return (x - lo) / (hi - lo)


def make_domain(
    name: str,
    df: pd.DataFrame,
    support_map: dict[str, str],
    quality: np.ndarray,
    quality_label: str,
    notes: str,
) -> DomainData:
    supports = np.column_stack([df[support_map[s]].to_numpy(dtype=float) for s in SECTORS])
    supports = np.clip(supports, EPS, 1.0)
    defects = -np.log(EPS + supports)
    quality = normalize01(np.asarray(quality, dtype=float))
    return DomainData(
        name=name,
        df=df,
        supports=supports,
        defects=defects,
        quality=quality,
        quality_label=quality_label,
        notes=notes,
    )


def load_domains(root: Path) -> list[DomainData]:
    domains: list[DomainData] = []

    sat50b = pd.read_csv(
        root
        / "experiments/sat_frustration_fields_paper50b/outputs_paper50b/paper50b_sat_instances.csv"
    )
    domains.append(
        make_domain(
            "SAT-Frustration",
            sat50b,
            {s: s for s in SECTORS},
            1.0 - normalize01(sat50b["log_nodes"].to_numpy(dtype=float)),
            "1 - minmax(log DPLL nodes)",
            "Paper 50b random-3-SAT frustration-field benchmark",
        )
    )

    quantum = pd.read_csv(
        root
        / "experiments/quantum_pointer_state_selection_paper51/outputs_paper51/paper51_pointer_instances.csv"
    )
    domains.append(
        make_domain(
            "Quantum-Pointer",
            quantum,
            {s: s for s in SECTORS},
            quantum["target"].to_numpy(dtype=float),
            "minmax(pointer robustness target)",
            "Paper 51 non-commuting Lindblad pointer-state benchmark",
        )
    )

    active = pd.read_csv(
        root / "experiments/active_respect_significance/outputs/active_respect_ensemble.csv"
    )
    domains.append(
        make_domain(
            "Active-Significance",
            active,
            {"H": "H", "B": "B", "S": "S_eff", "V": "V", "R_rob": "R_rob"},
            active["R_sig"].to_numpy(dtype=float),
            "minmax(active significance R_sig)",
            "Paper 44 active-structural-significance toy ensemble",
        )
    )

    societal = pd.read_csv(root / "experiments/societal_cci/outputs/societal_cci_ensemble.csv")
    domains.append(
        make_domain(
            "Societal-CCI",
            societal,
            {"H": "H", "B": "B", "S": "S_eff", "V": "V", "R_rob": "R_rob_soc"},
            societal["ASI_soc"].to_numpy(dtype=float),
            "minmax(societal constructive-impact ASI_soc)",
            "Societal CCI toy ensemble",
        )
    )

    return domains


def normalize_weights(w: np.ndarray) -> np.ndarray:
    w = np.asarray(w, dtype=float)
    w = np.where(np.isfinite(w), w, 0.0)
    w = np.clip(w, 0.0, None)
    if w.sum() <= EPS:
        return np.ones(len(SECTORS), dtype=float) / len(SECTORS)
    return w / w.sum()


def fit_weights(source: DomainData, rule: str) -> np.ndarray:
    D = source.defects
    cost = 1.0 - normalize01(source.quality)
    n = len(SECTORS)

    if rule == "equal":
        return np.ones(n, dtype=float) / n

    if rule == "v_only":
        w = np.zeros(n, dtype=float)
        w[SECTORS.index("V")] = 1.0
        return w

    if rule == "nnls":
        w, _ = nnls(D, cost)
        return normalize_weights(w)

    if rule == "ridge_pos":
        alpha = 0.05 * (np.trace(D.T @ D) / max(n, 1)) / max(len(D), 1)
        w = np.linalg.solve(D.T @ D + alpha * np.eye(n), D.T @ cost)
        return normalize_weights(w)

    if rule == "lasso_pos":
        # Positive sparse regression.  A small alpha is used because all
        # features are support-defects on comparable log scales.
        model = Lasso(alpha=0.001, positive=True, fit_intercept=True, max_iter=20000)
        model.fit(D, cost)
        return normalize_weights(model.coef_)

    if rule == "response_top20":
        # Source-only covariance/response proxy: compare the full ensemble to
        # the top-quality quintile and regularise the covariance inverse.
        threshold = np.quantile(source.quality, 0.80)
        top = D[source.quality >= threshold]
        delta = D.mean(axis=0) - top.mean(axis=0)
        C = np.cov(D, rowvar=False)
        eta = 0.05
        reg = eta * np.trace(C) / max(n, 1)
        w = np.linalg.solve(C + reg * np.eye(n), delta)
        return normalize_weights(w)

    raise ValueError(f"Unknown fit rule: {rule}")


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    value = spearmanr(x, y).correlation
    return float(0.0 if value is None or not np.isfinite(value) else value)


def safe_pearson(x: np.ndarray, y: np.ndarray) -> float:
    try:
        value = pearsonr(x, y).statistic
    except Exception:
        value = 0.0
    return float(0.0 if not np.isfinite(value) else value)


def safe_mi(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float).reshape(-1, 1)
    y = np.asarray(y, dtype=float)
    if np.nanstd(x) < EPS or np.nanstd(y) < EPS:
        return 0.0
    value = mutual_info_regression(x, y, random_state=SEED, n_neighbors=5)[0]
    return float(0.0 if not np.isfinite(value) else value)


def score_domain(target: DomainData, weights: np.ndarray) -> np.ndarray:
    return -target.defects @ weights


def evaluate(
    domains: list[DomainData],
    n_null: int,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows: list[dict] = []
    weight_rows: list[dict] = []
    null_rows: list[dict] = []

    for source in domains:
        for rule in FIT_RULES:
            w = fit_weights(source, rule)
            dominant_sector = SECTORS[int(np.argmax(w))]
            weight_rows.append(
                {
                    "source_domain": source.name,
                    "fit_rule": rule,
                    "dominant_sector": dominant_sector,
                    "v_weight": float(w[SECTORS.index("V")]),
                    **{f"lambda_{sector}": float(value) for sector, value in zip(SECTORS, w)},
                }
            )

            for target in domains:
                score = score_domain(target, w)
                transfer = f"{source.name}->{target.name}"
                is_cross = source.name != target.name
                rho = safe_spearman(score, target.quality)
                mi = safe_mi(score, target.quality)
                rows.append(
                    {
                        "source_domain": source.name,
                        "target_domain": target.name,
                        "transfer": transfer,
                        "fit_rule": rule,
                        "is_cross_domain": bool(is_cross),
                        "spearman": rho,
                        "pearson": safe_pearson(score, target.quality),
                        "mutual_info": mi,
                        "score_mean": float(score.mean()),
                        "score_std": float(score.std()),
                        "quality_mean": float(target.quality.mean()),
                        "quality_std": float(target.quality.std()),
                        "source_quality_label": source.quality_label,
                        "target_quality_label": target.quality_label,
                        "dominant_sector": dominant_sector,
                        "v_weight": float(w[SECTORS.index("V")]),
                    }
                )

                if is_cross and rule in {"nnls", "ridge_pos", "response_top20", "v_only"}:
                    null_rhos = []
                    null_mis = []
                    for _ in range(n_null):
                        shuffled = target.defects.copy()
                        for col in range(shuffled.shape[1]):
                            rng.shuffle(shuffled[:, col])
                        null_score = -shuffled @ w
                        null_rhos.append(safe_spearman(null_score, target.quality))
                        null_mis.append(safe_mi(null_score, target.quality))
                    null_rhos = np.asarray(null_rhos, dtype=float)
                    null_mis = np.asarray(null_mis, dtype=float)
                    null_rows.append(
                        {
                            "transfer": transfer,
                            "source_domain": source.name,
                            "target_domain": target.name,
                            "fit_rule": rule,
                            "n_null": n_null,
                            "real_spearman": rho,
                            "shuffle_spearman_mean": float(null_rhos.mean()),
                            "shuffle_spearman_p95": float(np.quantile(null_rhos, 0.95)),
                            "shuffle_spearman_p_value_greater": float(
                                (1 + np.sum(null_rhos >= rho)) / (n_null + 1)
                            ),
                            "real_mutual_info": mi,
                            "shuffle_mi_mean": float(null_mis.mean()),
                            "shuffle_mi_p95": float(np.quantile(null_mis, 0.95)),
                            "shuffle_mi_p_value_greater": float(
                                (1 + np.sum(null_mis >= mi)) / (n_null + 1)
                            ),
                        }
                    )

    results = pd.DataFrame(rows)
    weights = pd.DataFrame(weight_rows)
    nulls = pd.DataFrame(null_rows)

    cross = results[results["is_cross_domain"]].copy()
    summary_rows = []
    for rule, subset in cross.groupby("fit_rule"):
        summary_rows.append(
            {
                "fit_rule": rule,
                "cross_mean_spearman": float(subset["spearman"].mean()),
                "cross_median_spearman": float(subset["spearman"].median()),
                "cross_positive_fraction": float((subset["spearman"] > 0).mean()),
                "cross_mean_mutual_info": float(subset["mutual_info"].mean()),
                "mean_v_weight": float(weights[weights["fit_rule"] == rule]["v_weight"].mean()),
                "v_dominant_fraction": float(
                    (weights[weights["fit_rule"] == rule]["dominant_sector"] == "V").mean()
                ),
            }
        )
    rule_summary = pd.DataFrame(summary_rows).sort_values("cross_mean_spearman", ascending=False)

    return results, weights, nulls, rule_summary


def plot_rule_summary(rule_summary: pd.DataFrame, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    ordered = rule_summary.sort_values("cross_mean_spearman", ascending=False)
    x = np.arange(len(ordered))
    ax.bar(x, ordered["cross_mean_spearman"], color="#2f9df5", alpha=0.85)
    ax.axhline(0, color="#222222", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ordered["fit_rule"], rotation=30, ha="right")
    ax.set_ylabel("mean cross-domain Spearman")
    ax.set_title("Paper 53: transfer robustness across fitting rules")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "fig1_fit_rule_cross_domain_summary.png", dpi=180)
    plt.close(fig)


def plot_v_weight(weights: pd.DataFrame, outdir: Path) -> None:
    pivot = weights.pivot(index="source_domain", columns="fit_rule", values="v_weight").loc[
        :, FIT_RULES
    ]
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    image = ax.imshow(pivot.to_numpy(dtype=float), cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_title("Paper 53: V/connectivity weight by source and fit rule")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            value = float(pivot.iloc[i, j])
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", color="white" if value > 0.45 else "black")
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("lambda_V")
    fig.tight_layout()
    fig.savefig(outdir / "fig2_v_weight_heatmap.png", dpi=180)
    plt.close(fig)


def plot_transfer_heatmaps(results: pd.DataFrame, outdir: Path) -> None:
    domains = sorted(results["source_domain"].unique())
    rules = ["nnls", "ridge_pos", "response_top20", "v_only"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 9.5), sharex=True, sharey=True)
    vmax = max(0.75, float(np.nanmax(np.abs(results["spearman"]))) * 1.05)

    for ax, rule in zip(axes.flat, rules):
        subset = results[results["fit_rule"] == rule]
        matrix = np.full((len(domains), len(domains)), np.nan)
        for i, source in enumerate(domains):
            for j, target in enumerate(domains):
                row = subset[(subset["source_domain"] == source) & (subset["target_domain"] == target)]
                if not row.empty:
                    matrix[i, j] = float(row.iloc[0]["spearman"])
        image = ax.imshow(matrix, cmap="coolwarm", vmin=-vmax, vmax=vmax)
        ax.set_title(rule)
        ax.set_xticks(np.arange(len(domains)))
        ax.set_yticks(np.arange(len(domains)))
        ax.set_xticklabels(domains, rotation=35, ha="right", fontsize=8)
        ax.set_yticklabels(domains, fontsize=8)
        for i in range(len(domains)):
            for j in range(len(domains)):
                value = matrix[i, j]
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=8)
    fig.suptitle("Paper 53: frozen transfer matrices by fit rule", y=1.01)
    fig.colorbar(image, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02, label="Spearman rho")
    fig.tight_layout()
    fig.savefig(outdir / "fig3_transfer_heatmaps_by_rule.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_null_margin(nulls: pd.DataFrame, outdir: Path) -> None:
    if nulls.empty:
        return
    df = nulls.copy()
    df["margin_over_shuffle_p95"] = df["real_spearman"] - df["shuffle_spearman_p95"]
    df = df.sort_values("margin_over_shuffle_p95", ascending=False).head(14)
    labels = df["source_domain"] + "→" + df["target_domain"] + "\n" + df["fit_rule"]
    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    x = np.arange(len(df))
    colors = ["#35c48b" if v > 0 else "#d94a4a" for v in df["margin_over_shuffle_p95"]]
    ax.bar(x, df["margin_over_shuffle_p95"], color=colors, alpha=0.85)
    ax.axhline(0, color="#222222", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("real Spearman - shuffled p95")
    ax.set_title("Paper 53: strongest margins over shuffled-defect null")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "fig4_shuffle_null_margins.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help="Root of structural-selection-principle repository.",
    )
    parser.add_argument("--n-null", type=int, default=200)
    args = parser.parse_args()

    root = args.repo_root.resolve()
    outdir = Path(__file__).resolve().parent / "outputs_paper53"
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    domains = load_domains(root)
    results, weights, nulls, rule_summary = evaluate(domains, args.n_null, rng)

    results.to_csv(outdir / "paper53_transfer_results.csv", index=False)
    weights.to_csv(outdir / "paper53_frozen_weights.csv", index=False)
    nulls.to_csv(outdir / "paper53_shuffle_nulls.csv", index=False)
    rule_summary.to_csv(outdir / "paper53_fit_rule_summary.csv", index=False)

    plot_rule_summary(rule_summary, outdir)
    plot_v_weight(weights, outdir)
    plot_transfer_heatmaps(results, outdir)
    plot_null_margin(nulls, outdir)

    cross = results[results["is_cross_domain"]].copy()
    best_rule = rule_summary.iloc[0].to_dict() if len(rule_summary) else {}
    v_dominant_overall = float((weights["dominant_sector"] == "V").mean())
    strong_cross = cross[cross["spearman"] > 0.5]
    summary = {
        "seed": SEED,
        "n_null": int(args.n_null),
        "domains": [
            {
                "name": d.name,
                "n": int(len(d.df)),
                "quality_label": d.quality_label,
                "notes": d.notes,
            }
            for d in domains
        ],
        "fit_rules": FIT_RULES,
        "best_fit_rule_by_mean_cross_spearman": best_rule,
        "overall_v_dominant_fraction": v_dominant_overall,
        "max_cross_spearman": float(cross["spearman"].max()),
        "min_cross_spearman": float(cross["spearman"].min()),
        "mean_cross_spearman": float(cross["spearman"].mean()),
        "positive_cross_fraction": float((cross["spearman"] > 0).mean()),
        "strong_cross_transfers": strong_cross[
            ["source_domain", "target_domain", "fit_rule", "spearman", "v_weight", "dominant_sector"]
        ].to_dict(orient="records"),
    }
    (outdir / "paper53_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print("\nFit-rule summary:")
    print(rule_summary.to_string(index=False))
    print("\nWeights:")
    print(weights.to_string(index=False))


if __name__ == "__main__":
    main()
