#!/usr/bin/env python3
"""Paper 54: External ML sample-hardness benchmark for MAAT structural modes.

This script is intentionally conservative:

* It uses public scikit-learn datasets, not repository-internal toy targets.
* Sample hardness is predeclared as repeated-CV misclassification frequency.
* MAAT supports are computed from raw feature geometry and labels only.
* Baselines and shuffled-feature nulls are reported alongside MAAT features.

The goal is not to prove universality.  The goal is to break the internal
artifact loop of Papers 52--53 and ask whether MAAT-style structural features
carry signal on public external benchmarks.
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

os.environ.setdefault("LOKY_MAX_CPU_COUNT", "4")
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).parent / ".mplconfig"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import entropy, pearsonr, spearmanr
from sklearn.base import clone
from sklearn.datasets import load_breast_cancer, load_digits, load_wine
from sklearn.decomposition import PCA
from sklearn.ensemble import IsolationForest, RandomForestClassifier
from sklearn.linear_model import Ridge
from sklearn.metrics import (
    mean_squared_error,
    r2_score,
    roc_auc_score,
)
from sklearn.model_selection import RepeatedKFold, RepeatedStratifiedKFold
from sklearn.neighbors import NearestNeighbors
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


EPS = 1.0e-9
SEED = 4081
SECTORS = ["H", "B", "S", "V", "R_rob"]
DATASETS: dict[str, Callable] = {
    "breast_cancer": load_breast_cancer,
    "wine": load_wine,
    "digits": load_digits,
}


@dataclass
class DatasetResult:
    name: str
    sample_df: pd.DataFrame
    scalar_df: pd.DataFrame
    performance_df: pd.DataFrame
    null_df: pd.DataFrame
    diagnostics: dict


def normalize01(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    lo = np.nanmin(x)
    hi = np.nanmax(x)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo < EPS:
        return np.zeros_like(x, dtype=float)
    return (x - lo) / (hi - lo)


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    value = spearmanr(x, y).correlation
    return float(0.0 if value is None or not np.isfinite(value) else value)


def safe_pearson(x: np.ndarray, y: np.ndarray) -> float:
    try:
        value = pearsonr(x, y).statistic
    except Exception:
        value = 0.0
    return float(0.0 if not np.isfinite(value) else value)


def safe_auc(y_true: np.ndarray, score: np.ndarray) -> float:
    y_true = np.asarray(y_true, dtype=int)
    if len(np.unique(y_true)) < 2:
        return float("nan")
    try:
        value = roc_auc_score(y_true, score)
    except Exception:
        value = float("nan")
    return float(value)


def entropy_normalized(counts: np.ndarray) -> float:
    total = counts.sum()
    if total <= 0:
        return 0.0
    probs = counts[counts > 0] / total
    if len(probs) <= 1:
        return 0.0
    return float(entropy(probs, base=np.e) / np.log(len(counts)))


def compute_sample_hardness(
    X: np.ndarray,
    y: np.ndarray,
    *,
    repeats: int,
    folds: int,
    seed: int,
) -> pd.DataFrame:
    """Repeated-CV sample hardness using a fixed Random Forest classifier."""

    n = len(y)
    errors = np.zeros(n, dtype=float)
    counts = np.zeros(n, dtype=float)
    entropy_sum = np.zeros(n, dtype=float)
    true_conf_sum = np.zeros(n, dtype=float)
    fold_accs: list[float] = []

    cv = RepeatedStratifiedKFold(n_splits=folds, n_repeats=repeats, random_state=seed)
    for fold_id, (train_idx, test_idx) in enumerate(cv.split(X, y)):
        clf = RandomForestClassifier(
            n_estimators=100,
            max_depth=5,
            min_samples_leaf=2,
            class_weight="balanced",
            random_state=seed + fold_id,
            n_jobs=1,
        )
        clf.fit(X[train_idx], y[train_idx])
        pred = clf.predict(X[test_idx])
        proba = clf.predict_proba(X[test_idx])
        classes = clf.classes_
        class_index = np.array([np.where(classes == yi)[0][0] for yi in y[test_idx]])

        errors[test_idx] += pred != y[test_idx]
        counts[test_idx] += 1.0
        entropy_sum[test_idx] += entropy(proba.T, base=np.e)
        true_conf_sum[test_idx] += proba[np.arange(len(test_idx)), class_index]
        fold_accs.append(float(np.mean(pred == y[test_idx])))

    misclassification_rate = errors / np.maximum(counts, 1.0)
    boundary_entropy = entropy_sum / np.maximum(counts, 1.0)
    true_confidence = true_conf_sum / np.maximum(counts, 1.0)

    # Primary target: repeated-CV error rate.
    # Secondary target: entropy-smoothed boundary hardness, useful when a small
    # benchmark is nearly solved by the classifier and error rates are sparse.
    boundary_hardness = normalize01(boundary_entropy)
    confidence_loss = normalize01(1.0 - true_confidence)
    hybrid_hardness = normalize01(
        0.70 * normalize01(misclassification_rate)
        + 0.20 * boundary_hardness
        + 0.10 * confidence_loss
    )

    return pd.DataFrame(
        {
            "sample_id": np.arange(n, dtype=int),
            "label": y,
            "cv_count": counts,
            "misclassification_rate": misclassification_rate,
            "boundary_entropy": boundary_entropy,
            "true_class_confidence": true_confidence,
            "boundary_hardness": boundary_hardness,
            "hybrid_hardness": hybrid_hardness,
            "mean_fold_accuracy": np.mean(fold_accs),
        }
    )


def local_geometry_features(X: np.ndarray, y: np.ndarray, *, k: int, seed: int) -> pd.DataFrame:
    n, p = X.shape
    classes = np.unique(y)
    n_classes = len(classes)
    k_eff = int(min(k, n - 1))

    nn = NearestNeighbors(n_neighbors=k_eff + 1, metric="euclidean")
    nn.fit(X)
    distances, indices = nn.kneighbors(X)
    distances = distances[:, 1:]
    indices = indices[:, 1:]

    same_frac = np.zeros(n, dtype=float)
    label_entropy = np.zeros(n, dtype=float)
    label_margin = np.zeros(n, dtype=float)
    local_spread = distances.mean(axis=1)

    for i in range(n):
        neigh_labels = y[indices[i]]
        counts = np.array([np.sum(neigh_labels == c) for c in classes], dtype=float)
        same = np.sum(neigh_labels == y[i])
        same_frac[i] = same / k_eff
        label_entropy[i] = entropy_normalized(counts)
        other_counts = counts.copy()
        own_pos = int(np.where(classes == y[i])[0][0])
        own = other_counts[own_pos] / k_eff
        other_counts[own_pos] = 0.0
        strongest_other = other_counts.max() / k_eff
        label_margin[i] = np.clip(own - strongest_other, 0.0, 1.0)

    knn_density = normalize01(1.0 / (local_spread + EPS))

    centroids = np.vstack([X[y == c].mean(axis=0) for c in classes])
    own_centroid = np.vstack([centroids[int(np.where(classes == yi)[0][0])] for yi in y])
    centroid_dist = np.linalg.norm(X - own_centroid, axis=1)
    centroid_dist_norm = normalize01(centroid_dist)

    pca_components = max(1, min(p - 1, n - 1, max(2, int(np.ceil(0.30 * p)))))
    pca = PCA(n_components=pca_components, random_state=seed)
    Z = pca.fit_transform(X)
    X_hat = pca.inverse_transform(Z)
    pca_error = np.mean((X - X_hat) ** 2, axis=1)
    pca_error_norm = normalize01(pca_error)

    iso = IsolationForest(n_estimators=200, contamination="auto", random_state=seed)
    iso.fit(X)
    # sklearn returns higher decision_function for inliers, so invert it.
    isolation_anomaly = normalize01(-iso.decision_function(X))

    feature_norm = normalize01(np.linalg.norm(X, axis=1))
    local_spread_norm = normalize01(local_spread)

    # MAAT support language.  These are deliberately simple and predeclared:
    # H: local label coherence,
    # B: local class-margin balance,
    # S: controlled activity in a normal local-spread window,
    # V: local kNN connectivity/density.
    log_spread = np.log(local_spread + EPS)
    spread_med = np.median(log_spread)
    spread_mad = np.median(np.abs(log_spread - spread_med)) + EPS
    S_eff = np.exp(-0.5 * ((log_spread - spread_med) / (1.4826 * spread_mad + EPS)) ** 2)

    H = np.clip(1.0 - label_entropy, EPS, 1.0)
    B = np.clip(label_margin, EPS, 1.0)
    S = np.clip(S_eff, EPS, 1.0)
    V = np.clip(knn_density, EPS, 1.0)
    R_resp = np.power(H * B * V, 1.0 / 3.0)
    R_rob = np.minimum(R_resp, np.power(H * B * S * V, 1.0 / 4.0))
    R_rob = np.clip(R_rob, EPS, 1.0)

    supports = np.column_stack([H, B, S, V, R_rob])
    defects = -np.log(np.clip(supports, EPS, 1.0))
    defect_mean = defects.mean(axis=1)
    defect_max = defects.max(axis=1)
    defect_std = defects.std(axis=1)
    local_defect_mean = defect_mean[indices].mean(axis=1)
    local_defect_max = defect_mean[indices].max(axis=1)
    local_defect_std = defect_mean[indices].std(axis=1)
    defect_tail_rank = pd.Series(defect_mean).rank(pct=True).to_numpy()
    hotspot_threshold = np.quantile(defect_mean, 0.90)
    hotspot_neighbor_fraction = (defect_mean[indices] >= hotspot_threshold).mean(axis=1)

    return pd.DataFrame(
        {
            "same_label_fraction": same_frac,
            "label_entropy": label_entropy,
            "label_margin": label_margin,
            "knn_density": knn_density,
            "local_spread": local_spread_norm,
            "centroid_distance": centroid_dist_norm,
            "pca_reconstruction_error": pca_error_norm,
            "isolation_anomaly": isolation_anomaly,
            "feature_norm": feature_norm,
            "H": H,
            "B": B,
            "S": S,
            "V": V,
            "R_resp": R_resp,
            "R_rob": R_rob,
            "D_H": defects[:, 0],
            "D_B": defects[:, 1],
            "D_S": defects[:, 2],
            "D_V": defects[:, 3],
            "D_R_rob": defects[:, 4],
            "defect_mean": defect_mean,
            "defect_max": defect_max,
            "defect_std": defect_std,
            "local_defect_mean": local_defect_mean,
            "local_defect_max": local_defect_max,
            "local_defect_std": local_defect_std,
            "defect_tail_rank": defect_tail_rank,
            "hotspot_neighbor_fraction": hotspot_neighbor_fraction,
            "maat_scalar_cost": defect_mean,
            "maat_robustness_loss": 1.0 - R_rob,
            "label_disagreement": 1.0 - same_frac,
            "density_anomaly": 1.0 - knn_density,
        }
    )


FEATURE_SETS: dict[str, list[str]] = {
    "MAAT_supports": ["H", "B", "S", "V", "R_rob"],
    "MAAT_defects": ["D_H", "D_B", "D_S", "D_V", "D_R_rob"],
    "MAAT_v14_field": [
        "H",
        "B",
        "S",
        "V",
        "R_rob",
        "D_H",
        "D_B",
        "D_S",
        "D_V",
        "D_R_rob",
        "defect_mean",
        "defect_max",
        "defect_std",
        "local_defect_mean",
        "local_defect_max",
        "local_defect_std",
        "defect_tail_rank",
        "hotspot_neighbor_fraction",
    ],
    "external_geometry": [
        "knn_density",
        "local_spread",
        "centroid_distance",
        "pca_reconstruction_error",
        "isolation_anomaly",
        "feature_norm",
    ],
    "label_geometry": [
        "same_label_fraction",
        "label_entropy",
        "label_margin",
        "label_disagreement",
    ],
    "combined_geometry": [
        "knn_density",
        "local_spread",
        "centroid_distance",
        "pca_reconstruction_error",
        "isolation_anomaly",
        "feature_norm",
        "same_label_fraction",
        "label_entropy",
        "label_margin",
        "label_disagreement",
    ],
}


SCALAR_SCORES = {
    "MAAT_scalar_cost": "maat_scalar_cost",
    "MAAT_robustness_loss": "maat_robustness_loss",
    "label_disagreement": "label_disagreement",
    "density_anomaly": "density_anomaly",
    "centroid_distance": "centroid_distance",
    "pca_reconstruction_error": "pca_reconstruction_error",
    "isolation_anomaly": "isolation_anomaly",
}


def cross_val_predict_ridge(
    X_feat: np.ndarray,
    target: np.ndarray,
    *,
    seed: int,
    repeats: int,
) -> np.ndarray:
    pred_sum = np.zeros(len(target), dtype=float)
    pred_count = np.zeros(len(target), dtype=float)
    cv = RepeatedKFold(n_splits=5, n_repeats=repeats, random_state=seed)
    model = make_pipeline(StandardScaler(), Ridge(alpha=1.0))
    for train_idx, test_idx in cv.split(X_feat):
        estimator = clone(model)
        estimator.fit(X_feat[train_idx], target[train_idx])
        pred = estimator.predict(X_feat[test_idx])
        pred_sum[test_idx] += pred
        pred_count[test_idx] += 1.0
    return pred_sum / np.maximum(pred_count, 1.0)


def evaluate_predictions(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    *,
    dataset: str,
    model_name: str,
    target_name: str,
) -> dict:
    # Rank-based top-hard split avoids degenerate AUROC labels when many easy
    # samples have exactly zero misclassification rate.
    top_k = max(1, int(np.ceil(0.25 * len(y_true))))
    hard_binary = np.zeros(len(y_true), dtype=bool)
    hard_binary[np.argsort(y_true)[-top_k:]] = True
    hard_threshold = np.min(y_true[hard_binary])
    rmse = float(np.sqrt(mean_squared_error(y_true, y_pred)))
    return {
        "dataset": dataset,
        "model": model_name,
        "target": target_name,
        "spearman": safe_spearman(y_pred, y_true),
        "pearson": safe_pearson(y_pred, y_true),
        "r2": float(r2_score(y_true, y_pred)),
        "rmse": rmse,
        "auroc_top_hard": safe_auc(hard_binary.astype(int), y_pred),
        "hard_threshold": float(hard_threshold),
    }


def evaluate_dataset(
    name: str,
    loader: Callable,
    *,
    output_dir: Path,
    hardness_repeats: int,
    eval_repeats: int,
    n_null: int,
    seed: int,
) -> DatasetResult:
    raw = loader()
    X = np.asarray(raw.data, dtype=float)
    y = np.asarray(raw.target)
    X_scaled = StandardScaler().fit_transform(X)
    k = min(15, max(5, int(np.sqrt(len(y)))))

    hardness_df = compute_sample_hardness(
        X_scaled,
        y,
        repeats=hardness_repeats,
        folds=5,
        seed=seed,
    )
    geometry_df = local_geometry_features(X_scaled, y, k=k, seed=seed)
    sample_df = pd.concat([hardness_df, geometry_df], axis=1)
    sample_df.insert(0, "dataset", name)

    primary = sample_df["misclassification_rate"].to_numpy(dtype=float)
    secondary = sample_df["hybrid_hardness"].to_numpy(dtype=float)

    scalar_rows: list[dict] = []
    for score_name, col in SCALAR_SCORES.items():
        score = sample_df[col].to_numpy(dtype=float)
        scalar_rows.append(
            evaluate_predictions(
                primary,
                score,
                dataset=name,
                model_name=score_name,
                target_name="misclassification_rate",
            )
        )
        scalar_rows.append(
            evaluate_predictions(
                secondary,
                score,
                dataset=name,
                model_name=score_name,
                target_name="hybrid_hardness",
            )
        )

    perf_rows: list[dict] = []
    rng = np.random.default_rng(seed + 991)
    null_rows: list[dict] = []

    for feature_set, columns in FEATURE_SETS.items():
        X_feat = sample_df[columns].to_numpy(dtype=float)
        for target_name, target in [
            ("misclassification_rate", primary),
            ("hybrid_hardness", secondary),
        ]:
            pred = cross_val_predict_ridge(
                X_feat,
                target,
                seed=seed,
                repeats=eval_repeats,
            )
            perf_rows.append(
                evaluate_predictions(
                    target,
                    pred,
                    dataset=name,
                    model_name=feature_set,
                    target_name=target_name,
                )
            )

            if feature_set == "MAAT_v14_field":
                for null_id in range(n_null):
                    X_null = X_feat.copy()
                    for j in range(X_null.shape[1]):
                        rng.shuffle(X_null[:, j])
                    null_pred = cross_val_predict_ridge(
                        X_null,
                        target,
                        seed=seed + 1000 + null_id,
                        repeats=max(2, eval_repeats // 2),
                    )
                    row = evaluate_predictions(
                        target,
                        null_pred,
                        dataset=name,
                        model_name="MAAT_v14_field_shuffled_null",
                        target_name=target_name,
                    )
                    row["null_id"] = null_id
                    null_rows.append(row)

    scalar_df = pd.DataFrame(scalar_rows)
    performance_df = pd.DataFrame(perf_rows)
    null_df = pd.DataFrame(null_rows)

    diagnostics = {
        "dataset": name,
        "n_samples": int(len(y)),
        "n_features": int(X.shape[1]),
        "n_classes": int(len(np.unique(y))),
        "knn_k": int(k),
        "hardness_repeats": int(hardness_repeats),
        "eval_repeats": int(eval_repeats),
        "mean_cv_accuracy": float(sample_df["mean_fold_accuracy"].iloc[0]),
        "mean_misclassification_rate": float(primary.mean()),
        "fraction_nonzero_hardness": float(np.mean(primary > 0)),
        "max_misclassification_rate": float(primary.max()),
    }

    sample_df.to_csv(output_dir / f"paper54_{name}_samples.csv", index=False)
    scalar_df.to_csv(output_dir / f"paper54_{name}_scalar_scores.csv", index=False)
    performance_df.to_csv(output_dir / f"paper54_{name}_performance.csv", index=False)
    null_df.to_csv(output_dir / f"paper54_{name}_nulls.csv", index=False)

    return DatasetResult(name, sample_df, scalar_df, performance_df, null_df, diagnostics)


def make_figures(
    all_samples: pd.DataFrame,
    scalar_scores: pd.DataFrame,
    performance: pd.DataFrame,
    nulls: pd.DataFrame,
    output_dir: Path,
) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.8), sharey=False)
    for ax, (dataset, group) in zip(axes, all_samples.groupby("dataset")):
        ax.hist(group["misclassification_rate"], bins=20, color="#1f77b4", alpha=0.80)
        ax.set_title(dataset.replace("_", " "))
        ax.set_xlabel("CV misclassification rate")
        ax.set_ylabel("samples")
    fig.suptitle("External sample-hardness targets")
    fig.tight_layout()
    fig.savefig(output_dir / "fig1_hardness_distributions.png", dpi=180)
    plt.close(fig)

    primary = performance[performance["target"] == "misclassification_rate"]
    pivot = primary.pivot(index="model", columns="dataset", values="spearman")
    ordered_models = [
        "MAAT_supports",
        "MAAT_defects",
        "MAAT_v14_field",
        "external_geometry",
        "label_geometry",
        "combined_geometry",
    ]
    pivot = pivot.reindex(ordered_models)
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    im = ax.imshow(pivot.values, cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_xticks(np.arange(len(pivot.columns)), labels=[c.replace("_", "\n") for c in pivot.columns])
    ax.set_yticks(np.arange(len(pivot.index)), labels=pivot.index)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            value = pivot.iloc[i, j]
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=9)
    ax.set_title("Repeated-CV Ridge prediction of external sample hardness")
    fig.colorbar(im, ax=ax, label="Spearman correlation")
    fig.tight_layout()
    fig.savefig(output_dir / "fig2_feature_set_spearman_heatmap.png", dpi=180)
    plt.close(fig)

    summary = (
        primary.groupby("model", as_index=False)
        .agg(mean_spearman=("spearman", "mean"), mean_auroc=("auroc_top_hard", "mean"))
        .sort_values("mean_spearman", ascending=False)
    )
    fig, ax = plt.subplots(figsize=(9.0, 4.6))
    x = np.arange(len(summary))
    ax.bar(x, summary["mean_spearman"], color="#2ca02c", alpha=0.82)
    ax.set_ylabel("Mean Spearman across datasets")
    ax.set_title("Aggregate external benchmark performance")
    ax.set_xticks(x, labels=summary["model"], rotation=35, ha="right")
    fig.tight_layout()
    fig.savefig(output_dir / "fig3_aggregate_performance.png", dpi=180)
    plt.close(fig)

    scalar_primary = scalar_scores[scalar_scores["target"] == "misclassification_rate"]
    scalar_pivot = scalar_primary.pivot(index="model", columns="dataset", values="spearman")
    fig, ax = plt.subplots(figsize=(8.4, 5.1))
    im = ax.imshow(scalar_pivot.values, cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_xticks(
        np.arange(len(scalar_pivot.columns)),
        labels=[c.replace("_", "\n") for c in scalar_pivot.columns],
    )
    ax.set_yticks(np.arange(len(scalar_pivot.index)), labels=scalar_pivot.index)
    for i in range(scalar_pivot.shape[0]):
        for j in range(scalar_pivot.shape[1]):
            value = scalar_pivot.iloc[i, j]
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=9)
    ax.set_title("Predeclared scalar scores versus external hardness")
    fig.colorbar(im, ax=ax, label="Spearman correlation")
    fig.tight_layout()
    fig.savefig(output_dir / "fig4_scalar_score_heatmap.png", dpi=180)
    plt.close(fig)

    rows = []
    actual = performance[
        (performance["model"] == "MAAT_v14_field")
        & (performance["target"] == "misclassification_rate")
    ]
    null_primary = nulls[nulls["target"] == "misclassification_rate"]
    for dataset in actual["dataset"].unique():
        a = actual[actual["dataset"] == dataset]["spearman"].iloc[0]
        ns = null_primary[null_primary["dataset"] == dataset]["spearman"].to_numpy(dtype=float)
        rows.append(
            {
                "dataset": dataset,
                "actual": float(a),
                "null_mean": float(np.mean(ns)),
                "null_p95": float(np.quantile(ns, 0.95)),
                "margin_vs_p95": float(a - np.quantile(ns, 0.95)),
            }
        )
    null_summary = pd.DataFrame(rows)
    x = np.arange(len(null_summary))
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    ax.bar(x - 0.18, null_summary["actual"], width=0.36, label="MAAT v1.4 field")
    ax.bar(x + 0.18, null_summary["null_p95"], width=0.36, label="shuffled-null p95")
    ax.set_xticks(x, labels=[d.replace("_", "\n") for d in null_summary["dataset"]])
    ax.set_ylabel("Spearman correlation")
    ax.set_title("MAAT field signal against shuffled-defect null")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "fig5_shuffled_null_margin.png", dpi=180)
    plt.close(fig)

    null_summary.to_csv(output_dir / "paper54_null_summary.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent / "outputs_paper54")
    parser.add_argument("--hardness-repeats", type=int, default=25)
    parser.add_argument("--eval-repeats", type=int, default=10)
    parser.add_argument("--n-null", type=int, default=40)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    results: list[DatasetResult] = []
    for name, loader in DATASETS.items():
        print(f"[paper54] processing {name}...")
        results.append(
            evaluate_dataset(
                name,
                loader,
                output_dir=args.output_dir,
                hardness_repeats=args.hardness_repeats,
                eval_repeats=args.eval_repeats,
                n_null=args.n_null,
                seed=args.seed,
            )
        )

    all_samples = pd.concat([r.sample_df for r in results], ignore_index=True)
    scalar_scores = pd.concat([r.scalar_df for r in results], ignore_index=True)
    performance = pd.concat([r.performance_df for r in results], ignore_index=True)
    nulls = pd.concat([r.null_df for r in results], ignore_index=True)

    all_samples.to_csv(args.output_dir / "paper54_all_samples.csv", index=False)
    scalar_scores.to_csv(args.output_dir / "paper54_scalar_scores.csv", index=False)
    performance.to_csv(args.output_dir / "paper54_feature_set_performance.csv", index=False)
    nulls.to_csv(args.output_dir / "paper54_shuffled_nulls.csv", index=False)

    primary = performance[performance["target"] == "misclassification_rate"].copy()
    aggregate = (
        primary.groupby("model", as_index=False)
        .agg(
            mean_spearman=("spearman", "mean"),
            median_spearman=("spearman", "median"),
            mean_auroc=("auroc_top_hard", "mean"),
            mean_r2=("r2", "mean"),
            positive_fraction=("spearman", lambda x: float(np.mean(np.asarray(x) > 0))),
        )
        .sort_values("mean_spearman", ascending=False)
    )
    aggregate.to_csv(args.output_dir / "paper54_aggregate_model_summary.csv", index=False)

    scalar_primary = scalar_scores[scalar_scores["target"] == "misclassification_rate"].copy()
    scalar_aggregate = (
        scalar_primary.groupby("model", as_index=False)
        .agg(
            mean_spearman=("spearman", "mean"),
            median_spearman=("spearman", "median"),
            mean_auroc=("auroc_top_hard", "mean"),
            positive_fraction=("spearman", lambda x: float(np.mean(np.asarray(x) > 0))),
        )
        .sort_values("mean_spearman", ascending=False)
    )
    scalar_aggregate.to_csv(args.output_dir / "paper54_aggregate_scalar_summary.csv", index=False)

    make_figures(all_samples, scalar_scores, performance, nulls, args.output_dir)

    summary = {
        "seed": args.seed,
        "hardness_repeats": args.hardness_repeats,
        "eval_repeats": args.eval_repeats,
        "n_null": args.n_null,
        "datasets": {r.name: r.diagnostics for r in results},
        "aggregate_model_summary": aggregate.to_dict(orient="records"),
        "aggregate_scalar_summary": scalar_aggregate.to_dict(orient="records"),
    }
    with open(args.output_dir / "paper54_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\n[paper54] aggregate model summary")
    print(aggregate.to_string(index=False))
    print("\n[paper54] aggregate scalar summary")
    print(scalar_aggregate.to_string(index=False))
    print(f"\n[paper54] outputs written to {args.output_dir}")


if __name__ == "__main__":
    main()
