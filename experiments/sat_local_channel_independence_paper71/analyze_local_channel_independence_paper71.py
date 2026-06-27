#!/usr/bin/env python3
"""Paper 71: local SAT channel independence diagnostics.

This is a post-hoc hypothesis-generation analysis. It does not modify the
Paper-69 gate, does not define a new SAT gate, and does not re-run the solver.

Question:
    Are local SAT channels L merely reparameterizations of existing MAAT
    supports H,B,S,V,R_rob, or do they contain residual information about the
    Gate-vs-MOMS gap?

Input:
    Paper-70 derived CSV outputs. No raw CNF files are required.

Output:
    Derived CSV/JSON/PNG diagnostics only.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


EPS = 1.0e-12
SEED = 71071


SUPPORTS = ["H", "B", "S", "V", "R_rob", "G_gate"]
PRIMARY_SUPPORTS = ["H", "B", "S", "V", "R_rob"]
LOCAL_CHANNELS = [
    "short_variable_entropy",
    "paper69_degree_cv",
    "variable_occurrence_cv",
    "short_clause_pressure",
    "moms_local_signal",
    "jw_variable_concentration",
    "clause_length_entropy",
    "weighted_clause_pressure",
]
TARGET = "regret_gate_to_moms"


def read_csv_dict(path: Path) -> List[dict]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv_dict(path: Path, rows: Sequence[dict], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def mean(vals: Sequence[float]) -> float:
    return sum(vals) / len(vals) if vals else 0.0


def rankdata(values: Sequence[float]) -> List[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        avg_rank = (i + j - 1) / 2.0 + 1.0
        for k in range(i, j):
            ranks[order[k]] = avg_rank
        i = j
    return ranks


def pearson(x: Sequence[float], y: Sequence[float]) -> float:
    if len(x) != len(y) or len(x) < 2:
        return 0.0
    mx, my = mean(x), mean(y)
    dx = [v - mx for v in x]
    dy = [v - my for v in y]
    sx = math.sqrt(sum(v * v for v in dx))
    sy = math.sqrt(sum(v * v for v in dy))
    if sx <= EPS or sy <= EPS:
        return 0.0
    return sum(a * b for a, b in zip(dx, dy)) / (sx * sy)


def spearman(x: Sequence[float], y: Sequence[float]) -> float:
    return pearson(rankdata(x), rankdata(y))


def matrix_from_rows(rows: Sequence[dict], columns: Sequence[str]) -> np.ndarray:
    return np.array([[float(r[c]) for c in columns] for r in rows], dtype=float)


def vector_from_rows(rows: Sequence[dict], column: str) -> np.ndarray:
    return np.array([float(r[column]) for r in rows], dtype=float)


def zscore_matrix(x: np.ndarray) -> np.ndarray:
    mu = x.mean(axis=0)
    sig = x.std(axis=0)
    sig[sig < EPS] = 1.0
    return (x - mu) / sig


def zscore_vector(y: np.ndarray) -> np.ndarray:
    sig = y.std()
    if sig < EPS:
        return y * 0.0
    return (y - y.mean()) / sig


def fit_ols(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """Return coefficients, predictions, R^2. Adds intercept automatically."""
    xz = zscore_matrix(x)
    yz = zscore_vector(y)
    x_design = np.column_stack([np.ones(len(xz)), xz])
    beta, *_ = np.linalg.lstsq(x_design, yz, rcond=None)
    pred = x_design @ beta
    ss_res = float(np.sum((yz - pred) ** 2))
    ss_tot = float(np.sum((yz - yz.mean()) ** 2))
    r2 = 1.0 - ss_res / (ss_tot + EPS)
    return beta, pred, r2


def lofo_r2(rows: Sequence[dict], x_cols: Sequence[str], y_col: str) -> float:
    """Leave-one-family-out predictive R^2 on z-scored train/test folds."""
    families = sorted({r["family"] for r in rows})
    y_all = vector_from_rows(rows, y_col)
    preds: List[float] = []
    actuals: List[float] = []
    for fam in families:
        train = [r for r in rows if r["family"] != fam]
        test = [r for r in rows if r["family"] == fam]
        x_train = matrix_from_rows(train, x_cols)
        y_train = vector_from_rows(train, y_col)
        x_test = matrix_from_rows(test, x_cols)
        y_test = vector_from_rows(test, y_col)
        mu_x = x_train.mean(axis=0)
        sig_x = x_train.std(axis=0)
        sig_x[sig_x < EPS] = 1.0
        mu_y = y_train.mean()
        sig_y = y_train.std()
        if sig_y < EPS:
            sig_y = 1.0
        xz_train = (x_train - mu_x) / sig_x
        yz_train = (y_train - mu_y) / sig_y
        beta, *_ = np.linalg.lstsq(np.column_stack([np.ones(len(xz_train)), xz_train]), yz_train, rcond=None)
        xz_test = (x_test - mu_x) / sig_x
        pred_z = np.column_stack([np.ones(len(xz_test)), xz_test]) @ beta
        pred = pred_z * sig_y + mu_y
        preds.extend(pred.tolist())
        actuals.extend(y_test.tolist())
    ss_res = float(np.sum((np.array(actuals) - np.array(preds)) ** 2))
    ss_tot = float(np.sum((y_all - y_all.mean()) ** 2))
    return 1.0 - ss_res / (ss_tot + EPS)


def residualize(rows: Sequence[dict], x_cols: Sequence[str], y_col: str) -> List[float]:
    x = matrix_from_rows(rows, x_cols)
    y = vector_from_rows(rows, y_col)
    _, pred, _ = fit_ols(x, y)
    yz = zscore_vector(y)
    return (yz - pred).tolist()


def model_summary(rows: Sequence[dict], x_cols: Sequence[str], target: str) -> dict:
    x = matrix_from_rows(rows, x_cols)
    y = vector_from_rows(rows, target)
    _, pred, r2 = fit_ols(x, y)
    return {
        "model": "+".join(x_cols),
        "n_features": len(x_cols),
        "in_sample_r2": r2,
        "lofo_r2": lofo_r2(rows, x_cols, target),
        "spearman_pred_target": spearman(pred.tolist(), zscore_vector(y).tolist()),
        "pearson_pred_target": pearson(pred.tolist(), zscore_vector(y).tolist()),
    }


def shuffle_null_for_residual(
    rows: Sequence[dict], x_cols: Sequence[str], local_col: str, target: str, n: int = 2000
) -> dict:
    rng = random.Random(SEED)
    base_resid = residualize(rows, x_cols, local_col)
    y = [float(r[target]) for r in rows]
    observed = spearman(base_resid, y)
    local_vals = [float(r[local_col]) for r in rows]
    null_vals: List[float] = []
    for _ in range(n):
        shuffled = local_vals[:]
        rng.shuffle(shuffled)
        tmp_rows = [dict(r, **{local_col: v}) for r, v in zip(rows, shuffled)]
        resid = residualize(tmp_rows, x_cols, local_col)
        null_vals.append(spearman(resid, y))
    abs_obs = abs(observed)
    p_two_sided = (1 + sum(abs(v) >= abs_obs for v in null_vals)) / (n + 1)
    return {
        "local_channel": local_col,
        "observed_residual_spearman": observed,
        "shuffle_abs_p": p_two_sided,
        "null_mean": mean(null_vals),
        "null_q025": float(np.quantile(null_vals, 0.025)),
        "null_q975": float(np.quantile(null_vals, 0.975)),
        "n_shuffles": n,
    }


def plot_residual_bars(rows: Sequence[dict], output: Path) -> None:
    data = []
    y = [float(r[TARGET]) for r in rows]
    for col in LOCAL_CHANNELS:
        resid = residualize(rows, SUPPORTS, col)
        data.append((col, spearman(resid, y)))
    data.sort(key=lambda item: abs(item[1]), reverse=True)
    names = [d[0] for d in data]
    vals = [d[1] for d in data]
    colors = ["#b84a3a" if v > 0 else "#386fa4" for v in vals]
    plt.figure(figsize=(8.5, 4.8))
    plt.barh(range(len(names)), vals, color=colors)
    plt.axvline(0, color="black", linewidth=0.8)
    plt.yticks(range(len(names)), [n.replace("_", " ") for n in names], fontsize=8)
    plt.xlabel("Spearman corr. of L residual with gate-to-MOMS regret")
    plt.title("Paper 71: residual local-channel information after H,B,S,V,R")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(output, dpi=180)
    plt.close()


def plot_model_r2(rows: Sequence[dict], output: Path) -> None:
    models = [
        ("supports", PRIMARY_SUPPORTS),
        ("supports+gate", SUPPORTS),
        ("local-only", LOCAL_CHANNELS),
        ("supports+local", PRIMARY_SUPPORTS + LOCAL_CHANNELS),
    ]
    summaries = [model_summary(rows, cols, TARGET) | {"label": label} for label, cols in models]
    labels = [s["label"] for s in summaries]
    in_sample = [s["in_sample_r2"] for s in summaries]
    lofo = [s["lofo_r2"] for s in summaries]
    lofo_clipped = [max(v, -5.0) for v in lofo]
    x = np.arange(len(labels))
    width = 0.36
    plt.figure(figsize=(8.2, 4.8))
    plt.bar(x - width / 2, in_sample, width, label="in-sample R2", color="#4c78a8")
    plt.bar(x + width / 2, lofo_clipped, width, label="leave-family-out R2 (clipped)", color="#f58518")
    plt.axhline(0, color="black", linewidth=0.8)
    for i, true_val in enumerate(lofo):
        label = f"{true_val:.1f}" if true_val < -5 else f"{true_val:.2f}"
        plt.text(i + width / 2, lofo_clipped[i] - 0.12, label, ha="center", va="top", fontsize=8)
    plt.xticks(x, labels, rotation=15, ha="right")
    plt.ylim(-5.8, 1.2)
    plt.ylabel("R2 against gate-to-MOMS regret")
    plt.title("Paper 71: explanatory value of supports vs local channels")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output, dpi=180)
    plt.close()


def plot_scatter(rows: Sequence[dict], output: Path, local_col: str = "short_variable_entropy") -> None:
    resid = residualize(rows, SUPPORTS, local_col)
    regrets = [float(r[TARGET]) for r in rows]
    families = sorted({r["family"] for r in rows})
    cmap = plt.get_cmap("tab10")
    colors = {fam: cmap(i % 10) for i, fam in enumerate(families)}
    plt.figure(figsize=(7.2, 5.0))
    for fam in families:
        idx = [i for i, r in enumerate(rows) if r["family"] == fam]
        plt.scatter(
            [resid[i] for i in idx],
            [regrets[i] for i in idx],
            label=fam.replace("satlib_", ""),
            s=44,
            alpha=0.88,
            edgecolor="black",
            linewidth=0.4,
            color=colors[fam],
        )
    plt.xlabel(f"residual {local_col} after supports")
    plt.ylabel("gate-to-MOMS regret")
    plt.title("Paper 71: local-channel residual vs missing MOMS advantage")
    plt.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(output, dpi=180)
    plt.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--paper70-output-dir",
        type=Path,
        default=Path("../sat_local_structural_channels_paper70/outputs_paper70_local_sat_channels"),
    )
    ap.add_argument("--output-dir", type=Path, default=Path("outputs_paper71_local_channel_independence"))
    args = ap.parse_args()

    rows = read_csv_dict(args.paper70_output_dir / "paper70_local_channel_features.csv")
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)

    support_corr_rows: List[dict] = []
    for l_col in LOCAL_CHANNELS:
        l_vals = [float(r[l_col]) for r in rows]
        for s_col in SUPPORTS:
            s_vals = [float(r[s_col]) for r in rows]
            support_corr_rows.append(
                {
                    "local_channel": l_col,
                    "support": s_col,
                    "spearman": spearman(l_vals, s_vals),
                    "pearson": pearson(l_vals, s_vals),
                }
            )
    write_csv_dict(
        output / "paper71_local_support_correlations.csv",
        support_corr_rows,
        ["local_channel", "support", "spearman", "pearson"],
    )

    independence_rows: List[dict] = []
    y = [float(r[TARGET]) for r in rows]
    for l_col in LOCAL_CHANNELS:
        l_vals = [float(r[l_col]) for r in rows]
        _, _, r2_from_supports = fit_ols(matrix_from_rows(rows, SUPPORTS), vector_from_rows(rows, l_col))
        resid = residualize(rows, SUPPORTS, l_col)
        independence_rows.append(
            {
                "local_channel": l_col,
                "spearman_L_regret": spearman(l_vals, y),
                "pearson_L_regret": pearson(l_vals, y),
                "r2_L_from_supports": r2_from_supports,
                "residual_spearman_regret": spearman(resid, y),
                "residual_pearson_regret": pearson(resid, y),
                "interpretation": "candidate_independent" if r2_from_supports < 0.5 and abs(spearman(resid, y)) > 0.25 else "weak_or_redundant",
            }
        )
    write_csv_dict(
        output / "paper71_local_independence_summary.csv",
        independence_rows,
        [
            "local_channel",
            "spearman_L_regret",
            "pearson_L_regret",
            "r2_L_from_supports",
            "residual_spearman_regret",
            "residual_pearson_regret",
            "interpretation",
        ],
    )

    model_rows = [
        model_summary(rows, PRIMARY_SUPPORTS, TARGET) | {"label": "supports"},
        model_summary(rows, SUPPORTS, TARGET) | {"label": "supports+gate"},
        model_summary(rows, LOCAL_CHANNELS, TARGET) | {"label": "local-only"},
        model_summary(rows, PRIMARY_SUPPORTS + LOCAL_CHANNELS, TARGET) | {"label": "supports+local"},
    ]
    write_csv_dict(
        output / "paper71_regret_model_comparison.csv",
        model_rows,
        ["label", "model", "n_features", "in_sample_r2", "lofo_r2", "spearman_pred_target", "pearson_pred_target"],
    )

    shuffle_rows = [
        shuffle_null_for_residual(rows, SUPPORTS, "short_variable_entropy", TARGET),
        shuffle_null_for_residual(rows, SUPPORTS, "paper69_degree_cv", TARGET),
        shuffle_null_for_residual(rows, SUPPORTS, "variable_occurrence_cv", TARGET),
    ]
    write_csv_dict(
        output / "paper71_residual_shuffle_nulls.csv",
        shuffle_rows,
        [
            "local_channel",
            "observed_residual_spearman",
            "shuffle_abs_p",
            "null_mean",
            "null_q025",
            "null_q975",
            "n_shuffles",
        ],
    )

    plot_residual_bars(rows, output / "fig1_residual_local_channels.png")
    plot_model_r2(rows, output / "fig2_model_comparison.png")
    plot_scatter(rows, output / "fig3_short_entropy_residual_vs_regret.png")

    top_residual = sorted(
        independence_rows,
        key=lambda r: abs(float(r["residual_spearman_regret"])),
        reverse=True,
    )[:5]
    summary = {
        "paper": 71,
        "title": "Decision-Local Structural Selection in SAT",
        "analysis_type": "post_hoc_hypothesis_generation_not_gate_update",
        "source": "Paper 70 derived local-channel features",
        "n_instances": len(rows),
        "families": sorted({r["family"] for r in rows}),
        "support_basis": SUPPORTS,
        "local_channels": LOCAL_CHANNELS,
        "target": TARGET,
        "top_residual_channels": top_residual,
        "model_comparison": model_rows,
        "shuffle_nulls": shuffle_rows,
        "interpretation_rule": (
            "Candidate local channels may motivate a future preregistered test, "
            "but do not modify Paper 69 or define a new gate here."
        ),
    }
    (output / "paper71_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
