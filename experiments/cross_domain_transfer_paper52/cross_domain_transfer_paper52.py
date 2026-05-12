#!/usr/bin/env python3
"""Paper 52: Cross-domain transfer for frozen MAAT defect architectures.

The experiment intentionally avoids retuning on the target domain.

Protocol
--------
1. Load one source domain and one target domain.
2. Convert the common supports H,B,S,V,R_rob into support-defects
   D_a = -log(epsilon + Gamma_a).
3. Fit non-negative sector weights on the source domain only by regressing
   source defects against the predeclared source cost 1-quality.
4. Freeze the weights and evaluate the score
       S_frozen(X) = - sum_a lambda_a D_a(X)
   on the target domain.
5. Compare against equal-weight MAAT, scalar MAAT, R_rob, shuffled-defect
   nulls, and lambda-permutation nulls.

The score convention is predeclared: higher score means higher structural
quality.  No sign flips are allowed after observing target-domain results.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import pearsonr, spearmanr
from sklearn.feature_selection import mutual_info_regression


EPS = 1.0e-9
SEED = 4081
SECTORS = ["H", "B", "S", "V", "R_rob"]


@dataclass
class DomainData:
    name: str
    df: pd.DataFrame
    quality: np.ndarray
    defects: np.ndarray
    scalar_cost: np.ndarray
    notes: str


def normalize01(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    lo = np.nanmin(x)
    hi = np.nanmax(x)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo < EPS:
        return np.zeros_like(x, dtype=float)
    return (x - lo) / (hi - lo)


def support_defects(df: pd.DataFrame) -> np.ndarray:
    supports = df[SECTORS].to_numpy(dtype=float)
    supports = np.clip(supports, EPS, None)
    return -np.log(EPS + supports)


def load_sat(root: Path) -> DomainData:
    path = root / "experiments/sat_frustration_fields_paper50b/outputs_paper50b/paper50b_sat_instances.csv"
    df = pd.read_csv(path)
    # Predeclared quality: easier instances are structurally higher quality.
    quality = 1.0 - normalize01(df["log_nodes"].to_numpy(dtype=float))
    scalar_cost = df["F_MAAT"].to_numpy(dtype=float)
    return DomainData(
        name="SAT",
        df=df,
        quality=quality,
        defects=support_defects(df),
        scalar_cost=scalar_cost,
        notes="quality = 1 - minmax(log DPLL nodes)",
    )


def load_quantum(root: Path) -> DomainData:
    path = root / "experiments/quantum_pointer_state_selection_paper51/outputs_paper51/paper51_pointer_instances.csv"
    df = pd.read_csv(path)
    # Predeclared quality: higher pointer-state fidelity is better.
    quality = normalize01(df["target"].to_numpy(dtype=float))
    scalar_cost = df["F_maat"].to_numpy(dtype=float)
    return DomainData(
        name="Quantum",
        df=df,
        quality=quality,
        defects=support_defects(df),
        scalar_cost=scalar_cost,
        notes="quality = minmax(pointer robustness / fidelity target)",
    )


def fit_frozen_weights(source: DomainData) -> np.ndarray:
    """Fit non-negative source-only weights and normalise them to sum one."""
    cost = 1.0 - normalize01(source.quality)
    weights, _ = nnls(source.defects, cost)
    if weights.sum() <= EPS:
        weights = np.ones(len(SECTORS), dtype=float)
    weights = weights / weights.sum()
    return weights


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


def metric_row(
    transfer: str,
    model: str,
    score: np.ndarray,
    quality: np.ndarray,
    details: str = "",
) -> dict:
    score = np.asarray(score, dtype=float)
    quality = np.asarray(quality, dtype=float)
    return {
        "transfer": transfer,
        "model": model,
        "spearman": safe_spearman(score, quality),
        "pearson": safe_pearson(score, quality),
        "mutual_info": safe_mi(score, quality),
        "score_mean": float(np.mean(score)),
        "score_std": float(np.std(score)),
        "quality_mean": float(np.mean(quality)),
        "quality_std": float(np.std(quality)),
        "details": details,
    }


def evaluate_transfer(
    source: DomainData,
    target: DomainData,
    n_null: int,
    rng: np.random.Generator,
) -> tuple[list[dict], list[dict], np.ndarray]:
    transfer_name = f"{source.name}->{target.name}"
    weights = fit_frozen_weights(source)
    equal = np.ones(len(SECTORS), dtype=float) / len(SECTORS)

    frozen_score = -target.defects @ weights
    equal_score = -target.defects @ equal
    scalar_score = -target.scalar_cost
    rrob_score = target.df["R_rob"].to_numpy(dtype=float)

    rows = [
        metric_row(
            transfer_name,
            "MAAT frozen source weights",
            frozen_score,
            target.quality,
            "weights fitted on source only",
        ),
        metric_row(
            transfer_name,
            "MAAT equal weights",
            equal_score,
            target.quality,
            "predeclared equal-sector architecture",
        ),
        metric_row(
            transfer_name,
            "Scalar MAAT cost",
            scalar_score,
            target.quality,
            "negative scalar F_MAAT/F_maat",
        ),
        metric_row(
            transfer_name,
            "R_rob only",
            rrob_score,
            target.quality,
            "emergent robustness support only",
        ),
    ]

    real_spearman = rows[0]["spearman"]
    real_mi = rows[0]["mutual_info"]

    shuffled_spearman = []
    shuffled_mi = []
    lambda_spearman = []
    lambda_mi = []

    for _ in range(n_null):
        shuffled_defects = target.defects.copy()
        for col in range(shuffled_defects.shape[1]):
            rng.shuffle(shuffled_defects[:, col])
        shuffled_score = -shuffled_defects @ weights
        shuffled_spearman.append(safe_spearman(shuffled_score, target.quality))
        shuffled_mi.append(safe_mi(shuffled_score, target.quality))

        permuted_weights = rng.permutation(weights)
        permuted_score = -target.defects @ permuted_weights
        lambda_spearman.append(safe_spearman(permuted_score, target.quality))
        lambda_mi.append(safe_mi(permuted_score, target.quality))

    null_rows = []
    for null_name, rho_vals, mi_vals in [
        ("shuffled defects", np.asarray(shuffled_spearman), np.asarray(shuffled_mi)),
        ("lambda permutation", np.asarray(lambda_spearman), np.asarray(lambda_mi)),
    ]:
        null_rows.append(
            {
                "transfer": transfer_name,
                "null_model": null_name,
                "n_null": n_null,
                "spearman_mean": float(np.mean(rho_vals)),
                "spearman_std": float(np.std(rho_vals)),
                "spearman_p95": float(np.quantile(rho_vals, 0.95)),
                "spearman_p_value_greater": float((1 + np.sum(rho_vals >= real_spearman)) / (n_null + 1)),
                "mutual_info_mean": float(np.mean(mi_vals)),
                "mutual_info_std": float(np.std(mi_vals)),
                "mutual_info_p95": float(np.quantile(mi_vals, 0.95)),
                "mutual_info_p_value_greater": float((1 + np.sum(mi_vals >= real_mi)) / (n_null + 1)),
            }
        )

    return rows, null_rows, weights


def plot_weights(weights_df: pd.DataFrame, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    x = np.arange(len(SECTORS))
    width = 0.36
    for idx, (_, row) in enumerate(weights_df.iterrows()):
        vals = [row[f"lambda_{s}"] for s in SECTORS]
        ax.bar(x + (idx - 0.5) * width, vals, width=width, label=row["source_domain"])
    ax.set_xticks(x)
    ax.set_xticklabels(SECTORS)
    ax.set_ylim(0, max(0.55, float(weights_df[[f"lambda_{s}" for s in SECTORS]].to_numpy().max()) * 1.15))
    ax.set_ylabel("normalised frozen weight")
    ax.set_title("Paper 52: source-only frozen sector weights")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "fig1_frozen_weights.png", dpi=180)
    plt.close(fig)


def plot_transfer_results(results: pd.DataFrame, nulls: pd.DataFrame, outdir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    transfers = ["SAT->Quantum", "Quantum->SAT"]
    colors = {
        "MAAT frozen source weights": "#24c6dc",
        "MAAT equal weights": "#72e06a",
        "Scalar MAAT cost": "#f5a623",
        "R_rob only": "#ff5a78",
    }
    for ax, transfer in zip(axes, transfers):
        subset = results[results["transfer"] == transfer]
        models = subset["model"].tolist()
        vals = subset["spearman"].to_numpy(dtype=float)
        xpos = np.arange(len(models))
        ax.bar(xpos, vals, color=[colors.get(m, "#999999") for m in models], alpha=0.9)
        nsub = nulls[nulls["transfer"] == transfer]
        for _, row in nsub.iterrows():
            style = "--" if row["null_model"] == "shuffled defects" else ":"
            ax.axhline(row["spearman_p95"], color="#cccccc", linestyle=style, linewidth=1.2)
        ax.set_xticks(xpos)
        ax.set_xticklabels(["frozen", "equal", "scalar", "Rrob"], rotation=0)
        ax.set_title(transfer)
        ax.grid(axis="y", alpha=0.25)
        ax.axhline(0, color="#222222", linewidth=0.8)
    axes[0].set_ylabel("Spearman rho with predeclared quality")
    fig.suptitle("Paper 52: cross-domain transfer without target retuning", y=1.02)
    fig.tight_layout()
    fig.savefig(outdir / "fig2_transfer_spearman.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_scatter(
    source: DomainData,
    target: DomainData,
    weights: np.ndarray,
    outdir: Path,
) -> None:
    score = -target.defects @ weights
    rho = safe_spearman(score, target.quality)
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    ax.scatter(score, target.quality, s=18, alpha=0.55, edgecolor="none")
    ax.set_xlabel("frozen MAAT score")
    ax.set_ylabel("target quality")
    ax.set_title(f"{source.name}->{target.name}: frozen transfer (rho={rho:.3f})")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fname = f"fig3_scatter_{source.name.lower()}_to_{target.name.lower()}.png"
    fig.savefig(outdir / fname, dpi=180)
    plt.close(fig)


def plot_transfer_matrix(results: pd.DataFrame, outdir: Path) -> None:
    """Plot the frozen-weight source/target transfer matrix."""
    sources = ["SAT", "Quantum"]
    targets = ["SAT", "Quantum"]
    matrix = np.full((len(sources), len(targets)), np.nan)

    frozen = results[results["model"] == "MAAT frozen source weights"]
    for i, source in enumerate(sources):
        for j, target in enumerate(targets):
            transfer = f"{source}->{target}"
            row = frozen[frozen["transfer"] == transfer]
            if not row.empty:
                matrix[i, j] = float(row.iloc[0]["spearman"])

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    vmax = max(0.65, float(np.nanmax(np.abs(matrix))) * 1.05)
    image = ax.imshow(matrix, cmap="coolwarm", vmin=-vmax, vmax=vmax)

    ax.set_xticks(np.arange(len(targets)))
    ax.set_yticks(np.arange(len(sources)))
    ax.set_xticklabels(targets)
    ax.set_yticklabels(sources)
    ax.set_xlabel("Test domain")
    ax.set_ylabel("Frozen source domain")
    ax.set_title("Paper 52: frozen-weight transfer matrix")

    for i in range(len(sources)):
        for j in range(len(targets)):
            value = matrix[i, j]
            ax.text(
                j,
                i,
                f"{value:.3f}",
                ha="center",
                va="center",
                color="white" if abs(value) > 0.32 else "#111111",
                fontsize=12,
                fontweight="bold",
            )

    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Spearman rho")
    fig.tight_layout()
    fig.savefig(outdir / "fig4_transfer_matrix.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help="Root of structural-selection-principle repository.",
    )
    parser.add_argument("--n-null", type=int, default=500)
    args = parser.parse_args()

    root = args.repo_root.resolve()
    outdir = Path(__file__).resolve().parent / "outputs_paper52"
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    sat = load_sat(root)
    quantum = load_quantum(root)
    domains = [sat, quantum]

    all_rows: list[dict] = []
    all_nulls: list[dict] = []
    weight_rows: list[dict] = []
    transfer_weights: dict[str, np.ndarray] = {}

    for source in domains:
        weights = fit_frozen_weights(source)
        row = {
            "source_domain": source.name,
            "source_notes": source.notes,
        }
        row.update({f"lambda_{sector}": float(value) for sector, value in zip(SECTORS, weights)})
        weight_rows.append(row)

    for source, target in [(sat, sat), (quantum, quantum), (sat, quantum), (quantum, sat)]:
        rows, null_rows, weights = evaluate_transfer(source, target, args.n_null, rng)
        all_rows.extend(rows)
        all_nulls.extend(null_rows)
        transfer_weights[f"{source.name}->{target.name}"] = weights

    results = pd.DataFrame(all_rows)
    nulls = pd.DataFrame(all_nulls)
    weights_df = pd.DataFrame(weight_rows)

    results.to_csv(outdir / "paper52_transfer_results.csv", index=False)
    nulls.to_csv(outdir / "paper52_shuffle_nulls.csv", index=False)
    weights_df.to_csv(outdir / "paper52_frozen_weights.csv", index=False)

    # Main cross-domain rows only for headline summary.
    cross = results[
        (results["model"] == "MAAT frozen source weights")
        & (results["transfer"].isin(["SAT->Quantum", "Quantum->SAT"]))
    ].copy()
    equal_cross = results[
        (results["model"] == "MAAT equal weights")
        & (results["transfer"].isin(["SAT->Quantum", "Quantum->SAT"]))
    ].set_index("transfer")
    scalar_cross = results[
        (results["model"] == "Scalar MAAT cost")
        & (results["transfer"].isin(["SAT->Quantum", "Quantum->SAT"]))
    ].set_index("transfer")

    summary = {
        "seed": SEED,
        "sectors": SECTORS,
        "n_sat": int(len(sat.df)),
        "n_quantum": int(len(quantum.df)),
        "n_null": int(args.n_null),
        "score_convention": "higher score = higher structural quality; no target-domain sign retuning",
        "transfers": {},
    }
    for _, row in cross.iterrows():
        transfer = row["transfer"]
        summary["transfers"][transfer] = {
            "frozen_spearman": float(row["spearman"]),
            "frozen_mutual_info": float(row["mutual_info"]),
            "equal_spearman": float(equal_cross.loc[transfer, "spearman"]),
            "scalar_spearman": float(scalar_cross.loc[transfer, "spearman"]),
            "delta_spearman_vs_equal": float(row["spearman"] - equal_cross.loc[transfer, "spearman"]),
            "delta_spearman_vs_scalar": float(row["spearman"] - scalar_cross.loc[transfer, "spearman"]),
        }
    (outdir / "paper52_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    plot_weights(weights_df, outdir)
    plot_transfer_results(results, nulls, outdir)
    plot_scatter(sat, quantum, transfer_weights["SAT->Quantum"], outdir)
    plot_scatter(quantum, sat, transfer_weights["Quantum->SAT"], outdir)
    plot_transfer_matrix(results, outdir)

    print(json.dumps(summary, indent=2))
    print("\nTransfer results:")
    print(results.to_string(index=False))
    print("\nNull results:")
    print(nulls.to_string(index=False))


if __name__ == "__main__":
    main()
