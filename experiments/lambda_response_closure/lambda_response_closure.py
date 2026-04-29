#!/usr/bin/env python3
"""Response-theoretic closure for MAAT structural weights.

This script implements the linear-response / covariance closure

    delta <d_a> = - sum_b C_ab lambda_b

around the empirical reference measure mu_0.  Given a target selected
sub-ensemble, the sector weights are estimated as

    lambda = (C + eta tr(C)/A I)^(-1) (<d>_0 - <d>_target).

The result is not a claim that lambda_a are universal constants.  It is a
physically motivated effective derivation of lambda_a as response coefficients
of a concrete defect ensemble.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SECTORS = ["H", "B", "S", "V", "R"]
DEFECT_COLS = [f"d_{s}" for s in SECTORS]
RIDGE = 1.0e-3
EPS = 1.0e-12

HERE = Path(__file__).resolve().parent
INPUT_CSV = HERE.parent / "boundary_aware_lambda_calibration" / "maat_defects_fused.csv"
CLOSED_FIT_JSON = (
    HERE.parent
    / "boundary_aware_lambda_calibration"
    / "closed_maat_lambda_fit_results.json"
)
OUT_JSON = HERE / "lambda_response_closure_results.json"
PLOT_DIR = HERE / "plots"


def load_data() -> pd.DataFrame:
    if not INPUT_CSV.exists():
        raise FileNotFoundError(f"Missing input dataset: {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV)
    missing = [c for c in DEFECT_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Dataset is missing defect columns: {missing}")
    return df


def target_subensembles(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Define non-parametric target ensembles without using lambda weights."""
    score = df[DEFECT_COLS].mean(axis=1)
    n20 = max(1, int(0.2 * len(df)))
    targets = {
        "low_defect_20pct": df.assign(_score=score).nsmallest(n20, "_score"),
        "safe_and_core_safe": df[df["label"].isin(["safe", "core_safe"])],
        "safe_boundary_only": df[df["label"].eq("safe")],
        "not_violated": df[~df["label"].eq("violated")],
    }
    return {k: v.copy() for k, v in targets.items() if len(v) > 0}


def response_lambda(
    defects: np.ndarray,
    target_mean: np.ndarray,
    ridge: float = RIDGE,
) -> dict[str, np.ndarray | float]:
    """Compute raw and positivity-projected response weights."""
    ref_mean = defects.mean(axis=0)
    cov = np.cov(defects, rowvar=False, bias=False)
    scale = float(np.trace(cov) / cov.shape[0])
    cov_reg = cov + ridge * scale * np.eye(cov.shape[0])

    delta = ref_mean - target_mean
    raw = np.linalg.solve(cov_reg, delta)
    positive = np.clip(raw, 0.0, None)
    shares = positive / (positive.sum() + EPS)

    pred_linear = ref_mean - cov @ positive
    residual = pred_linear - target_mean
    residual_norm = float(
        np.linalg.norm(residual) / (np.linalg.norm(delta) + EPS)
    )

    return {
        "reference_mean": ref_mean,
        "target_mean": target_mean,
        "delta": delta,
        "covariance": cov,
        "covariance_regularized": cov_reg,
        "ridge": float(ridge),
        "ridge_scale": scale,
        "lambda_raw": raw,
        "lambda_positive": positive,
        "lambda_shares": shares,
        "linear_response_predicted_mean": pred_linear,
        "linear_response_residual": residual,
        "relative_residual_norm": residual_norm,
    }


def softmax_tilt(defects: np.ndarray, lambdas: np.ndarray) -> np.ndarray:
    energy = defects @ lambdas
    logw = -energy
    logw = logw - np.max(logw)
    weights = np.exp(logw)
    return weights / (weights.sum() + EPS)


def tilted_means(defects: np.ndarray, lambdas: np.ndarray) -> np.ndarray:
    return softmax_tilt(defects, lambdas) @ defects


def as_named(values: np.ndarray) -> dict[str, float]:
    return {s: float(v) for s, v in zip(SECTORS, values)}


def load_closed_fit() -> dict[str, float] | None:
    if not CLOSED_FIT_JSON.exists():
        return None
    data = json.loads(CLOSED_FIT_JSON.read_text())
    lambdas = data.get("lambdas", {})
    if not lambdas:
        return None
    return {s: float(lambdas[s]) for s in SECTORS}


def make_plots(
    results: dict[str, dict],
    corr: np.ndarray,
    closed_fit: dict[str, float] | None,
) -> None:
    PLOT_DIR.mkdir(exist_ok=True)

    target_names = list(results.keys())
    x = np.arange(len(SECTORS))

    # 1) Share comparison across target definitions.
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    width = 0.16
    offsets = np.linspace(-width * 1.5, width * 1.5, len(target_names))
    for offset, name in zip(offsets, target_names):
        shares = np.array([results[name]["lambda_shares"][s] for s in SECTORS])
        ax.bar(x + offset, shares, width=width, label=name.replace("_", " "))
    if closed_fit:
        closed = np.array([closed_fit[s] for s in SECTORS])
        closed_shares = closed / closed.sum()
        ax.plot(x, closed_shares, color="black", marker="o", linewidth=2.2,
                label="previous closed MaxEnt fit")
    ax.set_xticks(x)
    ax.set_xticklabels(SECTORS)
    ax.set_ylabel("normalised share")
    ax.set_title("Response-derived MAAT sector weights")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "lambda_response_shares.png", dpi=220)
    plt.close(fig)

    # 2) Correlation matrix of primitive defects under mu_0.
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    im = ax.imshow(corr, cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_xticks(x)
    ax.set_xticklabels(SECTORS)
    ax.set_yticks(x)
    ax.set_yticklabels(SECTORS)
    for i in range(len(SECTORS)):
        for j in range(len(SECTORS)):
            ax.text(j, i, f"{corr[i, j]:.2f}", ha="center", va="center", fontsize=9)
    ax.set_title("Defect correlation matrix under reference measure $\\mu_0$")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "lambda_response_correlation.png", dpi=220)
    plt.close(fig)

    # 3) Target matching for the safest target.
    default = "safe_and_core_safe" if "safe_and_core_safe" in results else target_names[0]
    ref = np.array([results[default]["reference_mean"][s] for s in SECTORS])
    target = np.array([results[default]["target_mean"][s] for s in SECTORS])
    pred = np.array([
        results[default]["linear_response_predicted_mean"][s] for s in SECTORS
    ])
    tilt = np.array([results[default]["tilted_mean"][s] for s in SECTORS])

    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    w = 0.20
    ax.bar(x - 1.5 * w, ref, width=w, label="reference $\\mu_0$")
    ax.bar(x - 0.5 * w, target, width=w, label="target")
    ax.bar(x + 0.5 * w, pred, width=w, label="linear response")
    ax.bar(x + 1.5 * w, tilt, width=w, label="exponential tilt")
    ax.set_xticks(x)
    ax.set_xticklabels(SECTORS)
    ax.set_ylabel("mean defect")
    ax.set_title(f"Target matching for {default.replace('_', ' ')}")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "lambda_response_target_match.png", dpi=220)
    plt.close(fig)

    # 4) Absolute lambda comparison to the previous closed fit.
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    selected = [
        name for name in ["low_defect_20pct", "safe_and_core_safe", "safe_boundary_only"]
        if name in results
    ]
    width = 0.18
    offsets = np.linspace(-width, width, len(selected))
    for offset, name in zip(offsets, selected):
        vals = np.array([results[name]["lambda_positive"][s] for s in SECTORS])
        ax.bar(x + offset, vals, width=width, label=name.replace("_", " "))
    if closed_fit:
        vals = np.array([closed_fit[s] for s in SECTORS])
        ax.plot(x, vals, color="black", marker="o", linewidth=2.2,
                label="previous closed MaxEnt fit")
    ax.set_xticks(x)
    ax.set_xticklabels(SECTORS)
    ax.set_ylabel("$\\lambda_a$")
    ax.set_title("Response closure versus previous closed fit")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "lambda_response_vs_closed_fit.png", dpi=220)
    plt.close(fig)


def main() -> None:
    df = load_data()
    defects = df[DEFECT_COLS].to_numpy(dtype=float)
    cov = np.cov(defects, rowvar=False, bias=False)
    std = np.sqrt(np.diag(cov))
    corr = cov / np.outer(std + EPS, std + EPS)

    closed_fit = load_closed_fit()
    targets = target_subensembles(df)

    results: dict[str, dict] = {}
    for name, sub in targets.items():
        target_mean = sub[DEFECT_COLS].mean(axis=0).to_numpy(dtype=float)
        rr = response_lambda(defects, target_mean)
        tilt_mean = tilted_means(defects, rr["lambda_positive"])

        results[name] = {
            "n_target": int(len(sub)),
            "reference_mean": as_named(rr["reference_mean"]),
            "target_mean": as_named(rr["target_mean"]),
            "delta": as_named(rr["delta"]),
            "lambda_raw": as_named(rr["lambda_raw"]),
            "lambda_positive": as_named(rr["lambda_positive"]),
            "lambda_shares": as_named(rr["lambda_shares"]),
            "linear_response_predicted_mean": as_named(
                rr["linear_response_predicted_mean"]
            ),
            "tilted_mean": as_named(tilt_mean),
            "linear_response_residual": as_named(rr["linear_response_residual"]),
            "relative_residual_norm": rr["relative_residual_norm"],
        }

    output = {
        "status": "effective_response_closure_not_universal_constants",
        "method": {
            "reference_measure": "empirical uniform measure over fused defect ensemble",
            "formula": "lambda = (Cov_mu0[d] + eta tr(C)/A I)^(-1)(<d>_mu0 - <d>_target)",
            "ridge": RIDGE,
            "positive_projection": "lambda_a < 0 is projected to 0 because lambda_a is interpreted as a non-negative penalty weight",
            "sectors": SECTORS,
            "input_csv": str(INPUT_CSV),
        },
        "dataset": {
            "n_samples": int(len(df)),
            "source_counts": df["source"].value_counts().to_dict(),
            "label_counts": df["label"].value_counts().to_dict(),
        },
        "reference_mean": as_named(defects.mean(axis=0)),
        "covariance": {
            SECTORS[i]: {SECTORS[j]: float(cov[i, j]) for j in range(len(SECTORS))}
            for i in range(len(SECTORS))
        },
        "correlation": {
            SECTORS[i]: {SECTORS[j]: float(corr[i, j]) for j in range(len(SECTORS))}
            for i in range(len(SECTORS))
        },
        "targets": results,
        "previous_closed_maxent_fit": closed_fit,
    }

    OUT_JSON.write_text(json.dumps(output, indent=2))
    make_plots(results, corr, closed_fit)

    print("\n=== MAAT Lambda Response Closure ===")
    print("Input:", INPUT_CSV)
    print("Samples:", len(df))
    print("Ridge:", RIDGE)
    print("\nReference mean defects:")
    for s, value in output["reference_mean"].items():
        print(f"  d_{s}: {value:.6f}")

    for name, data in results.items():
        print(f"\n--- Target: {name} (n={data['n_target']}) ---")
        print("lambda_positive:")
        for s in SECTORS:
            print(f"  lambda_{s}: {data['lambda_positive'][s]:.6f}"
                  f"  share={data['lambda_shares'][s]:.4f}")
        print("relative linear-response residual:",
              f"{data['relative_residual_norm']:.6f}")

    print("\nSaved:", OUT_JSON)
    print("Plots:", PLOT_DIR)


if __name__ == "__main__":
    main()
