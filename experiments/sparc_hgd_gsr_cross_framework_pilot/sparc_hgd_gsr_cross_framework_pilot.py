#!/usr/bin/env python3
"""
SPARC MAAT x HGD-GSR cross-framework pilot.

Purpose
-------
This script joins galaxy-level MAAT structural-residual summaries from Paper 65
with an externally supplied HGD-GSR Q_bar table. It tests whether MAAT's
D_struct and HGD-GSR's Q_bar are related, complementary, redundant, or
independently predictive of the same mismatch channels.

The script does not redistribute HGD-GSR data. If the Q_bar file is missing, it
creates a template CSV listing the SPARC galaxy names expected by the MAAT
summary and exits with a clear "waiting for data" status.

Expected HGD-GSR CSV columns
----------------------------
Required:
    galaxy, Q_bar

Accepted aliases:
    galaxy: galaxy, Galaxy, name, Name, SPARC, galaxy_name
    Q_bar: Q_bar, q_bar, Qbar, qbar, hgd_qbar, HGD_Q_bar

Run
---
    python3 sparc_hgd_gsr_cross_framework_pilot.py

Then place a filled CSV at:
    data/hgd_gsr_qbar.csv
and rerun the script.
"""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

_CACHE = Path(tempfile.gettempdir()) / "maat_sparc_hgd_cross_framework_cache"
_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

GALAXY_ALIASES = ["galaxy", "Galaxy", "name", "Name", "SPARC", "galaxy_name"]
QBAR_ALIASES = ["Q_bar", "q_bar", "Qbar", "qbar", "hgd_qbar", "HGD_Q_bar"]
EPS = 1.0e-12


def find_first(columns: Iterable[str], aliases: list[str]) -> str | None:
    cols = set(columns)
    for alias in aliases:
        if alias in cols:
            return alias
    return None


def safe_corr(x: Iterable[float], y: Iterable[float], method: str = "spearman") -> float:
    s = pd.DataFrame({"x": list(x), "y": list(y)}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(s) < 4 or s["x"].nunique() < 2 or s["y"].nunique() < 2:
        return float("nan")
    return float(s["x"].corr(s["y"], method=method))


def standardize_frame(df: pd.DataFrame) -> pd.DataFrame:
    out = df.astype(float).copy()
    for col in out.columns:
        sd = out[col].std(ddof=0)
        if not np.isfinite(sd) or sd < EPS:
            out[col] = 0.0
        else:
            out[col] = (out[col] - out[col].mean()) / sd
    return out


def r2_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    denom = float(np.sum((y_true - y_true.mean()) ** 2))
    if denom < EPS:
        return float("nan")
    return float(1.0 - np.sum((y_true - y_pred) ** 2) / denom)


def make_folds(n: int, k: int, seed: int) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    return [fold for fold in np.array_split(idx, k) if len(fold)]


def cv_ols(df: pd.DataFrame, features: list[str], target: str, k: int, seed: int) -> dict[str, float]:
    cols = features + [target]
    work = df[cols].replace([np.inf, -np.inf], np.nan).dropna()
    if len(work) < max(12, 2 * len(features) + 4):
        return {"n": len(work), "cv_r2": float("nan"), "cv_spearman": float("nan")}

    x_all = standardize_frame(work[features]).to_numpy(dtype=float)
    y_all = work[target].astype(float).to_numpy()
    folds = make_folds(len(work), min(k, len(work)), seed)
    pred = np.full(len(work), np.nan, dtype=float)

    for test_idx in folds:
        train_mask = np.ones(len(work), dtype=bool)
        train_mask[test_idx] = False
        x_train = x_all[train_mask]
        y_train = y_all[train_mask]
        x_test = x_all[test_idx]
        x_train_i = np.column_stack([np.ones(len(x_train)), x_train])
        x_test_i = np.column_stack([np.ones(len(x_test)), x_test])
        beta, *_ = np.linalg.lstsq(x_train_i, y_train, rcond=None)
        pred[test_idx] = x_test_i @ beta

    return {
        "n": int(len(work)),
        "cv_r2": r2_score(y_all, pred),
        "cv_spearman": safe_corr(y_all, pred, "spearman"),
    }


def bootstrap_correlations(
    df: pd.DataFrame,
    pairs: list[tuple[str, str]],
    n_bootstrap: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    n = len(df)
    for x, y in pairs:
        vals = []
        sub = df[[x, y]].replace([np.inf, -np.inf], np.nan).dropna()
        if len(sub) < 8:
            continue
        arr = sub.to_numpy()
        for _ in range(n_bootstrap):
            sample = arr[rng.integers(0, len(arr), len(arr))]
            vals.append(safe_corr(sample[:, 0], sample[:, 1], "spearman"))
        vals = np.asarray(vals, dtype=float)
        rows.append(
            {
                "x": x,
                "y": y,
                "n": int(len(sub)),
                "spearman": safe_corr(sub[x], sub[y], "spearman"),
                "boot_p05": float(np.nanpercentile(vals, 5)),
                "boot_p50": float(np.nanpercentile(vals, 50)),
                "boot_p95": float(np.nanpercentile(vals, 95)),
            }
        )
    return pd.DataFrame(rows)


def write_template(maat: pd.DataFrame, template_path: Path, output_dir: Path) -> None:
    template_path.parent.mkdir(parents=True, exist_ok=True)
    if not template_path.exists():
        pd.DataFrame({"galaxy": sorted(maat["galaxy"].astype(str).unique()), "Q_bar": ""}).to_csv(
            template_path, index=False
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    contract = output_dir / "DATA_CONTRACT.md"
    contract.write_text(
        "# HGD-GSR Q_bar Data Contract\n\n"
        "Place a collaborator-provided CSV at `data/hgd_gsr_qbar.csv` with columns:\n\n"
        "```text\n"
        "galaxy,Q_bar\n"
        "NGC2403,0.123\n"
        "...\n"
        "```\n\n"
        "Accepted galaxy aliases: galaxy, Galaxy, name, Name, SPARC, galaxy_name.\n"
        "Accepted Q_bar aliases: Q_bar, q_bar, Qbar, qbar, hgd_qbar, HGD_Q_bar.\n\n"
        "The pipeline joins by exact galaxy name after stripping whitespace.\n"
        "No HGD-GSR data are redistributed by this repository unless explicitly\n"
        "provided under clear permission/licence terms by the data owner.\n",
        encoding="utf-8",
    )
    summary = {
        "status": "waiting_for_hgd_gsr_qbar",
        "message": "Fill data/hgd_gsr_qbar.csv using the generated template and rerun.",
        "template_path": str(template_path),
        "n_template_galaxies": int(maat["galaxy"].nunique()),
    }
    (output_dir / "cross_framework_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


def load_qbar(path: Path) -> pd.DataFrame:
    q = pd.read_csv(path)
    galaxy_col = find_first(q.columns, GALAXY_ALIASES)
    qbar_col = find_first(q.columns, QBAR_ALIASES)
    if galaxy_col is None or qbar_col is None:
        raise ValueError(
            f"Q_bar CSV needs galaxy and Q_bar columns. Found columns: {list(q.columns)}"
        )
    out = q[[galaxy_col, qbar_col]].rename(columns={galaxy_col: "galaxy", qbar_col: "Q_bar"})
    out["galaxy"] = out["galaxy"].astype(str).str.strip()
    out["Q_bar"] = pd.to_numeric(out["Q_bar"], errors="coerce")
    out = out.dropna(subset=["galaxy", "Q_bar"]).drop_duplicates("galaxy")
    return out


def plot_scatter(joined: pd.DataFrame, output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    sc = ax.scatter(
        joined["mean_D_struct"],
        joined["Q_bar"],
        c=joined["nfw_like_chi2"],
        cmap="viridis",
        s=42,
        alpha=0.82,
        edgecolor="black",
        linewidth=0.25,
    )
    ax.set_xlabel("MAAT galaxy-level mean D_struct")
    ax.set_ylabel("HGD-GSR Q_bar")
    ax.set_title("Cross-framework structural coordinate comparison")
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("NFW-like chi2 mismatch")
    fig.tight_layout()
    fig.savefig(output_dir / "fig1_dstruct_vs_qbar.png", dpi=180)
    plt.close(fig)


def plot_cv(cv: pd.DataFrame, output_dir: Path) -> None:
    if cv.empty:
        return
    targets = list(cv["target"].drop_duplicates())
    models = list(cv["model"].drop_duplicates())
    pivot = cv.pivot(index="model", columns="target", values="cv_spearman").reindex(models)
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    im = ax.imshow(pivot.to_numpy(dtype=float), aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_xticks(range(len(targets)), targets, rotation=25, ha="right")
    ax.set_yticks(range(len(models)), models)
    ax.set_title("Cross-validated Spearman prediction by target")
    for i, model in enumerate(models):
        for j, target in enumerate(targets):
            val = pivot.loc[model, target]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, label="CV Spearman")
    fig.tight_layout()
    fig.savefig(output_dir / "fig2_cv_model_comparison.png", dpi=180)
    plt.close(fig)


def plot_bootstrap(boot: pd.DataFrame, output_dir: Path) -> None:
    if boot.empty:
        return
    labels = [f"{r.x}\nvs\n{r.y}" for r in boot.itertuples()]
    y = boot["spearman"].to_numpy(dtype=float)
    low = y - boot["boot_p05"].to_numpy(dtype=float)
    high = boot["boot_p95"].to_numpy(dtype=float) - y
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    ax.errorbar(range(len(boot)), y, yerr=[low, high], fmt="o", capsize=4)
    ax.axhline(0, color="black", lw=1)
    ax.set_xticks(range(len(boot)), labels, rotation=35, ha="right")
    ax.set_ylabel("Spearman rho with 90% bootstrap interval")
    ax.set_title("Bootstrap stability of cross-framework signals")
    fig.tight_layout()
    fig.savefig(output_dir / "fig3_bootstrap_correlations.png", dpi=180)
    plt.close(fig)


def main() -> None:
    here = Path(__file__).resolve().parent
    default_maat = here.parent / "sparc_ii_paper65" / "outputs_paper65" / (
        "paper65_galaxy_summary_with_proxy_morphology.csv"
    )
    parser = argparse.ArgumentParser()
    parser.add_argument("--maat-summary", type=Path, default=default_maat)
    parser.add_argument("--hgd-qbar", type=Path, default=here / "data" / "hgd_gsr_qbar.csv")
    parser.add_argument("--output-dir", type=Path, default=here / "outputs_cross_framework")
    parser.add_argument(
        "--targets",
        default="nfw_like_chi2,nfw_like_rmse_v2,rar_mean_sigma_residual,mean_residual_fraction",
    )
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=66066)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    maat = pd.read_csv(args.maat_summary)
    maat["galaxy"] = maat["galaxy"].astype(str).str.strip()

    if not args.hgd_qbar.exists():
        write_template(maat, here / "data" / "hgd_gsr_qbar_template.csv", args.output_dir)
        return

    qbar = load_qbar(args.hgd_qbar)
    joined = maat.merge(qbar, on="galaxy", how="inner")
    joined.to_csv(args.output_dir / "cross_framework_joined.csv", index=False)

    targets = [t.strip() for t in args.targets.split(",") if t.strip() and t.strip() in joined.columns]
    correlation_pairs = [("mean_D_struct", "Q_bar"), ("peak_D_struct", "Q_bar"), ("mean_R_rob", "Q_bar")]
    for target in targets:
        correlation_pairs.extend(
            [
                ("mean_D_struct", target),
                ("Q_bar", target),
                ("mean_R_rob", target),
            ]
        )

    corr_rows = []
    for x, y in correlation_pairs:
        if x in joined.columns and y in joined.columns:
            corr_rows.append(
                {
                    "x": x,
                    "y": y,
                    "n": int(joined[[x, y]].dropna().shape[0]),
                    "spearman": safe_corr(joined[x], joined[y], "spearman"),
                    "pearson": safe_corr(joined[x], joined[y], "pearson"),
                }
            )
    corr = pd.DataFrame(corr_rows)
    corr.to_csv(args.output_dir / "cross_framework_correlations.csv", index=False)

    boot = bootstrap_correlations(joined, correlation_pairs, args.bootstrap, args.seed)
    boot.to_csv(args.output_dir / "cross_framework_bootstrap.csv", index=False)

    base_features = [
        "vobs_max",
        "median_gas_fraction_proxy",
        "median_bulge_fraction_proxy",
        "n_points",
        "outer_residual_fraction_proxy",
    ]
    model_features = {
        "baseline_baryonic_proxy": base_features,
        "baseline_plus_Dstruct": base_features + ["mean_D_struct"],
        "baseline_plus_Qbar": base_features + ["Q_bar"],
        "baseline_plus_both": base_features + ["mean_D_struct", "Q_bar"],
        "Dstruct_only": ["mean_D_struct"],
        "Qbar_only": ["Q_bar"],
        "Dstruct_Qbar_only": ["mean_D_struct", "Q_bar"],
    }
    cv_rows = []
    for target in targets:
        for model, features in model_features.items():
            available = [f for f in features if f in joined.columns and f != target]
            if not available:
                continue
            metrics = cv_ols(joined, available, target, args.folds, args.seed)
            cv_rows.append({"target": target, "model": model, "features": ",".join(available), **metrics})
    cv = pd.DataFrame(cv_rows)
    cv.to_csv(args.output_dir / "cross_framework_cv_results.csv", index=False)

    plot_scatter(joined, args.output_dir)
    plot_cv(cv, args.output_dir)
    plot_bootstrap(boot, args.output_dir)

    summary = {
        "status": "complete",
        "n_maat_galaxies": int(maat["galaxy"].nunique()),
        "n_qbar_galaxies": int(qbar["galaxy"].nunique()),
        "n_joined_galaxies": int(joined["galaxy"].nunique()),
        "targets": targets,
        "spearman_Dstruct_Qbar": safe_corr(joined["mean_D_struct"], joined["Q_bar"], "spearman"),
        "spearman_peakDstruct_Qbar": safe_corr(joined["peak_D_struct"], joined["Q_bar"], "spearman"),
        "license_note": (
            "MAAT/SPARC-derived inputs follow Paper 65 attribution. HGD-GSR Q_bar "
            "values are not redistributed unless provided with explicit permission "
            "and citation/licence terms by the HGD-GSR author."
        ),
    }
    (args.output_dir / "cross_framework_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
