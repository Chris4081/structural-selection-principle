#!/usr/bin/env python3
"""
External Q_bar comparison for the SPARC MAAT structural residual pipeline.

This script compares the externally supplied HGD-GSR Q_bar values against the
already generated Paper-65 SPARC galaxy-level diagnostics. It is deliberately
standalone and non-destructive: it reads existing Paper-65 outputs, joins the
external Q_bar table by galaxy name, and writes a separate output folder.

Run:
    python3 sparc_qbar_external_comparison.py
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

_CACHE = Path(tempfile.gettempdir()) / "maat_sparc_qbar_matplotlib_cache"
_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

EPS = 1.0e-12


def safe_corr(x: Iterable[float], y: Iterable[float], method: str = "spearman") -> float:
    df = pd.DataFrame({"x": list(x), "y": list(y)}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(df) < 4 or df["x"].nunique() < 2 or df["y"].nunique() < 2:
        return float("nan")
    return float(df["x"].corr(df["y"], method=method))


def rank_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").rank(method="average")


def residualize(y: pd.Series, controls: pd.DataFrame) -> pd.Series:
    data = pd.concat([rank_series(y).rename("y"), controls.apply(rank_series)], axis=1)
    data = data.replace([np.inf, -np.inf], np.nan).dropna()
    if len(data) < 8:
        return pd.Series(dtype=float)
    yy = data["y"].to_numpy(dtype=float)
    xx = data.drop(columns=["y"]).to_numpy(dtype=float)
    xx = np.column_stack([np.ones(len(xx)), xx])
    beta, *_ = np.linalg.lstsq(xx, yy, rcond=None)
    resid = yy - xx @ beta
    return pd.Series(resid, index=data.index)


def partial_rank_corr(df: pd.DataFrame, x_col: str, y_col: str, control_cols: list[str]) -> float:
    cols = [x_col, y_col, *control_cols]
    work = df[cols].replace([np.inf, -np.inf], np.nan).dropna()
    if len(work) < 12:
        return float("nan")
    x_resid = residualize(work[x_col], work[control_cols])
    y_resid = residualize(work[y_col], work[control_cols])
    common = x_resid.index.intersection(y_resid.index)
    if len(common) < 8:
        return float("nan")
    return safe_corr(x_resid.loc[common], y_resid.loc[common], method="pearson")


def permutation_pvalue(x: pd.Series, y: pd.Series, observed: float, n_perm: int, seed: int) -> float:
    if not np.isfinite(observed):
        return float("nan")
    work = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(work) < 8:
        return float("nan")
    rng = np.random.default_rng(seed)
    xr = rank_series(work["x"]).to_numpy(dtype=float)
    yr = rank_series(work["y"]).to_numpy(dtype=float)
    count = 0
    for _ in range(n_perm):
        val = np.corrcoef(rng.permutation(xr), yr)[0, 1]
        if abs(val) >= abs(observed):
            count += 1
    return float((count + 1) / (n_perm + 1))


def load_inputs(summary_path: Path, qbar_path: Path, radial_corr_path: Path | None) -> tuple[pd.DataFrame, dict]:
    summary = pd.read_csv(summary_path)
    qbar = pd.read_csv(qbar_path)
    qbar["galaxy"] = qbar["galaxy"].astype(str)
    qbar["Q_bar"] = pd.to_numeric(qbar["Q_bar"], errors="coerce")
    merged = summary.merge(qbar, on="galaxy", how="inner")
    unmatched_summary = sorted(set(summary["galaxy"].astype(str)) - set(qbar["galaxy"].astype(str)))
    unmatched_qbar = sorted(set(qbar["galaxy"].astype(str)) - set(summary["galaxy"].astype(str)))

    if radial_corr_path and radial_corr_path.exists():
        radial = pd.read_csv(radial_corr_path)
        keep = [
            "galaxy",
            "rho_Dstruct_nfw_radius",
            "rho_Dstruct_rar_radius",
            "rho_Dstruct_radius",
        ]
        merged = merged.merge(radial[[c for c in keep if c in radial.columns]], on="galaxy", how="left")

    meta = {
        "summary_rows": int(len(summary)),
        "qbar_rows": int(len(qbar)),
        "matched_rows": int(len(merged)),
        "unmatched_summary": unmatched_summary,
        "unmatched_qbar": unmatched_qbar,
    }
    return merged, meta


def compute_correlations(df: pd.DataFrame, n_perm: int, seed: int) -> pd.DataFrame:
    targets = [
        "mean_D_struct",
        "peak_D_struct",
        "mean_R_rob",
        "mean_V",
        "mean_residual_fraction",
        "outer_residual_fraction",
        "nfw_like_chi2",
        "nfw_like_rmse_v2",
        "rar_mean_sigma_residual",
        "vobs_max",
        "median_gas_fraction_proxy",
        "median_bulge_fraction_proxy",
        "outer_residual_fraction_proxy",
        "rho_Dstruct_nfw_radius",
        "rho_Dstruct_rar_radius",
        "rho_Dstruct_radius",
    ]
    rows = []
    for target in targets:
        if target not in df.columns:
            continue
        pearson = safe_corr(df["Q_bar"], df[target], "pearson")
        spearman = safe_corr(df["Q_bar"], df[target], "spearman")
        p_perm = permutation_pvalue(df["Q_bar"], df[target], spearman, n_perm=n_perm, seed=seed)
        rows.append(
            {
                "target": target,
                "n": int(pd.DataFrame({"x": df["Q_bar"], "y": df[target]}).dropna().shape[0]),
                "pearson_Qbar": pearson,
                "spearman_Qbar": spearman,
                "spearman_perm_p_two_sided": p_perm,
            }
        )
    return pd.DataFrame(rows).sort_values("spearman_Qbar", key=lambda s: s.abs(), ascending=False)


def compute_partial_checks(df: pd.DataFrame) -> pd.DataFrame:
    control_sets = {
        "vobs": ["vobs_max"],
        "gas_bulge": ["median_gas_fraction_proxy", "median_bulge_fraction_proxy"],
        "vobs_gas_bulge": ["vobs_max", "median_gas_fraction_proxy", "median_bulge_fraction_proxy"],
    }
    targets = [
        "mean_D_struct",
        "peak_D_struct",
        "mean_residual_fraction",
        "outer_residual_fraction",
        "nfw_like_chi2",
        "rar_mean_sigma_residual",
    ]
    rows = []
    for target in targets:
        if target not in df.columns:
            continue
        for label, controls in control_sets.items():
            if not all(c in df.columns for c in controls):
                continue
            rows.append(
                {
                    "target": target,
                    "controls": label,
                    "partial_rank_corr_Qbar": partial_rank_corr(df, "Q_bar", target, controls),
                }
            )
    return pd.DataFrame(rows)


def compute_group_summaries(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    morph_rows = []
    if "morphology_proxy" in df.columns:
        for morph, g in df.groupby("morphology_proxy"):
            morph_rows.append(
                {
                    "morphology_proxy": morph,
                    "n": int(len(g)),
                    "median_Q_bar": float(g["Q_bar"].median()),
                    "mean_Q_bar": float(g["Q_bar"].mean()),
                    "mean_D_struct": float(g["mean_D_struct"].mean()),
                    "outer_residual_fraction": float(g["outer_residual_fraction"].mean()),
                    "nfw_like_chi2": float(g["nfw_like_chi2"].mean()),
                }
            )
    morph_summary = pd.DataFrame(morph_rows).sort_values("median_Q_bar") if morph_rows else pd.DataFrame()

    bins = pd.qcut(df["Q_bar"], q=4, labels=["Q1_low", "Q2", "Q3", "Q4_high"], duplicates="drop")
    binned = df.assign(Q_bar_quartile=bins)
    quartile_summary = (
        binned.groupby("Q_bar_quartile", observed=True)
        .agg(
            n=("galaxy", "count"),
            Q_bar_min=("Q_bar", "min"),
            Q_bar_max=("Q_bar", "max"),
            Q_bar_median=("Q_bar", "median"),
            mean_D_struct=("mean_D_struct", "mean"),
            peak_D_struct=("peak_D_struct", "mean"),
            outer_residual_fraction=("outer_residual_fraction", "mean"),
            nfw_like_chi2=("nfw_like_chi2", "mean"),
            rar_mean_sigma_residual=("rar_mean_sigma_residual", "mean"),
        )
        .reset_index()
    )
    return morph_summary, quartile_summary


def plot_scatter(df: pd.DataFrame, x: str, y: str, out: Path, title: str, ylabel: str) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    if "morphology_proxy" in df.columns:
        for morph, g in df.groupby("morphology_proxy"):
            ax.scatter(g[x], g[y], label=morph, alpha=0.75, s=34)
        ax.legend(fontsize=7, frameon=False)
    else:
        ax.scatter(df[x], df[y], alpha=0.75, s=34)
    rho = safe_corr(df[x], df[y], "spearman")
    ax.set_title(f"{title}\nSpearman rho={rho:.3f}")
    ax.set_xlabel("External Q_bar")
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_correlation_bars(corr: pd.DataFrame, out: Path) -> None:
    top = corr.copy()
    top = top[np.isfinite(top["spearman_Qbar"])].head(12)
    fig, ax = plt.subplots(figsize=(8.0, 5.8))
    labels = top["target"].tolist()
    vals = top["spearman_Qbar"].to_numpy(dtype=float)
    y = np.arange(len(labels))
    ax.barh(y, vals, color=["#0f766e" if v >= 0 else "#9f1239" for v in vals])
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0.0, color="black", linewidth=0.8)
    ax.set_xlabel("Spearman rho with Q_bar")
    ax.set_title("External Q_bar: strongest rank associations")
    ax.grid(axis="x", alpha=0.25)
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_morphology(df: pd.DataFrame, out: Path) -> None:
    if "morphology_proxy" not in df.columns:
        return
    order = df.groupby("morphology_proxy")["Q_bar"].median().sort_values().index.tolist()
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    data = [df.loc[df["morphology_proxy"] == m, "Q_bar"].dropna().to_numpy() for m in order]
    ax.boxplot(data, tick_labels=order, showfliers=False)
    ax.set_ylabel("External Q_bar")
    ax.set_title("Q_bar distribution by Paper-65 morphology proxy")
    ax.tick_params(axis="x", labelrotation=25)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def write_readme(outdir: Path, summary: dict, corr: pd.DataFrame) -> None:
    strongest = corr.head(6)
    lines = [
        "# External Q_bar comparison",
        "",
        "This folder contains a neutral comparison between the externally supplied HGD-GSR Q_bar table and the existing Paper-65 SPARC galaxy-level diagnostics.",
        "",
        "## Coverage",
        "",
        f"- Paper-65 summary rows: `{summary['coverage']['summary_rows']}`",
        f"- Q_bar rows: `{summary['coverage']['qbar_rows']}`",
        f"- Matched galaxies: `{summary['coverage']['matched_rows']}`",
        "",
        "## Strongest rank associations",
        "",
        "| target | Spearman rho | permutation p |",
        "|---|---:|---:|",
    ]
    for _, row in strongest.iterrows():
        lines.append(
            f"| `{row['target']}` | {row['spearman_Qbar']:.4f} | {row['spearman_perm_p_two_sided']:.4f} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation note",
            "",
            "These results are diagnostic correlations only. They do not establish a physical model or causal relation. Q_bar is treated as an external input and compared against existing MAAT/SPARC structural summaries without changing the Paper-61 or Paper-65 pipelines.",
            "",
            "No endorsement by SPARC, Zenodo, VizieR, CDS, the original data providers, or the external Q_bar author is implied.",
        ]
    )
    outdir.joinpath("README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, default=Path("outputs_paper65/paper65_galaxy_summary_with_proxy_morphology.csv"))
    parser.add_argument("--radial-corr", type=Path, default=Path("outputs_paper65/paper65_galaxywise_radial_correlations.csv"))
    parser.add_argument("--qbar", type=Path, default=Path("external_qbar_alhawarat.csv"))
    parser.add_argument("--output", type=Path, default=Path("outputs_qbar_external"))
    parser.add_argument("--n-perm", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=65066)
    args = parser.parse_args()

    outdir = args.output
    outdir.mkdir(parents=True, exist_ok=True)

    df, coverage = load_inputs(args.summary, args.qbar, args.radial_corr)
    if df.empty:
        raise SystemExit("No galaxy overlap between Q_bar table and Paper-65 summary.")

    df["log_Q_bar"] = np.log(np.maximum(df["Q_bar"], EPS))
    df.to_csv(outdir / "qbar_joined_summary.csv", index=False)

    corr = compute_correlations(df, n_perm=args.n_perm, seed=args.seed)
    corr.to_csv(outdir / "qbar_correlations.csv", index=False)

    partial = compute_partial_checks(df)
    partial.to_csv(outdir / "qbar_partial_rank_checks.csv", index=False)

    morph_summary, quartile_summary = compute_group_summaries(df)
    morph_summary.to_csv(outdir / "qbar_morphology_proxy_summary.csv", index=False)
    quartile_summary.to_csv(outdir / "qbar_quartile_summary.csv", index=False)

    plot_scatter(df, "Q_bar", "mean_D_struct", outdir / "fig1_qbar_vs_mean_dstruct.png", "Q_bar vs mean D_struct", "mean D_struct")
    plot_scatter(df, "Q_bar", "outer_residual_fraction", outdir / "fig2_qbar_vs_outer_residual.png", "Q_bar vs outer residual fraction", "outer residual fraction")
    plot_correlation_bars(corr, outdir / "fig3_qbar_rank_associations.png")
    plot_morphology(df, outdir / "fig4_qbar_by_morphology_proxy.png")

    summary = {
        "coverage": coverage,
        "qbar": {
            "min": float(df["Q_bar"].min()),
            "median": float(df["Q_bar"].median()),
            "mean": float(df["Q_bar"].mean()),
            "max": float(df["Q_bar"].max()),
        },
        "strongest_rank_associations": corr.head(10).to_dict(orient="records"),
        "partial_rank_checks": partial.to_dict(orient="records"),
        "notes": [
            "Q_bar is treated as a neutral external input.",
            "Correlations are diagnostic and not causal claims.",
            "Existing Paper-61 and Paper-65 outputs are not modified.",
        ],
    }
    (outdir / "qbar_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    write_readme(outdir, summary, corr)

    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
