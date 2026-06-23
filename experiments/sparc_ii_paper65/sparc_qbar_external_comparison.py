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
MASS_PROXY_FACTOR = 2.325e5  # M_sun for v^2[km/s] * r[kpc] / G.


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


def derive_baryonic_size_controls(radius_path: Path | None) -> pd.DataFrame:
    """Return simple galaxy-level mass/size controls from SPARC radius rows.

    The mass term is a baryonic dynamical mass-scale proxy,
    M_bar(<R_max>) ~= V_bar(R_max)^2 R_max / G, not a photometric stellar mass.
    It is used only as a rank-control variable for the independence test.
    """
    if radius_path is None or not radius_path.exists():
        return pd.DataFrame()

    profile = pd.read_csv(radius_path)
    required = {"galaxy", "r_kpc", "v_baryon_km_s"}
    if not required.issubset(profile.columns):
        return pd.DataFrame()

    rows = []
    for galaxy, g in profile.groupby("galaxy"):
        work = g[["r_kpc", "v_baryon_km_s"]].replace([np.inf, -np.inf], np.nan).dropna()
        if work.empty:
            continue
        work = work.sort_values("r_kpc")
        rmax = float(work["r_kpc"].max())
        outer = work.iloc[-1]
        vbar_outer = float(max(outer["v_baryon_km_s"], 0.0))
        mbar_proxy = MASS_PROXY_FACTOR * vbar_outer * vbar_outer * max(rmax, EPS)
        rows.append(
            {
                "galaxy": str(galaxy),
                "r_max_kpc": rmax,
                "vbar_outer_km_s": vbar_outer,
                "mbar_proxy_msun": mbar_proxy,
                "log_r_max_kpc": float(np.log10(max(rmax, EPS))),
                "log_mbar_proxy_msun": float(np.log10(max(mbar_proxy, EPS))),
            }
        )
    return pd.DataFrame(rows)


def load_inputs(
    summary_path: Path,
    qbar_path: Path,
    radial_corr_path: Path | None,
    radius_diagnostics_path: Path | None,
) -> tuple[pd.DataFrame, dict]:
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

    controls = derive_baryonic_size_controls(radius_diagnostics_path)
    if not controls.empty:
        merged = merged.merge(controls, on="galaxy", how="left")

    meta = {
        "summary_rows": int(len(summary)),
        "qbar_rows": int(len(qbar)),
        "matched_rows": int(len(merged)),
        "mass_size_control_rows": int(len(controls)) if not controls.empty else 0,
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
        "baryonic_mass_size": ["log_mbar_proxy_msun", "log_r_max_kpc"],
        "baryonic_mass_size_vobs": ["log_mbar_proxy_msun", "log_r_max_kpc", "vobs_max"],
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


def zscore_train_test(x_train: np.ndarray, x_test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mu = np.nanmean(x_train, axis=0)
    sd = np.nanstd(x_train, axis=0)
    sd = np.where(sd < EPS, 1.0, sd)
    return (x_train - mu) / sd, (x_test - mu) / sd


def cv_linear_metrics(df: pd.DataFrame, target: str, features: list[str], seed: int, n_fold: int = 5) -> dict:
    work = df[[target, *features]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(work) < 20:
        return {
            "target": target,
            "model": "",
            "features": ",".join(features),
            "n": int(len(work)),
            "cv_r2": float("nan"),
            "cv_spearman": float("nan"),
        }
    rng = np.random.default_rng(seed)
    indices = np.arange(len(work))
    rng.shuffle(indices)
    folds = np.array_split(indices, n_fold)
    pred = np.full(len(work), np.nan)
    y = work[target].to_numpy(dtype=float)
    x = work[features].to_numpy(dtype=float)
    for fold in folds:
        train = np.setdiff1d(indices, fold)
        x_train, x_test = zscore_train_test(x[train], x[fold])
        design_train = np.column_stack([np.ones(len(train)), x_train])
        design_test = np.column_stack([np.ones(len(fold)), x_test])
        beta, *_ = np.linalg.lstsq(design_train, y[train], rcond=None)
        pred[fold] = design_test @ beta
    good = np.isfinite(pred) & np.isfinite(y)
    ss_res = float(np.sum((y[good] - pred[good]) ** 2))
    ss_tot = float(np.sum((y[good] - np.mean(y[good])) ** 2))
    r2 = 1.0 - ss_res / max(ss_tot, EPS)
    rho = safe_corr(pd.Series(pred[good]), pd.Series(y[good]), method="spearman")
    return {
        "target": target,
        "model": "",
        "features": ",".join(features),
        "n": int(np.sum(good)),
        "cv_r2": float(r2),
        "cv_spearman": float(rho),
    }


def compute_independence_model_comparison(df: pd.DataFrame, seed: int) -> pd.DataFrame:
    targets = [
        "nfw_like_chi2",
        "nfw_like_rmse_v2",
        "rar_mean_sigma_residual",
    ]
    models = {
        "Dstruct_only": ["mean_D_struct"],
        "Qbar_only": ["Q_bar"],
        "Dstruct_Qbar": ["mean_D_struct", "Q_bar"],
    }
    rows = []
    for target in targets:
        if target not in df.columns:
            continue
        for model, features in models.items():
            if not all(f in df.columns for f in features):
                continue
            row = cv_linear_metrics(df, target, features, seed=seed)
            row["model"] = model
            rows.append(row)
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


def plot_model_comparison(model_df: pd.DataFrame, out: Path) -> None:
    if model_df.empty:
        return
    targets = model_df["target"].drop_duplicates().tolist()
    models = ["Dstruct_only", "Qbar_only", "Dstruct_Qbar"]
    x = np.arange(len(targets))
    width = 0.24
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    for i, model in enumerate(models):
        vals = []
        for target in targets:
            row = model_df[(model_df["target"] == target) & (model_df["model"] == model)]
            vals.append(float(row["cv_spearman"].iloc[0]) if not row.empty else np.nan)
        ax.bar(x + (i - 1) * width, vals, width, label=model)
    ax.set_xticks(x)
    ax.set_xticklabels(targets, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("5-fold CV Spearman")
    ax.set_title("Q_bar independence test: mismatch prediction")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
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
            "## Independence-test highlights",
            "",
            "- Direct `Q_bar`--`mean_D_struct` Spearman correlation is reported in `qbar_correlations.csv`.",
            "- Partial rank checks include controls for `log_mbar_proxy_msun` and `log_r_max_kpc` in `qbar_partial_rank_checks.csv`.",
            "- Mismatch prediction comparisons for `D_struct only`, `Q_bar only`, and `D_struct + Q_bar` are reported in `qbar_independence_model_comparison.csv`.",
            "",
            "The baryonic mass control is a simple enclosed baryonic mass-scale proxy derived from `V_bar(R_max)^2 R_max / G`; it is not a full photometric stellar-mass estimate.",
            "",
            "## Interpretation note",
            "",
            "These results are diagnostic correlations only. They do not establish a physical model or causal relation. Q_bar is treated as an external input and compared against existing MAAT/SPARC structural summaries without changing the Paper-61 or Paper-65 pipelines.",
            "",
            "## Data attribution and license note",
            "",
            "SPARC-derived MAAT inputs follow the Paper-65 SPARC attribution and CC-BY-4.0 notes. SPARC rotation-curve data should be cited to Lelli, McGaugh, and Schombert and to the Zenodo-hosted SPARC record when reused.",
            "",
            "The `Q_bar` table is a collaborator-supplied derived HGD-GSR descriptor from Ali Alhawarat. It is not an original SPARC measurement and not a MAAT-derived quantity. Unless broader redistribution terms are separately supplied by the HGD-GSR author, treat the table as a neutral comparison input for this pilot and cite/acknowledge Ali Alhawarat and HGD-GSR when discussing or reusing the comparison.",
            "",
            "No endorsement by SPARC, Zenodo, VizieR, CDS, HGD-GSR, Ali Alhawarat, the original SPARC data providers, or any catalogue authors is implied.",
        ]
    )
    outdir.joinpath("README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, default=Path("outputs_paper65/paper65_galaxy_summary_with_proxy_morphology.csv"))
    parser.add_argument("--radial-corr", type=Path, default=Path("outputs_paper65/paper65_galaxywise_radial_correlations.csv"))
    parser.add_argument("--radius-diagnostics", type=Path, default=Path("../sparc_structural_residual_pilot_paper61/outputs_sparc_real/paper61_radius_diagnostics.csv"))
    parser.add_argument("--qbar", type=Path, default=Path("external_qbar_alhawarat.csv"))
    parser.add_argument("--output", type=Path, default=Path("outputs_qbar_external"))
    parser.add_argument("--n-perm", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=65066)
    args = parser.parse_args()

    outdir = args.output
    outdir.mkdir(parents=True, exist_ok=True)

    df, coverage = load_inputs(args.summary, args.qbar, args.radial_corr, args.radius_diagnostics)
    if df.empty:
        raise SystemExit("No galaxy overlap between Q_bar table and Paper-65 summary.")

    df["log_Q_bar"] = np.log(np.maximum(df["Q_bar"], EPS))
    df.to_csv(outdir / "qbar_joined_summary.csv", index=False)

    corr = compute_correlations(df, n_perm=args.n_perm, seed=args.seed)
    corr.to_csv(outdir / "qbar_correlations.csv", index=False)

    partial = compute_partial_checks(df)
    partial.to_csv(outdir / "qbar_partial_rank_checks.csv", index=False)

    model_comparison = compute_independence_model_comparison(df, seed=args.seed)
    model_comparison.to_csv(outdir / "qbar_independence_model_comparison.csv", index=False)

    morph_summary, quartile_summary = compute_group_summaries(df)
    morph_summary.to_csv(outdir / "qbar_morphology_proxy_summary.csv", index=False)
    quartile_summary.to_csv(outdir / "qbar_quartile_summary.csv", index=False)

    plot_scatter(df, "Q_bar", "mean_D_struct", outdir / "fig1_qbar_vs_mean_dstruct.png", "Q_bar vs mean D_struct", "mean D_struct")
    plot_scatter(df, "Q_bar", "outer_residual_fraction", outdir / "fig2_qbar_vs_outer_residual.png", "Q_bar vs outer residual fraction", "outer residual fraction")
    plot_correlation_bars(corr, outdir / "fig3_qbar_rank_associations.png")
    plot_morphology(df, outdir / "fig4_qbar_by_morphology_proxy.png")
    plot_model_comparison(model_comparison, outdir / "fig5_qbar_independence_model_comparison.png")

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
        "independence_model_comparison": model_comparison.to_dict(orient="records"),
        "notes": [
            "Q_bar is treated as a neutral external input.",
            "Correlations are diagnostic and not causal claims.",
            "Existing Paper-61 and Paper-65 outputs are not modified.",
            "Baryonic mass and galaxy size controls are approximate proxies derived from SPARC rotation-curve rows.",
        ],
    }
    (outdir / "qbar_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    write_readme(outdir, summary, corr)

    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
