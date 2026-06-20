#!/usr/bin/env python3
"""
Paper 65: SPARC II robustness checks for structural dark-matter residuals.

This script extends Paper 61 with:

1. shuffled-radius nulls that destroy radial pairing within each galaxy,
2. morphology/proxy splits,
3. an optional published halo-catalogue cross-check.

The script never requires external network access.  It consumes the prepared
SPARC-derived Paper-61 outputs.  If a published halo-catalogue CSV is provided
locally, it is joined by galaxy name; otherwise the halo-catalogue check is
reported as not available and the internal NFW-like channel is used only as a
lightweight comparison proxy.

Run:
    python3 sparc_ii_paper65.py
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

_CACHE = Path(tempfile.gettempdir()) / "maat_sparc_paper65_matplotlib_cache"
_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

EPS = 1.0e-12

GALAXY_ALIASES = ["galaxy", "Galaxy", "name", "Name", "ID", "SPARC", "galaxy_name"]
HALO_MODEL_ALIASES = ["model", "Model", "halo_model", "profile", "Profile"]
CHI2_ALIASES = [
    "chi2",
    "chi2_red",
    "red_chi2",
    "chisq",
    "chisq_red",
    "Chi2",
    "Chi2_red",
    "best_chi2",
    "reduced_chi2",
]


def safe_corr(x: Iterable[float], y: Iterable[float], method: str = "spearman") -> float:
    s = pd.DataFrame({"x": list(x), "y": list(y)}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(s) < 4 or s["x"].nunique() < 2 or s["y"].nunique() < 2:
        return float("nan")
    return float(s["x"].corr(s["y"], method=method))


def zscore_within_group(df: pd.DataFrame, column: str, group_col: str = "galaxy") -> pd.Series:
    def _z(s: pd.Series) -> pd.Series:
        vals = s.astype(float)
        sd = vals.std(ddof=0)
        if not np.isfinite(sd) or sd < EPS:
            return pd.Series(np.zeros(len(vals)), index=s.index)
        return (vals - vals.mean()) / sd

    return df.groupby(group_col, group_keys=False)[column].apply(_z)


def find_first(columns: Iterable[str], aliases: list[str]) -> str | None:
    cols = set(columns)
    for alias in aliases:
        if alias in cols:
            return alias
    return None


def classify_proxy_morphology(profile: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for galaxy, g in profile.groupby("galaxy"):
        vobs_max = float(g["v_obs_km_s"].max())
        vgas2 = np.square(g["v_gas_km_s"].astype(float).to_numpy())
        vdisk2 = np.square(g["v_disk_km_s"].astype(float).to_numpy())
        vbul2 = np.square(g["v_bulge_km_s"].astype(float).to_numpy())
        denom = np.maximum(vgas2 + 0.5 * vdisk2 + 0.7 * vbul2, EPS)
        gas_frac = float(np.nanmedian(vgas2 / denom))
        bulge_frac = float(np.nanmedian(0.7 * vbul2 / denom))
        residual_outer = float(
            g.loc[g["r_kpc"] >= g["r_kpc"].median(), "residual_fraction"].mean()
        )
        if bulge_frac > 0.12:
            proxy = "bulge_proxy"
        elif vobs_max < 85.0:
            proxy = "dwarf_proxy"
        elif gas_frac > 0.42:
            proxy = "gas_rich_late_proxy"
        elif residual_outer > 0.72:
            proxy = "dm_dominated_disk_proxy"
        else:
            proxy = "disk_proxy"
        rows.append(
            {
                "galaxy": galaxy,
                "morphology_proxy": proxy,
                "vobs_max": vobs_max,
                "median_gas_fraction_proxy": gas_frac,
                "median_bulge_fraction_proxy": bulge_frac,
                "outer_residual_fraction_proxy": residual_outer,
            }
        )
    proxy_df = pd.DataFrame(rows)
    out = summary.merge(proxy_df, on="galaxy", how="left")
    return out


def load_morphology_catalog(path: Path) -> pd.DataFrame | None:
    if not path or not path.exists():
        return None
    df = pd.read_csv(path)
    gcol = find_first(df.columns, GALAXY_ALIASES)
    tcol = find_first(df.columns, ["type", "Type", "morphology", "Morphology", "HubbleType", "T", "Ttype"])
    if gcol is None or tcol is None:
        raise ValueError("Morphology catalog must contain galaxy/name and type/morphology/T columns.")
    out = df[[gcol, tcol]].copy()
    out.columns = ["galaxy", "morphology_catalogue"]
    return out


def load_halo_catalog(path: Path) -> tuple[pd.DataFrame | None, dict]:
    if not path or not path.exists():
        return None, {
            "halo_catalogue_status": "not_available",
            "halo_catalogue_note": "No local published halo-catalogue CSV was provided.",
        }
    df = pd.read_csv(path)
    gcol = find_first(df.columns, GALAXY_ALIASES)
    mcol = find_first(df.columns, HALO_MODEL_ALIASES)
    ccol = find_first(df.columns, CHI2_ALIASES)
    if gcol is None or ccol is None:
        raise ValueError("Halo catalog must contain a galaxy/name column and a chi2-like column.")
    out = pd.DataFrame({"galaxy": df[gcol].astype(str), "halo_catalogue_chi2": pd.to_numeric(df[ccol], errors="coerce")})
    if mcol is not None:
        out["halo_catalogue_model"] = df[mcol].astype(str)
    else:
        out["halo_catalogue_model"] = "unspecified"
    out = out.dropna(subset=["galaxy", "halo_catalogue_chi2"])
    # Multi-model catalogues are collapsed to the best available fit per galaxy.
    best = out.sort_values("halo_catalogue_chi2").groupby("galaxy", as_index=False).first()
    meta = {
        "halo_catalogue_status": "loaded",
        "halo_catalogue_path": str(path),
        "halo_catalogue_rows": int(len(out)),
        "halo_catalogue_matched_galaxies_before_join": int(best["galaxy"].nunique()),
        "halo_catalogue_chi2_column": ccol,
        "halo_catalogue_model_column": mcol or "not provided",
    }
    return best, meta


def pooled_radial_correlations(profile: pd.DataFrame, d_column: str = "D_struct") -> dict:
    df = profile.copy()
    df["D_z"] = zscore_within_group(df, d_column)
    df["nfw_z"] = zscore_within_group(df, "nfw_fractional_residual")
    df["rar_z"] = zscore_within_group(df, "rar_sigma_residual")
    df["radius_z"] = zscore_within_group(df, "r_kpc")
    return {
        "pooled_radius_spearman_Dstruct_vs_nfw": safe_corr(df["D_z"], df["nfw_z"], "spearman"),
        "pooled_radius_spearman_Dstruct_vs_rar": safe_corr(df["D_z"], df["rar_z"], "spearman"),
        "pooled_radius_spearman_Dstruct_vs_radius": safe_corr(df["D_z"], df["radius_z"], "spearman"),
    }


def galaxywise_radial_correlations(profile: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for galaxy, g in profile.groupby("galaxy"):
        if len(g) < 5:
            continue
        rows.append(
            {
                "galaxy": galaxy,
                "n_points": int(len(g)),
                "rho_Dstruct_nfw_radius": safe_corr(g["D_struct"], g["nfw_fractional_residual"], "spearman"),
                "rho_Dstruct_rar_radius": safe_corr(g["D_struct"], g["rar_sigma_residual"], "spearman"),
                "rho_Dstruct_radius": safe_corr(g["D_struct"], g["r_kpc"], "spearman"),
            }
        )
    return pd.DataFrame(rows)


def radius_shuffle_null(profile: pd.DataFrame, seed: int = 65065, n_null: int = 500) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(seed)
    real = pooled_radial_correlations(profile)
    null_rows = []
    base = profile.copy()
    for i in range(n_null):
        shuffled = base.copy()
        shuffled["D_struct"] = shuffled.groupby("galaxy", group_keys=False)["D_struct"].transform(
            lambda s: pd.Series(rng.permutation(s.to_numpy()), index=s.index)
        )
        vals = pooled_radial_correlations(shuffled)
        vals["null_id"] = i
        null_rows.append(vals)
    null_df = pd.DataFrame(null_rows)
    summary: dict[str, float | int] = {"n_null": int(n_null), **real}
    for key in [
        "pooled_radius_spearman_Dstruct_vs_nfw",
        "pooled_radius_spearman_Dstruct_vs_rar",
        "pooled_radius_spearman_Dstruct_vs_radius",
    ]:
        vals = null_df[key].dropna().to_numpy()
        summary[f"{key}_null_mean"] = float(np.mean(vals))
        summary[f"{key}_null_std"] = float(np.std(vals))
        summary[f"{key}_null_abs95"] = float(np.quantile(np.abs(vals), 0.95))
        summary[f"{key}_exceeds_abs95"] = bool(abs(summary[key]) > summary[f"{key}_null_abs95"])
    return null_df, summary


def morphology_splits(profile: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    merged = profile.merge(summary[["galaxy", "morphology_proxy"]], on="galaxy", how="left")
    rows = []
    for morph, g in merged.groupby("morphology_proxy"):
        if g["galaxy"].nunique() < 3:
            continue
        corr = pooled_radial_correlations(g)
        rows.append(
            {
                "morphology_proxy": morph,
                "n_galaxies": int(g["galaxy"].nunique()),
                "n_rows": int(len(g)),
                "mean_D_struct": float(g["D_struct"].mean()),
                "mean_residual_fraction": float(g["residual_fraction"].mean()),
                "median_R_rob": float(g["R_rob"].median()),
                **corr,
            }
        )
    return pd.DataFrame(rows).sort_values("morphology_proxy")


def halo_catalogue_crosscheck(summary: pd.DataFrame, halo: pd.DataFrame | None) -> tuple[pd.DataFrame, dict]:
    if halo is None:
        return pd.DataFrame(), {
            "halo_catalogue_matched_galaxies": 0,
            "spearman_Dstruct_vs_halo_catalogue_chi2": float("nan"),
            "spearman_Dstruct_vs_internal_nfw_chi2": safe_corr(summary["mean_D_struct"], summary["nfw_like_chi2"], "spearman"),
            "halo_catalogue_crosscheck_status": "not_run_no_local_catalogue",
        }
    joined = summary.merge(halo, on="galaxy", how="inner")
    return joined, {
        "halo_catalogue_matched_galaxies": int(joined["galaxy"].nunique()),
        "spearman_Dstruct_vs_halo_catalogue_chi2": safe_corr(joined["mean_D_struct"], joined["halo_catalogue_chi2"], "spearman"),
        "spearman_Dstruct_vs_internal_nfw_chi2_on_matched": safe_corr(joined["mean_D_struct"], joined["nfw_like_chi2"], "spearman"),
        "halo_catalogue_crosscheck_status": "run" if len(joined) >= 5 else "insufficient_matches",
    }


def plot_radius_null(null_df: pd.DataFrame, summary: dict, outdir: Path) -> None:
    labels = ["NFW-like", "RAR", "radius"]
    keys = [
        "pooled_radius_spearman_Dstruct_vs_nfw",
        "pooled_radius_spearman_Dstruct_vs_rar",
        "pooled_radius_spearman_Dstruct_vs_radius",
    ]
    real = [abs(summary[k]) for k in keys]
    null95 = [summary[f"{k}_null_abs95"] for k in keys]
    x = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.bar(x - width / 2, real, width, color="#006d77", label="real radial pairing")
    ax.bar(x + width / 2, null95, width, color="#adb5bd", label="radius-shuffle |rho| 95%")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("|within-galaxy pooled Spearman|")
    ax.set_title("Shuffled-radius null test")
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "fig1_shuffled_radius_null.png", dpi=180)
    plt.close(fig)


def plot_morphology_splits(splits: pd.DataFrame, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 5.0))
    order = splits.sort_values("mean_D_struct", ascending=False)
    ax.bar(order["morphology_proxy"], order["mean_D_struct"], color="#2a9d8f")
    ax.set_ylabel("mean D_struct")
    ax.set_title("Structural residual by baryonic morphology proxy")
    ax.tick_params(axis="x", rotation=25)
    for idx, row in enumerate(order.itertuples(index=False)):
        ax.text(idx, row.mean_D_struct, f"n={row.n_galaxies}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "fig2_morphology_proxy_splits.png", dpi=180)
    plt.close(fig)


def plot_halo_crosscheck(summary: pd.DataFrame, joined: pd.DataFrame, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    if not joined.empty and "halo_catalogue_chi2" in joined.columns:
        ax.scatter(joined["mean_D_struct"], joined["halo_catalogue_chi2"], s=45, alpha=0.78, color="#5e548e")
        ax.set_ylabel("published halo-catalogue chi2")
        ax.set_title("Optional halo-catalogue cross-check")
    else:
        ax.scatter(summary["mean_D_struct"], summary["nfw_like_chi2"], s=45, alpha=0.72, color="#006d77")
        ax.text(
            0.03,
            0.95,
            "No local published halo catalogue supplied;\nshowing internal NFW-like proxy.",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
        )
        ax.set_ylabel("internal NFW-like chi2")
        ax.set_title("Halo-catalogue slot and internal proxy")
    ax.set_xlabel("mean D_struct")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "fig3_halo_catalogue_crosscheck.png", dpi=180)
    plt.close(fig)


def plot_galaxywise_correlations(gal_corr: pd.DataFrame, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    vals = [
        gal_corr["rho_Dstruct_nfw_radius"].dropna().to_numpy(),
        gal_corr["rho_Dstruct_rar_radius"].dropna().to_numpy(),
        gal_corr["rho_Dstruct_radius"].dropna().to_numpy(),
    ]
    ax.boxplot(vals, tick_labels=["NFW-like", "RAR", "radius"], showmeans=True)
    ax.axhline(0.0, color="#333333", lw=1)
    ax.set_ylabel("within-galaxy Spearman")
    ax.set_title("Galaxywise radial correlation distribution")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "fig4_galaxywise_radial_correlations.png", dpi=180)
    plt.close(fig)


def write_readme(outdir: Path, payload: dict) -> None:
    text = f"""# Paper 65: SPARC II

SPARC II extends Paper 61 with shuffled-radius nulls, morphology/proxy splits,
and an optional published halo-catalogue cross-check.

Headline radius-null result:

- D_struct vs NFW-like pooled radial Spearman: `{payload['radius_null_summary']['pooled_radius_spearman_Dstruct_vs_nfw']:.4f}`
- D_struct vs RAR pooled radial Spearman: `{payload['radius_null_summary']['pooled_radius_spearman_Dstruct_vs_rar']:.4f}`
- NFW-like exceeds radius-shuffle |rho|95: `{payload['radius_null_summary']['pooled_radius_spearman_Dstruct_vs_nfw_exceeds_abs95']}`
- RAR exceeds radius-shuffle |rho|95: `{payload['radius_null_summary']['pooled_radius_spearman_Dstruct_vs_rar_exceeds_abs95']}`

Halo-catalogue status: `{payload['halo_catalogue']['halo_catalogue_crosscheck_status']}`.

Outputs:

- `paper65_radius_shuffle_nulls.csv`
- `paper65_radius_null_summary.json`
- `paper65_galaxywise_radial_correlations.csv`
- `paper65_morphology_proxy_splits.csv`
- `paper65_halo_catalogue_join.csv` if a local halo catalogue is supplied
- `paper65_summary.json`
- four PNG figures

No endorsement by SPARC, Li et al., Zenodo, VizieR, CDS, or any original data
provider is implied.
"""
    (outdir / "README.md").write_text(text, encoding="utf-8")


def run(args: argparse.Namespace) -> dict:
    profile = pd.read_csv(args.radius_diagnostics)
    summary = pd.read_csv(args.galaxy_summary)
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    summary = classify_proxy_morphology(profile, summary)
    morph_cat = load_morphology_catalog(Path(args.morphology_catalog)) if args.morphology_catalog else None
    if morph_cat is not None:
        summary = summary.merge(morph_cat, on="galaxy", how="left")

    halo, halo_meta = load_halo_catalog(Path(args.halo_catalog)) if args.halo_catalog else (None, {
        "halo_catalogue_status": "not_available",
        "halo_catalogue_note": "No local published halo-catalogue CSV was provided.",
    })

    null_df, null_summary = radius_shuffle_null(profile, seed=args.seed, n_null=args.n_null)
    gal_corr = galaxywise_radial_correlations(profile)
    splits = morphology_splits(profile, summary)
    halo_join, halo_summary = halo_catalogue_crosscheck(summary, halo)

    summary.to_csv(outdir / "paper65_galaxy_summary_with_proxy_morphology.csv", index=False)
    null_df.to_csv(outdir / "paper65_radius_shuffle_nulls.csv", index=False)
    gal_corr.to_csv(outdir / "paper65_galaxywise_radial_correlations.csv", index=False)
    splits.to_csv(outdir / "paper65_morphology_proxy_splits.csv", index=False)
    if not halo_join.empty:
        halo_join.to_csv(outdir / "paper65_halo_catalogue_join.csv", index=False)

    plot_radius_null(null_df, null_summary, outdir)
    plot_morphology_splits(splits, outdir)
    plot_halo_crosscheck(summary, halo_join, outdir)
    plot_galaxywise_correlations(gal_corr, outdir)

    payload = {
        "paper": 65,
        "title": "SPARC II: Shuffled-Radius Nulls, Halo-Catalogue Cross-Check, and Morphology Splits",
        "seed": args.seed,
        "n_null": args.n_null,
        "n_galaxies": int(summary["galaxy"].nunique()),
        "n_radius_rows": int(len(profile)),
        "radius_null_summary": null_summary,
        "galaxywise_radial_correlations": {
            "median_rho_Dstruct_nfw_radius": float(gal_corr["rho_Dstruct_nfw_radius"].median()),
            "median_rho_Dstruct_rar_radius": float(gal_corr["rho_Dstruct_rar_radius"].median()),
            "median_rho_Dstruct_radius": float(gal_corr["rho_Dstruct_radius"].median()),
            "n_galaxies_with_ge5_points": int(len(gal_corr)),
        },
        "morphology_proxy_counts": summary["morphology_proxy"].value_counts().to_dict(),
        "morphology_catalogue_status": "loaded" if morph_cat is not None else "not_available_proxy_used",
        "halo_catalogue_metadata": halo_meta,
        "halo_catalogue": halo_summary,
        "license_note": (
            "SPARC-derived inputs follow the Paper 61 SPARC attribution note. "
            "If a published halo catalogue is supplied locally, cite the original "
            "catalogue source, e.g. Li et al. for SPARC halo fits, and do not imply "
            "endorsement by catalogue authors."
        ),
    }
    (outdir / "paper65_radius_null_summary.json").write_text(json.dumps(null_summary, indent=2), encoding="utf-8")
    (outdir / "paper65_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_readme(outdir, payload)
    print(json.dumps(payload, indent=2))
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Paper 65 SPARC II robustness checks")
    parser.add_argument(
        "--radius-diagnostics",
        default="../sparc_structural_residual_pilot_paper61/outputs_sparc_real/paper61_radius_diagnostics.csv",
    )
    parser.add_argument(
        "--galaxy-summary",
        default="../sparc_structural_residual_pilot_paper61/outputs_sparc_real/paper61_galaxy_summary.csv",
    )
    parser.add_argument("--output", default="outputs_paper65")
    parser.add_argument("--halo-catalog", default="", help="Optional local published halo-catalogue CSV")
    parser.add_argument("--morphology-catalog", default="", help="Optional local morphology catalogue CSV")
    parser.add_argument("--seed", type=int, default=65065)
    parser.add_argument("--n-null", type=int, default=500)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
