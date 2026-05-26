#!/usr/bin/env python3
"""
Paper 61 - SPARC structural dark-matter residual pilot.

The script is SPARC-ready but does not redistribute SPARC data.  If a local
CSV is provided, it computes MAAT v1.6 structural residual diagnostics for real
rotation-curve rows.  If no CSV is present, it runs a clearly marked synthetic
schema test so the pipeline, plots, and output tables remain reproducible.

Expected CSV columns, with common aliases accepted:

    galaxy, r_kpc, v_obs_km_s, e_v_obs_km_s,
    v_gas_km_s, v_disk_km_s, v_bulge_km_s

Optional:

    upsilon_disk, upsilon_bulge

For SPARC-like mass models, the baryonic velocity is computed as

    v_bar^2 = Vgas |Vgas| + Upsilon_disk Vdisk |Vdisk|
              + Upsilon_bulge Vbulge |Vbulge|

using fixed default mass-to-light factors unless per-row values are supplied.
"""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path
from typing import Iterable

_cache_dir = Path(tempfile.gettempdir()) / "maat_sparc_paper61_matplotlib_cache"
_cache_dir.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_cache_dir))
os.environ.setdefault("XDG_CACHE_HOME", str(_cache_dir))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


EPS = 1e-12
KPC_M = 3.085677581491367e19
G_DAG_M_S2 = 1.2e-10


COLUMN_ALIASES = {
    "galaxy": ["galaxy", "name", "id", "Galaxy", "Name", "ID"],
    "r_kpc": ["r_kpc", "R", "Rad", "radius", "Radius", "R_kpc"],
    "v_obs_km_s": ["v_obs_km_s", "Vobs", "V_obs", "vobs", "Vflat_obs"],
    "e_v_obs_km_s": ["e_v_obs_km_s", "eVobs", "e_Vobs", "Vobs_err", "err"],
    "v_gas_km_s": ["v_gas_km_s", "Vgas", "V_gas", "vgas"],
    "v_disk_km_s": ["v_disk_km_s", "Vdisk", "V_disk", "vdisk"],
    "v_bulge_km_s": ["v_bulge_km_s", "Vbul", "Vbulge", "V_bulge", "vbulge"],
    "upsilon_disk": ["upsilon_disk", "Ydisk", "MLdisk", "Upsilon_disk"],
    "upsilon_bulge": ["upsilon_bulge", "Ybulge", "MLbulge", "Upsilon_bulge"],
}


def signed_square(v: np.ndarray) -> np.ndarray:
    """Preserve SPARC-style signed velocity contributions if present."""
    return np.sign(v) * v * v


def first_present(columns: Iterable[str], aliases: list[str]) -> str | None:
    available = set(columns)
    for alias in aliases:
        if alias in available:
            return alias
    return None


def standardise_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame()
    for target, aliases in COLUMN_ALIASES.items():
        source = first_present(df.columns, aliases)
        if source is not None:
            out[target] = df[source]

    required = [
        "galaxy",
        "r_kpc",
        "v_obs_km_s",
        "v_gas_km_s",
        "v_disk_km_s",
        "v_bulge_km_s",
    ]
    missing = [col for col in required if col not in out.columns]
    if missing:
        raise ValueError(
            "Missing required columns after alias mapping: "
            + ", ".join(missing)
            + ". Expected at least galaxy, radius, observed velocity, gas, disk, and bulge velocities."
        )

    if "e_v_obs_km_s" not in out.columns:
        out["e_v_obs_km_s"] = np.maximum(0.05 * np.asarray(out["v_obs_km_s"], dtype=float), 3.0)
    if "upsilon_disk" not in out.columns:
        out["upsilon_disk"] = np.nan
    if "upsilon_bulge" not in out.columns:
        out["upsilon_bulge"] = np.nan

    numeric = [col for col in out.columns if col != "galaxy"]
    for col in numeric:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out = out.dropna(subset=["galaxy", "r_kpc", "v_obs_km_s"])
    return out.sort_values(["galaxy", "r_kpc"]).reset_index(drop=True)


def make_mock_sparc_like(seed: int = 61) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    classes = [
        ("LSB_01", "low_surface_brightness", 96.0, 7.0, 46.0, 0.18, 0.0),
        ("LSB_02", "low_surface_brightness", 112.0, 8.5, 54.0, 0.15, 0.0),
        ("HSB_01", "high_surface_brightness", 188.0, 3.0, 138.0, 0.42, 0.0),
        ("HSB_02", "high_surface_brightness", 214.0, 3.5, 158.0, 0.40, 0.0),
        ("DW_01", "dwarf", 72.0, 5.8, 38.0, 0.20, 0.0),
        ("DW_02", "dwarf", 82.0, 6.2, 42.0, 0.22, 0.0),
        ("BULGE_01", "bulge_dominated", 236.0, 2.5, 130.0, 0.48, 86.0),
        ("BULGE_02", "bulge_dominated", 252.0, 2.8, 142.0, 0.46, 96.0),
    ]
    radii = np.linspace(0.7, 24.0, 34)
    for name, family, vflat, rturn, disk_peak, gas_frac, bulge_peak in classes:
        for r in radii:
            v_obs = vflat * (1.0 - np.exp(-r / rturn))
            if family == "dwarf":
                v_obs *= 1.0 + 0.04 * np.tanh((r - 10.0) / 3.0)
            if family == "bulge_dominated":
                v_obs += 18.0 * np.exp(-0.5 * ((r - 2.0) / 1.2) ** 2)
            noise = rng.normal(0.0, 0.012 * max(vflat, 50.0))
            v_obs = max(v_obs + noise, 5.0)

            disk = disk_peak * (r / 3.0) * np.exp(0.5 * (1.0 - r / 3.0))
            gas = gas_frac * vflat * (1.0 - np.exp(-r / 9.0))
            bulge = bulge_peak * np.exp(-r / 2.2)
            rows.append(
                {
                    "galaxy": name,
                    "family": family,
                    "r_kpc": r,
                    "v_obs_km_s": v_obs,
                    "e_v_obs_km_s": max(3.0, 0.04 * v_obs),
                    "v_gas_km_s": gas,
                    "v_disk_km_s": disk,
                    "v_bulge_km_s": bulge,
                    "upsilon_disk": 0.5,
                    "upsilon_bulge": 0.7,
                }
            )
    return pd.DataFrame(rows)


def gradient_log(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    if len(y) < 3:
        return np.zeros_like(y, dtype=float)
    return np.gradient(np.log(np.maximum(y, EPS)), np.log(np.maximum(x, EPS)))


def nfw_shape(r: np.ndarray, rs: float) -> np.ndarray:
    x = np.maximum(r / max(rs, EPS), EPS)
    shape = np.log1p(x) - x / (1.0 + x)
    return np.maximum(shape / np.maximum(r, EPS), 0.0)


def fit_nfw_like(r: np.ndarray, residual_v2: np.ndarray, err_v: np.ndarray) -> tuple[np.ndarray, dict]:
    weights = 1.0 / np.maximum(2.0 * np.sqrt(np.maximum(residual_v2, 1.0)) * err_v, 5.0) ** 2
    rs_grid = np.geomspace(0.6, 40.0, 96)
    best = None
    for rs in rs_grid:
        shape = nfw_shape(r, rs)
        denom = np.sum(weights * shape * shape) + EPS
        amp = max(float(np.sum(weights * residual_v2 * shape) / denom), 0.0)
        pred = amp * shape
        rmse = float(np.sqrt(np.mean((pred - residual_v2) ** 2)))
        chi2 = float(np.mean(weights * (pred - residual_v2) ** 2))
        if best is None or chi2 < best["chi2"]:
            best = {"rs": float(rs), "amp": float(amp), "rmse_v2": rmse, "chi2": chi2, "pred": pred}
    assert best is not None
    pred = np.asarray(best.pop("pred"))
    return pred, best


def rar_velocity(r: np.ndarray, vbar2: np.ndarray) -> np.ndarray:
    gbar = np.maximum(vbar2, 0.0) * 1.0e6 / (np.maximum(r, EPS) * KPC_M)
    gobs_pred = gbar / np.maximum(1.0 - np.exp(-np.sqrt(gbar / G_DAG_M_S2)), EPS)
    return np.sqrt(np.maximum(gobs_pred * r * KPC_M, 0.0)) / 1000.0


def compute_galaxy(df_gal: pd.DataFrame, default_ups_disk: float, default_ups_bulge: float) -> tuple[pd.DataFrame, dict]:
    g = df_gal.sort_values("r_kpc").copy()
    r = g["r_kpc"].to_numpy(float)
    vobs = g["v_obs_km_s"].to_numpy(float)
    err = g["e_v_obs_km_s"].to_numpy(float)
    vgas = g["v_gas_km_s"].to_numpy(float)
    vdisk = g["v_disk_km_s"].to_numpy(float)
    vbulge = g["v_bulge_km_s"].to_numpy(float)
    ups_disk = g["upsilon_disk"].fillna(default_ups_disk).to_numpy(float)
    ups_bulge = g["upsilon_bulge"].fillna(default_ups_bulge).to_numpy(float)

    vbar2 = signed_square(vgas) + ups_disk * signed_square(vdisk) + ups_bulge * signed_square(vbulge)
    vbar2 = np.maximum(vbar2, 0.0)
    vobs2 = np.maximum(vobs * vobs, EPS)
    residual_v2 = np.maximum(vobs2 - vbar2, 0.0)
    vbar = np.sqrt(vbar2)
    vres = np.sqrt(residual_v2)
    residual_fraction = residual_v2 / vobs2

    overshoot = np.maximum(vbar2 - vobs2, 0.0) / vobs2
    H = 1.0 / (1.0 + 10.0 * overshoot)
    slope = gradient_log(vobs, r)
    curvature = np.abs(np.gradient(slope, np.log(np.maximum(r, EPS)))) if len(r) >= 3 else np.zeros_like(r)
    B = 1.0 / (1.0 + np.abs(slope) + 0.35 * curvature)

    activity = np.abs(np.gradient(residual_fraction, r)) if len(r) >= 3 else np.zeros_like(r)
    activity_norm = activity / (np.nanmedian(activity[activity > 0]) + EPS)
    S_eff = np.exp(-0.5 * (np.log(np.maximum(activity_norm, EPS)) / 1.35) ** 2)

    m_res = r * residual_v2
    rho_res = np.maximum(np.gradient(m_res, r) / (4.0 * np.pi * r * r + EPS), 0.0) if len(r) >= 3 else np.zeros_like(r)
    rho_res_norm = rho_res / (np.nanmax(rho_res) + EPS)
    m_bar = r * vbar2
    rho_bar = np.maximum(np.gradient(m_bar, r) / (4.0 * np.pi * r * r + EPS), 0.0) if len(r) >= 3 else np.zeros_like(r)
    rho_bar_norm = rho_bar / (np.nanmax(rho_bar) + EPS)
    res_slope = np.abs(gradient_log(rho_res_norm + 0.03, r))
    bar_slope = np.abs(gradient_log(rho_bar_norm + 0.03, r))
    shape_mismatch = np.abs(res_slope - bar_slope) / (res_slope + bar_slope + EPS)
    radial_overlap = np.sqrt(np.maximum(rho_res_norm, 0.0) * np.maximum(rho_bar_norm, 0.0))
    V = np.clip(0.65 / (1.0 + shape_mismatch) + 0.35 * radial_overlap, 0.0, 1.0)

    R_resp = np.power(np.maximum(H * B * V, 0.0), 1.0 / 3.0)
    R_rob = np.minimum(R_resp, np.power(np.maximum(H * B * S_eff * V, 0.0), 0.25))
    D_struct = residual_fraction * R_rob * V

    nfw_v2, nfw_info = fit_nfw_like(r, residual_v2, err)
    v_rar = rar_velocity(r, vbar2)
    rar_residual = np.abs(vobs - v_rar) / np.maximum(err, 1.0)
    nfw_residual = np.abs(residual_v2 - nfw_v2) / np.maximum(vobs2, EPS)

    out = g.copy()
    out["v_baryon_km_s"] = vbar
    out["v_residual_km_s"] = vres
    out["residual_fraction"] = residual_fraction
    out["H"] = H
    out["B"] = B
    out["S_eff"] = S_eff
    out["V"] = V
    out["R_resp"] = R_resp
    out["R_rob"] = R_rob
    out["D_struct"] = D_struct
    out["v_nfw_like_km_s"] = np.sqrt(np.maximum(nfw_v2, 0.0))
    out["v_rar_km_s"] = v_rar
    out["nfw_fractional_residual"] = nfw_residual
    out["rar_sigma_residual"] = rar_residual

    summary = {
        "galaxy": str(g["galaxy"].iloc[0]),
        "family": str(g["family"].iloc[0]) if "family" in g.columns else "unknown",
        "n_points": int(len(g)),
        "mean_residual_fraction": float(np.mean(residual_fraction)),
        "outer_residual_fraction": float(np.mean(residual_fraction[r >= np.nanmedian(r)])),
        "mean_D_struct": float(np.mean(D_struct)),
        "peak_D_struct": float(np.max(D_struct)),
        "mean_R_rob": float(np.mean(R_rob)),
        "mean_V": float(np.mean(V)),
        "nfw_like_rmse_v2": float(nfw_info["rmse_v2"]),
        "nfw_like_chi2": float(nfw_info["chi2"]),
        "nfw_like_rs_kpc": float(nfw_info["rs"]),
        "rar_mean_sigma_residual": float(np.mean(rar_residual)),
    }
    return out, summary


def safe_corr(x: Iterable[float], y: Iterable[float], method: str = "pearson") -> float:
    s = pd.DataFrame({"x": list(x), "y": list(y)}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(s) < 3 or s["x"].nunique() < 2 or s["y"].nunique() < 2:
        return float("nan")
    return float(s["x"].corr(s["y"], method=method))


def make_nulls(summary: pd.DataFrame, seed: int = 61, n_null: int = 300) -> dict:
    rng = np.random.default_rng(seed)
    real_nfw = safe_corr(summary["mean_D_struct"], summary["nfw_like_chi2"], "spearman")
    real_rar = safe_corr(summary["mean_D_struct"], summary["rar_mean_sigma_residual"], "spearman")
    null_nfw = []
    null_rar = []
    dvals = summary["mean_D_struct"].to_numpy(float)
    for _ in range(n_null):
        shuffled = rng.permutation(dvals)
        null_nfw.append(safe_corr(shuffled, summary["nfw_like_chi2"], "spearman"))
        null_rar.append(safe_corr(shuffled, summary["rar_mean_sigma_residual"], "spearman"))
    return {
        "real_spearman_Dstruct_vs_nfw_chi2": real_nfw,
        "real_spearman_Dstruct_vs_rar_residual": real_rar,
        "null_mean_abs_spearman_nfw": float(np.nanmean(np.abs(null_nfw))),
        "null_mean_abs_spearman_rar": float(np.nanmean(np.abs(null_rar))),
        "n_null": int(n_null),
    }


def plot_examples(profile: pd.DataFrame, outdir: Path) -> None:
    galaxies = list(profile["galaxy"].drop_duplicates())[:4]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=False, sharey=False)
    for ax, galaxy in zip(axes.ravel(), galaxies):
        g = profile[profile["galaxy"] == galaxy]
        ax.errorbar(g["r_kpc"], g["v_obs_km_s"], yerr=g["e_v_obs_km_s"], fmt="o", ms=3, label="obs", color="#1b4965")
        ax.plot(g["r_kpc"], g["v_baryon_km_s"], lw=2, label="baryon", color="#ca6702")
        ax.plot(g["r_kpc"], g["v_residual_km_s"], lw=2, label="residual", color="#5e548e")
        ax.plot(g["r_kpc"], g["v_nfw_like_km_s"], lw=1.6, ls="--", label="NFW-like residual", color="#006d77")
        ax.set_title(str(galaxy))
        ax.set_xlabel("r [kpc]")
        ax.set_ylabel("v [km/s]")
        ax.grid(alpha=0.25)
    handles, labels = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(outdir / "fig1_rotation_decomposition_examples.png", dpi=180)
    plt.close(fig)


def plot_population(summary: pd.DataFrame, outdir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    family = summary["family"].astype("category")
    colors = family.cat.codes
    sc = axes[0].scatter(summary["mean_D_struct"], summary["nfw_like_chi2"], c=colors, cmap="tab10", s=70, alpha=0.85)
    axes[0].set_xlabel("mean structural dark support")
    axes[0].set_ylabel("NFW-like chi2")
    axes[0].set_title("Structural support vs NFW-like mismatch")
    axes[0].grid(alpha=0.25)
    axes[1].scatter(summary["mean_D_struct"], summary["rar_mean_sigma_residual"], c=colors, cmap="tab10", s=70, alpha=0.85)
    axes[1].set_xlabel("mean structural dark support")
    axes[1].set_ylabel("RAR mean sigma residual")
    axes[1].set_title("Structural support vs RAR mismatch")
    axes[1].grid(alpha=0.25)
    labels = list(family.cat.categories)
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", label=label, markerfacecolor=sc.cmap(sc.norm(i)), markersize=8)
        for i, label in enumerate(labels)
    ]
    fig.legend(handles=handles, loc="upper center", ncol=min(4, len(handles)), frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(outdir / "fig2_population_structural_vs_baselines.png", dpi=180)
    plt.close(fig)


def plot_nulls(nulls: dict, outdir: Path) -> None:
    labels = ["Dstruct vs NFW", "Dstruct vs RAR"]
    real = [abs(nulls["real_spearman_Dstruct_vs_nfw_chi2"]), abs(nulls["real_spearman_Dstruct_vs_rar_residual"])]
    null = [nulls["null_mean_abs_spearman_nfw"], nulls["null_mean_abs_spearman_rar"]]
    x = np.arange(len(labels))
    width = 0.35
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.bar(x - width / 2, real, width, label="real pairing", color="#006d77")
    ax.bar(x + width / 2, null, width, label="shuffled galaxy null", color="#adb5bd")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("|Spearman rho|")
    ax.set_title("Predeclared galaxy-level null comparison")
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "fig3_shuffled_galaxy_nulls.png", dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict:
    in_path = Path(args.input)
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    if in_path.exists() and args.mode != "mock":
        raw = pd.read_csv(in_path)
        data = standardise_columns(raw)
        mode = "real_input"
        data_note = f"Loaded local input file: {in_path}"
    else:
        data = make_mock_sparc_like(args.seed)
        mode = "mock_schema_test"
        data_note = "No local SPARC CSV was used; generated synthetic SPARC-like schema test data."

    profiles = []
    summaries = []
    for _, group in data.groupby("galaxy", sort=False):
        profile, summary = compute_galaxy(group, args.upsilon_disk, args.upsilon_bulge)
        profiles.append(profile)
        summaries.append(summary)

    profile_df = pd.concat(profiles, ignore_index=True)
    summary_df = pd.DataFrame(summaries)
    nulls = make_nulls(summary_df, args.seed, args.n_null)

    profile_df.to_csv(outdir / "paper61_radius_diagnostics.csv", index=False)
    summary_df.to_csv(outdir / "paper61_galaxy_summary.csv", index=False)
    plot_examples(profile_df, outdir)
    plot_population(summary_df, outdir)
    plot_nulls(nulls, outdir)

    result = {
        "paper": 61,
        "title": "SPARC Pilot Protocol for Structural Dark-Matter Residuals",
        "mode": mode,
        "data_note": data_note,
        "n_galaxies": int(summary_df["galaxy"].nunique()),
        "n_radius_points": int(len(profile_df)),
        "mean_residual_fraction": float(profile_df["residual_fraction"].mean()),
        "mean_D_struct": float(profile_df["D_struct"].mean()),
        "mean_R_rob": float(profile_df["R_rob"].mean()),
        "mean_V": float(profile_df["V"].mean()),
        **nulls,
        "license_note": (
            "If SPARC data are used or redistributed, cite Lelli, McGaugh & "
            "Schombert (2016), cite DOI 10.5281/zenodo.16284118 for the Zenodo "
            "record when applicable, respect CC-BY-4.0 attribution, and include "
            "VizieR/CDS acknowledgment if using VizieR. Derived MAAT artifacts "
            "must be distinguished from original SPARC measurements."
        ),
    }
    with open(outdir / "paper61_summary.json", "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Paper 61 SPARC structural residual pilot")
    parser.add_argument("--input", default="data/sparc_rotation_curves.csv", help="Local SPARC-like CSV path")
    parser.add_argument("--output", default="outputs_paper61", help="Output directory")
    parser.add_argument("--mode", choices=["auto", "mock"], default="auto", help="auto loads input if present; mock forces synthetic schema test")
    parser.add_argument("--upsilon-disk", type=float, default=0.5, help="Default disk mass-to-light factor")
    parser.add_argument("--upsilon-bulge", type=float, default=0.7, help="Default bulge mass-to-light factor")
    parser.add_argument("--seed", type=int, default=61)
    parser.add_argument("--n-null", type=int, default=300)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
