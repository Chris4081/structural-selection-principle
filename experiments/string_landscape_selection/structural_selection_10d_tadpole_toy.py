#!/usr/bin/env python3
"""10D-near tadpole-closure toy test for structural selection.

The script builds a small type-IIB-inspired flux ensemble with an explicit
symplectic D3 tadpole proxy

    N_flux = F^T J H,
    Delta Q_D3 = N_flux + N_D3 - L.

It is not a full compactification scan. Its purpose is narrower: operationalise
the B_10 sector and test whether a structural support score carries ranking
information that is invisible to energy-only ordering.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/codex-mpl-cache")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/codex-mpl-cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


@dataclass(frozen=True)
class ScanConfig:
    seed: int = 4081
    n_candidates: int = 320
    flux_dim: int = 8
    flux_abs_max: int = 4
    tadpole_bound: int = 24
    lambda_scale: float = 1.0e-4
    lambda_obs: float = 0.0
    lambda_sigma: float = 6.0e-5
    epsilon: float = 1.0e-9
    z_st: float = 1.0


def symplectic_tadpole(f_flux: np.ndarray, h_flux: np.ndarray) -> int:
    """Return the integer proxy F^T J H for an even-dimensional flux vector."""
    half = f_flux.size // 2
    return int(np.dot(f_flux[:half], h_flux[half:]) - np.dot(f_flux[half:], h_flux[:half]))


def support(defect: float) -> float:
    return 1.0 / (1.0 + max(0.0, defect))


def rankdata(values: Iterable[float]) -> np.ndarray:
    """Simple one-based ranks without tie averaging, sufficient for this scan."""
    values_arr = np.asarray(list(values), dtype=float)
    order = np.argsort(values_arr, kind="mergesort")
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(values_arr) + 1, dtype=float)
    return ranks


def spearman(x_values: Iterable[float], y_values: Iterable[float]) -> float:
    x_rank = rankdata(x_values)
    y_rank = rankdata(y_values)
    if np.std(x_rank) == 0.0 or np.std(y_rank) == 0.0:
        return float("nan")
    return float(np.corrcoef(x_rank, y_rank)[0, 1])


def candidate_from_fluxes(
    idx: int,
    f_flux: np.ndarray,
    h_flux: np.ndarray,
    rng: np.random.Generator,
    config: ScanConfig,
) -> dict:
    n_flux = symplectic_tadpole(f_flux, h_flux)
    flux_norm = float(np.dot(f_flux, f_flux) + np.dot(h_flux, h_flux))
    max_flux_norm = float(2 * config.flux_dim * config.flux_abs_max**2)
    flux_activity = flux_norm / max_flux_norm

    brane_noise = int(rng.integers(-5, 6))
    n_d3_ideal = config.tadpole_bound - n_flux
    n_d3 = max(0, int(round(n_d3_ideal + brane_noise)))
    delta_q = int(n_flux + n_d3 - config.tadpole_bound)

    volume = (
        18.0
        + 1.8 * max(0.0, config.tadpole_bound - abs(n_flux))
        + rng.gamma(shape=2.2, scale=3.0)
        - 0.18 * flux_norm
    )
    volume = float(max(4.0, volume))
    compact_length = volume ** (1.0 / 6.0)

    alignment_den = float(np.linalg.norm(f_flux) * np.linalg.norm(h_flux))
    alignment = 0.0 if alignment_den == 0.0 else float(abs(np.dot(f_flux, h_flux)) / alignment_den)

    lambda_eff = float(
        rng.normal(config.lambda_obs, config.lambda_scale)
        + 2.5e-5 * (flux_activity - 0.45)
        + rng.normal(0.0, 1.4e-5)
    )

    stationarity_residual = abs(rng.normal(0.0, 0.035)) + 0.035 * flux_activity
    activity_target = 0.48

    d_h = stationarity_residual**2
    d_b = 18.0 * (delta_q / (config.tadpole_bound + 1.0)) ** 2
    d_s = 3.0 * (flux_activity - activity_target) ** 2
    d_v = 2.5 * (alignment - 0.35) ** 2
    d_r = max(0.0, (16.0 - volume) / 16.0) ** 2 + 0.15 * math.exp(-volume / 24.0)
    d_sigma = (1.0 / compact_length) ** 2 + 0.08 * alignment**2

    gamma_h = support(d_h)
    gamma_b = support(d_b)
    gamma_s = support(d_s)
    gamma_v = support(d_v)
    gamma_r = support(d_r)
    gamma_sigma = support(d_sigma)
    gammas = {
        "H": gamma_h,
        "B": gamma_b,
        "S": gamma_s,
        "V": gamma_v,
        "R": gamma_r,
        "Sigma": gamma_sigma,
    }

    f_sel_product = -float(np.mean([math.log(config.epsilon + val) for val in gammas.values()]))

    delta_e_10 = abs(lambda_eff - config.lambda_obs) / config.lambda_scale
    delta_q_10 = abs(delta_q) / (config.tadpole_bound + 1.0)
    delta_dim_10 = (1.0 / compact_length) ** 2 + (1.0 / max(compact_length - 1.0, 0.25)) ** 2
    maat_numerator = float(np.mean([gamma_h, gamma_b, gamma_s, gamma_v, gamma_r]))
    m10_support = (
        config.z_st
        * maat_numerator
        / (config.epsilon + delta_e_10 + delta_q_10 + delta_dim_10)
    )
    f_sel_10 = -math.log(config.epsilon + m10_support)

    sw_penalty = 12.0 * delta_q_10**2 + 3.0 * max(0.0, 1.45 - compact_length) ** 2
    lambda_penalty = 0.5 * ((lambda_eff - config.lambda_obs) / config.lambda_sigma) ** 2
    f_pheno = f_sel_10 + sw_penalty + lambda_penalty

    return {
        "id": idx,
        "F_flux": [int(x) for x in f_flux.tolist()],
        "H_flux": [int(x) for x in h_flux.tolist()],
        "N_flux": n_flux,
        "N_D3": n_d3,
        "Delta_Q_D3": delta_q,
        "abs_Delta_Q_D3": abs(delta_q),
        "flux_norm": flux_norm,
        "flux_activity": flux_activity,
        "alignment": alignment,
        "volume": volume,
        "compact_length": compact_length,
        "lambda_eff": lambda_eff,
        "abs_lambda_eff": abs(lambda_eff - config.lambda_obs),
        "defects": {
            "d_H": d_h,
            "d_B": d_b,
            "d_S": d_s,
            "d_V": d_v,
            "d_R": d_r,
            "d_Sigma": d_sigma,
            "Delta_E_10": delta_e_10,
            "Delta_Q_10": delta_q_10,
            "Delta_dim_10": delta_dim_10,
        },
        "supports": gammas,
        "M10_MAAT": m10_support,
        "F_sel_product": f_sel_product,
        "F_sel_10": f_sel_10,
        "F_SW": sw_penalty,
        "lambda_penalty": lambda_penalty,
        "F_pheno": f_pheno,
    }


def generate_scan(config: ScanConfig) -> list[dict]:
    rng = np.random.default_rng(config.seed)
    candidates: list[dict] = []
    seen: set[tuple[int, ...]] = set()

    attempts = 0
    while len(candidates) < config.n_candidates and attempts < config.n_candidates * 200:
        attempts += 1
        f_flux = rng.integers(-config.flux_abs_max, config.flux_abs_max + 1, size=config.flux_dim)
        h_flux = rng.integers(-config.flux_abs_max, config.flux_abs_max + 1, size=config.flux_dim)
        if not np.any(f_flux) or not np.any(h_flux):
            continue
        key = tuple(int(x) for x in np.concatenate([f_flux, h_flux]))
        if key in seen:
            continue
        seen.add(key)
        n_flux = symplectic_tadpole(f_flux, h_flux)
        if n_flux < -config.tadpole_bound or n_flux > 2 * config.tadpole_bound:
            continue
        candidates.append(candidate_from_fluxes(len(candidates), f_flux, h_flux, rng, config))

    if len(candidates) < config.n_candidates:
        raise RuntimeError(f"Only generated {len(candidates)} candidates after {attempts} attempts.")

    energy_ranks = rankdata(item["abs_lambda_eff"] for item in candidates)
    structural_ranks = rankdata(item["F_sel_10"] for item in candidates)
    pheno_ranks = rankdata(item["F_pheno"] for item in candidates)
    for item, erank, srank, prank in zip(candidates, energy_ranks, structural_ranks, pheno_ranks):
        item["rank_energy"] = int(erank)
        item["rank_structural"] = int(srank)
        item["rank_pheno"] = int(prank)
    return candidates


def top_summary(candidates: list[dict], key: str, top_n: int = 20) -> dict:
    top = sorted(candidates, key=lambda item: item[key])[:top_n]
    return {
        "score": key,
        "top_n": top_n,
        "mean_abs_delta_q": float(np.mean([item["abs_Delta_Q_D3"] for item in top])),
        "median_abs_delta_q": float(np.median([item["abs_Delta_Q_D3"] for item in top])),
        "mean_abs_lambda": float(np.mean([item["abs_lambda_eff"] for item in top])),
        "mean_volume": float(np.mean([item["volume"] for item in top])),
        "mean_M10_MAAT": float(np.mean([item["M10_MAAT"] for item in top])),
        "ids": [item["id"] for item in top],
    }


def summarize(candidates: list[dict], config: ScanConfig) -> dict:
    top_energy = set(top_summary(candidates, "abs_lambda_eff")["ids"])
    top_struct = set(top_summary(candidates, "F_sel_10")["ids"])
    top_pheno = set(top_summary(candidates, "F_pheno")["ids"])
    return {
        "config": asdict(config),
        "n_candidates": len(candidates),
        "spearman_abs_lambda_vs_F_sel_10": spearman(
            [item["abs_lambda_eff"] for item in candidates],
            [item["F_sel_10"] for item in candidates],
        ),
        "spearman_abs_delta_q_vs_F_sel_10": spearman(
            [item["abs_Delta_Q_D3"] for item in candidates],
            [item["F_sel_10"] for item in candidates],
        ),
        "spearman_abs_delta_q_vs_abs_lambda": spearman(
            [item["abs_Delta_Q_D3"] for item in candidates],
            [item["abs_lambda_eff"] for item in candidates],
        ),
        "top_energy": top_summary(candidates, "abs_lambda_eff"),
        "top_structural": top_summary(candidates, "F_sel_10"),
        "top_pheno": top_summary(candidates, "F_pheno"),
        "top20_overlap_energy_structural": len(top_energy & top_struct),
        "top20_overlap_energy_pheno": len(top_energy & top_pheno),
        "best_energy_candidate": min(candidates, key=lambda item: item["abs_lambda_eff"]),
        "best_structural_candidate": min(candidates, key=lambda item: item["F_sel_10"]),
        "best_pheno_candidate": min(candidates, key=lambda item: item["F_pheno"]),
    }


def plot_energy_vs_structural(candidates: list[dict], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    scatter = ax.scatter(
        [item["abs_lambda_eff"] for item in candidates],
        [item["F_sel_10"] for item in candidates],
        c=[item["abs_Delta_Q_D3"] for item in candidates],
        cmap="viridis",
        s=38,
        alpha=0.78,
        edgecolors="none",
    )
    best_struct = sorted(candidates, key=lambda item: item["F_sel_10"])[:5]
    for item in best_struct:
        ax.annotate(
            f"id={item['id']}\n|dQ|={item['abs_Delta_Q_D3']}",
            (item["abs_lambda_eff"], item["F_sel_10"]),
            textcoords="offset points",
            xytext=(6, 5),
            fontsize=7,
            bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.7},
        )
    ax.set_xscale("log")
    ax.set_xlabel(r"$|\Lambda_{\mathrm{eff}}-\Lambda_{\mathrm{obs}}|$")
    ax.set_ylabel(r"$F_{\mathrm{sel}}^{(10)}$")
    ax.set_title("Energy-Only Proximity Versus 10D Structural Selection")
    ax.grid(True, alpha=0.25)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(r"$|\Delta Q_{D3}|$")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_rank_comparison(candidates: list[dict], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 6.0))
    scatter = ax.scatter(
        [item["rank_energy"] for item in candidates],
        [item["rank_structural"] for item in candidates],
        c=[item["abs_Delta_Q_D3"] for item in candidates],
        cmap="plasma",
        s=38,
        alpha=0.8,
    )
    max_rank = len(candidates)
    ax.plot([1, max_rank], [1, max_rank], color="gray", linestyle="--", linewidth=1.0)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.set_xlabel("Rank by energy proximity")
    ax.set_ylabel(r"Rank by $F_{\mathrm{sel}}^{(10)}$")
    ax.set_title("Structural Ranking Is Not Energy Ranking")
    ax.grid(True, alpha=0.24)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(r"$|\Delta Q_{D3}|$")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_tadpole_bins(candidates: list[dict], output: Path) -> None:
    bins = [(0, 0), (1, 2), (3, 5), (6, 99)]
    labels = ["0", "1-2", "3-5", ">=6"]
    grouped = []
    for lo, hi in bins:
        grouped.append([item["F_sel_10"] for item in candidates if lo <= item["abs_Delta_Q_D3"] <= hi])

    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    box = ax.boxplot(grouped, tick_labels=labels, patch_artist=True, showfliers=False)
    colors = plt.get_cmap("Set3")(np.linspace(0.15, 0.85, len(grouped)))
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.85)
    means = [np.mean(values) if values else np.nan for values in grouped]
    ax.scatter(range(1, len(means) + 1), means, marker="D", color="black", s=34, label="Mean")
    ax.set_xlabel(r"Tadpole mismatch bin $|\Delta Q_{D3}|$")
    ax.set_ylabel(r"$F_{\mathrm{sel}}^{(10)}$")
    ax.set_title(r"Operational $B_{10}$ Sector Penalises Tadpole Mismatch")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_top_decomposition(candidates: list[dict], output: Path, top_n: int = 8) -> None:
    top = sorted(candidates, key=lambda item: item["F_sel_10"])[:top_n]
    parts = [
        ("Delta_E_10", r"$\Delta_E^{(10)}$"),
        ("Delta_Q_10", r"$\Delta_Q^{(10)}$"),
        ("Delta_dim_10", r"$\Delta_{\mathrm{dim}}^{(10)}$"),
    ]
    labels = [f"id={item['id']}\n|dQ|={item['abs_Delta_Q_D3']}" for item in top]
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    bottom = np.zeros(len(top))
    colors = ["#3b6ea8", "#d95f02", "#2ca25f"]
    for color, (key, label) in zip(colors, parts):
        values = np.array([item["defects"][key] for item in top], dtype=float)
        ax.bar(labels, values, bottom=bottom, label=label, color=color, alpha=0.88)
        bottom += values
    ax.set_ylabel("Denominator contribution")
    ax.set_title(r"Top 10D Structural Candidates: Obstruction Decomposition")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def make_plots(candidates: list[dict], output_dir: Path) -> list[str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    plots = {
        "energy_vs_structural": output_dir / "energy_vs_structural.png",
        "rank_comparison": output_dir / "rank_comparison.png",
        "tadpole_bins": output_dir / "tadpole_bins.png",
        "top_obstruction_decomposition": output_dir / "top_obstruction_decomposition.png",
    }
    plot_energy_vs_structural(candidates, plots["energy_vs_structural"])
    plot_rank_comparison(candidates, plots["rank_comparison"])
    plot_tadpole_bins(candidates, plots["tadpole_bins"])
    plot_top_decomposition(candidates, plots["top_obstruction_decomposition"])
    return [str(path) for path in plots.values()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=4081)
    parser.add_argument("--n-candidates", type=int, default=320)
    parser.add_argument("--output", type=Path, default=Path("structural_selection_10d_tadpole_toy_results.json"))
    parser.add_argument("--plot-dir", type=Path, default=Path("structural_selection_10d_tadpole_toy_plots"))
    args = parser.parse_args()

    config = ScanConfig(seed=args.seed, n_candidates=args.n_candidates)
    candidates = generate_scan(config)
    summary = summarize(candidates, config)
    plots = make_plots(candidates, args.plot_dir)

    payload = {
        "description": "10D-near type-IIB-inspired tadpole-closure toy test for structural selection.",
        "summary": summary,
        "plots": plots,
        "candidates": candidates,
    }
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print("plots:")
    for path in plots:
        print(f"  {path}")


if __name__ == "__main__":
    main()
