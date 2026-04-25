#!/usr/bin/env python3
"""IIB/KKLT bridge test for structural selection.

This script bridges the 10D-near flux/tadpole toy model and the reduced KKLT
benchmark. Integer flux vectors generate a D3-tadpole residual and a small
effective flux superpotential W0. The resulting KKLT model is then solved for
stationary points, and candidates are ranked by energy proximity, ordinary KKLT
structural score, and a bridge score combining B_10 tadpole closure with
stationarity, metastability, and EFT-control diagnostics.
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

from structural_selection_iib_kklt_scan import (
    DynamicParams,
    StructuralParams,
    candidate_record,
    solve_stationary_points,
)


@dataclass(frozen=True)
class BridgeConfig:
    seed: int = 2604
    n_backgrounds: int = 90
    flux_dim: int = 8
    flux_abs_max: int = 4
    tadpole_bound: int = 24
    w0_min: float = 5.0e-5
    w0_max: float = 3.5e-4
    A: float = 1.0
    a: float = 0.1
    n_uplift: float = 2.0
    d_uplift_min: float = 0.0
    d_uplift_max: float = 1.2e-9
    tau_seed_min: float = 18.0
    tau_seed_max: float = 180.0
    num_tau_seeds: int = 18
    tau_min_stationary: float = 5.0
    tau_max_stationary: float = 180.0
    tau_barrier_max: float = 240.0
    barrier_samples: int = 260
    grad_tol: float = 1.0e-8
    epsilon: float = 1.0e-9
    lambda_obs: float = 0.0
    lambda_sigma: float = 2.5e-10


def symplectic_tadpole(f_flux: np.ndarray, h_flux: np.ndarray) -> int:
    half = f_flux.size // 2
    return int(np.dot(f_flux[:half], h_flux[half:]) - np.dot(f_flux[half:], h_flux[:half]))


def support(defect: float) -> float:
    return 1.0 / (1.0 + max(0.0, defect))


def rankdata(values: Iterable[float]) -> np.ndarray:
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


def flux_to_w0(f_flux: np.ndarray, h_flux: np.ndarray, config: BridgeConfig) -> float:
    """Map integer fluxes to a controlled small negative W0 proxy."""
    coeff_f = np.array([0.31, -0.19, 0.23, -0.11, 0.17, -0.29, 0.13, 0.07])
    coeff_h = np.array([-0.08, 0.21, -0.15, 0.27, -0.18, 0.12, -0.24, 0.16])
    raw = float(np.dot(coeff_f, f_flux) + np.dot(coeff_h, h_flux))
    squashed = 0.5 * (1.0 + math.tanh(raw / 4.0))
    magnitude = config.w0_min + (config.w0_max - config.w0_min) * squashed
    return -magnitude


def flux_to_uplift(f_flux: np.ndarray, h_flux: np.ndarray, config: BridgeConfig) -> float:
    flux_norm = float(np.dot(f_flux, f_flux) + np.dot(h_flux, h_flux))
    max_norm = float(2 * config.flux_dim * config.flux_abs_max**2)
    return config.d_uplift_min + (config.d_uplift_max - config.d_uplift_min) * (flux_norm / max_norm)


def make_background(idx: int, rng: np.random.Generator, config: BridgeConfig) -> dict:
    f_flux = rng.integers(-config.flux_abs_max, config.flux_abs_max + 1, size=config.flux_dim)
    h_flux = rng.integers(-config.flux_abs_max, config.flux_abs_max + 1, size=config.flux_dim)
    while not np.any(f_flux) or not np.any(h_flux):
        f_flux = rng.integers(-config.flux_abs_max, config.flux_abs_max + 1, size=config.flux_dim)
        h_flux = rng.integers(-config.flux_abs_max, config.flux_abs_max + 1, size=config.flux_dim)

    n_flux = symplectic_tadpole(f_flux, h_flux)
    n_d3_ideal = config.tadpole_bound - n_flux
    n_d3 = max(0, int(round(n_d3_ideal + int(rng.integers(-5, 6)))))
    delta_q = int(n_flux + n_d3 - config.tadpole_bound)
    flux_norm = float(np.dot(f_flux, f_flux) + np.dot(h_flux, h_flux))
    alignment_den = float(np.linalg.norm(f_flux) * np.linalg.norm(h_flux))
    alignment = 0.0 if alignment_den == 0.0 else float(abs(np.dot(f_flux, h_flux)) / alignment_den)

    return {
        "background_id": idx,
        "F_flux": [int(x) for x in f_flux.tolist()],
        "H_flux": [int(x) for x in h_flux.tolist()],
        "N_flux": n_flux,
        "N_D3": n_d3,
        "Delta_Q_D3": delta_q,
        "abs_Delta_Q_D3": abs(delta_q),
        "flux_norm": flux_norm,
        "flux_alignment": alignment,
        "W0": flux_to_w0(f_flux, h_flux, config),
        "D": flux_to_uplift(f_flux, h_flux, config),
    }


def bridge_record(
    background: dict,
    kklt_record: dict,
    config: BridgeConfig,
) -> dict:
    residuals = kklt_record["residuals"]
    scores = kklt_record["scores"]
    tau = float(kklt_record["tau"])
    energy = float(kklt_record["energy"])
    lambda_min = float(kklt_record["lambda_min"])
    abs_delta_q = float(background["abs_Delta_Q_D3"])
    delta_q_10 = abs_delta_q / (config.tadpole_bound + 1.0)
    d_b10 = 18.0 * delta_q_10**2
    gamma_b10 = support(d_b10)

    energy_scale = max(abs(background["W0"]) ** 2, abs(energy), 1.0e-18)
    delta_e_kklt = abs(energy - config.lambda_obs) / energy_scale
    delta_lambda_obs = abs(energy - config.lambda_obs) / config.lambda_sigma
    delta_dim_kklt = (10.0 / tau) ** 2 + math.exp(-config.a * tau)
    tachyon_defect = max(0.0, -lambda_min / max(abs(lambda_min), max(abs(x) for x in kklt_record["hessian_eigenvalues"]), 1.0e-18))

    gamma_h = scores["H"]
    gamma_s = scores["S"]
    gamma_v = scores["V"]
    gamma_r = scores["R"]
    gamma_sigma = scores["Sigma"]
    numerator = float(np.mean([gamma_h, gamma_b10, gamma_s, gamma_v, gamma_r]))

    bridge_support = numerator / (
        config.epsilon
        + delta_e_kklt
        + delta_q_10
        + delta_dim_kklt
        + 0.2 * tachyon_defect
    )
    f_bridge_10 = -math.log(config.epsilon + bridge_support)
    f_sw_bridge = 12.0 * delta_q_10**2 + 4.0 * max(0.0, 12.0 - tau) ** 2 / 144.0 + 2.0 * tachyon_defect
    lambda_penalty = 0.5 * delta_lambda_obs**2
    f_pheno_bridge = f_bridge_10 + f_sw_bridge + lambda_penalty

    merged = {
        **background,
        "tau": tau,
        "theta": kklt_record["theta"],
        "energy": energy,
        "abs_energy": abs(energy - config.lambda_obs),
        "vacuum_type": kklt_record["vacuum_type"],
        "lambda_min": lambda_min,
        "hessian_eigenvalues": kklt_record["hessian_eigenvalues"],
        "kklt_f_master": kklt_record["f_master"],
        "kklt_residuals": residuals,
        "kklt_scores": scores,
        "bridge_defects": {
            "d_B10": d_b10,
            "Delta_Q_10": delta_q_10,
            "Delta_E_KKLT": delta_e_kklt,
            "Delta_dim_KKLT": delta_dim_kklt,
            "tachyon_defect": tachyon_defect,
        },
        "bridge_support": bridge_support,
        "F_bridge_10": f_bridge_10,
        "F_SW_bridge": f_sw_bridge,
        "lambda_penalty": lambda_penalty,
        "F_pheno_bridge": f_pheno_bridge,
    }
    return merged


def generate_bridge_scan(config: BridgeConfig) -> tuple[list[dict], list[dict]]:
    rng = np.random.default_rng(config.seed)
    structural = StructuralParams()
    tau_seeds = np.linspace(config.tau_seed_min, config.tau_seed_max, config.num_tau_seeds)
    theta_seeds = [0.0, math.pi / config.a]

    backgrounds: list[dict] = []
    candidates: list[dict] = []
    seen_flux: set[tuple[int, ...]] = set()
    attempts = 0

    while len(backgrounds) < config.n_backgrounds and attempts < config.n_backgrounds * 100:
        attempts += 1
        background = make_background(len(backgrounds), rng, config)
        key = tuple(background["F_flux"] + background["H_flux"])
        if key in seen_flux:
            continue
        seen_flux.add(key)
        if background["N_flux"] < -config.tadpole_bound or background["N_flux"] > 2 * config.tadpole_bound:
            continue

        params = DynamicParams(
            W0=background["W0"],
            A=config.A,
            a=config.a,
            D=background["D"],
            n=config.n_uplift,
        )
        points = solve_stationary_points(
            params=params,
            tau_seeds=tau_seeds,
            theta_seeds=theta_seeds,
            tau_min=config.tau_min_stationary,
            tau_max=config.tau_max_stationary,
            grad_tol=config.grad_tol,
        )
        background = {**background, "stationary_points": [{"tau": tau, "theta": theta} for tau, theta in points]}
        backgrounds.append(background)

        for tau, theta in points:
            kklt = candidate_record(
                params=params,
                structural=structural,
                delta_q_d3=background["Delta_Q_D3"],
                tau=tau,
                theta=theta,
                tau_barrier_max=config.tau_barrier_max,
                barrier_samples=config.barrier_samples,
            )
            candidates.append(bridge_record(background, kklt, config))

    if not candidates:
        raise RuntimeError("No KKLT stationary candidates found. Adjust scan parameters.")

    energy_ranks = rankdata(item["abs_energy"] for item in candidates)
    kklt_ranks = rankdata(item["kklt_f_master"] for item in candidates)
    bridge_ranks = rankdata(item["F_bridge_10"] for item in candidates)
    pheno_ranks = rankdata(item["F_pheno_bridge"] for item in candidates)
    for item, erank, krank, brank, prank in zip(candidates, energy_ranks, kklt_ranks, bridge_ranks, pheno_ranks):
        item["rank_energy"] = int(erank)
        item["rank_kklt"] = int(krank)
        item["rank_bridge"] = int(brank)
        item["rank_pheno"] = int(prank)

    return backgrounds, candidates


def top_summary(candidates: list[dict], key: str, top_n: int = 12) -> dict:
    top = sorted(candidates, key=lambda item: item[key])[:top_n]
    return {
        "score": key,
        "top_n": top_n,
        "mean_abs_delta_q": float(np.mean([item["abs_Delta_Q_D3"] for item in top])),
        "median_abs_delta_q": float(np.median([item["abs_Delta_Q_D3"] for item in top])),
        "mean_abs_energy": float(np.mean([item["abs_energy"] for item in top])),
        "mean_tau": float(np.mean([item["tau"] for item in top])),
        "mean_lambda_min": float(np.mean([item["lambda_min"] for item in top])),
        "ids": [item["background_id"] for item in top],
    }


def summarize(backgrounds: list[dict], candidates: list[dict], config: BridgeConfig) -> dict:
    top_energy = set(top_summary(candidates, "abs_energy")["ids"])
    top_bridge = set(top_summary(candidates, "F_bridge_10")["ids"])
    top_pheno = set(top_summary(candidates, "F_pheno_bridge")["ids"])
    class_counts: dict[str, int] = {}
    for item in candidates:
        class_counts[item["vacuum_type"]] = class_counts.get(item["vacuum_type"], 0) + 1
    return {
        "config": asdict(config),
        "n_backgrounds": len(backgrounds),
        "n_stationary_candidates": len(candidates),
        "class_counts": class_counts,
        "spearman_abs_energy_vs_F_bridge_10": spearman(
            [item["abs_energy"] for item in candidates],
            [item["F_bridge_10"] for item in candidates],
        ),
        "spearman_abs_delta_q_vs_abs_energy": spearman(
            [item["abs_Delta_Q_D3"] for item in candidates],
            [item["abs_energy"] for item in candidates],
        ),
        "spearman_abs_delta_q_vs_F_bridge_10": spearman(
            [item["abs_Delta_Q_D3"] for item in candidates],
            [item["F_bridge_10"] for item in candidates],
        ),
        "top_energy": top_summary(candidates, "abs_energy"),
        "top_kklt": top_summary(candidates, "kklt_f_master"),
        "top_bridge": top_summary(candidates, "F_bridge_10"),
        "top_pheno": top_summary(candidates, "F_pheno_bridge"),
        "top12_overlap_energy_bridge": len(top_energy & top_bridge),
        "top12_overlap_energy_pheno": len(top_energy & top_pheno),
        "best_energy_candidate": min(candidates, key=lambda item: item["abs_energy"]),
        "best_bridge_candidate": min(candidates, key=lambda item: item["F_bridge_10"]),
        "best_pheno_candidate": min(candidates, key=lambda item: item["F_pheno_bridge"]),
    }


def plot_energy_bridge(candidates: list[dict], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    scatter = ax.scatter(
        [item["abs_energy"] for item in candidates],
        [item["F_bridge_10"] for item in candidates],
        c=[item["abs_Delta_Q_D3"] for item in candidates],
        cmap="viridis",
        s=52,
        alpha=0.82,
        edgecolors="black",
        linewidths=0.25,
    )
    for item in sorted(candidates, key=lambda obj: obj["F_bridge_10"])[:5]:
        ax.annotate(
            f"id={item['background_id']}\n|dQ|={item['abs_Delta_Q_D3']}\ntau={item['tau']:.1f}",
            (item["abs_energy"], item["F_bridge_10"]),
            textcoords="offset points",
            xytext=(6, 5),
            fontsize=7,
            bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.75},
        )
    ax.set_xscale("log")
    ax.set_xlabel(r"$|V_{\mathrm{eff}}-\Lambda_{\mathrm{obs}}|$")
    ax.set_ylabel(r"$F_{\mathrm{bridge}}^{(10)}$")
    ax.set_title("IIB/KKLT Bridge: Energy Proximity Versus 10D Structural Score")
    ax.grid(True, alpha=0.25)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(r"$|\Delta Q_{D3}|$")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_rank_bridge(candidates: list[dict], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 6.0))
    scatter = ax.scatter(
        [item["rank_energy"] for item in candidates],
        [item["rank_bridge"] for item in candidates],
        c=[item["tau"] for item in candidates],
        cmap="magma",
        s=54,
        alpha=0.82,
    )
    max_rank = len(candidates)
    ax.plot([1, max_rank], [1, max_rank], color="gray", linestyle="--", linewidth=1.0)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.set_xlabel("Rank by energy proximity")
    ax.set_ylabel(r"Rank by $F_{\mathrm{bridge}}^{(10)}$")
    ax.set_title("Bridge Ranking Versus Energy Ranking")
    ax.grid(True, alpha=0.25)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(r"KKLT modulus $\tau$")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_tadpole_tau(candidates: list[dict], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    scatter = ax.scatter(
        [item["abs_Delta_Q_D3"] for item in candidates],
        [item["tau"] for item in candidates],
        c=[item["F_bridge_10"] for item in candidates],
        cmap="cividis_r",
        s=55,
        alpha=0.84,
        edgecolors="black",
        linewidths=0.2,
    )
    ax.set_xlabel(r"$|\Delta Q_{D3}|$")
    ax.set_ylabel(r"KKLT volume modulus $\tau$")
    ax.set_title(r"Tadpole Closure and KKLT Volume Control")
    ax.grid(True, alpha=0.25)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(r"$F_{\mathrm{bridge}}^{(10)}$")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_bridge_decomposition(candidates: list[dict], output: Path, top_n: int = 8) -> None:
    top = sorted(candidates, key=lambda item: item["F_bridge_10"])[:top_n]
    labels = [f"id={item['background_id']}\n|dQ|={item['abs_Delta_Q_D3']}" for item in top]
    parts = [
        ("Delta_E_KKLT", r"$\Delta_E^{\mathrm{KKLT}}$"),
        ("Delta_Q_10", r"$\Delta_Q^{(10)}$"),
        ("Delta_dim_KKLT", r"$\Delta_{\mathrm{dim}}^{\mathrm{KKLT}}$"),
        ("tachyon_defect", r"$d_{\mathrm{tach}}$"),
    ]
    fig, ax = plt.subplots(figsize=(9.4, 5.2))
    bottom = np.zeros(len(top))
    colors = ["#386cb0", "#fdb462", "#7fc97f", "#ef3b2c"]
    for color, (key, label) in zip(colors, parts):
        values = np.array([item["bridge_defects"][key] for item in top], dtype=float)
        ax.bar(labels, values, bottom=bottom, color=color, label=label, alpha=0.9)
        bottom += values
    ax.set_ylabel("Bridge denominator contribution")
    ax.set_title("Top IIB/KKLT Bridge Candidates: Obstruction Decomposition")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="best", ncol=2)
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def make_plots(candidates: list[dict], output_dir: Path) -> list[str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    plots = {
        "energy_vs_bridge": output_dir / "energy_vs_bridge.png",
        "rank_bridge": output_dir / "rank_bridge.png",
        "tadpole_tau": output_dir / "tadpole_tau.png",
        "bridge_obstruction_decomposition": output_dir / "bridge_obstruction_decomposition.png",
    }
    plot_energy_bridge(candidates, plots["energy_vs_bridge"])
    plot_rank_bridge(candidates, plots["rank_bridge"])
    plot_tadpole_tau(candidates, plots["tadpole_tau"])
    plot_bridge_decomposition(candidates, plots["bridge_obstruction_decomposition"])
    return [str(path) for path in plots.values()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=2604)
    parser.add_argument("--n-backgrounds", type=int, default=90)
    parser.add_argument("--output", type=Path, default=Path("structural_selection_iib_kklt_bridge_results.json"))
    parser.add_argument("--plot-dir", type=Path, default=Path("structural_selection_iib_kklt_bridge_plots"))
    args = parser.parse_args()

    config = BridgeConfig(seed=args.seed, n_backgrounds=args.n_backgrounds)
    backgrounds, candidates = generate_bridge_scan(config)
    summary = summarize(backgrounds, candidates, config)
    plots = make_plots(candidates, args.plot_dir)
    payload = {
        "description": "Bridge scan linking flux/tadpole data to reduced KKLT stationary points.",
        "summary": summary,
        "plots": plots,
        "backgrounds": backgrounds,
        "candidates": candidates,
    }
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print("plots:")
    for path in plots:
        print(f"  {path}")


if __name__ == "__main__":
    main()
