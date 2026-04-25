#!/usr/bin/env python3
"""Backreaction and Standard-Model-sector bridge for structural selection.

This script extends the Picard-Fuchs exact-period-to-KKLT benchmark with two
additional phenomenological diagnostic layers:

1. A 10D-near backreaction proxy. Localized D3/O-like source data are placed on
   a small internal lattice and a Poisson equation for a warp-potential proxy is
   solved spectrally. The resulting warp energy, gradient scale, and net-source
   mismatch define an operational backreaction defect.

2. A Standard-Model-sector proxy. Integer topological data are mapped to gauge
   rank, chirality/generation count, hypercharge-integrity, and Yukawa-rank
   diagnostics. These are not a real SM compactification, but they provide the
   first explicit likelihood layer for SM-like data.
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

from structural_selection_iib_exact_period_kklt_bridge import (
    ExactPeriodConfig,
    exact_bridge_record,
    make_background,
    mirror_quintic_period_vector,
    rankdata,
    spearman,
)
from structural_selection_iib_kklt_scan import (
    DynamicParams,
    StructuralParams,
    candidate_record,
    solve_stationary_points,
)


@dataclass(frozen=True)
class BackreactionSMConfig:
    seed: int = 2904
    n_backgrounds: int = 90
    internal_grid: int = 12
    sm_sigma: float = 1.0
    backreaction_weight: float = 0.85
    sm_weight: float = 0.65
    exact_period: ExactPeriodConfig = ExactPeriodConfig(seed=2904, n_backgrounds=90)


def fft_poisson_defect(
    rng: np.random.Generator,
    n_d3: int,
    delta_q: int,
    flux_norm: float,
    config: BackreactionSMConfig,
) -> dict:
    """Solve a small toroidal Poisson problem for a warp-potential proxy."""
    n = config.internal_grid
    density = np.zeros((n, n, n), dtype=float)
    num_positive = int(min(max(n_d3, 1), 48))
    for _ in range(num_positive):
        idx = tuple(int(x) for x in rng.integers(0, n, size=3))
        density[idx] += 1.0

    negative_total = num_positive + delta_q
    num_negative = int(min(max(abs(negative_total), 1), 56))
    negative_charge = negative_total / num_negative if num_negative else 0.0
    for _ in range(num_negative):
        idx = tuple(int(x) for x in rng.integers(0, n, size=3))
        density[idx] -= negative_charge

    flux_density = 0.015 * flux_norm / (n**3)
    density += flux_density
    density -= float(np.mean(density))

    rho_hat = np.fft.fftn(density)
    k = np.fft.fftfreq(n) * 2.0 * math.pi
    kx, ky, kz = np.meshgrid(k, k, k, indexing="ij")
    k2 = kx**2 + ky**2 + kz**2
    phi_hat = np.zeros_like(rho_hat)
    mask = k2 > 0.0
    phi_hat[mask] = -rho_hat[mask] / k2[mask]
    phi = np.real(np.fft.ifftn(phi_hat))

    grad = np.gradient(phi)
    grad_sq = sum(component**2 for component in grad)
    warp_energy = float(np.mean(phi**2))
    grad_energy = float(np.mean(grad_sq))
    max_warp = float(np.max(np.abs(phi)))
    source_clump = float(np.std(density))
    net_mismatch = abs(delta_q) / 33.0
    defect = (
        0.55 * math.tanh(10.0 * warp_energy)
        + 0.30 * math.tanh(8.0 * grad_energy)
        + 0.45 * net_mismatch**2
        + 0.12 * math.tanh(source_clump / 2.0)
    )
    return {
        "warp_energy": warp_energy,
        "warp_gradient_energy": grad_energy,
        "max_warp_proxy": max_warp,
        "source_clumping": source_clump,
        "net_mismatch": net_mismatch,
        "backreaction_defect": float(defect),
    }


def sm_sector_proxy(background: dict, rng: np.random.Generator, config: BackreactionSMConfig) -> dict:
    """Generate deterministic-ish SM-like proxy diagnostics from flux/topology data."""
    f_flux = np.array(background["F_flux"], dtype=int)
    h_flux = np.array(background["H_flux"], dtype=int)
    flux_sum = int(np.sum(f_flux - h_flux))
    flux_pairing = int(np.dot(f_flux, h_flux))
    chirality_index = int((flux_sum + background["N_flux"]) % 9 - 4)
    generations = abs(chirality_index)
    gauge_rank = int(8 + ((abs(flux_pairing) + background["N_D3"]) % 10))
    hidden_rank = int(max(0, gauge_rank - 12))
    hypercharge_mass = abs((f_flux[0] + h_flux[1] - f_flux[2]) % 5 - 2) / 2.0
    yukawa_rank = int(min(3, abs(f_flux[1] - h_flux[2]) % 5))
    susy_break_proxy = float(abs(background["W0"]) * 1.0e4 + background["g_s"])

    d_rank = ((gauge_rank - 12) / 6.0) ** 2
    d_gen = ((generations - 3) / 3.0) ** 2
    d_hyper = hypercharge_mass**2
    d_yukawa = ((3 - yukawa_rank) / 3.0) ** 2
    d_hidden = (hidden_rank / 8.0) ** 2
    d_susy = max(0.0, susy_break_proxy - 3.2) ** 2 / 4.0
    sm_defect = d_rank + d_gen + d_hyper + d_yukawa + 0.35 * d_hidden + 0.15 * d_susy
    sm_likelihood = math.exp(-sm_defect / config.sm_sigma)
    return {
        "gauge_rank_proxy": gauge_rank,
        "hidden_rank_proxy": hidden_rank,
        "chirality_index_proxy": chirality_index,
        "generation_count_proxy": generations,
        "hypercharge_mass_proxy": float(hypercharge_mass),
        "yukawa_rank_proxy": yukawa_rank,
        "susy_break_proxy": susy_break_proxy,
        "sm_defect": float(sm_defect),
        "sm_likelihood": float(sm_likelihood),
    }


def enhanced_record(
    base_record: dict,
    backreaction: dict,
    sm_sector: dict,
    config: BackreactionSMConfig,
) -> dict:
    f_exact = float(base_record["F_exact_period"])
    lambda_penalty = float(base_record["lambda_penalty"])
    f_sw = float(base_record["F_SW_exact_period"])
    backreaction_defect = float(backreaction["backreaction_defect"])
    sm_defect = float(sm_sector["sm_defect"])
    f_backreaction_sm = f_exact + config.backreaction_weight * backreaction_defect + config.sm_weight * sm_defect
    f_pheno_full = f_backreaction_sm + f_sw + lambda_penalty
    return {
        **base_record,
        "backreaction": backreaction,
        "sm_sector": sm_sector,
        "F_backreaction_sm": f_backreaction_sm,
        "F_pheno_full": f_pheno_full,
    }


def generate_scan(config: BackreactionSMConfig) -> tuple[list[dict], list[dict]]:
    rng = np.random.default_rng(config.seed)
    period_config = config.exact_period
    structural = StructuralParams()
    tau_seeds = np.linspace(period_config.tau_seed_min, period_config.tau_seed_max, period_config.num_tau_seeds)
    theta_seeds = [0.0, math.pi / period_config.a]
    backgrounds: list[dict] = []
    candidates: list[dict] = []
    seen: set[tuple[int, ...]] = set()
    attempts = 0

    while len(backgrounds) < config.n_backgrounds and attempts < config.n_backgrounds * 160:
        attempts += 1
        background = make_background(len(backgrounds), rng, period_config)
        key = tuple(background["F_flux"] + background["H_flux"] + [round(background["z_abs"] * 1.0e7)])
        if key in seen:
            continue
        seen.add(key)
        if background["N_flux"] < -period_config.tadpole_bound or background["N_flux"] > 2 * period_config.tadpole_bound:
            continue

        params = DynamicParams(
            W0=background["W0"],
            A=period_config.A,
            a=period_config.a,
            D=background["D"],
            n=period_config.n_uplift,
        )
        points = solve_stationary_points(
            params=params,
            tau_seeds=tau_seeds,
            theta_seeds=theta_seeds,
            tau_min=period_config.tau_min_stationary,
            tau_max=period_config.tau_max_stationary,
            grad_tol=period_config.grad_tol,
        )
        background = {**background, "stationary_points": [{"tau": tau, "theta": theta} for tau, theta in points]}
        backgrounds.append(background)

        backreaction = fft_poisson_defect(
            rng=rng,
            n_d3=background["N_D3"],
            delta_q=background["Delta_Q_D3"],
            flux_norm=background["flux_norm"],
            config=config,
        )
        sm_sector = sm_sector_proxy(background, rng, config)

        for tau, theta in points:
            kklt = candidate_record(
                params=params,
                structural=structural,
                delta_q_d3=background["Delta_Q_D3"],
                tau=tau,
                theta=theta,
                tau_barrier_max=period_config.tau_barrier_max,
                barrier_samples=period_config.barrier_samples,
            )
            base = exact_bridge_record(background, kklt, period_config)
            candidates.append(enhanced_record(base, backreaction, sm_sector, config))

    if not candidates:
        raise RuntimeError("No backreaction/SM bridge candidates found.")

    ranks = {
        "rank_energy": rankdata(item["abs_energy"] for item in candidates),
        "rank_exact_period": rankdata(item["F_exact_period"] for item in candidates),
        "rank_backreaction_sm": rankdata(item["F_backreaction_sm"] for item in candidates),
        "rank_pheno_full": rankdata(item["F_pheno_full"] for item in candidates),
    }
    for idx, item in enumerate(candidates):
        for name, values in ranks.items():
            item[name] = int(values[idx])
    return backgrounds, candidates


def top_summary(candidates: list[dict], key: str, top_n: int = 12) -> dict:
    top = sorted(candidates, key=lambda item: item[key])[:top_n]
    return {
        "score": key,
        "top_n": top_n,
        "mean_abs_delta_q": float(np.mean([item["abs_Delta_Q_D3"] for item in top])),
        "mean_abs_energy": float(np.mean([item["abs_energy"] for item in top])),
        "mean_backreaction_defect": float(np.mean([item["backreaction"]["backreaction_defect"] for item in top])),
        "mean_sm_defect": float(np.mean([item["sm_sector"]["sm_defect"] for item in top])),
        "mean_generation_count": float(np.mean([item["sm_sector"]["generation_count_proxy"] for item in top])),
        "mean_hypercharge_mass": float(np.mean([item["sm_sector"]["hypercharge_mass_proxy"] for item in top])),
        "mean_tau": float(np.mean([item["tau"] for item in top])),
        "ids": [item["background_id"] for item in top],
    }


def summarize(backgrounds: list[dict], candidates: list[dict], config: BackreactionSMConfig) -> dict:
    top_energy = set(top_summary(candidates, "abs_energy")["ids"])
    top_full = set(top_summary(candidates, "F_pheno_full")["ids"])
    class_counts: dict[str, int] = {}
    for item in candidates:
        class_counts[item["vacuum_type"]] = class_counts.get(item["vacuum_type"], 0) + 1
    return {
        "config": {
            **asdict(config),
            "exact_period": asdict(config.exact_period),
        },
        "n_backgrounds": len(backgrounds),
        "n_stationary_candidates": len(candidates),
        "class_counts": class_counts,
        "spearman_abs_energy_vs_F_pheno_full": spearman(
            [item["abs_energy"] for item in candidates],
            [item["F_pheno_full"] for item in candidates],
        ),
        "spearman_backreaction_vs_F_pheno_full": spearman(
            [item["backreaction"]["backreaction_defect"] for item in candidates],
            [item["F_pheno_full"] for item in candidates],
        ),
        "spearman_sm_defect_vs_F_pheno_full": spearman(
            [item["sm_sector"]["sm_defect"] for item in candidates],
            [item["F_pheno_full"] for item in candidates],
        ),
        "spearman_abs_delta_q_vs_F_pheno_full": spearman(
            [item["abs_Delta_Q_D3"] for item in candidates],
            [item["F_pheno_full"] for item in candidates],
        ),
        "top_energy": top_summary(candidates, "abs_energy"),
        "top_exact_period": top_summary(candidates, "F_exact_period"),
        "top_backreaction_sm": top_summary(candidates, "F_backreaction_sm"),
        "top_pheno_full": top_summary(candidates, "F_pheno_full"),
        "top12_overlap_energy_full": len(top_energy & top_full),
        "best_energy_candidate": min(candidates, key=lambda item: item["abs_energy"]),
        "best_full_candidate": min(candidates, key=lambda item: item["F_pheno_full"]),
    }


def plot_full_energy(candidates: list[dict], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    scatter = ax.scatter(
        [item["abs_energy"] for item in candidates],
        [item["F_pheno_full"] for item in candidates],
        c=[item["sm_sector"]["sm_defect"] for item in candidates],
        cmap="viridis",
        s=50,
        alpha=0.82,
        edgecolors="black",
        linewidths=0.2,
    )
    ax.set_xscale("log")
    ax.set_xlabel(r"$|V_{\mathrm{eff}}-\Lambda_{\mathrm{obs}}|$")
    ax.set_ylabel(r"$F_{\mathrm{full}}^{(10+\mathrm{SM})}$")
    ax.set_title("Backreaction+SM Full Score Versus Energy Proximity")
    ax.grid(True, alpha=0.25)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("SM-sector defect")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_backreaction_sm(candidates: list[dict], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 5.3))
    scatter = ax.scatter(
        [item["backreaction"]["backreaction_defect"] for item in candidates],
        [item["sm_sector"]["sm_defect"] for item in candidates],
        c=[item["F_pheno_full"] for item in candidates],
        cmap="cividis_r",
        s=50,
        alpha=0.84,
    )
    ax.set_xlabel("Backreaction defect")
    ax.set_ylabel("SM-sector defect")
    ax.set_title("Backreaction and SM-Likeness Define Additional Selection Axes")
    ax.grid(True, alpha=0.25)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(r"$F_{\mathrm{full}}^{(10+\mathrm{SM})}$")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_rank_full(candidates: list[dict], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 6.0))
    scatter = ax.scatter(
        [item["rank_energy"] for item in candidates],
        [item["rank_pheno_full"] for item in candidates],
        c=[item["backreaction"]["backreaction_defect"] for item in candidates],
        cmap="magma",
        s=52,
        alpha=0.82,
    )
    max_rank = len(candidates)
    ax.plot([1, max_rank], [1, max_rank], color="gray", linestyle="--", linewidth=1.0)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.set_xlabel("Rank by energy proximity")
    ax.set_ylabel(r"Rank by full score")
    ax.set_title("Full Backreaction+SM Ranking Versus Energy Ranking")
    ax.grid(True, alpha=0.25)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("Backreaction defect")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_full_decomposition(candidates: list[dict], output: Path, top_n: int = 8) -> None:
    top = sorted(candidates, key=lambda item: item["F_pheno_full"])[:top_n]
    labels = [f"id={item['background_id']}\ng={item['sm_sector']['generation_count_proxy']}" for item in top]
    parts = [
        ("exact", r"$F_{\mathrm{exact}}$"),
        ("backreaction", r"$d_{\mathrm{back}}$"),
        ("sm", r"$d_{\mathrm{SM}}$"),
        ("lambda", r"$\Lambda$ penalty"),
    ]
    values_by_part = {
        "exact": [item["F_exact_period"] for item in top],
        "backreaction": [item["backreaction"]["backreaction_defect"] for item in top],
        "sm": [item["sm_sector"]["sm_defect"] for item in top],
        "lambda": [item["lambda_penalty"] for item in top],
    }
    fig, ax = plt.subplots(figsize=(9.4, 5.3))
    x = np.arange(len(top))
    width = 0.2
    colors = ["#386cb0", "#fdb462", "#7fc97f", "#ef3b2c"]
    for offset, color, (key, label) in zip(np.linspace(-0.3, 0.3, len(parts)), colors, parts):
        ax.bar(x + offset, values_by_part[key], width=width, color=color, label=label, alpha=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Component value")
    ax.set_title("Top Full Candidates: Component Diagnostics")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="best", ncol=2)
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def make_plots(candidates: list[dict], output_dir: Path) -> list[str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    plots = {
        "full_energy": output_dir / "full_energy.png",
        "backreaction_sm_plane": output_dir / "backreaction_sm_plane.png",
        "rank_full": output_dir / "rank_full.png",
        "full_component_diagnostics": output_dir / "full_component_diagnostics.png",
    }
    plot_full_energy(candidates, plots["full_energy"])
    plot_backreaction_sm(candidates, plots["backreaction_sm_plane"])
    plot_rank_full(candidates, plots["rank_full"])
    plot_full_decomposition(candidates, plots["full_component_diagnostics"])
    return [str(path) for path in plots.values()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=2904)
    parser.add_argument("--n-backgrounds", type=int, default=90)
    parser.add_argument("--output", type=Path, default=Path("structural_selection_iib_backreaction_sm_bridge_results.json"))
    parser.add_argument("--plot-dir", type=Path, default=Path("structural_selection_iib_backreaction_sm_bridge_plots"))
    args = parser.parse_args()

    exact_config = ExactPeriodConfig(seed=args.seed, n_backgrounds=args.n_backgrounds)
    config = BackreactionSMConfig(seed=args.seed, n_backgrounds=args.n_backgrounds, exact_period=exact_config)
    backgrounds, candidates = generate_scan(config)
    summary = summarize(backgrounds, candidates, config)
    plots = make_plots(candidates, args.plot_dir)
    payload = {
        "description": "Exact-period-to-KKLT bridge augmented with backreaction and SM-sector proxies.",
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
