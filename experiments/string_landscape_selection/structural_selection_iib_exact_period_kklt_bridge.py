#!/usr/bin/env python3
"""Mirror-quintic Picard-Fuchs period bridge for structural selection.

This benchmark replaces the polynomial LCS period approximation by numerical
Frobenius periods of the mirror-quintic Picard-Fuchs system around z=0:

    theta^4 Pi - 5 z (5 theta + 1)(5 theta + 2)(5 theta + 3)(5 theta + 4) Pi = 0.

The periods are evaluated from the truncated hypergeometric/Frobenius series.
They are then used to compute W0_raw = (F - tau_ax H).Pi(z), which is passed to
the reduced KKLT stationary-point solver and structural-selection pipeline.
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
import mpmath as mp
import numpy as np

from structural_selection_iib_kklt_scan import (
    DynamicParams,
    StructuralParams,
    candidate_record,
    solve_stationary_points,
)


@dataclass(frozen=True)
class ExactPeriodConfig:
    seed: int = 2804
    n_backgrounds: int = 110
    flux_dim: int = 4
    flux_abs_max: int = 6
    tadpole_bound: int = 32
    z_abs_min: float = 1.0e-6
    z_abs_max: float = 1.8e-4
    z_arg_max: float = 0.75
    gs_min: float = 0.08
    gs_max: float = 0.34
    period_terms: int = 42
    mp_dps: int = 50
    w0_min: float = 5.0e-5
    w0_max: float = 3.5e-4
    w0_raw_scale: float = 65.0
    A: float = 1.0
    a: float = 0.1
    n_uplift: float = 2.0
    d_uplift_max: float = 1.2e-9
    tau_seed_min: float = 18.0
    tau_seed_max: float = 185.0
    num_tau_seeds: int = 20
    tau_min_stationary: float = 5.0
    tau_max_stationary: float = 185.0
    tau_barrier_max: float = 245.0
    barrier_samples: int = 260
    grad_tol: float = 1.0e-8
    epsilon: float = 1.0e-9
    lambda_obs: float = 0.0
    lambda_sigma: float = 5.0e-15


def mirror_quintic_period_vector(z_value: complex, config: ExactPeriodConfig) -> np.ndarray:
    """Return four Frobenius periods near the mirror-quintic LCS point.

    The basis is local and not symplectically integral. It is sufficient for a
    controlled bridge test because the same basis is used for all candidates.
    """
    mp.mp.dps = config.mp_dps
    z = mp.mpc(z_value)
    log_z = mp.log(z)

    sums = [mp.mpc(0), mp.mpc(0), mp.mpc(0), mp.mpc(0)]
    for n in range(config.period_terms):
        n_mp = mp.mpf(n)
        coeff = mp.factorial(5 * n) / (mp.factorial(n) ** 5)
        term = coeff * (z**n)

        l1 = 5 * mp.digamma(5 * n_mp + 1) - 5 * mp.digamma(n_mp + 1) + log_z
        l2 = 25 * mp.polygamma(1, 5 * n_mp + 1) - 5 * mp.polygamma(1, n_mp + 1)
        l3 = 125 * mp.polygamma(2, 5 * n_mp + 1) - 5 * mp.polygamma(2, n_mp + 1)

        sums[0] += term
        sums[1] += term * l1
        sums[2] += term * (l1**2 + l2)
        sums[3] += term * (l1**3 + 3 * l1 * l2 + l3)

    two_pi_i = 2 * mp.pi * 1j
    periods = [
        sums[0],
        sums[1] / two_pi_i,
        sums[2] / (2 * two_pi_i**2),
        sums[3] / (6 * two_pi_i**3),
    ]
    return np.array([complex(item) for item in periods], dtype=complex)


def symplectic_tadpole(f_flux: np.ndarray, h_flux: np.ndarray) -> int:
    return int(f_flux[0] * h_flux[2] + f_flux[1] * h_flux[3] - f_flux[2] * h_flux[0] - f_flux[3] * h_flux[1])


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


def w0_from_exact_periods(
    f_flux: np.ndarray,
    h_flux: np.ndarray,
    periods: np.ndarray,
    tau_ax: complex,
    config: ExactPeriodConfig,
) -> tuple[complex, float]:
    w_raw = complex(np.dot(f_flux.astype(complex) - tau_ax * h_flux.astype(complex), periods))
    fraction = abs(w_raw) / (abs(w_raw) + config.w0_raw_scale)
    w0_abs = config.w0_min + (config.w0_max - config.w0_min) * fraction
    return w_raw, -float(w0_abs)


def make_background(idx: int, rng: np.random.Generator, config: ExactPeriodConfig) -> dict:
    f_flux = rng.integers(-config.flux_abs_max, config.flux_abs_max + 1, size=config.flux_dim)
    h_flux = rng.integers(-config.flux_abs_max, config.flux_abs_max + 1, size=config.flux_dim)
    while not np.any(f_flux) or not np.any(h_flux):
        f_flux = rng.integers(-config.flux_abs_max, config.flux_abs_max + 1, size=config.flux_dim)
        h_flux = rng.integers(-config.flux_abs_max, config.flux_abs_max + 1, size=config.flux_dim)

    log_abs_z = rng.uniform(math.log(config.z_abs_min), math.log(config.z_abs_max))
    z_abs = float(math.exp(log_abs_z))
    z_arg = float(rng.uniform(-config.z_arg_max, config.z_arg_max))
    z_value = z_abs * complex(math.cos(z_arg), math.sin(z_arg))
    periods = mirror_quintic_period_vector(z_value, config)

    c0 = float(rng.uniform(-0.5, 0.5))
    gs = float(rng.uniform(config.gs_min, config.gs_max))
    tau_ax = complex(c0, 1.0 / gs)
    w_raw, w0 = w0_from_exact_periods(f_flux, h_flux, periods, tau_ax, config)

    n_flux = symplectic_tadpole(f_flux, h_flux)
    n_d3_ideal = config.tadpole_bound - n_flux
    n_d3 = max(0, int(round(n_d3_ideal + int(rng.integers(-5, 6)))))
    delta_q = int(n_flux + n_d3 - config.tadpole_bound)

    flux_norm = float(np.dot(f_flux, f_flux) + np.dot(h_flux, h_flux))
    max_norm = float(2 * config.flux_dim * config.flux_abs_max**2)
    d_uplift = config.d_uplift_max * (0.32 + 0.68 * flux_norm / max_norm) * (gs / config.gs_max)

    lcs_distance = float(abs(5**5 * z_abs))
    period_control = lcs_distance**2 + gs**2 + (0.18 / max(abs(periods[0]), 1.0e-8)) ** 2

    return {
        "background_id": idx,
        "F_flux": [int(x) for x in f_flux.tolist()],
        "H_flux": [int(x) for x in h_flux.tolist()],
        "N_flux": n_flux,
        "N_D3": n_d3,
        "Delta_Q_D3": delta_q,
        "abs_Delta_Q_D3": abs(delta_q),
        "flux_norm": flux_norm,
        "z_re": float(z_value.real),
        "z_im": float(z_value.imag),
        "z_abs": z_abs,
        "z_arg": z_arg,
        "lcs_distance_5p5z": lcs_distance,
        "g_s": gs,
        "C0": c0,
        "tau_ax_real": c0,
        "tau_ax_imag": 1.0 / gs,
        "periods_re": [float(x.real) for x in periods],
        "periods_im": [float(x.imag) for x in periods],
        "W0_raw_re": float(w_raw.real),
        "W0_raw_im": float(w_raw.imag),
        "W0_raw_abs": float(abs(w_raw)),
        "W0": w0,
        "D": float(d_uplift),
        "period_control_defect": period_control,
    }


def exact_bridge_record(background: dict, kklt_record: dict, config: ExactPeriodConfig) -> dict:
    tau = float(kklt_record["tau"])
    energy = float(kklt_record["energy"])
    lambda_min = float(kklt_record["lambda_min"])
    scores = kklt_record["scores"]

    delta_q_10 = background["abs_Delta_Q_D3"] / (config.tadpole_bound + 1.0)
    d_b10 = 20.0 * delta_q_10**2
    gamma_b10 = support(d_b10)
    energy_scale = max(abs(background["W0"]) ** 2, abs(energy), 1.0e-18)
    delta_e = abs(energy - config.lambda_obs) / energy_scale
    delta_dim = (10.0 / tau) ** 2 + math.exp(-config.a * tau) + 0.4 * background["period_control_defect"]
    tachyon_defect = max(0.0, -lambda_min / max(abs(lambda_min), max(abs(x) for x in kklt_record["hessian_eigenvalues"]), 1.0e-18))

    numerator = float(np.mean([scores["H"], gamma_b10, scores["S"], scores["V"], scores["R"]]))
    exact_support = numerator / (
        config.epsilon + delta_e + delta_q_10 + delta_dim + 0.2 * tachyon_defect
    )
    f_exact = -math.log(config.epsilon + exact_support)
    f_sw = 12.0 * delta_q_10**2 + 4.0 * max(0.0, background["lcs_distance_5p5z"] - 0.85) ** 2 + 2.0 * tachyon_defect
    lambda_penalty = 0.5 * ((energy - config.lambda_obs) / config.lambda_sigma) ** 2
    f_pheno = f_exact + f_sw + lambda_penalty

    return {
        **background,
        "tau": tau,
        "theta": kklt_record["theta"],
        "energy": energy,
        "abs_energy": abs(energy - config.lambda_obs),
        "vacuum_type": kklt_record["vacuum_type"],
        "lambda_min": lambda_min,
        "hessian_eigenvalues": kklt_record["hessian_eigenvalues"],
        "kklt_f_master": kklt_record["f_master"],
        "kklt_scores": scores,
        "kklt_residuals": kklt_record["residuals"],
        "exact_period_defects": {
            "d_B10": d_b10,
            "Delta_Q_10": delta_q_10,
            "Delta_E_KKLT": delta_e,
            "Delta_dim_exact_period_KKLT": delta_dim,
            "period_control_defect": background["period_control_defect"],
            "tachyon_defect": tachyon_defect,
        },
        "exact_period_support": exact_support,
        "F_exact_period": f_exact,
        "F_SW_exact_period": f_sw,
        "lambda_penalty": lambda_penalty,
        "F_pheno_exact_period": f_pheno,
    }


def generate_scan(config: ExactPeriodConfig) -> tuple[list[dict], list[dict]]:
    rng = np.random.default_rng(config.seed)
    structural = StructuralParams()
    tau_seeds = np.linspace(config.tau_seed_min, config.tau_seed_max, config.num_tau_seeds)
    theta_seeds = [0.0, math.pi / config.a]
    backgrounds: list[dict] = []
    candidates: list[dict] = []
    seen: set[tuple[int, ...]] = set()
    attempts = 0

    while len(backgrounds) < config.n_backgrounds and attempts < config.n_backgrounds * 140:
        attempts += 1
        background = make_background(len(backgrounds), rng, config)
        key = tuple(background["F_flux"] + background["H_flux"] + [round(background["z_abs"] * 1.0e7)])
        if key in seen:
            continue
        seen.add(key)
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
            candidates.append(exact_bridge_record(background, kklt, config))

    if not candidates:
        raise RuntimeError("No exact-period-to-KKLT stationary candidates found.")

    ranks = {
        "rank_energy": rankdata(item["abs_energy"] for item in candidates),
        "rank_exact_period": rankdata(item["F_exact_period"] for item in candidates),
        "rank_pheno_exact_period": rankdata(item["F_pheno_exact_period"] for item in candidates),
        "rank_kklt": rankdata(item["kklt_f_master"] for item in candidates),
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
        "median_abs_delta_q": float(np.median([item["abs_Delta_Q_D3"] for item in top])),
        "mean_abs_energy": float(np.mean([item["abs_energy"] for item in top])),
        "mean_tau": float(np.mean([item["tau"] for item in top])),
        "mean_z_abs": float(np.mean([item["z_abs"] for item in top])),
        "mean_g_s": float(np.mean([item["g_s"] for item in top])),
        "mean_W0_abs": float(np.mean([abs(item["W0"]) for item in top])),
        "ids": [item["background_id"] for item in top],
    }


def summarize(backgrounds: list[dict], candidates: list[dict], config: ExactPeriodConfig) -> dict:
    top_energy = set(top_summary(candidates, "abs_energy")["ids"])
    top_exact = set(top_summary(candidates, "F_exact_period")["ids"])
    top_pheno = set(top_summary(candidates, "F_pheno_exact_period")["ids"])
    class_counts: dict[str, int] = {}
    for item in candidates:
        class_counts[item["vacuum_type"]] = class_counts.get(item["vacuum_type"], 0) + 1
    return {
        "config": asdict(config),
        "n_backgrounds": len(backgrounds),
        "n_stationary_candidates": len(candidates),
        "class_counts": class_counts,
        "spearman_abs_energy_vs_F_exact_period": spearman(
            [item["abs_energy"] for item in candidates],
            [item["F_exact_period"] for item in candidates],
        ),
        "spearman_abs_delta_q_vs_abs_energy": spearman(
            [item["abs_Delta_Q_D3"] for item in candidates],
            [item["abs_energy"] for item in candidates],
        ),
        "spearman_abs_delta_q_vs_F_exact_period": spearman(
            [item["abs_Delta_Q_D3"] for item in candidates],
            [item["F_exact_period"] for item in candidates],
        ),
        "spearman_period_control_vs_F_exact_period": spearman(
            [item["period_control_defect"] for item in candidates],
            [item["F_exact_period"] for item in candidates],
        ),
        "top_energy": top_summary(candidates, "abs_energy"),
        "top_kklt": top_summary(candidates, "kklt_f_master"),
        "top_exact_period": top_summary(candidates, "F_exact_period"),
        "top_pheno_exact_period": top_summary(candidates, "F_pheno_exact_period"),
        "top12_overlap_energy_exact_period": len(top_energy & top_exact),
        "top12_overlap_energy_pheno_exact": len(top_energy & top_pheno),
        "best_energy_candidate": min(candidates, key=lambda item: item["abs_energy"]),
        "best_exact_period_candidate": min(candidates, key=lambda item: item["F_exact_period"]),
        "best_pheno_exact_period_candidate": min(candidates, key=lambda item: item["F_pheno_exact_period"]),
    }


def plot_energy_exact(candidates: list[dict], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    scatter = ax.scatter(
        [item["abs_energy"] for item in candidates],
        [item["F_exact_period"] for item in candidates],
        c=[item["abs_Delta_Q_D3"] for item in candidates],
        cmap="viridis",
        s=50,
        alpha=0.82,
        edgecolors="black",
        linewidths=0.2,
    )
    for item in sorted(candidates, key=lambda obj: obj["F_exact_period"])[:5]:
        ax.annotate(
            f"id={item['background_id']}\n|dQ|={item['abs_Delta_Q_D3']}\n|z|={item['z_abs']:.1e}",
            (item["abs_energy"], item["F_exact_period"]),
            textcoords="offset points",
            xytext=(6, 5),
            fontsize=7,
            bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.75},
        )
    ax.set_xscale("log")
    ax.set_xlabel(r"$|V_{\mathrm{eff}}-\Lambda_{\mathrm{obs}}|$")
    ax.set_ylabel(r"$F_{\mathrm{exact}}^{(10)}$")
    ax.set_title("Exact-Period-to-KKLT Bridge: Energy Versus Structural Score")
    ax.grid(True, alpha=0.25)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(r"$|\Delta Q_{D3}|$")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_rank_exact(candidates: list[dict], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 6.0))
    scatter = ax.scatter(
        [item["rank_energy"] for item in candidates],
        [item["rank_exact_period"] for item in candidates],
        c=[item["z_abs"] for item in candidates],
        cmap="magma",
        s=50,
        alpha=0.82,
    )
    max_rank = len(candidates)
    ax.plot([1, max_rank], [1, max_rank], color="gray", linestyle="--", linewidth=1.0)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.set_xlabel("Rank by energy proximity")
    ax.set_ylabel(r"Rank by $F_{\mathrm{exact}}^{(10)}$")
    ax.set_title("Exact-Period Ranking Versus Energy Ranking")
    ax.grid(True, alpha=0.25)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(r"$|z|$")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_w0_exact(candidates: list[dict], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    scatter = ax.scatter(
        [item["W0_raw_abs"] for item in candidates],
        [abs(item["W0"]) for item in candidates],
        c=[item["lcs_distance_5p5z"] for item in candidates],
        cmap="cividis",
        s=48,
        alpha=0.84,
    )
    ax.set_xlabel(r"$|W_{0,\mathrm{raw}}|=|(F-\tau_{\mathrm{ax}}H)\cdot\Pi_{\mathrm{PF}}(z)|$")
    ax.set_ylabel(r"Reduced KKLT $|W_0|$")
    ax.set_title("Picard-Fuchs Periods Generate the Effective Flux Superpotential")
    ax.grid(True, alpha=0.25)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(r"$|5^5 z|$")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_exact_decomposition(candidates: list[dict], output: Path, top_n: int = 8) -> None:
    top = sorted(candidates, key=lambda item: item["F_exact_period"])[:top_n]
    labels = [f"id={item['background_id']}\n|dQ|={item['abs_Delta_Q_D3']}" for item in top]
    parts = [
        ("Delta_E_KKLT", r"$\Delta_E^{\mathrm{KKLT}}$"),
        ("Delta_Q_10", r"$\Delta_Q^{(10)}$"),
        ("Delta_dim_exact_period_KKLT", r"$\Delta_{\mathrm{dim+PF}}$"),
        ("tachyon_defect", r"$d_{\mathrm{tach}}$"),
    ]
    fig, ax = plt.subplots(figsize=(9.4, 5.2))
    bottom = np.zeros(len(top))
    colors = ["#386cb0", "#fdb462", "#7fc97f", "#ef3b2c"]
    for color, (key, label) in zip(colors, parts):
        values = np.array([item["exact_period_defects"][key] for item in top], dtype=float)
        ax.bar(labels, values, bottom=bottom, color=color, label=label, alpha=0.9)
        bottom += values
    ax.set_ylabel("Exact-period denominator contribution")
    ax.set_title("Top Exact-Period-to-KKLT Candidates: Obstruction Decomposition")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="best", ncol=2)
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def make_plots(candidates: list[dict], output_dir: Path) -> list[str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    plots = {
        "energy_vs_exact_period": output_dir / "energy_vs_exact_period.png",
        "rank_exact_period": output_dir / "rank_exact_period.png",
        "w0_exact_period_map": output_dir / "w0_exact_period_map.png",
        "exact_period_obstruction_decomposition": output_dir / "exact_period_obstruction_decomposition.png",
    }
    plot_energy_exact(candidates, plots["energy_vs_exact_period"])
    plot_rank_exact(candidates, plots["rank_exact_period"])
    plot_w0_exact(candidates, plots["w0_exact_period_map"])
    plot_exact_decomposition(candidates, plots["exact_period_obstruction_decomposition"])
    return [str(path) for path in plots.values()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=2804)
    parser.add_argument("--n-backgrounds", type=int, default=110)
    parser.add_argument("--output", type=Path, default=Path("structural_selection_iib_exact_period_kklt_bridge_results.json"))
    parser.add_argument("--plot-dir", type=Path, default=Path("structural_selection_iib_exact_period_kklt_bridge_plots"))
    args = parser.parse_args()

    config = ExactPeriodConfig(seed=args.seed, n_backgrounds=args.n_backgrounds)
    backgrounds, candidates = generate_scan(config)
    summary = summarize(backgrounds, candidates, config)
    plots = make_plots(candidates, args.plot_dir)
    payload = {
        "description": "Mirror-quintic Picard-Fuchs exact-period-to-KKLT bridge scan.",
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
