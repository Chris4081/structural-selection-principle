#!/usr/bin/env python3
"""MAAT structural-selection toy benchmark for FLRW scalar-field histories.

The script implements a compact flat-FLRW scalar-field toy model and evaluates
a five-sector structural score on several histories:

1. slow-roll coherent expansion,
2. kinetic relaxation,
3. weakly generative static expansion,
4. collapse from a negative potential,
5. ghost-like accelerated expansion.

The diagnostics are operational proxies, not a complete cosmological model.
The purpose is to test whether the same structural-selection architecture used
in static field benchmarks can rank dynamical cosmological histories.
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/codex-mpl-cache")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/codex-mpl-cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


EPS = 1.0e-8


@dataclass(frozen=True)
class CosmoConfig:
    dt: float = 0.005
    total_time: float = 60.0
    v0: float = 1.0
    delta: float = 3.0


def v_plateau(phi: np.ndarray | float, config: CosmoConfig) -> np.ndarray | float:
    return config.v0 * (1.0 - np.tanh(phi / config.delta) ** 2) + 1.0e-6


def dv_plateau(phi: np.ndarray | float, config: CosmoConfig) -> np.ndarray | float:
    z = phi / config.delta
    return config.v0 * (-2.0 / config.delta) * np.tanh(z) * (1.0 - np.tanh(z) ** 2)


def v_negative(phi: np.ndarray | float, _: CosmoConfig) -> np.ndarray | float:
    return 0.2 - 0.08 * phi**2


def dv_negative(phi: np.ndarray | float, _: CosmoConfig) -> np.ndarray | float:
    return -0.16 * phi


def evolve(
    phi0: float,
    pi0: float,
    label: str,
    config: CosmoConfig,
    potential: str = "plateau",
    ghost: bool = False,
) -> dict:
    steps = int(config.total_time / config.dt)
    phi = np.zeros(steps)
    pi = np.zeros(steps)
    hubble = np.zeros(steps)
    rho = np.zeros(steps)
    w = np.zeros(steps)
    eps_h = np.zeros(steps)
    scale_factor = np.ones(steps)
    valid = np.ones(steps)

    phi[0] = phi0
    pi[0] = pi0

    if potential == "plateau":
        pot = lambda p: v_plateau(p, config)
        dpot = lambda p: dv_plateau(p, config)
    elif potential == "negative":
        pot = lambda p: v_negative(p, config)
        dpot = lambda p: dv_negative(p, config)
    else:
        raise ValueError(f"unknown potential: {potential}")

    def rhs(p: float, q: float) -> tuple[float, float, float, float, float]:
        kinetic_sign = -1.0 if ghost else 1.0
        density = 0.5 * kinetic_sign * q**2 + pot(p)
        if density <= 0.0:
            h_val = 0.0
            ok = 0.0
        else:
            h_val = float(np.sqrt(density / 3.0))
            ok = 1.0
        dq = -3.0 * h_val * q - kinetic_sign * dpot(p)
        dp = q
        return dp, dq, h_val, float(density), ok

    for idx in range(steps - 1):
        p = float(phi[idx])
        q = float(pi[idx])

        k1p, k1q, _, _, _ = rhs(p, q)
        k2p, k2q, _, _, _ = rhs(p + 0.5 * config.dt * k1p, q + 0.5 * config.dt * k1q)
        k3p, k3q, _, _, _ = rhs(p + 0.5 * config.dt * k2p, q + 0.5 * config.dt * k2q)
        k4p, k4q, _, _, _ = rhs(p + config.dt * k3p, q + config.dt * k3q)

        phi[idx + 1] = p + (config.dt / 6.0) * (k1p + 2.0 * k2p + 2.0 * k3p + k4p)
        pi[idx + 1] = q + (config.dt / 6.0) * (k1q + 2.0 * k2q + 2.0 * k3q + k4q)

        _, _, hubble[idx], rho[idx], valid[idx] = rhs(float(phi[idx]), float(pi[idx]))

        kinetic_sign = -1.0 if ghost else 1.0
        kinetic = 0.5 * kinetic_sign * pi[idx] ** 2
        potential_energy = pot(float(phi[idx]))
        if abs(kinetic + potential_energy) < EPS:
            w[idx] = 999.0
        else:
            w[idx] = (kinetic - potential_energy) / (kinetic + potential_energy + EPS)
        eps_h[idx] = 1.5 * (1.0 + w[idx])

        if valid[idx] > 0.0:
            scale_factor[idx + 1] = scale_factor[idx] * np.exp(hubble[idx] * config.dt)
        else:
            scale_factor[idx + 1] = scale_factor[idx]

    _, _, hubble[-1], rho[-1], valid[-1] = rhs(float(phi[-1]), float(pi[-1]))
    kinetic_sign = -1.0 if ghost else 1.0
    kinetic = 0.5 * kinetic_sign * pi[-1] ** 2
    potential_energy = pot(float(phi[-1]))
    if abs(kinetic + potential_energy) < EPS:
        w[-1] = 999.0
    else:
        w[-1] = (kinetic - potential_energy) / (kinetic + potential_energy + EPS)
    eps_h[-1] = 1.5 * (1.0 + w[-1])

    return {
        "label": label,
        "phi": phi,
        "pi": pi,
        "H": hubble,
        "rho": rho,
        "w": w,
        "epsH": eps_h,
        "a": scale_factor,
        "valid": valid,
        "ghost": ghost,
        "potential": potential,
    }


def maat_cosmo_score(run: dict) -> dict:
    hubble = run["H"]
    rho = run["rho"]
    w = run["w"]
    eps_h = run["epsH"]
    valid = run["valid"]
    phi = run["phi"]
    pi = run["pi"]

    friedmann_res = hubble**2 - np.maximum(rho, 0.0) / 3.0
    invalid_fraction = float(np.mean(valid < 0.5))
    d_h = float(np.mean(friedmann_res**2) / (np.mean(hubble**4) + EPS) + 100.0 * invalid_fraction)
    h_score = 1.0 / (1.0 + d_h)

    negative_density = float(np.mean(np.maximum(0.0, -rho) ** 2))
    w_bad = float(np.mean(np.maximum(0.0, np.abs(w) - 3.0) ** 2))
    d_b = negative_density + 0.01 * w_bad
    b_score = 1.0 / (1.0 + d_b)

    inflation_fraction = float(np.mean((eps_h < 1.0) & (valid > 0.5)))
    efolds = float(np.log((run["a"][-1] + EPS) / (run["a"][0] + EPS)))
    efold_support = efolds / (20.0 + efolds) if efolds > 0.0 else 0.0
    activity = float(np.var(phi) / (1.0 + np.var(phi)))
    s_score = 0.45 * inflation_fraction + 0.35 * efold_support + 0.20 * activity

    d_hubble = np.diff(hubble)
    roughness_h = float(np.mean(d_hubble**2) / (np.mean(hubble**2) + EPS))
    d_pi = np.diff(pi)
    roughness_pi = float(np.mean(d_pi**2) / (np.mean(pi**2) + EPS))
    d_v = roughness_h + 0.2 * roughness_pi
    v_score = 1.0 / (1.0 + d_v)

    kinetic_bad = float(np.mean(np.maximum(0.0, w - 0.0) ** 2))
    eps_bad = float(np.mean(np.maximum(0.0, eps_h - 1.0) ** 2))
    ghost_penalty = 10.0 if run["ghost"] else 0.0
    collapse_penalty = 20.0 * invalid_fraction
    d_r = kinetic_bad + eps_bad + ghost_penalty + collapse_penalty
    r_score = 1.0 / (1.0 + d_r)

    scores = np.array([h_score, b_score, max(s_score, EPS), max(v_score, EPS), max(r_score, EPS)])
    f_maat = float(-np.sum(np.log(scores + EPS)))
    stability = float(min(r_score, (h_score * b_score * max(s_score, EPS) * max(v_score, EPS)) ** 0.25))

    return {
        "F_MAAT": f_maat,
        "Stability": stability,
        "H": h_score,
        "B": b_score,
        "S": s_score,
        "V": v_score,
        "R": r_score,
        "mean_w": float(np.mean(np.clip(w, -10.0, 10.0))),
        "inflation_fraction": inflation_fraction,
        "invalid_fraction": invalid_fraction,
        "efolds": efolds,
    }


def make_runs(config: CosmoConfig) -> list[dict]:
    return [
        evolve(0.2, 0.00, "slow_roll_coherent", config, potential="plateau"),
        evolve(0.2, 3.00, "kinetic_relaxing", config, potential="plateau"),
        evolve(8.0, 0.00, "static_dead", config, potential="plateau"),
        evolve(2.0, 0.50, "collapse_negative_potential", config, potential="negative"),
        evolve(0.2, 1.20, "ghost_like", config, potential="plateau", ghost=True),
    ]


def plot_histories(runs: list[dict], results: dict[str, dict], config: CosmoConfig, plot_dir: Path) -> list[str]:
    plot_dir.mkdir(parents=True, exist_ok=True)
    t = np.linspace(0.0, config.total_time, int(config.total_time / config.dt))

    plt.figure(figsize=(10, 6))
    for run in runs:
        plt.plot(t, np.clip(run["w"], -2.0, 2.0), label=run["label"])
    plt.axhline(-1.0 / 3.0, linestyle="--", color="black", alpha=0.5, label="acceleration threshold")
    plt.xlabel("t")
    plt.ylabel("w(t), clipped")
    plt.title("Cosmology v2: equation-of-state evolution")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.25)
    plt.tight_layout()
    w_plot = plot_dir / "fig_cosmo_v2_w.png"
    plt.savefig(w_plot, dpi=300)
    plt.close()

    plt.figure(figsize=(10, 6))
    for run in runs:
        plt.plot(t, np.log(run["a"] + EPS), label=run["label"])
    plt.xlabel("t")
    plt.ylabel("ln a(t)")
    plt.title("Cosmology v2: expansion history")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.25)
    plt.tight_layout()
    expansion_plot = plot_dir / "fig_cosmo_v2_expansion.png"
    plt.savefig(expansion_plot, dpi=300)
    plt.close()

    labels = list(results.keys())
    f_maat = [results[label]["F_MAAT"] for label in labels]
    plt.figure(figsize=(10, 6))
    plt.bar(labels, f_maat)
    plt.ylabel("F_MAAT")
    plt.title("Cosmology v2: MAAT ranking")
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    ranking_plot = plot_dir / "fig_cosmo_v2_maat_ranking.png"
    plt.savefig(ranking_plot, dpi=300)
    plt.close()

    return [w_plot.name, expansion_plot.name, ranking_plot.name]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("maat_cosmology_toy_v2_results.json"))
    parser.add_argument("--plot-dir", type=Path, default=Path("maat_cosmology_toy_v2_plots"))
    args = parser.parse_args()

    config = CosmoConfig()
    runs = make_runs(config)
    results = {run["label"]: maat_cosmo_score(run) for run in runs}
    ranking = [label for label, _ in sorted(results.items(), key=lambda item: item[1]["F_MAAT"])]
    plots = plot_histories(runs, results, config, args.plot_dir)

    payload = {
        "config": asdict(config),
        "ranking_by_F_MAAT": ranking,
        "results": results,
        "plots": plots,
        "status": "Operational FLRW scalar-field toy benchmark; not a complete cosmological model.",
    }
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("MAAT cosmology toy v2 complete.")
    print(f"Results: {args.output}")
    print(f"Plots:   {args.plot_dir}")
    print("Ranking:", " < ".join(ranking))


if __name__ == "__main__":
    main()

