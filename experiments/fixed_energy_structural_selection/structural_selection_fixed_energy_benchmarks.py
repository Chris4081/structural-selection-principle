#!/usr/bin/env python3
"""Fixed-energy and cross-model structural-selection benchmarks.

This script consolidates a small reproducibility bundle for testing whether a
structural MAAT/selection functional contains information beyond ordinary
energy minimisation.

The benchmarks are deliberately modest:

1. Static 1D phi^4 kink ranking.
2. Static 2D phi^4 domain-wall ranking.
3. Equal-energy 1D phi^4 perturbation test.
4. Static 1D Sine-Gordon soliton ranking.

All diagnostics are operational proxies. The script is not a proof of a
fundamental theory; it is a compact numerical stress test for the claim that
coherent structures can be distinguished from rough or unstable perturbations
even when total energy is approximately fixed.
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/codex-mpl-cache")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/codex-mpl-cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


EPS = 1.0e-8


@dataclass(frozen=True)
class BenchmarkConfig:
    seed: int = 42
    phi4_1d_n: int = 512
    phi4_2d_n: int = 160
    sine_gordon_n: int = 512
    length_phi4: float = 20.0
    length_sine_gordon: float = 30.0
    equal_energy_offset: float = 0.30


def open_gradient_1d(phi: np.ndarray, dx: float) -> np.ndarray:
    grad = np.zeros_like(phi)
    grad[1:-1] = (phi[2:] - phi[:-2]) / (2.0 * dx)
    grad[0] = (phi[1] - phi[0]) / dx
    grad[-1] = (phi[-1] - phi[-2]) / dx
    return grad


def open_laplacian_1d(phi: np.ndarray, dx: float) -> np.ndarray:
    lap = np.zeros_like(phi)
    lap[1:-1] = (phi[2:] - 2.0 * phi[1:-1] + phi[:-2]) / dx**2
    lap[0] = lap[1]
    lap[-1] = lap[-2]
    return lap


def open_gradient_2d(phi: np.ndarray, dx: float, dy: float) -> tuple[np.ndarray, np.ndarray]:
    gx = np.zeros_like(phi)
    gy = np.zeros_like(phi)
    gx[1:-1, :] = (phi[2:, :] - phi[:-2, :]) / (2.0 * dx)
    gx[0, :] = (phi[1, :] - phi[0, :]) / dx
    gx[-1, :] = (phi[-1, :] - phi[-2, :]) / dx
    gy[:, 1:-1] = (phi[:, 2:] - phi[:, :-2]) / (2.0 * dy)
    gy[:, 0] = (phi[:, 1] - phi[:, 0]) / dy
    gy[:, -1] = (phi[:, -1] - phi[:, -2]) / dy
    return gx, gy


def open_laplacian_2d(phi: np.ndarray, dx: float, dy: float) -> np.ndarray:
    lap = np.zeros_like(phi)
    lap[1:-1, 1:-1] = (
        (phi[2:, 1:-1] - 2.0 * phi[1:-1, 1:-1] + phi[:-2, 1:-1]) / dx**2
        + (phi[1:-1, 2:] - 2.0 * phi[1:-1, 1:-1] + phi[1:-1, :-2]) / dy**2
    )
    lap[0, :] = lap[1, :]
    lap[-1, :] = lap[-2, :]
    lap[:, 0] = lap[:, 1]
    lap[:, -1] = lap[:, -2]
    return lap


def nearest_neighbor_information_1d(phi: np.ndarray) -> float:
    a = phi[:-1]
    b = phi[1:]
    if np.std(a) < EPS or np.std(b) < EPS:
        return 0.0
    corr = float(np.corrcoef(a, b)[0, 1])
    if np.isnan(corr):
        return 0.0
    return max(0.0, corr)


def nearest_neighbor_information_2d(phi: np.ndarray) -> float:
    pairs = [
        (phi[:-1, :].ravel(), phi[1:, :].ravel()),
        (phi[:, :-1].ravel(), phi[:, 1:].ravel()),
    ]
    corrs: list[float] = []
    for a, b in pairs:
        if np.std(a) < EPS or np.std(b) < EPS:
            continue
        corr = float(np.corrcoef(a, b)[0, 1])
        if not np.isnan(corr):
            corrs.append(max(0.0, corr))
    if not corrs:
        return 0.0
    return float(np.mean(corrs))


def activity_support(phi: np.ndarray) -> float:
    var = float(np.var(phi))
    return var / (1.0 + var)


def phi4_energy_1d(phi: np.ndarray, dx: float) -> float:
    grad = open_gradient_1d(phi, dx)
    potential = 0.125 * (phi**2 - 1.0) ** 2
    return float(np.sum(0.5 * grad**2 + potential) * dx)


def phi4_residual_1d(phi: np.ndarray, dx: float) -> np.ndarray:
    return -open_laplacian_1d(phi, dx) + 0.5 * phi * (phi**2 - 1.0)


def phi4_energy_2d(phi: np.ndarray, dx: float, dy: float) -> float:
    gx, gy = open_gradient_2d(phi, dx, dy)
    potential = 0.125 * (phi**2 - 1.0) ** 2
    return float(np.sum(0.5 * (gx**2 + gy**2) + potential) * dx * dy)


def phi4_residual_2d(phi: np.ndarray, dx: float, dy: float) -> np.ndarray:
    return -open_laplacian_2d(phi, dx, dy) + 0.5 * phi * (phi**2 - 1.0)


def sine_gordon_energy(phi: np.ndarray, dx: float) -> float:
    grad = open_gradient_1d(phi, dx)
    potential = 1.0 - np.cos(phi)
    return float(np.sum(0.5 * grad**2 + potential) * dx)


def sine_gordon_residual(phi: np.ndarray, dx: float) -> np.ndarray:
    return -open_laplacian_1d(phi, dx) + np.sin(phi)


def structural_score_1d(
    phi: np.ndarray,
    dx: float,
    energy_fn: Callable[[np.ndarray, float], float],
    residual_fn: Callable[[np.ndarray, float], np.ndarray],
) -> dict:
    residual = residual_fn(phi, dx)
    lap = open_laplacian_1d(phi, dx)
    d_h = float(np.mean(residual**2))
    d_r = float(np.mean(lap**2))
    h = 1.0 / (1.0 + d_h)
    b = 1.0
    s = activity_support(phi)
    info = nearest_neighbor_information_1d(phi)
    v = info / (1.0 + info)
    r = 1.0 / (1.0 + d_r)
    supports = np.array([h, b, max(s, EPS), max(v, EPS), r], dtype=float)
    f_maat = float(-np.sum(np.log(supports + EPS)))
    stability = float(min(r, (h * b * max(s, EPS) * max(v, EPS)) ** 0.25))
    return {
        "Energy": energy_fn(phi, dx),
        "F_MAAT": f_maat,
        "Stability": stability,
        "H": h,
        "B": b,
        "S": s,
        "V": v,
        "R": r,
        "d_H": d_h,
        "d_R": d_r,
        "Info": info,
    }


def structural_score_2d(phi: np.ndarray, dx: float, dy: float) -> dict:
    residual = phi4_residual_2d(phi, dx, dy)
    lap = open_laplacian_2d(phi, dx, dy)
    d_h = float(np.mean(residual**2))
    d_r = float(np.mean(lap**2))
    h = 1.0 / (1.0 + d_h)
    b = 1.0
    s = activity_support(phi)
    info = nearest_neighbor_information_2d(phi)
    v = info / (1.0 + info)
    r = 1.0 / (1.0 + d_r)
    supports = np.array([h, b, max(s, EPS), max(v, EPS), r], dtype=float)
    f_maat = float(-np.sum(np.log(supports + EPS)))
    stability = float(min(r, (h * b * max(s, EPS) * max(v, EPS)) ** 0.25))
    return {
        "Energy": phi4_energy_2d(phi, dx, dy),
        "F_MAAT": f_maat,
        "Stability": stability,
        "H": h,
        "B": b,
        "S": s,
        "V": v,
        "R": r,
        "d_H": d_h,
        "d_R": d_r,
        "Info": info,
    }


def serialise_results(results: dict[str, dict]) -> dict[str, dict]:
    return {
        name: {
            key: float(value) if isinstance(value, (np.floating, np.integer)) else value
            for key, value in result.items()
            if key != "field"
        }
        for name, result in results.items()
    }


def ranking(results: dict[str, dict], key: str = "F_MAAT") -> list[str]:
    return [name for name, _ in sorted(results.items(), key=lambda item: item[1][key])]


def make_phi4_1d(config: BenchmarkConfig) -> tuple[np.ndarray, float, dict[str, dict]]:
    rng = np.random.default_rng(config.seed)
    x = np.linspace(-config.length_phi4 / 2.0, config.length_phi4 / 2.0, config.phi4_1d_n)
    dx = float(x[1] - x[0])
    kink = np.tanh(x / 2.0)
    states = {
        "vacuum_plus": np.ones_like(x),
        "vacuum_minus": -np.ones_like(x),
        "kink": kink,
        "distorted_kink": kink + 0.15 * np.sin(3.0 * x) * np.exp(-0.05 * x**2),
        "strong_distorted_kink": kink + 0.35 * np.sin(5.0 * x) * np.exp(-0.04 * x**2),
        "chaos": rng.uniform(-1.0, 1.0, size=config.phi4_1d_n),
    }
    results = {
        name: {"field": field, **structural_score_1d(field, dx, phi4_energy_1d, phi4_residual_1d)}
        for name, field in states.items()
    }
    return x, dx, results


def make_phi4_2d(config: BenchmarkConfig) -> tuple[np.ndarray, np.ndarray, float, float, dict[str, dict]]:
    rng = np.random.default_rng(config.seed)
    x = np.linspace(-config.length_phi4 / 2.0, config.length_phi4 / 2.0, config.phi4_2d_n)
    y = np.linspace(-config.length_phi4 / 2.0, config.length_phi4 / 2.0, config.phi4_2d_n)
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])
    x_grid, y_grid = np.meshgrid(x, y, indexing="ij")
    domain_wall = np.tanh(x_grid / 2.0)
    shift = 0.75 * np.sin(2.0 * np.pi * y_grid / config.length_phi4)
    strong_shift = 1.50 * np.sin(2.0 * np.pi * y_grid / config.length_phi4) + 0.70 * np.sin(
        6.0 * np.pi * y_grid / config.length_phi4
    )
    states = {
        "vacuum_plus": np.ones_like(x_grid),
        "vacuum_minus": -np.ones_like(x_grid),
        "domain_wall": domain_wall,
        "distorted_domain_wall": np.tanh((x_grid - shift) / 2.0),
        "strong_distorted_domain_wall": np.tanh((x_grid - strong_shift) / 2.0),
        "chaos": rng.uniform(-1.0, 1.0, size=(config.phi4_2d_n, config.phi4_2d_n)),
    }
    results = {
        name: {"field": field, **structural_score_2d(field, dx, dy)}
        for name, field in states.items()
    }
    return x, y, dx, dy, results


def find_amplitude(
    generator: Callable[[float], np.ndarray],
    target_energy: float,
    dx: float,
    lo: float = 0.0,
    hi: float = 2.0,
    tol: float = 1.0e-5,
) -> float:
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        energy_mid = phi4_energy_1d(generator(mid), dx)
        if energy_mid < target_energy:
            lo = mid
        else:
            hi = mid
        if abs(energy_mid - target_energy) < tol:
            break
    return 0.5 * (lo + hi)


def make_equal_energy(config: BenchmarkConfig) -> tuple[np.ndarray, float, float, dict[str, dict]]:
    x = np.linspace(-config.length_phi4 / 2.0, config.length_phi4 / 2.0, config.phi4_1d_n)
    dx = float(x[1] - x[0])
    kink = np.tanh(x / 2.0)
    base_energy = phi4_energy_1d(kink, dx)
    target_energy = base_energy + config.equal_energy_offset

    def smooth_wave(amplitude: float) -> np.ndarray:
        return kink + amplitude * np.sin(3.0 * x) * np.exp(-0.05 * x**2)

    def localized_bump(amplitude: float) -> np.ndarray:
        return kink + amplitude * np.exp(-(x**2) / 2.0)

    def rough_noise(amplitude: float) -> np.ndarray:
        rng = np.random.default_rng(7)
        return kink + amplitude * rng.normal(0.0, 1.0, size=config.phi4_1d_n) * np.exp(-0.04 * x**2)

    def two_bump(amplitude: float) -> np.ndarray:
        return kink + amplitude * np.exp(-((x - 3.0) ** 2) / (2.0 * 0.8**2)) - amplitude * np.exp(
            -((x + 3.0) ** 2) / (2.0 * 0.8**2)
        )

    generators = {
        "smooth_wave": smooth_wave,
        "localized_bump": localized_bump,
        "rough_noise": rough_noise,
        "two_bump": two_bump,
    }
    results: dict[str, dict] = {}
    for name, generator in generators.items():
        amplitude = find_amplitude(generator, target_energy, dx)
        field = generator(amplitude)
        results[name] = {
            "field": field,
            "amplitude": float(amplitude),
            **structural_score_1d(field, dx, phi4_energy_1d, phi4_residual_1d),
        }
    results["kink_reference"] = {
        "field": kink,
        "amplitude": 0.0,
        **structural_score_1d(kink, dx, phi4_energy_1d, phi4_residual_1d),
    }
    return x, base_energy, target_energy, results


def make_sine_gordon(config: BenchmarkConfig) -> tuple[np.ndarray, float, dict[str, dict]]:
    rng = np.random.default_rng(config.seed)
    x = np.linspace(
        -config.length_sine_gordon / 2.0,
        config.length_sine_gordon / 2.0,
        config.sine_gordon_n,
    )
    dx = float(x[1] - x[0])
    soliton = 4.0 * np.arctan(np.exp(x))
    states = {
        "vacuum_0": np.zeros_like(x),
        "vacuum_2pi": 2.0 * np.pi * np.ones_like(x),
        "soliton": soliton,
        "distorted_soliton": soliton + 0.20 * np.sin(3.0 * x) * np.exp(-0.04 * x**2),
        "strong_distorted_soliton": soliton + 0.50 * np.sin(5.0 * x) * np.exp(-0.035 * x**2),
        "chaos": rng.uniform(0.0, 2.0 * np.pi, size=config.sine_gordon_n),
    }
    results = {
        name: {"field": field, **structural_score_1d(field, dx, sine_gordon_energy, sine_gordon_residual)}
        for name, field in states.items()
    }
    return x, dx, results


def robustness_phi4_1d(config: BenchmarkConfig) -> dict:
    def run_once(n: int, seed: int, distortion: float, strong_distortion: float) -> tuple[bool, list[str]]:
        rng = np.random.default_rng(seed)
        x = np.linspace(-config.length_phi4 / 2.0, config.length_phi4 / 2.0, n)
        dx = float(x[1] - x[0])
        kink = np.tanh(x / 2.0)
        states = {
            "kink": kink,
            "distorted": kink + distortion * np.sin(3.0 * x) * np.exp(-0.05 * x**2),
            "strong_distorted": kink + strong_distortion * np.sin(5.0 * x) * np.exp(-0.04 * x**2),
            "chaos": rng.uniform(-1.0, 1.0, size=n),
        }
        scores = {
            name: structural_score_1d(field, dx, phi4_energy_1d, phi4_residual_1d)["F_MAAT"]
            for name, field in states.items()
        }
        order = [name for name, _ in sorted(scores.items(), key=lambda item: item[1])]
        return order == ["kink", "distorted", "strong_distorted", "chaos"], order

    total = 0
    passed = 0
    failures = []
    for n in [256, 512, 1024]:
        for seed in range(20):
            for distortion in [0.05, 0.10, 0.15, 0.20]:
                for strong_distortion in [0.25, 0.35, 0.50]:
                    total += 1
                    ok, order = run_once(n, seed, distortion, strong_distortion)
                    passed += int(ok)
                    if not ok and len(failures) < 10:
                        failures.append(
                            {
                                "N": n,
                                "seed": seed,
                                "distortion": distortion,
                                "strong_distortion": strong_distortion,
                                "observed_order": order,
                            }
                        )
    return {"total": total, "passed": passed, "pass_rate": passed / total, "first_failures": failures}


def robustness_phi4_2d(config: BenchmarkConfig) -> dict:
    def run_once(n: int, seed: int, distortion: float, strong_distortion: float) -> tuple[bool, list[str]]:
        rng = np.random.default_rng(seed)
        x = np.linspace(-config.length_phi4 / 2.0, config.length_phi4 / 2.0, n)
        y = np.linspace(-config.length_phi4 / 2.0, config.length_phi4 / 2.0, n)
        dx = float(x[1] - x[0])
        dy = float(y[1] - y[0])
        x_grid, y_grid = np.meshgrid(x, y, indexing="ij")
        domain_wall = np.tanh(x_grid / 2.0)
        shift = distortion * np.sin(2.0 * np.pi * y_grid / config.length_phi4)
        strong_shift = strong_distortion * np.sin(2.0 * np.pi * y_grid / config.length_phi4) + 0.70 * np.sin(
            6.0 * np.pi * y_grid / config.length_phi4
        )
        states = {
            "domain_wall": domain_wall,
            "distorted": np.tanh((x_grid - shift) / 2.0),
            "strong_distorted": np.tanh((x_grid - strong_shift) / 2.0),
            "chaos": rng.uniform(-1.0, 1.0, size=(n, n)),
        }
        scores = {name: structural_score_2d(field, dx, dy)["F_MAAT"] for name, field in states.items()}
        order = [name for name, _ in sorted(scores.items(), key=lambda item: item[1])]
        return order == ["domain_wall", "distorted", "strong_distorted", "chaos"], order

    total = 0
    passed = 0
    failures = []
    for n in [64, 96, 128]:
        for seed in range(10):
            for distortion in [0.25, 0.50, 0.75, 1.00]:
                for strong_distortion in [1.25, 1.50, 2.00]:
                    total += 1
                    ok, order = run_once(n, seed, distortion, strong_distortion)
                    passed += int(ok)
                    if not ok and len(failures) < 10:
                        failures.append(
                            {
                                "N": n,
                                "seed": seed,
                                "distortion": distortion,
                                "strong_distortion": strong_distortion,
                                "observed_order": order,
                            }
                        )
    return {"total": total, "passed": passed, "pass_rate": passed / total, "first_failures": failures}


def robustness_sine_gordon(config: BenchmarkConfig) -> dict:
    def run_once(n: int, seed: int, distortion: float, strong_distortion: float) -> tuple[bool, list[str]]:
        rng = np.random.default_rng(seed)
        x = np.linspace(-config.length_sine_gordon / 2.0, config.length_sine_gordon / 2.0, n)
        dx = float(x[1] - x[0])
        soliton = 4.0 * np.arctan(np.exp(x))
        states = {
            "soliton": soliton,
            "distorted": soliton + distortion * np.sin(3.0 * x) * np.exp(-0.04 * x**2),
            "strong_distorted": soliton + strong_distortion * np.sin(5.0 * x) * np.exp(-0.035 * x**2),
            "chaos": rng.uniform(0.0, 2.0 * np.pi, size=n),
        }
        scores = {
            name: structural_score_1d(field, dx, sine_gordon_energy, sine_gordon_residual)["F_MAAT"]
            for name, field in states.items()
        }
        order = [name for name, _ in sorted(scores.items(), key=lambda item: item[1])]
        return order == ["soliton", "distorted", "strong_distorted", "chaos"], order

    total = 0
    passed = 0
    failures = []
    for n in [256, 512, 1024]:
        for seed in range(20):
            for distortion in [0.05, 0.10, 0.20, 0.30]:
                for strong_distortion in [0.40, 0.50, 0.70]:
                    total += 1
                    ok, order = run_once(n, seed, distortion, strong_distortion)
                    passed += int(ok)
                    if not ok and len(failures) < 10:
                        failures.append(
                            {
                                "N": n,
                                "seed": seed,
                                "distortion": distortion,
                                "strong_distortion": strong_distortion,
                                "observed_order": order,
                            }
                        )
    return {"total": total, "passed": passed, "pass_rate": passed / total, "first_failures": failures}


def annotate_points(ax: plt.Axes, labels: list[str], x_values: list[float], y_values: list[float], offsets: dict[str, tuple[int, int]]) -> None:
    for label, x_val, y_val in zip(labels, x_values, y_values):
        dx, dy = offsets.get(label, (6, 6))
        ax.annotate(
            label.replace("_", " "),
            (x_val, y_val),
            textcoords="offset points",
            xytext=(dx, dy),
            fontsize=8,
            bbox={"boxstyle": "round,pad=0.2", "fc": "white", "ec": "0.75", "alpha": 0.9},
            arrowprops={"arrowstyle": "-", "lw": 0.5, "alpha": 0.45},
        )


def annotate_point_data(
    ax: plt.Axes,
    label: str,
    xy: tuple[float, float],
    xytext: tuple[float, float],
    fontsize: int = 8,
) -> None:
    ax.annotate(
        label.replace("_", " "),
        xy,
        xycoords="data",
        xytext=xytext,
        textcoords="data",
        fontsize=fontsize,
        bbox={"boxstyle": "round,pad=0.2", "fc": "white", "ec": "0.75", "alpha": 0.92},
        arrowprops={"arrowstyle": "-", "lw": 0.55, "alpha": 0.5},
    )


def plot_phi4_1d_fields(x: np.ndarray, results: dict[str, dict], plot_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 6.0), constrained_layout=True)
    for name in ["kink", "distorted_kink", "strong_distorted_kink", "chaos"]:
        linewidth = 2.8 if name == "kink" else 1.8
        alpha = 0.45 if name == "chaos" else 1.0
        ax.plot(x, results[name]["field"], linewidth=linewidth, alpha=alpha, label=name.replace("_", " "))
    ax.set_xlabel("x")
    ax.set_ylabel("phi(x)")
    ax.set_title("1D phi^4 structural candidates")
    ax.grid(alpha=0.25)
    ax.legend(loc="upper left", frameon=True, fontsize=9)
    fig.savefig(plot_dir / "phi4_1d_fields.png", dpi=260)
    plt.close(fig)


def plot_phi4_2d_fields(x: np.ndarray, y: np.ndarray, results: dict[str, dict], plot_dir: Path) -> None:
    fields = [
        ("domain wall", results["domain_wall"]["field"]),
        ("distorted wall", results["distorted_domain_wall"]["field"]),
        ("strong distorted wall", results["strong_distorted_domain_wall"]["field"]),
        ("chaos", results["chaos"]["field"]),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.6), constrained_layout=True, sharex=True, sharey=True)
    axes_flat = axes.ravel()
    last_image = None
    for ax, (title, field) in zip(axes_flat, fields):
        last_image = ax.imshow(
            field.T,
            origin="lower",
            extent=[float(x[0]), float(x[-1]), float(y[0]), float(y[-1])],
            cmap="coolwarm",
            vmin=-1.0,
            vmax=1.0,
            interpolation="nearest",
        )
        ax.set_title(title, fontsize=10, pad=6)
        ax.set_aspect("equal")
        ax.grid(False)
    for ax in axes[:, 0]:
        ax.set_ylabel("y")
    for ax in axes[-1, :]:
        ax.set_xlabel("x")
    fig.suptitle("2D phi^4 domain-wall configurations", fontsize=13)
    fig.colorbar(last_image, ax=axes_flat.tolist(), shrink=0.82, pad=0.02, label="phi")
    fig.savefig(plot_dir / "phi4_2d_fields.png", dpi=260)
    plt.close(fig)


def plot_energy_vs_structural(phi4_1d: dict[str, dict], phi4_2d: dict[str, dict], plot_dir: Path) -> None:
    labels_1d = ["vacuum_plus", "kink", "distorted_kink", "strong_distorted_kink", "chaos"]
    labels_2d = ["vacuum_plus", "domain_wall", "distorted_domain_wall", "strong_distorted_domain_wall", "chaos"]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.5), constrained_layout=True)

    x1 = [phi4_1d[label]["Energy"] for label in labels_1d]
    y1 = [phi4_1d[label]["F_MAAT"] for label in labels_1d]
    axes[0].scatter(x1, y1, s=80, color="#1f77b4")
    annotate_points(
        axes[0],
        labels_1d,
        x1,
        y1,
        {
            "vacuum_plus": (8, -20),
            "kink": (8, 10),
            "distorted_kink": (8, 16),
            "strong_distorted_kink": (-88, 12),
            "chaos": (-60, -22),
        },
    )
    axes[0].set_xscale("symlog", linthresh=1.0)
    axes[0].set_xlabel("Energy")
    axes[0].set_ylabel("F_MAAT")
    axes[0].set_title("1D phi^4")
    axes[0].grid(alpha=0.25)

    x2 = [phi4_2d[label]["Energy"] for label in labels_2d]
    y2 = [phi4_2d[label]["F_MAAT"] for label in labels_2d]
    axes[1].scatter(x2, y2, s=80, color="#d62728")
    # The three domain-wall points are intentionally close in both energy and
    # F_MAAT. Use fixed data-coordinate labels to keep the plot readable.
    point_map_2d = {label: (phi4_2d[label]["Energy"], phi4_2d[label]["F_MAAT"]) for label in labels_2d}
    annotate_point_data(axes[1], "vacuum_plus", point_map_2d["vacuum_plus"], (0.25, 33.2))
    annotate_point_data(axes[1], "domain_wall", point_map_2d["domain_wall"], (28.0, 3.35))
    annotate_point_data(axes[1], "distorted_domain_wall", point_map_2d["distorted_domain_wall"], (28.0, 0.55))
    annotate_point_data(axes[1], "strong_distorted_domain_wall", point_map_2d["strong_distorted_domain_wall"], (1.0, 3.95))
    annotate_point_data(axes[1], "chaos", point_map_2d["chaos"], (80.0, 26.5))
    axes[1].set_xscale("symlog", linthresh=1.0)
    axes[1].set_ylim(-1.0, max(y2) + 2.4)
    axes[1].set_xlabel("Energy")
    axes[1].set_ylabel("F_MAAT")
    axes[1].set_title("2D phi^4")
    axes[1].grid(alpha=0.25)

    fig.suptitle("Energy versus structural ranking")
    fig.savefig(plot_dir / "energy_vs_structural.png", dpi=260)
    plt.close(fig)


def plot_equal_energy(x: np.ndarray, target_energy: float, results: dict[str, dict], plot_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 6.0), constrained_layout=True)
    for name, result in results.items():
        linewidth = 2.8 if name == "kink_reference" else 1.8
        alpha = 0.70 if name == "rough_noise" else 1.0
        ax.plot(x, result["field"], linewidth=linewidth, alpha=alpha, label=name.replace("_", " "))
    ax.set_xlabel("x")
    ax.set_ylabel("phi(x)")
    ax.set_title("Equal-energy 1D phi^4 candidates")
    ax.legend(loc="upper left", fontsize=8, ncol=2, frameon=True)
    ax.grid(alpha=0.25)
    fig.savefig(plot_dir / "equal_energy_fields.png", dpi=260)
    plt.close(fig)

    labels = ["two_bump", "localized_bump", "smooth_wave", "rough_noise", "kink_reference"]
    energies = [results[label]["Energy"] for label in labels]
    f_values = [results[label]["F_MAAT"] for label in labels]
    fig, ax = plt.subplots(figsize=(8.5, 6.2), constrained_layout=True)
    ax.scatter(energies, f_values, s=90, color="#2ca02c")
    ax.axvline(target_energy, linestyle="--", alpha=0.45, color="black", label="target energy")
    annotate_points(
        ax,
        labels,
        energies,
        f_values,
        {
            "two_bump": (10, 14),
            "localized_bump": (10, -20),
            "smooth_wave": (-88, 12),
            "rough_noise": (-82, -18),
            "kink_reference": (10, 10),
        },
    )
    ax.set_xlabel("Energy")
    ax.set_ylabel("F_MAAT")
    ax.set_title("Fixed energy, different structural quality")
    ax.grid(alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    fig.savefig(plot_dir / "equal_energy_structural_score.png", dpi=260)
    plt.close(fig)


def plot_sine_gordon(x: np.ndarray, results: dict[str, dict], plot_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 6.0), constrained_layout=True)
    for name in ["soliton", "distorted_soliton", "strong_distorted_soliton", "chaos"]:
        linewidth = 2.8 if name == "soliton" else 1.8
        alpha = 0.45 if name == "chaos" else 1.0
        ax.plot(x, results[name]["field"], linewidth=linewidth, alpha=alpha, label=name.replace("_", " "))
    ax.set_xlabel("x")
    ax.set_ylabel("phi(x)")
    ax.set_title("Sine-Gordon structural candidates")
    ax.legend(loc="upper left", fontsize=9, frameon=True)
    ax.grid(alpha=0.25)
    fig.savefig(plot_dir / "sine_gordon_fields.png", dpi=260)
    plt.close(fig)

    labels = ["vacuum_0", "soliton", "distorted_soliton", "strong_distorted_soliton", "chaos"]
    energies = [results[label]["Energy"] for label in labels]
    f_values = [results[label]["F_MAAT"] for label in labels]
    fig, ax = plt.subplots(figsize=(8.5, 6.0), constrained_layout=True)
    ax.scatter(energies, f_values, s=90, color="#9467bd")
    annotate_points(
        ax,
        labels,
        energies,
        f_values,
        {
            "vacuum_0": (8, -20),
            "soliton": (10, 12),
            "distorted_soliton": (10, -22),
            "strong_distorted_soliton": (-102, 14),
            "chaos": (-58, -22),
        },
    )
    ax.set_xscale("symlog", linthresh=1.0)
    ax.set_xlabel("Energy")
    ax.set_ylabel("F_MAAT")
    ax.set_title("Sine-Gordon: energy versus structural ranking")
    ax.grid(alpha=0.25)
    fig.savefig(plot_dir / "sine_gordon_energy_vs_structural.png", dpi=260)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("fixed_energy_structural_selection_results.json"))
    parser.add_argument("--plot-dir", type=Path, default=Path("fixed_energy_structural_selection_plots"))
    args = parser.parse_args()

    config = BenchmarkConfig()
    args.plot_dir.mkdir(parents=True, exist_ok=True)

    phi4_1d_x, _, phi4_1d_results = make_phi4_1d(config)
    phi4_2d_x, phi4_2d_y, _, _, phi4_2d_results = make_phi4_2d(config)
    equal_x, base_energy, target_energy, equal_results = make_equal_energy(config)
    sine_x, _, sine_results = make_sine_gordon(config)

    robustness = {
        "phi4_1d": robustness_phi4_1d(config),
        "phi4_2d": robustness_phi4_2d(config),
        "sine_gordon": robustness_sine_gordon(config),
    }
    robustness["combined"] = {
        "total": sum(item["total"] for item in robustness.values()),
        "passed": sum(item["passed"] for item in robustness.values()),
    }
    robustness["combined"]["pass_rate"] = robustness["combined"]["passed"] / robustness["combined"]["total"]

    plot_phi4_1d_fields(phi4_1d_x, phi4_1d_results, args.plot_dir)
    plot_phi4_2d_fields(phi4_2d_x, phi4_2d_y, phi4_2d_results, args.plot_dir)
    plot_energy_vs_structural(phi4_1d_results, phi4_2d_results, args.plot_dir)
    plot_equal_energy(equal_x, target_energy, equal_results, args.plot_dir)
    plot_sine_gordon(sine_x, sine_results, args.plot_dir)

    payload = {
        "config": asdict(config),
        "phi4_1d": {
            "ranking_by_F_MAAT": ranking(phi4_1d_results),
            "results": serialise_results(phi4_1d_results),
        },
        "phi4_2d": {
            "ranking_by_F_MAAT": ranking(phi4_2d_results),
            "results": serialise_results(phi4_2d_results),
        },
        "equal_energy_phi4_1d": {
            "base_kink_energy": base_energy,
            "target_energy": target_energy,
            "ranking_by_F_MAAT": ranking(equal_results),
            "ranking_by_energy": ranking(equal_results, key="Energy"),
            "results": serialise_results(equal_results),
        },
        "sine_gordon": {
            "ranking_by_F_MAAT": ranking(sine_results),
            "results": serialise_results(sine_results),
        },
        "robustness": robustness,
        "plots": sorted(path.name for path in args.plot_dir.glob("*.png")),
    }
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("Fixed-energy structural-selection benchmark complete.")
    print(f"Results: {args.output}")
    print(f"Plots:   {args.plot_dir}")
    print(
        "Robustness:",
        f"{robustness['combined']['passed']}/{robustness['combined']['total']}",
        f"({100.0 * robustness['combined']['pass_rate']:.2f}%)",
    )
    print("2D plot layout uses a 2x2 constrained grid to avoid title/label/colorbar overlap.")


if __name__ == "__main__":
    main()
