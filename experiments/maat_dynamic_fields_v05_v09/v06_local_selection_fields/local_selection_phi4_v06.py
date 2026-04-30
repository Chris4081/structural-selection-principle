#!/usr/bin/env python3
"""MAAT v0.6 local selection-pressure fields in a 1D phi^4 kink.

This experiment promotes global response weights lambda_a to local fields
lambda_a(x,t).  It is an effective benchmark, not a microscopic derivation.

The field evolves by an overdamped 1D phi^4 equation with local MAAT feedback:

    phi_t = M(x,t) * (phi_xx + phi - phi^3) - g_S lambda_S highpass(phi)

where M(x,t) = 1 + g_H lambda_H + g_R lambda_R.

Local response fields obey:

    tau lambda_t = D_lambda lambda_xx - lambda + lambda_star(x,t)

with lambda_star computed from local coarse-grained defect covariance.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "outputs"
FIELDS = ["H", "S", "R"]
A = len(FIELDS)


def gaussian_kernel(sigma_pts: float, radius: int | None = None) -> np.ndarray:
    if radius is None:
        radius = max(3, int(4 * sigma_pts))
    grid = np.arange(-radius, radius + 1)
    kernel = np.exp(-0.5 * (grid / sigma_pts) ** 2)
    return kernel / kernel.sum()


def smooth_edge(values: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    radius = len(kernel) // 2
    padded = np.pad(values, radius, mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def dx(values: np.ndarray, h: float) -> np.ndarray:
    out = np.zeros_like(values)
    out[1:-1] = (values[2:] - values[:-2]) / (2 * h)
    out[0] = out[1]
    out[-1] = out[-2]
    return out


def dxx(values: np.ndarray, h: float) -> np.ndarray:
    out = np.zeros_like(values)
    out[1:-1] = (values[2:] - 2 * values[1:-1] + values[:-2]) / (h * h)
    out[0] = out[1]
    out[-1] = out[-2]
    return out


def kink_profile(x: np.ndarray) -> np.ndarray:
    return np.tanh(x / np.sqrt(2.0))


def phi4_residual(phi: np.ndarray, h: float) -> np.ndarray:
    return dxx(phi, h) + phi - phi**3


def static_energy(phi: np.ndarray, h: float) -> float:
    grad = dx(phi, h)
    density = 0.5 * grad**2 + 0.25 * (phi**2 - 1.0) ** 2
    return float(np.trapezoid(density, dx=h))


def robust_scale(values: np.ndarray) -> float:
    return float(np.percentile(values, 95) + 1e-10)


def compute_raw_defects(phi: np.ndarray, h: float) -> dict[str, np.ndarray]:
    residual = phi4_residual(phi, h)
    curvature = dxx(phi, h)
    potential = 0.25 * (phi**2 - 1.0) ** 2
    return {
        "H": residual**2,
        "S": curvature**2,
        "R": potential,
    }


def normalise_defects(raw: dict[str, np.ndarray], scales: dict[str, float]) -> np.ndarray:
    defects = []
    for name in FIELDS:
        q = raw[name]
        defects.append(q / (scales[name] + q + 1e-12))
    return np.stack(defects, axis=0)


def local_mean_cov(defects: np.ndarray, kernel: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return local mean[A,N] and covariance[A,A,N]."""
    mean = np.stack([smooth_edge(defects[i], kernel) for i in range(A)], axis=0)
    cov = np.zeros((A, A, defects.shape[1]))
    for i in range(A):
        for j in range(A):
            prod = smooth_edge(defects[i] * defects[j], kernel)
            cov[i, j] = prod - mean[i] * mean[j]
    return mean, cov


def lambda_star_field(
    local_mean: np.ndarray,
    local_cov: np.ndarray,
    target: np.ndarray,
    eta: float,
    max_lambda: float,
) -> np.ndarray:
    n = local_mean.shape[1]
    star = np.zeros_like(local_mean)
    for idx in range(n):
        mat = local_cov[:, :, idx]
        reg = eta * np.trace(mat) / A + 1e-8
        rhs = local_mean[:, idx] - target[:, idx]
        try:
            sol = np.linalg.solve(mat + reg * np.eye(A), rhs)
        except np.linalg.LinAlgError:
            sol = np.linalg.lstsq(mat + reg * np.eye(A), rhs, rcond=None)[0]
        star[:, idx] = np.clip(sol, 0.0, max_lambda)
    return star


def apply_boundary(phi: np.ndarray) -> None:
    phi[0] = -1.0
    phi[-1] = 1.0


def highpass(phi: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    return phi - smooth_edge(phi, kernel)


def simulate() -> dict:
    rng = np.random.default_rng(42)

    # Spatial grid and integration parameters.
    n = 512
    length = 60.0
    x = np.linspace(-length / 2, length / 2, n)
    h = x[1] - x[0]
    dt = 0.002
    steps = 6000
    t_final = dt * (steps - 1)

    # Local response parameters.
    tau_lambda = 0.35
    d_lambda = 0.08
    eta = 2e-2
    max_lambda = 5.0
    g_h = 0.55
    g_s = 0.28
    g_r = 0.22
    max_mobility = 3.0

    coarse_kernel = gaussian_kernel(sigma_pts=5.0)
    smooth_kernel = gaussian_kernel(sigma_pts=2.0)

    phi_ref = kink_profile(x)
    raw_ref = compute_raw_defects(phi_ref, h)

    noise = 0.16 * rng.normal(size=n)
    noise += 0.06 * np.sin(7.0 * np.pi * (x - x.min()) / length)
    envelope = np.exp(-(x / 18.0) ** 8)
    phi0 = phi_ref + envelope * noise
    apply_boundary(phi0)

    raw0 = compute_raw_defects(phi0, h)
    scales = {name: max(robust_scale(raw0[name]), robust_scale(raw_ref[name]), 1e-8) for name in FIELDS}

    target = normalise_defects(raw_ref, scales)
    # Give the target a small tolerance floor around the clean kink.
    target = np.minimum(0.95, target + 0.015)

    phi_base = phi0.copy()
    phi_maat = phi0.copy()
    lambdas = np.zeros((A, n))

    perturb_step = int(3.0 / dt)
    perturb_profile = 0.48 * np.exp(-0.5 * ((x - 8.0) / 1.1) ** 2)

    rows = []
    lambda_snapshots = []
    snapshot_times = [0.0, 1.5, 3.2, 4.5, 7.0, t_final]
    snapshot_steps = {min(steps - 1, int(t / dt)) for t in snapshot_times}
    field_snapshots = {}

    for step in range(steps):
        t = step * dt

        if step == perturb_step:
            phi_base += perturb_profile
            phi_maat += perturb_profile
            apply_boundary(phi_base)
            apply_boundary(phi_maat)

        # Baseline Allen-Cahn relaxation.
        r_base = phi4_residual(phi_base, h)
        phi_base += dt * r_base
        apply_boundary(phi_base)

        # MAAT local-response fields.
        raw = compute_raw_defects(phi_maat, h)
        defects = normalise_defects(raw, scales)
        local_mean, local_cov = local_mean_cov(defects, coarse_kernel)
        star = lambda_star_field(local_mean, local_cov, target, eta=eta, max_lambda=max_lambda)

        for i in range(A):
            lambdas[i] += dt * (
                d_lambda * dxx(lambdas[i], h)
                + (-lambdas[i] + star[i]) / tau_lambda
            )
        lambdas = np.clip(lambdas, 0.0, max_lambda)

        r_maat = phi4_residual(phi_maat, h)
        mobility = 1.0 + g_h * lambdas[0] + g_r * lambdas[2]
        mobility = np.clip(mobility, 1.0, max_mobility)
        phi_maat += dt * (
            mobility * r_maat
            - g_s * lambdas[1] * highpass(phi_maat, smooth_kernel)
        )
        apply_boundary(phi_maat)

        if step % 20 == 0 or step == steps - 1:
            base_residual = np.sqrt(np.mean(phi4_residual(phi_base, h) ** 2))
            maat_residual = np.sqrt(np.mean(phi4_residual(phi_maat, h) ** 2))
            base_dist = np.sqrt(np.mean((phi_base - phi_ref) ** 2))
            maat_dist = np.sqrt(np.mean((phi_maat - phi_ref) ** 2))
            base_rough = np.sqrt(np.mean(dxx(phi_base, h) ** 2))
            maat_rough = np.sqrt(np.mean(dxx(phi_maat, h) ** 2))
            support = 1.0 / (1.0 + np.mean(defects, axis=1))
            stability = min(float(support[2]), float((support[0] * support[1] * support[2]) ** (1 / 3)))
            rows.append(
                {
                    "step": step,
                    "t": t,
                    "baseline_energy": static_energy(phi_base, h),
                    "maat_energy": static_energy(phi_maat, h),
                    "baseline_residual_rms": float(base_residual),
                    "maat_residual_rms": float(maat_residual),
                    "baseline_distance_to_kink": float(base_dist),
                    "maat_distance_to_kink": float(maat_dist),
                    "baseline_roughness": float(base_rough),
                    "maat_roughness": float(maat_rough),
                    "mean_lambda_H": float(np.mean(lambdas[0])),
                    "mean_lambda_S": float(np.mean(lambdas[1])),
                    "mean_lambda_R": float(np.mean(lambdas[2])),
                    "max_lambda_H": float(np.max(lambdas[0])),
                    "max_lambda_S": float(np.max(lambdas[1])),
                    "max_lambda_R": float(np.max(lambdas[2])),
                    "maat_stability": stability,
                    "mean_defect_H": float(np.mean(defects[0])),
                    "mean_defect_S": float(np.mean(defects[1])),
                    "mean_defect_R": float(np.mean(defects[2])),
                }
            )

        if step in snapshot_steps:
            lambda_snapshots.append((t, lambdas.copy()))
            field_snapshots[f"{t:.2f}"] = {
                "baseline": phi_base.copy(),
                "maat": phi_maat.copy(),
            }

    return {
        "params": {
            "n": n,
            "length": length,
            "dx": h,
            "dt": dt,
            "steps": steps,
            "t_final": t_final,
            "tau_lambda": tau_lambda,
            "D_lambda": d_lambda,
            "eta": eta,
            "max_lambda": max_lambda,
            "g_H": g_h,
            "g_S": g_s,
            "g_R": g_r,
            "max_mobility": max_mobility,
            "perturb_time": perturb_step * dt,
            "seed": 42,
            "fields": FIELDS,
        },
        "x": x,
        "phi_ref": phi_ref,
        "phi0": phi0,
        "phi_base": phi_base,
        "phi_maat": phi_maat,
        "lambdas": lambdas,
        "lambda_snapshots": lambda_snapshots,
        "field_snapshots": field_snapshots,
        "target": target,
        "rows": rows,
    }


def plot_fields(result: dict) -> None:
    x = result["x"]
    plt.figure(figsize=(9.5, 5.7))
    plt.plot(x, result["phi_ref"], color="black", lw=2, label="clean kink")
    plt.plot(x, result["phi0"], color="gray", alpha=0.45, label="initial noisy kink")
    plt.plot(x, result["phi_base"], color="tab:orange", lw=1.8, label="baseline final")
    plt.plot(x, result["phi_maat"], color="tab:blue", lw=1.8, label="MAAT-local final")
    plt.xlabel("x")
    plt.ylabel(r"$\phi(x)$")
    plt.title(r"1D $\phi^4$ kink: baseline vs local MAAT selection fields")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUTDIR / "field_profiles.png", dpi=240)
    plt.close()


def plot_metrics(result: dict) -> None:
    rows = result["rows"]
    t = np.array([row["t"] for row in rows])
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.6), sharex=True)
    pairs = [
        ("distance_to_kink", "Distance to clean kink"),
        ("residual_rms", "EOM residual RMS"),
        ("roughness", "Curvature roughness"),
        ("energy", "Static energy"),
    ]
    for ax, (key, title) in zip(axes.ravel(), pairs):
        ax.plot(t, [row[f"baseline_{key}"] for row in rows], label="baseline", color="tab:orange")
        ax.plot(t, [row[f"maat_{key}"] for row in rows], label="MAAT local", color="tab:blue")
        ax.axvline(result["params"]["perturb_time"], color="tab:red", alpha=0.35, ls="--")
        ax.set_title(title)
        ax.grid(alpha=0.3)
    axes[1, 0].set_xlabel("time")
    axes[1, 1].set_xlabel("time")
    axes[0, 0].legend()
    fig.tight_layout()
    plt.savefig(OUTDIR / "metric_comparison.png", dpi=240)
    plt.close()


def plot_lambda_snapshots(result: dict) -> None:
    x = result["x"]
    snapshots = result["lambda_snapshots"]
    fig, axes = plt.subplots(A, 1, figsize=(9.5, 8.0), sharex=True)
    for i, name in enumerate(FIELDS):
        for t, lam in snapshots:
            axes[i].plot(x, lam[i], label=f"t={t:.1f}")
        axes[i].set_ylabel(rf"$\lambda_{name}(x)$")
        axes[i].grid(alpha=0.3)
    axes[-1].set_xlabel("x")
    axes[0].set_title("Local selection-pressure fields")
    axes[0].legend(ncol=3, fontsize=8)
    fig.tight_layout()
    plt.savefig(OUTDIR / "lambda_field_snapshots.png", dpi=240)
    plt.close()


def plot_lambda_heatmap(result: dict) -> None:
    # Reconstruct sparse heatmap from snapshots for a compact figure.
    x = result["x"]
    snapshots = result["lambda_snapshots"]
    fig, axes = plt.subplots(1, A, figsize=(13, 4.4), sharey=True)
    for i, name in enumerate(FIELDS):
        mat = np.stack([lam[i] for _, lam in snapshots], axis=0)
        times = [t for t, _ in snapshots]
        im = axes[i].imshow(
            mat,
            extent=[x.min(), x.max(), min(times), max(times)],
            aspect="auto",
            origin="lower",
            cmap="magma",
        )
        axes[i].set_title(rf"$\lambda_{name}(x,t)$")
        axes[i].set_xlabel("x")
        fig.colorbar(im, ax=axes[i], fraction=0.046, pad=0.04)
    axes[0].set_ylabel("time")
    fig.tight_layout()
    plt.savefig(OUTDIR / "lambda_field_heatmap.png", dpi=240)
    plt.close()


def plot_final_defects(result: dict) -> None:
    x = result["x"]
    raw_final = compute_raw_defects(result["phi_maat"], result["params"]["dx"])
    raw_ref = compute_raw_defects(result["phi_ref"], result["params"]["dx"])
    scales = {name: max(robust_scale(raw_final[name]), robust_scale(raw_ref[name]), 1e-8) for name in FIELDS}
    defects_final = normalise_defects(raw_final, scales)
    defects_ref = normalise_defects(raw_ref, scales)

    fig, axes = plt.subplots(A, 1, figsize=(9.5, 7.5), sharex=True)
    for i, name in enumerate(FIELDS):
        axes[i].plot(x, defects_ref[i], color="black", label="clean kink target")
        axes[i].plot(x, defects_final[i], color="tab:blue", alpha=0.85, label="MAAT final")
        axes[i].set_ylabel(rf"$d_{name}(x)$")
        axes[i].grid(alpha=0.3)
    axes[-1].set_xlabel("x")
    axes[0].set_title("Final local defects compared with clean-kink target")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    plt.savefig(OUTDIR / "final_defect_profiles.png", dpi=240)
    plt.close()


def write_outputs(result: dict) -> None:
    OUTDIR.mkdir(exist_ok=True)
    rows = result["rows"]

    csv_path = OUTDIR / "local_selection_phi4_timeseries.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    final = rows[-1]
    post_rows = [row for row in rows if row["t"] >= result["params"]["perturb_time"]]

    def mean_ratio(key: str) -> float:
        baseline = np.mean([row[f"baseline_{key}"] for row in post_rows])
        maat = np.mean([row[f"maat_{key}"] for row in post_rows])
        return float(maat / baseline)

    summary = {
        "params": result["params"],
        "final": final,
        "improvements": {
            "distance_to_kink_ratio_maat_over_baseline": final["maat_distance_to_kink"] / final["baseline_distance_to_kink"],
            "residual_ratio_maat_over_baseline": final["maat_residual_rms"] / final["baseline_residual_rms"],
            "roughness_ratio_maat_over_baseline": final["maat_roughness"] / final["baseline_roughness"],
            "energy_ratio_maat_over_baseline": final["maat_energy"] / final["baseline_energy"],
        },
        "post_perturbation_mean_ratios": {
            "distance_to_kink_ratio_maat_over_baseline": mean_ratio("distance_to_kink"),
            "residual_ratio_maat_over_baseline": mean_ratio("residual_rms"),
            "roughness_ratio_maat_over_baseline": mean_ratio("roughness"),
            "energy_ratio_maat_over_baseline": mean_ratio("energy"),
        },
        "max_lambdas": {
            "H": max(row["max_lambda_H"] for row in rows),
            "S": max(row["max_lambda_S"] for row in rows),
            "R": max(row["max_lambda_R"] for row in rows),
        },
    }
    (OUTDIR / "local_selection_phi4_summary.json").write_text(json.dumps(summary, indent=2))

    plot_fields(result)
    plot_metrics(result)
    plot_lambda_snapshots(result)
    plot_lambda_heatmap(result)
    plot_final_defects(result)

    readme = f"""# MAAT v0.6 Local Selection Fields

This folder contains a minimal 1D phi^4 benchmark for local MAAT
selection-pressure fields.

The experiment promotes global weights to local fields:

```text
lambda_a -> lambda_a(x,t)
```

and evolves them with:

```text
tau_lambda partial_t lambda_a =
    D_lambda partial_xx lambda_a - lambda_a + lambda_a_star(x,t)
```

The local fixed point `lambda_star(x,t)` is computed from coarse-grained local
defect covariance.

## Outputs

| File | Meaning |
|------|---------|
| `local_selection_phi4_timeseries.csv` | Time-series metrics |
| `local_selection_phi4_summary.json` | Parameters and final diagnostics |
| `field_profiles.png` | Initial, clean, baseline final, and MAAT final fields |
| `metric_comparison.png` | Baseline vs MAAT metrics |
| `lambda_field_snapshots.png` | Local lambda field snapshots |
| `lambda_field_heatmap.png` | Compact lambda heatmap |
| `final_defect_profiles.png` | Final defects vs clean-kink target |

## Final Ratios

| Metric | MAAT / baseline |
|--------|----------------:|
| distance to kink | {summary['improvements']['distance_to_kink_ratio_maat_over_baseline']:.4f} |
| residual RMS | {summary['improvements']['residual_ratio_maat_over_baseline']:.4f} |
| roughness | {summary['improvements']['roughness_ratio_maat_over_baseline']:.4f} |
| energy | {summary['improvements']['energy_ratio_maat_over_baseline']:.4f} |

## Post-Perturbation Mean Ratios

| Metric | MAAT / baseline |
|--------|----------------:|
| distance to kink | {summary['post_perturbation_mean_ratios']['distance_to_kink_ratio_maat_over_baseline']:.4f} |
| residual RMS | {summary['post_perturbation_mean_ratios']['residual_ratio_maat_over_baseline']:.4f} |
| roughness | {summary['post_perturbation_mean_ratios']['roughness_ratio_maat_over_baseline']:.4f} |
| energy | {summary['post_perturbation_mean_ratios']['energy_ratio_maat_over_baseline']:.4f} |

This is an effective local-field benchmark, not a first-principles microscopic
derivation.
"""
    (OUTDIR / "README.md").write_text(readme)


def main() -> None:
    result = simulate()
    write_outputs(result)
    summary = json.loads((OUTDIR / "local_selection_phi4_summary.json").read_text())
    print(json.dumps({
        "output_dir": str(OUTDIR),
        "final": summary["final"],
        "improvements": summary["improvements"],
        "post_perturbation_mean_ratios": summary["post_perturbation_mean_ratios"],
        "max_lambdas": summary["max_lambdas"],
    }, indent=2))


if __name__ == "__main__":
    main()
