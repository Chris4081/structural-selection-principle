#!/usr/bin/env python3
"""MAAT v0.5 dynamic selection-pressure flow.

This experiment turns the v0.4 static response closure

    lambda* = (C + eta tr(C)/A I)^(-1) (d_bar - d_target)

into a dynamical relaxation equation

    tau d lambda / dt = -lambda + lambda*(t).

The goal is not to model a specific physical system, but to test whether
selection pressures can be treated as stable response fields driven by
defect-covariance feedback.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


OUTDIR = Path(__file__).resolve().parent / "outputs"
FIELDS = ["H", "B", "S", "V", "R"]
A = len(FIELDS)


def softplus(x: np.ndarray) -> np.ndarray:
    return np.log1p(np.exp(-np.abs(x))) + np.maximum(x, 0.0)


def make_base_ensemble(
    n_samples: int = 1200,
    seed: int = 42,
) -> np.ndarray:
    """Create a correlated toy defect ensemble in [0, 1]^5."""
    rng = np.random.default_rng(seed)
    mean = np.array([0.22, 0.18, 0.25, 0.20, 0.16])
    cov = np.array(
        [
            [0.040, 0.020, 0.010, 0.018, 0.012],
            [0.020, 0.035, 0.012, 0.014, 0.018],
            [0.010, 0.012, 0.045, 0.020, 0.010],
            [0.018, 0.014, 0.020, 0.042, 0.022],
            [0.012, 0.018, 0.010, 0.022, 0.038],
        ]
    )
    raw = rng.multivariate_normal(mean, cov, size=n_samples)
    return np.clip(raw, 0.0, 1.0)


def target_schedule(t: float) -> np.ndarray:
    """Piecewise target defect schedule with a boundary-pressure event."""
    base = np.array([0.12, 0.10, 0.14, 0.11, 0.10])

    # From t=45 onward, robustness target tightens strongly.
    if 45 <= t < 95:
        base = base.copy()
        base[4] = 0.035
        base[1] = 0.075

    # After t=95, the system relaxes to a moderate target.
    if t >= 95:
        base = np.array([0.10, 0.085, 0.12, 0.095, 0.075])

    return base


def selected_ensemble(defects: np.ndarray, lam: np.ndarray, beta: float) -> np.ndarray:
    """Return weights for exponential selection over defect samples."""
    cost = defects @ np.maximum(lam, 0.0)
    score = -beta * (cost - cost.min())
    w = np.exp(score)
    return w / w.sum()


def weighted_mean_and_cov(defects: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = weights @ defects
    centered = defects - mean
    cov = (centered * weights[:, None]).T @ centered
    return mean, cov


def lambda_star(
    cov: np.ndarray,
    mean: np.ndarray,
    target: np.ndarray,
    eta: float,
) -> np.ndarray:
    reg = eta * np.trace(cov) / A
    mat = cov + reg * np.eye(A)
    response = np.linalg.solve(mat, mean - target)
    return np.maximum(response, 0.0)


def loss(mean: np.ndarray, target: np.ndarray, lam: np.ndarray, eta_loss: float) -> float:
    return float(0.5 * np.sum((mean - target) ** 2) + 0.5 * eta_loss * np.sum(lam**2))


def simulate(
    steps: int = 1500,
    dt: float = 0.08,
    tau: float = 4.0,
    eta: float = 1e-2,
    eta_loss: float = 2e-4,
    beta: float = 2.0,
    seed: int = 42,
) -> dict:
    defects = make_base_ensemble(seed=seed)
    lam = np.zeros(A)
    rows = []

    for step in range(steps):
        t = step * dt
        target = target_schedule(t)
        weights = selected_ensemble(defects, lam, beta)
        mean, cov = weighted_mean_and_cov(defects, weights)
        star = lambda_star(cov, mean, target, eta)
        dlam = (-lam + star) / tau
        lam = np.maximum(lam + dt * dlam, 0.0)
        support = 1.0 / (1.0 + mean)
        stability = min(
            support[4],
            float(np.prod(support[:4]) ** 0.25),
        )

        eig = np.linalg.eigvalsh(cov)
        eig_safe = np.maximum(eig, 1e-12)
        kappa = eig_safe.max() / eig_safe.min()

        row = {
            "step": step,
            "t": t,
            "loss": loss(mean, target, lam, eta_loss),
            "lambda_norm": float(np.linalg.norm(lam)),
            "lambda_star_norm": float(np.linalg.norm(star)),
            "tracking_error": float(np.linalg.norm(lam - star)),
            "mean_defect_norm": float(np.linalg.norm(mean)),
            "target_defect_norm": float(np.linalg.norm(target)),
            "stability": stability,
            "cov_trace": float(np.trace(cov)),
            "cov_log_kappa": float(np.log10(kappa + 1.0)),
        }
        for i, name in enumerate(FIELDS):
            row[f"lambda_{name}"] = float(lam[i])
            row[f"lambda_star_{name}"] = float(star[i])
            row[f"mean_d_{name}"] = float(mean[i])
            row[f"target_d_{name}"] = float(target[i])
        rows.append(row)

    return {
        "params": {
            "steps": steps,
            "dt": dt,
            "tau": tau,
            "eta": eta,
            "eta_loss": eta_loss,
            "beta": beta,
            "seed": seed,
            "fields": FIELDS,
        },
        "rows": rows,
    }


def as_arrays(rows: list[dict]) -> dict[str, np.ndarray]:
    keys = rows[0].keys()
    return {key: np.array([row[key] for row in rows]) for key in keys}


def plot_lambdas(data: dict[str, np.ndarray]) -> None:
    plt.figure(figsize=(10, 6))
    for name in FIELDS:
        plt.plot(data["t"], data[f"lambda_{name}"], label=rf"$\lambda_{name}$")
    plt.axvspan(45, 95, color="tab:red", alpha=0.08, label="boundary-pressure target")
    plt.xlabel("time")
    plt.ylabel("selection pressure")
    plt.title("MAAT v0.5 dynamic selection pressures")
    plt.grid(alpha=0.3)
    plt.legend(ncol=2)
    plt.tight_layout()
    plt.savefig(OUTDIR / "lambda_trajectory.png", dpi=240)
    plt.close()


def plot_loss(data: dict[str, np.ndarray]) -> None:
    plt.figure(figsize=(9, 5.4))
    plt.plot(data["t"], data["loss"], label=r"$\mathcal{L}(\lambda)$")
    plt.plot(data["t"], data["tracking_error"], label=r"$\|\lambda-\lambda^\star\|$")
    plt.yscale("log")
    plt.axvspan(45, 95, color="tab:red", alpha=0.08)
    plt.xlabel("time")
    plt.ylabel("log scale")
    plt.title("Tracking loss and response error")
    plt.grid(alpha=0.3, which="both")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUTDIR / "loss_and_tracking.png", dpi=240)
    plt.close()


def plot_defects(data: dict[str, np.ndarray]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharex=True)
    for name in FIELDS:
        axes[0].plot(data["t"], data[f"mean_d_{name}"], label=fr"$\bar d_{name}$")
        axes[1].plot(data["t"], data[f"target_d_{name}"], label=fr"$d^\ast_{name}$")
    for ax in axes:
        ax.axvspan(45, 95, color="tab:red", alpha=0.08)
        ax.grid(alpha=0.3)
        ax.set_xlabel("time")
    axes[0].set_ylabel("defect")
    axes[0].set_title("Selected mean defects")
    axes[1].set_title("Target defects")
    axes[0].legend(fontsize=8, ncol=2)
    plt.tight_layout()
    plt.savefig(OUTDIR / "mean_vs_target_defects.png", dpi=240)
    plt.close()


def plot_phase(data: dict[str, np.ndarray]) -> None:
    plt.figure(figsize=(7.4, 5.8))
    scatter = plt.scatter(
        data["lambda_R"],
        data["lambda_norm"],
        c=data["t"],
        s=18,
        cmap="viridis",
    )
    plt.colorbar(scatter, label="time")
    plt.xlabel(r"$\lambda_R$")
    plt.ylabel(r"$\|\lambda\|$")
    plt.title("Response phase portrait")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "lambda_phase_portrait.png", dpi=240)
    plt.close()


def plot_stability(data: dict[str, np.ndarray]) -> None:
    fig, ax1 = plt.subplots(figsize=(9, 5.4))
    ax1.plot(data["t"], data["stability"], color="tab:green", label="MAAT stability")
    ax1.set_xlabel("time")
    ax1.set_ylabel("stability", color="tab:green")
    ax1.tick_params(axis="y", labelcolor="tab:green")
    ax1.grid(alpha=0.3)

    ax2 = ax1.twinx()
    ax2.plot(data["t"], data["cov_log_kappa"], color="tab:purple", alpha=0.75, label=r"$\log_{10}(1+\kappa)$")
    ax2.set_ylabel("covariance log conditioning", color="tab:purple")
    ax2.tick_params(axis="y", labelcolor="tab:purple")

    ax1.axvspan(45, 95, color="tab:red", alpha=0.08)
    plt.title("Stability and covariance conditioning during dynamic selection")
    fig.tight_layout()
    plt.savefig(OUTDIR / "stability_and_conditioning.png", dpi=240)
    plt.close()


def write_outputs(result: dict) -> None:
    OUTDIR.mkdir(exist_ok=True)
    rows = result["rows"]
    data = as_arrays(rows)

    import csv

    csv_path = OUTDIR / "lambda_dynamic_flow_timeseries.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "params": result["params"],
        "final": rows[-1],
        "min_loss": min(row["loss"] for row in rows),
        "max_lambda_R": max(row["lambda_R"] for row in rows),
        "mean_tracking_error_last_100": float(np.mean(data["tracking_error"][-100:])),
        "mean_loss_last_100": float(np.mean(data["loss"][-100:])),
    }
    (OUTDIR / "lambda_dynamic_flow_summary.json").write_text(json.dumps(summary, indent=2))

    plot_lambdas(data)
    plot_loss(data)
    plot_defects(data)
    plot_phase(data)
    plot_stability(data)

    readme = f"""# MAAT v0.5 Dynamic Lambda Flow

This folder contains a toy simulation of dynamic structural selection
pressures.

The simulated equation is:

```text
tau d lambda / dt = -lambda + lambda_star(t)
lambda_star(t) = (C(t) + eta tr(C(t))/A I)^(-1) (mean(d)(t) - d_target(t))
```

The interval `45 <= t < 95` introduces a boundary-pressure event by tightening
the target robustness defect `d_R*`.

## Main Outputs

| File | Meaning |
|------|---------|
| `lambda_dynamic_flow_timeseries.csv` | Full time series |
| `lambda_dynamic_flow_summary.json` | Parameters and final diagnostics |
| `lambda_trajectory.png` | Evolution of all lambda sectors |
| `loss_and_tracking.png` | Loss and tracking error |
| `mean_vs_target_defects.png` | Selected mean defects vs target defects |
| `lambda_phase_portrait.png` | Phase portrait in lambda space |
| `stability_and_conditioning.png` | Stability and covariance conditioning |

## Final Diagnostics

| Quantity | Value |
|----------|------:|
| final loss | {summary['final']['loss']:.6g} |
| min loss | {summary['min_loss']:.6g} |
| final lambda norm | {summary['final']['lambda_norm']:.6g} |
| max lambda_R | {summary['max_lambda_R']:.6g} |
| mean tracking error last 100 steps | {summary['mean_tracking_error_last_100']:.6g} |

Interpretation: v0.4 response closure appears as the moving fixed point of a
v0.5 dynamical relaxation system. This is an effective toy model, not yet a
microscopic physical derivation.
"""
    (OUTDIR / "README.md").write_text(readme)


def main() -> None:
    result = simulate()
    write_outputs(result)
    print(json.dumps({
        "output_dir": str(OUTDIR),
        "final": result["rows"][-1],
    }, indent=2))


if __name__ == "__main__":
    main()
