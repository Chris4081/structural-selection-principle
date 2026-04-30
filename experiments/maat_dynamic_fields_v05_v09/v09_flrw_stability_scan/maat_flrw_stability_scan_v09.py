#!/usr/bin/env python3
"""MAAT v0.9 FLRW stability scan for the kinetic structural-selection branch.

This is a toy effective-cosmology benchmark for the v0.8 scalar worked
example.  It scans the sigma=+1 branch of

    P(X, lambda) = X + mu lambda U(X)

with U(X) the normalised MAAT support cost of the kinetic defect
d_X = X / Lambda_X^4.

The goal is not to fit cosmological data.  The goal is to test whether
homogeneous FLRW trajectories can remain ghost-free, gradient-stable, and
background-safe over a small parameter grid.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "outputs"


def support_functions(x: float, eps: float = 1e-6) -> tuple[float, float, float]:
    """Return U, U_X, U_XX for dimensionless X/Lambda_X^4 = x.

    U is normalised so U(0)=0:

        U = -log((eps + Gamma)/(eps + 1)), Gamma = 1/(1+x).

    The first and second derivatives are with respect to x.
    """
    x = max(float(x), 0.0)
    gamma = 1.0 / (1.0 + x)
    u = -np.log((eps + gamma) / (eps + 1.0))
    ux = gamma**2 / (eps + gamma)
    uxx = -(gamma**3 * (2.0 * eps + gamma)) / (eps + gamma) ** 2
    return float(u), float(ux), float(uxx)


def model_quantities(
    state: np.ndarray,
    params: dict,
) -> dict:
    """Compute H, stability, densities, and equation-of-state diagnostics."""
    a, phi, pi, lam, vlam = state
    if a <= 0:
        return {"valid": False}

    mu = params["mu"]
    m_lam = params["m_lambda"]
    z_lam = params["Z_lambda"]
    rho_m0 = params["rho_m0"]
    rho_r0 = params["rho_r0"]
    rho_l = params["rho_lambda_cosmo"]

    x = 0.5 * pi**2
    u, ux, uxx = support_functions(x, params["epsilon"])

    p_x = 1.0 + mu * lam * ux
    p_xx = mu * lam * uxx
    kinetic_denominator = p_x + 2.0 * x * p_xx
    c_s2 = p_x / kinetic_denominator if kinetic_denominator != 0 else np.nan

    v_lam = 0.5 * m_lam**2 * lam**2
    rho_phi = x
    p_phi = x
    rho_maat = 0.5 * z_lam * vlam**2 + v_lam + mu * lam * (2.0 * x * ux - u)
    p_maat = 0.5 * z_lam * vlam**2 - v_lam + mu * lam * u

    rho_m = rho_m0 / a**3
    rho_r = rho_r0 / a**4
    rho_total = rho_m + rho_r + rho_l + rho_phi + rho_maat
    if rho_total <= 0 or not np.isfinite(rho_total):
        return {"valid": False}

    h = np.sqrt(rho_total / 3.0)
    omega_maat = rho_maat / rho_total
    w_maat = p_maat / rho_maat if abs(rho_maat) > 1e-14 else np.nan

    return {
        "valid": True,
        "a": a,
        "H": h,
        "X": x,
        "U": u,
        "U_X": ux,
        "U_XX": uxx,
        "P_X": p_x,
        "P_XX": p_xx,
        "D": kinetic_denominator,
        "c_s2": c_s2,
        "rho_phi": rho_phi,
        "p_phi": p_phi,
        "rho_maat": rho_maat,
        "p_maat": p_maat,
        "rho_total": rho_total,
        "omega_maat": omega_maat,
        "w_maat": w_maat,
        "rho_m": rho_m,
        "rho_r": rho_r,
        "rho_lambda_cosmo": rho_l,
    }


def rhs(state: np.ndarray, params: dict) -> np.ndarray:
    q = model_quantities(state, params)
    if not q["valid"]:
        return np.full_like(state, np.nan)

    a, phi, pi, lam, vlam = state
    mu = params["mu"]
    m_lam = params["m_lambda"]
    z_lam = params["Z_lambda"]

    h = q["H"]
    p_x = q["P_X"]
    d = q["D"]
    ux = q["U_X"]
    u = q["U"]

    if d <= 0 or p_x <= 0:
        return np.full_like(state, np.nan)

    adot = a * h
    phidot = pi
    # d/dt(P_X phi_dot) + 3 H P_X phi_dot = 0.
    pidot = -(3.0 * h * p_x * pi + mu * ux * vlam * pi) / d
    lamdot = vlam
    # Z (lambda_ddot + 3H lambda_dot) + V_lambda,lambda - mu U = 0.
    vlamdot = (mu * u - m_lam**2 * lam - 3.0 * h * z_lam * vlam) / z_lam

    return np.array([adot, phidot, pidot, lamdot, vlamdot], dtype=float)


def rk4_step(state: np.ndarray, dt: float, params: dict) -> np.ndarray:
    k1 = rhs(state, params)
    k2 = rhs(state + 0.5 * dt * k1, params)
    k3 = rhs(state + 0.5 * dt * k2, params)
    k4 = rhs(state + dt * k3, params)
    if not all(np.all(np.isfinite(k)) for k in (k1, k2, k3, k4)):
        return np.full_like(state, np.nan)
    return state + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0


def simulate_case(mu: float, pi0: float, *, record: bool = False) -> dict:
    params = {
        "mu": float(mu),
        "m_lambda": 0.55,
        "Z_lambda": 1.0,
        "epsilon": 1e-6,
        # Units: 8 pi G = 1 and H0 ~= 1.  These are reference background
        # densities, not a data fit.
        "rho_m0": 0.90,
        "rho_r0": 3.0e-4,
        "rho_lambda_cosmo": 2.10,
    }
    dt = 0.01
    steps = 900
    state = np.array([0.30, 0.0, float(pi0), 0.05, 0.0], dtype=float)

    rows = []
    fail_reason = ""
    for step in range(steps):
        t = step * dt
        q = model_quantities(state, params)
        if not q["valid"]:
            fail_reason = "invalid_density"
            break

        stable_now = (
            q["P_X"] > 0.0
            and q["D"] > 0.0
            and q["c_s2"] > 0.0
            and np.isfinite(q["c_s2"])
            and q["rho_maat"] >= -1e-10
        )
        if not stable_now:
            fail_reason = "stability_violation"
            break

        if record or step % 5 == 0:
            rows.append(
                {
                    "step": step,
                    "t": t,
                    "a": state[0],
                    "phi": state[1],
                    "pi": state[2],
                    "lambda": state[3],
                    "lambda_dot": state[4],
                    **{k: float(v) for k, v in q.items() if k != "valid"},
                }
            )

        new_state = rk4_step(state, dt, params)
        if not np.all(np.isfinite(new_state)):
            fail_reason = "integration_failure"
            break
        state = new_state

    if not rows:
        return {
            "mu": float(mu),
            "pi0": float(pi0),
            "stable": False,
            "fail_reason": fail_reason or "no_rows",
        }

    min_px = min(row["P_X"] for row in rows)
    min_d = min(row["D"] for row in rows)
    min_cs2 = min(row["c_s2"] for row in rows)
    max_cs2 = max(row["c_s2"] for row in rows)
    max_omega = max(row["omega_maat"] for row in rows)
    final = rows[-1]
    stable = fail_reason == ""
    background_safe = stable and max_omega < 0.05 and min_cs2 > 0.01 and max_cs2 < 5.0

    return {
        "mu": float(mu),
        "pi0": float(pi0),
        "stable": bool(stable),
        "background_safe": bool(background_safe),
        "fail_reason": fail_reason,
        "min_P_X": float(min_px),
        "min_D": float(min_d),
        "min_c_s2": float(min_cs2),
        "max_c_s2": float(max_cs2),
        "max_omega_maat": float(max_omega),
        "final_a": float(final["a"]),
        "final_lambda": float(final["lambda"]),
        "final_X": float(final["X"]),
        "final_w_maat": float(final["w_maat"]),
        "final_omega_maat": float(final["omega_maat"]),
        "rows": rows if record else None,
    }


def run_scan() -> tuple[list[dict], dict]:
    # The grid intentionally extends beyond the perturbative region so that
    # ghost/gradient and background-safety boundaries become visible.
    mu_values = np.logspace(-3, 3.0, 28)
    pi_values = np.linspace(0.03, 3.20, 28)
    results = []
    for mu in mu_values:
        for pi0 in pi_values:
            results.append(simulate_case(float(mu), float(pi0), record=False))

    stable_results = [r for r in results if r.get("stable")]
    safe_results = [r for r in results if r.get("background_safe")]
    conservative_safe = [
        r for r in stable_results
        if r["max_omega_maat"] < 0.03 and r["max_c_s2"] < 2.0
    ]
    if conservative_safe:
        representative = max(conservative_safe, key=lambda r: r["max_omega_maat"])
    elif safe_results:
        representative = max(safe_results, key=lambda r: r["max_omega_maat"])
    elif stable_results:
        representative = min(stable_results, key=lambda r: abs(r["max_omega_maat"] - 0.05))
    else:
        representative = results[0]

    rep = simulate_case(representative["mu"], representative["pi0"], record=True)
    return results, rep


def write_scan_csv(results: list[dict]) -> None:
    keys = [
        "mu",
        "pi0",
        "stable",
        "background_safe",
        "fail_reason",
        "min_P_X",
        "min_D",
        "min_c_s2",
        "max_c_s2",
        "max_omega_maat",
        "final_a",
        "final_lambda",
        "final_X",
        "final_w_maat",
        "final_omega_maat",
    ]
    with (OUTDIR / "flrw_stability_scan_results.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for r in results:
            writer.writerow({k: r.get(k, "") for k in keys})


def grid_from_results(results: list[dict], key: str, default: float = np.nan) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mus = np.array(sorted({r["mu"] for r in results}))
    pis = np.array(sorted({r["pi0"] for r in results}))
    grid = np.full((len(pis), len(mus)), default, dtype=float)
    mu_index = {v: i for i, v in enumerate(mus)}
    pi_index = {v: i for i, v in enumerate(pis)}
    for r in results:
        val = r.get(key, default)
        if val == "":
            val = default
        grid[pi_index[r["pi0"]], mu_index[r["mu"]]] = float(val)
    return mus, pis, grid


def plot_scan(results: list[dict]) -> None:
    mus, pis, max_omega = grid_from_results(results, "max_omega_maat")
    _, _, max_cs2 = grid_from_results(results, "max_c_s2")
    _, _, final_w = grid_from_results(results, "final_w_maat")
    _, _, stable = grid_from_results(results, "stable", default=0.0)
    _, _, safe = grid_from_results(results, "background_safe", default=0.0)

    x_edges = np.geomspace(mus.min() / 1.2, mus.max() * 1.2, len(mus) + 1)
    y_step = pis[1] - pis[0]
    y_edges = np.linspace(pis.min() - y_step / 2, pis.max() + y_step / 2, len(pis) + 1)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8.8), sharex=True, sharey=True)
    plots = [
        (np.log10(np.maximum(max_omega, 1e-12)), r"$\log_{10}\max\Omega_{\rm MAAT}$", "viridis"),
        (np.clip(max_cs2, 0, 10), r"$\max c_s^2$ (clipped at 10)", "magma"),
        (np.clip(final_w, -2, 2), r"final $w_{\rm MAAT}$", "coolwarm"),
        (stable + safe, "stability mask (1 stable, 2 safe)", "cividis"),
    ]
    for ax, (data, title, cmap) in zip(axes.ravel(), plots):
        im = ax.pcolormesh(x_edges, y_edges, data, shading="auto", cmap=cmap)
        ax.set_xscale("log")
        ax.set_title(title)
        ax.set_ylabel(r"initial $\dot\phi$")
        ax.grid(alpha=0.25)
        fig.colorbar(im, ax=ax)
    axes[1, 0].set_xlabel(r"selection strength $\mu=M_\ast^4/\Lambda_X^4$")
    axes[1, 1].set_xlabel(r"selection strength $\mu=M_\ast^4/\Lambda_X^4$")
    fig.suptitle("MAAT v0.9 FLRW stability scan", y=0.995)
    fig.tight_layout()
    fig.savefig(OUTDIR / "flrw_stability_scan_heatmaps.png", dpi=240)
    plt.close(fig)


def plot_representative(rep: dict) -> None:
    rows = rep["rows"]
    t = np.array([row["t"] for row in rows])
    a = np.array([row["a"] for row in rows])

    fig, axes = plt.subplots(2, 2, figsize=(12, 8.4), sharex=True)
    axes[0, 0].plot(t, a, color="black")
    axes[0, 0].set_ylabel("scale factor a(t)")
    axes[0, 0].set_title("Representative stable trajectory")

    axes[0, 1].plot(t, [row["omega_maat"] for row in rows], color="tab:blue", label=r"$\Omega_{\rm MAAT}$")
    axes[0, 1].plot(t, [row["omega_maat"] + row["rho_phi"] / row["rho_total"] for row in rows],
                    color="tab:cyan", ls="--", label=r"$\Omega_{\rm MAAT}+\Omega_\phi$")
    axes[0, 1].set_ylabel("density fraction")
    axes[0, 1].legend()

    axes[1, 0].plot(t, [row["P_X"] for row in rows], label=r"$P_X$", color="tab:green")
    axes[1, 0].plot(t, [row["c_s2"] for row in rows], label=r"$c_s^2$", color="tab:purple")
    axes[1, 0].axhline(0, color="black", lw=0.8)
    axes[1, 0].set_ylabel("stability")
    axes[1, 0].set_xlabel("time")
    axes[1, 0].legend()

    axes[1, 1].plot(t, np.clip([row["w_maat"] for row in rows], -2, 2), color="tab:red")
    axes[1, 1].axhline(-1, color="black", lw=0.8, ls=":")
    axes[1, 1].axhline(0, color="black", lw=0.8, ls="--")
    axes[1, 1].set_ylabel(r"$w_{\rm MAAT}$ clipped to $[-2,2]$")
    axes[1, 1].set_xlabel("time")

    for ax in axes.ravel():
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUTDIR / "representative_flrw_trajectory.png", dpi=240)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    axes[0].plot(t, [row["lambda"] for row in rows], color="tab:orange")
    axes[0].set_title(r"Response field $\lambda(t)$")
    axes[0].set_xlabel("time")
    axes[0].set_ylabel(r"$\lambda$")

    axes[1].plot(t, [row["X"] for row in rows], color="tab:gray")
    axes[1].set_title(r"Kinetic invariant $X(t)$")
    axes[1].set_xlabel("time")
    axes[1].set_ylabel("X")

    for ax in axes:
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUTDIR / "lambda_and_kinetic_evolution.png", dpi=240)
    plt.close(fig)


def phase_space_plot() -> None:
    xs = np.linspace(0.0, 6.0, 360)
    qs = np.linspace(0.0, 60.0, 360)
    cs2 = np.zeros((len(xs), len(qs)))
    px = np.zeros_like(cs2)
    stable = np.zeros_like(cs2)
    for i, x in enumerate(xs):
        _, ux, uxx = support_functions(x)
        for j, q in enumerate(qs):
            p_x = 1.0 + q * ux
            d = p_x + 2.0 * x * q * uxx
            val = p_x / d if d != 0 else np.nan
            px[i, j] = p_x
            cs2[i, j] = val
            stable[i, j] = 1.0 if p_x > 0 and d > 0 and val > 0 else 0.0

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    im0 = axes[0].imshow(
        np.clip(cs2, 0, 3),
        extent=[qs.min(), qs.max(), xs.min(), xs.max()],
        aspect="auto",
        origin="lower",
        cmap="magma",
    )
    axes[0].set_title(r"Sound speed $c_s^2$ in $(X,\mu\lambda)$ space")
    axes[0].set_xlabel(r"$q=\mu\lambda$")
    axes[0].set_ylabel(r"$X/\Lambda_X^4$")
    fig.colorbar(im0, ax=axes[0])

    im1 = axes[1].imshow(
        stable,
        extent=[qs.min(), qs.max(), xs.min(), xs.max()],
        aspect="auto",
        origin="lower",
        cmap="Greens",
        vmin=0,
        vmax=1,
    )
    axes[1].set_title("No-ghost and gradient-stable region")
    axes[1].set_xlabel(r"$q=\mu\lambda$")
    fig.colorbar(im1, ax=axes[1])

    for ax in axes:
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUTDIR / "kinetic_branch_phase_space.png", dpi=240)
    plt.close(fig)


def write_outputs(results: list[dict], rep: dict) -> None:
    OUTDIR.mkdir(exist_ok=True)
    write_scan_csv(results)

    with (OUTDIR / "representative_flrw_trajectory.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rep["rows"][0].keys()))
        writer.writeheader()
        writer.writerows(rep["rows"])

    stable_count = sum(1 for r in results if r.get("stable"))
    safe_count = sum(1 for r in results if r.get("background_safe"))
    total = len(results)
    stable_results = [r for r in results if r.get("stable")]
    summary = {
        "total_cases": total,
        "stable_cases": stable_count,
        "background_safe_cases": safe_count,
        "stable_fraction": stable_count / total,
        "background_safe_fraction": safe_count / total,
        "representative": {k: v for k, v in rep.items() if k != "rows"},
        "stable_ranges": {
            "mu_min": min((r["mu"] for r in stable_results), default=None),
            "mu_max": max((r["mu"] for r in stable_results), default=None),
            "pi0_min": min((r["pi0"] for r in stable_results), default=None),
            "pi0_max": max((r["pi0"] for r in stable_results), default=None),
        },
    }
    (OUTDIR / "flrw_stability_scan_summary.json").write_text(json.dumps(summary, indent=2))

    plot_scan(results)
    plot_representative(rep)
    phase_space_plot()

    readme = f"""# MAAT v0.9 FLRW Stability Scan

This folder contains a toy FLRW stability scan for the stable sigma=+1 branch
of the v0.8 scalar worked example.

The scanned model is:

```text
P(X, lambda) = X + mu lambda U(X)
d_X = X / Lambda_X^4
U(X) = -log((epsilon + Gamma)/(epsilon + 1))
Gamma = 1/(1+d_X)
```

The scan checks:

- no ghost: `P_X > 0`
- gradient stability: `c_s^2 > 0`
- positive MAAT energy density
- background safety: `max Omega_MAAT < 0.05`

## Outputs

| File | Meaning |
|------|---------|
| `flrw_stability_scan_results.csv` | Full parameter scan |
| `flrw_stability_scan_summary.json` | Aggregate diagnostics |
| `representative_flrw_trajectory.csv` | Representative stable trajectory |
| `flrw_stability_scan_heatmaps.png` | Scan heatmaps |
| `representative_flrw_trajectory.png` | Stable trajectory observables |
| `lambda_and_kinetic_evolution.png` | Lambda and kinetic invariant |
| `kinetic_branch_phase_space.png` | Analytic stability map in `(X, mu lambda)` |

## Summary

| Quantity | Value |
|----------|------:|
| total cases | {summary['total_cases']} |
| stable cases | {summary['stable_cases']} |
| background-safe cases | {summary['background_safe_cases']} |
| stable fraction | {summary['stable_fraction']:.4f} |
| background-safe fraction | {summary['background_safe_fraction']:.4f} |

This is not a cosmological data fit. It is a stability and consistency scan
for the kinetic structural-selection branch.
"""
    (OUTDIR / "README.md").write_text(readme)


def main() -> None:
    results, rep = run_scan()
    write_outputs(results, rep)
    summary = json.loads((OUTDIR / "flrw_stability_scan_summary.json").read_text())
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
