#!/usr/bin/env python3
"""
MAAT Paper 45 — Effective Metric Response and Gravitational Slip
================================================================

This script extends the Paper-43 linear-growth benchmark from a growth-only
effective coupling

    mu(z) = G_eff/G = 1 + eta_g * C_hat_proj(z)

to a minimal metric-response benchmark with gravitational slip

    eta_slip(z) = Phi/Psi = 1 + beta_s * C_hat_proj(z)

and Weyl/lensing response

    Sigma(z) = mu(z) * [1 + eta_slip(z)] / 2.

The benchmark is intentionally conservative:

  * no Boltzmann code,
  * no weak-lensing likelihood,
  * no CMB calculation,
  * no claim of modified-gravity evidence.

It tests whether the existing MAAT projection kernel can be embedded into the
standard linear-perturbation metric-response parameterization while remaining
bounded and stable.

Run:
    python3 maat_metric_response_solver_v01.py
"""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "maat_paper45_matplotlib"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


BASE_DIR = Path(__file__).resolve().parent
OUT = BASE_DIR / "outputs_paper45"
OUT.mkdir(exist_ok=True)

EPS = 1e-12


@dataclass(frozen=True)
class Cosmology:
    H0: float = 67.4
    Omega_m0: float = 0.315
    Omega_L0: float = 0.685
    sigma8_0: float = 0.811


@dataclass(frozen=True)
class ProjectionParams:
    gamma_lambda: float = 1.3305071182030523
    Bstar_lambda: float = 4.5090578827501755
    alpha_lambda: float = 1.9936731161592685
    latent_zc: float = 1.1
    latent_sharpness: float = 3.5
    latent_floor: float = 0.20


COSMO = Cosmology()
PROJ = ProjectionParams()


# Compact f sigma_8 comparison set used by Papers 40, 42, and 43.
GROWTH_DATA = np.array([
    [0.020, 0.428, 0.046],
    [0.067, 0.423, 0.055],
    [0.100, 0.370, 0.130],
    [0.170, 0.510, 0.060],
    [0.220, 0.420, 0.070],
    [0.250, 0.351, 0.058],
    [0.370, 0.460, 0.038],
    [0.410, 0.450, 0.040],
    [0.570, 0.444, 0.038],
    [0.600, 0.430, 0.040],
    [0.780, 0.380, 0.040],
    [0.800, 0.470, 0.080],
    [1.400, 0.482, 0.116],
])


def E2_of_a(a: np.ndarray | float, c: Cosmology = COSMO) -> np.ndarray | float:
    a = np.asarray(a)
    return c.Omega_m0 * a ** (-3) + c.Omega_L0


def omega_m_a(a: np.ndarray | float, c: Cosmology = COSMO) -> np.ndarray | float:
    a = np.asarray(a)
    return c.Omega_m0 * a ** (-3) / np.maximum(E2_of_a(a, c), EPS)


def dlnH_dN(a: float, c: Cosmology = COSMO) -> float:
    return float(-1.5 * omega_m_a(a, c))


def E_lcdm_z(z: np.ndarray | float, c: Cosmology = COSMO) -> np.ndarray | float:
    z = np.asarray(z)
    return np.sqrt(c.Omega_m0 * (1 + z) ** 3 + c.Omega_L0)


def omega_m_z(z: np.ndarray | float, c: Cosmology = COSMO) -> np.ndarray | float:
    z = np.asarray(z)
    return c.Omega_m0 * (1 + z) ** 3 / np.maximum(E_lcdm_z(z, c) ** 2, EPS)


def growth_factor_approx(z: np.ndarray | float, c: Cosmology = COSMO) -> np.ndarray | float:
    """Carroll-Press-Turner style growth approximation, normalized to D(0)=1."""
    z = np.asarray(z)
    om = omega_m_z(z, c)
    ol = c.Omega_L0 / np.maximum(E_lcdm_z(z, c) ** 2, EPS)
    g = (5 * om / 2) / (
        om ** (4 / 7)
        - ol
        + (1 + om / 2) * (1 + ol / 70)
    )
    g0 = (5 * c.Omega_m0 / 2) / (
        c.Omega_m0 ** (4 / 7)
        - c.Omega_L0
        + (1 + c.Omega_m0 / 2) * (1 + c.Omega_L0 / 70)
    )
    return g / (g0 * (1 + z))


def latent_depth(z: np.ndarray | float, p: ProjectionParams = PROJ) -> np.ndarray | float:
    z = np.asarray(z)
    return p.latent_floor + (1 - p.latent_floor) / (
        1 + np.exp(p.latent_sharpness * (z - p.latent_zc))
    )


def projection_raw(z: np.ndarray | float, p: ProjectionParams = PROJ) -> np.ndarray | float:
    z = np.asarray(z)
    b_exp = E_lcdm_z(z) * (1 + z)
    d_acc = growth_factor_approx(z) * latent_depth(z, p)
    s_proj = b_exp * np.tanh((b_exp / (p.Bstar_lambda + EPS)) ** p.alpha_lambda)
    return s_proj / (1 + p.gamma_lambda * d_acc + EPS)


def projection_kernel(z: np.ndarray, zmax: float = 3.0) -> np.ndarray:
    c_raw = projection_raw(z)
    c0 = float(projection_raw(0.0))
    cmax = float(projection_raw(zmax))
    return np.clip((c_raw - c0) / (cmax - c0 + EPS), 0.0, 1.0)


def mu_of_z(z: np.ndarray | float, eta_g: float) -> np.ndarray | float:
    z = np.asarray(z)
    return 1.0 + eta_g * projection_kernel(z)


def slip_of_z(z: np.ndarray | float, beta_s: float) -> np.ndarray | float:
    z = np.asarray(z)
    return 1.0 + beta_s * projection_kernel(z)


def sigma_of_z(z: np.ndarray | float, eta_g: float, beta_s: float) -> np.ndarray | float:
    mu = mu_of_z(z, eta_g)
    slip = slip_of_z(z, beta_s)
    return mu * (1.0 + slip) / 2.0


def mu_of_a(a: float, eta_g: float) -> float:
    z = 1.0 / a - 1.0
    return float(mu_of_z(np.array([z]), eta_g)[0])


def growth_rhs(N: float, y: np.ndarray, eta_g: float) -> np.ndarray:
    a = float(np.exp(N))
    delta, delta_N = y
    friction = 2.0 + dlnH_dN(a)
    source = 1.5 * float(omega_m_a(a)) * mu_of_a(a, eta_g)
    return np.array([delta_N, -friction * delta_N + source * delta])


def integrate_growth(eta_g: float, a_ini: float = 0.01, n_steps: int = 2600) -> dict:
    N = np.linspace(np.log(a_ini), 0.0, n_steps)
    dN = float(N[1] - N[0])
    y = np.array([a_ini, a_ini], dtype=float)
    delta = np.empty(n_steps)
    delta_N = np.empty(n_steps)
    delta[0], delta_N[0] = y

    for i in range(1, n_steps):
        n0 = float(N[i - 1])
        k1 = growth_rhs(n0, y, eta_g)
        k2 = growth_rhs(n0 + 0.5 * dN, y + 0.5 * dN * k1, eta_g)
        k3 = growth_rhs(n0 + 0.5 * dN, y + 0.5 * dN * k2, eta_g)
        k4 = growth_rhs(n0 + dN, y + dN * k3, eta_g)
        y = y + (dN / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        delta[i], delta_N[i] = y

    a = np.exp(N)
    z = 1.0 / a - 1.0
    D = delta / np.maximum(delta[-1], EPS)
    f = delta_N / np.maximum(delta, EPS)
    fs8 = f * D * COSMO.sigma8_0
    chat = projection_kernel(z)
    mu = 1.0 + eta_g * chat

    return {
        "N": N,
        "a": a,
        "z": z,
        "delta": delta,
        "D": D,
        "f": f,
        "fsigma8": fs8,
        "mu": mu,
        "C_hat": chat,
    }


def interp_at_z(z_query: np.ndarray, z_grid_desc: np.ndarray, values: np.ndarray) -> np.ndarray:
    return np.interp(z_query, z_grid_desc[::-1], values[::-1])


def chi2_to_growth_data(model: dict) -> float:
    z = GROWTH_DATA[:, 0]
    obs = GROWTH_DATA[:, 1]
    err = GROWTH_DATA[:, 2]
    pred = interp_at_z(z, model["z"], model["fsigma8"])
    return float(np.sum(((obs - pred) / err) ** 2))


def main() -> None:
    print("=" * 72)
    print("MAAT Paper 45 — Effective Metric Response and Gravitational Slip")
    print("=" * 72)

    eta_values = np.linspace(0.0, 0.08, 41)
    beta_values = np.linspace(-0.06, 0.06, 61)
    representative_eta = 0.02
    representative_beta = 0.03

    lcdm = integrate_growth(0.0)
    growth_models = {float(eta): integrate_growth(float(eta)) for eta in eta_values}
    representative = growth_models[representative_eta]
    z_mask = lcdm["z"] <= 3.0
    z_curves = lcdm["z"][z_mask]

    dD_rep = (representative["D"] - lcdm["D"]) / np.maximum(lcdm["D"], EPS)
    dfs_rep = (
        (representative["fsigma8"] - lcdm["fsigma8"])
        / np.maximum(np.abs(lcdm["fsigma8"]), EPS)
    )

    slip_rep = slip_of_z(representative["z"], representative_beta)
    sigma_rep = sigma_of_z(representative["z"], representative_eta, representative_beta)
    weyl_rep = sigma_rep * representative["D"] / np.maximum(lcdm["D"], EPS)
    eg_rep = sigma_rep * lcdm["f"] / np.maximum(representative["f"], EPS)

    scan_rows = []
    for eta_g in eta_values:
        model = growth_models[float(eta_g)]
        chi2 = chi2_to_growth_data(model)
        for beta_s in beta_values:
            slip = slip_of_z(model["z"], float(beta_s))
            sigma = sigma_of_z(model["z"], float(eta_g), float(beta_s))
            weyl = sigma * model["D"] / np.maximum(lcdm["D"], EPS)
            eg_proxy = sigma * lcdm["f"] / np.maximum(model["f"], EPS)
            stable = bool(
                np.all(model["mu"] > 0.0)
                and np.all(slip > 0.0)
                and np.all(sigma > 0.0)
                and np.all(model["D"] > 0.0)
                and np.all(model["f"] > 0.0)
                and np.all(np.isfinite(weyl))
            )
            scan_rows.append([
                float(eta_g),
                float(beta_s),
                chi2,
                float(np.max(np.abs(model["mu"][z_mask] - 1.0)) * 100),
                float(np.max(np.abs(slip[z_mask] - 1.0)) * 100),
                float(np.max(np.abs(sigma[z_mask] - 1.0)) * 100),
                float(np.max(np.abs(weyl[z_mask] - 1.0)) * 100),
                float(np.max(np.abs(eg_proxy[z_mask] - 1.0)) * 100),
                stable,
            ])

    scan = np.array(scan_rows, dtype=object)
    stable_count = int(np.sum(scan[:, 8].astype(bool)))
    best_growth_idx = int(np.argmin(scan[:, 2].astype(float)))
    best_growth = scan[best_growth_idx]

    curves = np.column_stack([
        z_curves,
        lcdm["a"][z_mask],
        lcdm["D"][z_mask],
        representative["D"][z_mask],
        lcdm["f"][z_mask],
        representative["f"][z_mask],
        lcdm["fsigma8"][z_mask],
        representative["fsigma8"][z_mask],
        representative["C_hat"][z_mask],
        representative["mu"][z_mask],
        slip_rep[z_mask],
        sigma_rep[z_mask],
        weyl_rep[z_mask],
        eg_rep[z_mask],
        dD_rep[z_mask] * 100,
        dfs_rep[z_mask] * 100,
        (representative["mu"][z_mask] - 1.0) * 100,
        (slip_rep[z_mask] - 1.0) * 100,
        (sigma_rep[z_mask] - 1.0) * 100,
        (weyl_rep[z_mask] - 1.0) * 100,
        (eg_rep[z_mask] - 1.0) * 100,
    ])
    np.savetxt(
        OUT / "paper45_metric_curves.csv",
        curves,
        delimiter=",",
        header=(
            "z,a,D_LCDM,D_MAAT_eta002,f_LCDM,f_MAAT_eta002,"
            "fsigma8_LCDM,fsigma8_MAAT_eta002,C_hat,mu_eta002,"
            "eta_slip_beta003,Sigma_eta002_beta003,Weyl_proxy,EG_proxy,"
            "delta_D_pct,delta_fsigma8_pct,delta_mu_pct,delta_slip_pct,"
            "delta_Sigma_pct,delta_Weyl_pct,delta_EG_pct"
        ),
        comments="",
    )

    np.savetxt(
        OUT / "paper45_metric_scan.csv",
        scan,
        delimiter=",",
        fmt=["%.8g", "%.8g", "%.10g", "%.10g", "%.10g", "%.10g", "%.10g", "%.10g", "%s"],
        header=(
            "eta_g,beta_s,chi2_growth,max_abs_delta_mu_pct,"
            "max_abs_delta_slip_pct,max_abs_delta_Sigma_pct,"
            "max_abs_delta_Weyl_pct,max_abs_delta_EG_pct,stable"
        ),
        comments="",
    )

    z_data = GROWTH_DATA[:, 0]
    data_table = np.column_stack([
        z_data,
        GROWTH_DATA[:, 1],
        GROWTH_DATA[:, 2],
        interp_at_z(z_data, lcdm["z"], lcdm["fsigma8"]),
        interp_at_z(z_data, representative["z"], representative["fsigma8"]),
        interp_at_z(z_data, representative["z"], representative["mu"]),
        interp_at_z(z_data, representative["z"], slip_rep),
        interp_at_z(z_data, representative["z"], sigma_rep),
    ])
    np.savetxt(
        OUT / "paper45_growth_metric_comparison.csv",
        data_table,
        delimiter=",",
        header=(
            "z,fsigma8_obs,sigma,fsigma8_LCDM,fsigma8_MAAT_eta002,"
            "mu_eta002,eta_slip_beta003,Sigma_eta002_beta003"
        ),
        comments="",
    )

    plt.rcParams.update({"font.size": 11, "figure.dpi": 150, "savefig.dpi": 165})
    order = np.argsort(z_curves)

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    ax = axes[0, 0]
    ax.plot(z_curves[order], representative["mu"][z_mask][order], label=r"$\mu(z)$", lw=2)
    ax.plot(z_curves[order], slip_rep[z_mask][order], label=r"$\eta_{\rm slip}(z)$", lw=2)
    ax.plot(z_curves[order], sigma_rep[z_mask][order], label=r"$\Sigma(z)$", lw=2)
    ax.axhline(1.0, color="k", ls=":", lw=0.9)
    ax.set_xlabel("z")
    ax.set_ylabel("metric response")
    ax.set_title(r"Representative metric response: $\eta_g=0.02$, $\beta_s=0.03$")
    ax.grid(alpha=0.25)
    ax.legend()

    ax = axes[0, 1]
    ax.errorbar(GROWTH_DATA[:, 0], GROWTH_DATA[:, 1], yerr=GROWTH_DATA[:, 2], fmt="o", capsize=3, label="growth data")
    ax.plot(z_curves[order], lcdm["fsigma8"][z_mask][order], label="LCDM", lw=2)
    ax.plot(z_curves[order], representative["fsigma8"][z_mask][order], "--", label=r"MAAT $\eta_g=0.02$", lw=2)
    ax.set_xlabel("z")
    ax.set_ylabel(r"$f\sigma_8(z)$")
    ax.set_title("Growth data comparison")
    ax.grid(alpha=0.25)
    ax.legend()

    ax = axes[1, 0]
    ax.plot(z_curves[order], dD_rep[z_mask][order] * 100, label=r"$\Delta D/D$", lw=2)
    ax.plot(z_curves[order], dfs_rep[z_mask][order] * 100, label=r"$\Delta f\sigma_8/f\sigma_8$", lw=2)
    ax.plot(z_curves[order], (weyl_rep[z_mask][order] - 1.0) * 100, label=r"$\Delta$ Weyl proxy", lw=2)
    ax.plot(z_curves[order], (eg_rep[z_mask][order] - 1.0) * 100, label=r"$\Delta E_G$ proxy", lw=2)
    ax.axhline(0, color="k", ls=":", lw=0.9)
    ax.set_xlabel("z")
    ax.set_ylabel("relative deviation [%]")
    ax.set_title("Growth-lensing consistency proxies")
    ax.grid(alpha=0.25)
    ax.legend()

    ax = axes[1, 1]
    ax.plot(z_curves[order], representative["C_hat"][z_mask][order], lw=2.4, label=r"$\widehat{C}_{\rm proj}$")
    ax.plot(z_curves[order], representative["mu"][z_mask][order] - 1.0, lw=2, label=r"$\mu-1$")
    ax.plot(z_curves[order], slip_rep[z_mask][order] - 1.0, lw=2, label=r"$\eta_{\rm slip}-1$")
    ax.set_xlabel("z")
    ax.set_ylabel("bounded source")
    ax.set_title("Projection kernel sources both response channels")
    ax.grid(alpha=0.25)
    ax.legend()

    fig.tight_layout()
    fig.savefig(OUT / "fig1_metric_response_summary.png", bbox_inches="tight")
    plt.close(fig)

    # Heatmaps over the metric-response grid.
    eta_grid = eta_values
    beta_grid = beta_values
    sigma_grid = np.empty((len(beta_grid), len(eta_grid)))
    weyl_grid = np.empty_like(sigma_grid)
    eg_grid = np.empty_like(sigma_grid)
    chi_grid = np.empty_like(sigma_grid)

    for row in scan_rows:
        eta_g, beta_s = float(row[0]), float(row[1])
        i = int(np.where(np.isclose(beta_grid, beta_s))[0][0])
        j = int(np.where(np.isclose(eta_grid, eta_g))[0][0])
        chi_grid[i, j] = float(row[2])
        sigma_grid[i, j] = float(row[5])
        weyl_grid[i, j] = float(row[6])
        eg_grid[i, j] = float(row[7])

    fig2, axes2 = plt.subplots(1, 3, figsize=(16, 4.8), sharey=True)
    heatmaps = [
        (sigma_grid, r"max $|\Sigma-1|$ [%]"),
        (weyl_grid, r"max Weyl-proxy deviation [%]"),
        (eg_grid, r"max $E_G$-proxy deviation [%]"),
    ]
    for ax, (grid, title) in zip(axes2, heatmaps):
        im = ax.imshow(
            grid,
            origin="lower",
            aspect="auto",
            extent=[eta_grid.min(), eta_grid.max(), beta_grid.min(), beta_grid.max()],
            cmap="magma",
        )
        ax.scatter([representative_eta], [representative_beta], c="cyan", edgecolors="black", s=70, marker="*")
        ax.set_xlabel(r"$\eta_g$")
        ax.set_title(title)
        ax.grid(alpha=0.15)
        cbar = fig2.colorbar(im, ax=ax)
        cbar.ax.set_ylabel("%")
    axes2[0].set_ylabel(r"$\beta_s$")
    fig2.suptitle("Metric-response scan: bounded perturbative deviations")
    fig2.tight_layout()
    fig2.savefig(OUT / "fig2_metric_response_scan.png", bbox_inches="tight")
    plt.close(fig2)

    fig3, ax = plt.subplots(figsize=(9.0, 6.6))
    ax.scatter(
        scan[:, 3].astype(float),
        scan[:, 5].astype(float),
        c=scan[:, 6].astype(float),
        s=18,
        alpha=0.7,
        cmap="viridis",
        edgecolors="none",
    )
    ax.scatter(
        [np.max(np.abs(representative["mu"][z_mask] - 1.0)) * 100],
        [np.max(np.abs(sigma_rep[z_mask] - 1.0)) * 100],
        c="red",
        edgecolors="black",
        marker="*",
        s=160,
        label="representative",
    )
    ax.set_xlabel(r"max $|\mu-1|$ [%]")
    ax.set_ylabel(r"max $|\Sigma-1|$ [%]")
    ax.set_title("Growth response and metric response are distinct channels")
    ax.grid(alpha=0.25)
    ax.legend()
    cbar = fig3.colorbar(ax.collections[0], ax=ax)
    cbar.set_label("max Weyl-proxy deviation [%]")
    fig3.tight_layout()
    fig3.savefig(OUT / "fig3_growth_vs_metric_response.png", bbox_inches="tight")
    plt.close(fig3)

    summary = {
        "model": "MAAT Paper 45 Effective Metric Response Benchmark",
        "cosmology": asdict(COSMO),
        "projection_params": asdict(PROJ),
        "definitions": {
            "mu": "mu(z)=1+eta_g*C_hat_proj(z)",
            "eta_slip": "eta_slip(z)=Phi/Psi=1+beta_s*C_hat_proj(z)",
            "Sigma": "Sigma(z)=mu(z)*(1+eta_slip(z))/2",
            "Weyl_proxy": "Sigma(z)*D_MAAT(z)/D_LCDM(z)",
            "EG_proxy": "Sigma(z)*f_LCDM(z)/f_MAAT(z)",
        },
        "scan": {
            "eta_g_min": float(eta_values.min()),
            "eta_g_max": float(eta_values.max()),
            "eta_g_count": int(len(eta_values)),
            "beta_s_min": float(beta_values.min()),
            "beta_s_max": float(beta_values.max()),
            "beta_s_count": int(len(beta_values)),
            "total_branches": int(len(scan_rows)),
            "stable_count": stable_count,
            "best_growth_eta_g": float(best_growth[0]),
            "best_growth_beta_s_first_listed": float(best_growth[1]),
            "best_growth_beta_s_note": (
                "beta_s does not enter the growth-only chi2 and is therefore "
                "unconstrained by the compact f sigma_8 comparison."
            ),
            "best_growth_chi2": float(best_growth[2]),
        },
        "representative": {
            "eta_g": representative_eta,
            "beta_s": representative_beta,
            "chi2_growth": chi2_to_growth_data(representative),
            "max_abs_delta_D_pct_z_le_3": float(np.max(np.abs(dD_rep[z_mask])) * 100),
            "max_abs_delta_fsigma8_pct_z_le_3": float(np.max(np.abs(dfs_rep[z_mask])) * 100),
            "max_abs_delta_mu_pct_z_le_3": float(np.max(np.abs(representative["mu"][z_mask] - 1.0)) * 100),
            "max_abs_delta_slip_pct_z_le_3": float(np.max(np.abs(slip_rep[z_mask] - 1.0)) * 100),
            "max_abs_delta_Sigma_pct_z_le_3": float(np.max(np.abs(sigma_rep[z_mask] - 1.0)) * 100),
            "max_abs_delta_Weyl_proxy_pct_z_le_3": float(np.max(np.abs(weyl_rep[z_mask] - 1.0)) * 100),
            "max_abs_delta_EG_proxy_pct_z_le_3": float(np.max(np.abs(eg_rep[z_mask] - 1.0)) * 100),
            "mu_positive": bool(np.all(representative["mu"] > 0.0)),
            "slip_positive": bool(np.all(slip_rep > 0.0)),
            "Sigma_positive": bool(np.all(sigma_rep > 0.0)),
        },
        "growth_data": {
            "n_points": int(len(GROWTH_DATA)),
            "note": "The compact f sigma_8 set constrains eta_g through growth only and does not constrain beta_s without lensing data.",
        },
        "status": (
            "Effective metric-response benchmark only; no Boltzmann code, no weak-lensing likelihood, "
            "no CMB anisotropy calculation, and no evidence claim."
        ),
    }

    with open(OUT / "paper45_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\n--- Metric-response scan ---")
    print(f"Stable branches: {stable_count}/{len(scan_rows)}")
    print(f"Best growth-only eta_g: {float(best_growth[0]):.4f}")
    print(f"Growth-only beta_s is unconstrained; first listed best beta_s: {float(best_growth[1]):.4f}")
    print(f"Best growth chi2: {float(best_growth[2]):.6f}")

    print("\n--- Representative eta_g=0.02, beta_s=0.03 ---")
    rep = summary["representative"]
    print(f"Max |mu-1| z<=3: {rep['max_abs_delta_mu_pct_z_le_3']:.6f}%")
    print(f"Max |eta_slip-1| z<=3: {rep['max_abs_delta_slip_pct_z_le_3']:.6f}%")
    print(f"Max |Sigma-1| z<=3: {rep['max_abs_delta_Sigma_pct_z_le_3']:.6f}%")
    print(f"Max |Weyl proxy -1| z<=3: {rep['max_abs_delta_Weyl_proxy_pct_z_le_3']:.6f}%")
    print(f"Max |EG proxy -1| z<=3: {rep['max_abs_delta_EG_proxy_pct_z_le_3']:.6f}%")
    print(f"\nOutputs written to: {OUT.resolve()}")


if __name__ == "__main__":
    main()
