#!/usr/bin/env python3
"""
MAAT Paper 46 — Environmental Screening in Structural Cosmology
================================================================

This script extends the Paper-45 metric-response benchmark by making the
bounded projection kernel environment-dependent:

    C_env(z, Delta, Sigma_env) = C_hat_proj(z) * S_env(Delta, Sigma_env)

with

    S_env = [1 + alpha_rho Delta_+^n + alpha_sigma Sigma_env^m]^(-1),
    Delta_+ = max(Delta - 1, 0).

The resulting effective growth and metric-response channels are

    mu_env(z)        = 1 + eta_g  C_env(z),
    eta_slip_env(z)  = 1 + beta_s C_env(z),
    Sigma_lens_env   = mu_env (1 + eta_slip_env) / 2.

The purpose is not to derive a microscopic chameleon/Vainshtein/symmetron
mechanism. It is an effective screening benchmark showing how MAAT projection
response can be suppressed in dense or structurally stabilized environments.

Run:
    python3 maat_environmental_screening_solver_v01.py
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
    str(Path(tempfile.gettempdir()) / "maat_paper46_matplotlib"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


BASE_DIR = Path(__file__).resolve().parent
OUT = BASE_DIR / "outputs_paper46"
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


@dataclass(frozen=True)
class ScreeningParams:
    alpha_rho: float = 0.15
    alpha_sigma: float = 1.0
    n_rho: float = 0.75
    m_sigma: float = 2.0


@dataclass(frozen=True)
class Environment:
    name: str
    Delta: float
    Sigma_env: float
    description: str


COSMO = Cosmology()
PROJ = ProjectionParams()
SCREEN = ScreeningParams()

ETA_G = 0.02
BETA_S = 0.03


ENVIRONMENTS = [
    Environment(
        "void",
        Delta=0.20,
        Sigma_env=0.05,
        description="Underdense low-stabilization region; weakest screening.",
    ),
    Environment(
        "sheet",
        Delta=0.80,
        Sigma_env=0.15,
        description="Mildly underdense sheet-like environment.",
    ),
    Environment(
        "field",
        Delta=1.00,
        Sigma_env=0.20,
        description="Mean-density field environment.",
    ),
    Environment(
        "filament",
        Delta=5.00,
        Sigma_env=0.45,
        description="Moderately overdense connected large-scale structure.",
    ),
    Environment(
        "cluster",
        Delta=100.00,
        Sigma_env=0.85,
        description="Dense high-stability cluster-like environment.",
    ),
    Environment(
        "local_dense",
        Delta=1.0e6,
        Sigma_env=0.95,
        description="Very high-density local-test proxy; near-complete screening.",
    ),
]


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


def screening_factor(
    Delta: np.ndarray | float,
    Sigma_env: np.ndarray | float,
    s: ScreeningParams = SCREEN,
) -> np.ndarray | float:
    Delta = np.asarray(Delta)
    Sigma_env = np.asarray(Sigma_env)
    Delta_plus = np.maximum(Delta - 1.0, 0.0)
    denom = 1.0 + s.alpha_rho * Delta_plus ** s.n_rho + s.alpha_sigma * Sigma_env ** s.m_sigma
    return 1.0 / np.maximum(denom, EPS)


def c_env_z(z: np.ndarray, env: Environment) -> np.ndarray:
    return projection_kernel(z) * screening_factor(env.Delta, env.Sigma_env)


def mu_env_z(z: np.ndarray, env: Environment, eta_g: float = ETA_G) -> np.ndarray:
    return 1.0 + eta_g * c_env_z(z, env)


def slip_env_z(z: np.ndarray, env: Environment, beta_s: float = BETA_S) -> np.ndarray:
    return 1.0 + beta_s * c_env_z(z, env)


def sigma_lens_env_z(
    z: np.ndarray,
    env: Environment,
    eta_g: float = ETA_G,
    beta_s: float = BETA_S,
) -> np.ndarray:
    mu = mu_env_z(z, env, eta_g)
    slip = slip_env_z(z, env, beta_s)
    return mu * (1.0 + slip) / 2.0


def mu_env_a(a: float, env: Environment, eta_g: float = ETA_G) -> float:
    z = 1.0 / a - 1.0
    return float(mu_env_z(np.array([z]), env, eta_g)[0])


def growth_rhs(N: float, y: np.ndarray, env: Environment | None, eta_g: float) -> np.ndarray:
    a = float(np.exp(N))
    delta, delta_N = y
    friction = 2.0 + dlnH_dN(a)
    mu = 1.0 if env is None else mu_env_a(a, env, eta_g)
    source = 1.5 * float(omega_m_a(a)) * mu
    return np.array([delta_N, -friction * delta_N + source * delta])


def integrate_growth(
    env: Environment | None,
    eta_g: float = ETA_G,
    a_ini: float = 0.01,
    n_steps: int = 2600,
) -> dict:
    N = np.linspace(np.log(a_ini), 0.0, n_steps)
    dN = float(N[1] - N[0])
    y = np.array([a_ini, a_ini], dtype=float)
    delta = np.empty(n_steps)
    delta_N = np.empty(n_steps)
    delta[0], delta_N[0] = y

    for i in range(1, n_steps):
        n0 = float(N[i - 1])
        k1 = growth_rhs(n0, y, env, eta_g)
        k2 = growth_rhs(n0 + 0.5 * dN, y + 0.5 * dN * k1, env, eta_g)
        k3 = growth_rhs(n0 + 0.5 * dN, y + 0.5 * dN * k2, env, eta_g)
        k4 = growth_rhs(n0 + dN, y + dN * k3, env, eta_g)
        y = y + (dN / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        delta[i], delta_N[i] = y

    a = np.exp(N)
    z = 1.0 / a - 1.0
    D = delta / np.maximum(delta[-1], EPS)
    f = delta_N / np.maximum(delta, EPS)
    fs8 = f * D * COSMO.sigma8_0
    chat = projection_kernel(z)
    if env is None:
        cenv = np.zeros_like(z)
        mu = np.ones_like(z)
        slip = np.ones_like(z)
        sigma_lens = np.ones_like(z)
        senv = 0.0
    else:
        cenv = c_env_z(z, env)
        mu = mu_env_z(z, env, eta_g)
        slip = slip_env_z(z, env, BETA_S)
        sigma_lens = sigma_lens_env_z(z, env, eta_g, BETA_S)
        senv = float(screening_factor(env.Delta, env.Sigma_env))
    return {
        "N": N,
        "a": a,
        "z": z,
        "delta": delta,
        "D": D,
        "f": f,
        "fsigma8": fs8,
        "C_hat": chat,
        "C_env": cenv,
        "mu": mu,
        "eta_slip": slip,
        "Sigma_lens": sigma_lens,
        "S_env": senv,
    }


def interp_at_z(z_query: np.ndarray, z_grid_desc: np.ndarray, values: np.ndarray) -> np.ndarray:
    return np.interp(z_query, z_grid_desc[::-1], values[::-1])


def chi2_to_growth_data(model: dict) -> float:
    z = GROWTH_DATA[:, 0]
    obs = GROWTH_DATA[:, 1]
    err = GROWTH_DATA[:, 2]
    pred = interp_at_z(z, model["z"], model["fsigma8"])
    return float(np.sum(((obs - pred) / err) ** 2))


def env_metrics(env: Environment, model: dict, lcdm: dict, z_mask: np.ndarray) -> dict:
    dD = (model["D"] - lcdm["D"]) / np.maximum(lcdm["D"], EPS)
    dfs = (model["fsigma8"] - lcdm["fsigma8"]) / np.maximum(np.abs(lcdm["fsigma8"]), EPS)
    weyl_proxy = model["Sigma_lens"] * model["D"] / np.maximum(lcdm["D"], EPS)
    eg_proxy = model["Sigma_lens"] * lcdm["f"] / np.maximum(model["f"], EPS)
    return {
        "name": env.name,
        "Delta": env.Delta,
        "Sigma_env": env.Sigma_env,
        "S_env": model["S_env"],
        "chi2_growth": chi2_to_growth_data(model),
        "max_abs_C_env": float(np.max(np.abs(model["C_env"][z_mask]))),
        "max_abs_delta_D_pct": float(np.max(np.abs(dD[z_mask])) * 100),
        "max_abs_delta_fsigma8_pct": float(np.max(np.abs(dfs[z_mask])) * 100),
        "max_abs_delta_mu_pct": float(np.max(np.abs(model["mu"][z_mask] - 1.0)) * 100),
        "max_abs_delta_slip_pct": float(np.max(np.abs(model["eta_slip"][z_mask] - 1.0)) * 100),
        "max_abs_delta_Sigma_lens_pct": float(np.max(np.abs(model["Sigma_lens"][z_mask] - 1.0)) * 100),
        "max_abs_delta_Weyl_proxy_pct": float(np.max(np.abs(weyl_proxy[z_mask] - 1.0)) * 100),
        "max_abs_delta_EG_proxy_pct": float(np.max(np.abs(eg_proxy[z_mask] - 1.0)) * 100),
        "stable": bool(
            np.all(model["mu"] > 0.0)
            and np.all(model["eta_slip"] > 0.0)
            and np.all(model["Sigma_lens"] > 0.0)
            and np.all(model["D"] > 0.0)
            and np.all(model["f"] > 0.0)
        ),
        "description": env.description,
    }


def save_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(",".join(fields) + "\n")
        for row in rows:
            values = []
            for field in fields:
                value = row[field]
                if isinstance(value, str):
                    values.append('"' + value.replace('"', '""') + '"')
                elif isinstance(value, bool):
                    values.append(str(value))
                else:
                    values.append(f"{float(value):.10g}")
            f.write(",".join(values) + "\n")


def main() -> None:
    print("=" * 72)
    print("MAAT Paper 46 — Environmental Screening in Structural Cosmology")
    print("=" * 72)

    lcdm = integrate_growth(None)
    z_mask = lcdm["z"] <= 3.0
    z_curves = lcdm["z"][z_mask]
    order = np.argsort(z_curves)
    models = {env.name: integrate_growth(env) for env in ENVIRONMENTS}
    metrics = [env_metrics(env, models[env.name], lcdm, z_mask) for env in ENVIRONMENTS]

    fields = [
        "name",
        "Delta",
        "Sigma_env",
        "S_env",
        "chi2_growth",
        "max_abs_C_env",
        "max_abs_delta_D_pct",
        "max_abs_delta_fsigma8_pct",
        "max_abs_delta_mu_pct",
        "max_abs_delta_slip_pct",
        "max_abs_delta_Sigma_lens_pct",
        "max_abs_delta_Weyl_proxy_pct",
        "max_abs_delta_EG_proxy_pct",
        "stable",
        "description",
    ]
    save_csv(OUT / "paper46_environment_archetypes.csv", metrics, fields)

    # Wide curve table for reproducibility.
    base_cols = [z_curves, lcdm["a"][z_mask], lcdm["D"][z_mask], lcdm["f"][z_mask], lcdm["fsigma8"][z_mask], lcdm["C_hat"][z_mask]]
    headers = ["z", "a", "D_LCDM", "f_LCDM", "fsigma8_LCDM", "C_hat_proj"]
    for env in ENVIRONMENTS:
        m = models[env.name]
        base_cols.extend([
            m["C_env"][z_mask],
            m["mu"][z_mask],
            m["eta_slip"][z_mask],
            m["Sigma_lens"][z_mask],
            m["D"][z_mask],
            m["fsigma8"][z_mask],
        ])
        headers.extend([
            f"C_env_{env.name}",
            f"mu_{env.name}",
            f"eta_slip_{env.name}",
            f"Sigma_lens_{env.name}",
            f"D_{env.name}",
            f"fsigma8_{env.name}",
        ])
    np.savetxt(
        OUT / "paper46_environment_curves.csv",
        np.column_stack(base_cols),
        delimiter=",",
        header=",".join(headers),
        comments="",
    )

    # Environmental phase-space grid.
    Delta_grid = np.logspace(np.log10(0.05), 6.0, 120)
    Sigma_grid = np.linspace(0.0, 1.0, 121)
    DD, SS = np.meshgrid(Delta_grid, Sigma_grid)
    Sgrid = screening_factor(DD, SS)
    cmax = float(np.max(lcdm["C_hat"][z_mask]))
    max_mu = ETA_G * Sgrid * cmax * 100
    max_slip = BETA_S * Sgrid * cmax * 100
    max_sigma = ((1.0 + ETA_G * Sgrid * cmax) * (1.0 + 0.5 * BETA_S * Sgrid * cmax) - 1.0) * 100

    grid_rows = []
    for i in range(SS.shape[0]):
        for j in range(SS.shape[1]):
            grid_rows.append({
                "Delta": float(DD[i, j]),
                "Sigma_env": float(SS[i, j]),
                "S_env": float(Sgrid[i, j]),
                "max_delta_mu_pct": float(max_mu[i, j]),
                "max_delta_slip_pct": float(max_slip[i, j]),
                "max_delta_Sigma_lens_pct": float(max_sigma[i, j]),
            })
    save_csv(
        OUT / "paper46_screening_grid.csv",
        grid_rows,
        ["Delta", "Sigma_env", "S_env", "max_delta_mu_pct", "max_delta_slip_pct", "max_delta_Sigma_lens_pct"],
    )

    # Estimate transition density for S_env ~= 0.5 at Sigma_env=0.2.
    s_line = screening_factor(Delta_grid, 0.2)
    idx_transition = int(np.argmin(np.abs(s_line - 0.5)))
    delta_half = float(Delta_grid[idx_transition])
    s_half = float(s_line[idx_transition])

    # Plots.
    plt.rcParams.update({"font.size": 11, "figure.dpi": 150, "savefig.dpi": 165})

    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.0))
    selected_names = ["void", "field", "filament", "cluster", "local_dense"]
    colors = {
        "void": "tab:blue",
        "field": "tab:green",
        "filament": "tab:orange",
        "cluster": "tab:red",
        "local_dense": "black",
    }

    ax = axes[0, 0]
    ax.plot(z_curves[order], lcdm["C_hat"][z_mask][order], color="0.35", lw=2.4, label=r"$\widehat{C}_{proj}$")
    for name in selected_names:
        model = models[name]
        ax.plot(z_curves[order], model["C_env"][z_mask][order], lw=2, color=colors[name], label=name.replace("_", " "))
    ax.set_xlabel("z")
    ax.set_ylabel("kernel")
    ax.set_title("Environmental suppression of the projection kernel")
    ax.grid(alpha=0.25)
    ax.legend()

    ax = axes[0, 1]
    for name in selected_names:
        model = models[name]
        ax.plot(z_curves[order], (model["mu"][z_mask][order] - 1.0) * 100, lw=2, color=colors[name], label=name.replace("_", " "))
    ax.set_xlabel("z")
    ax.set_ylabel(r"$\mu_{env}-1$ [%]")
    ax.set_title("Growth response is screened in dense environments")
    ax.grid(alpha=0.25)

    ax = axes[1, 0]
    for name in selected_names:
        model = models[name]
        ax.plot(z_curves[order], (model["Sigma_lens"][z_mask][order] - 1.0) * 100, lw=2, color=colors[name], label=name.replace("_", " "))
    ax.set_xlabel("z")
    ax.set_ylabel(r"$\Sigma_{lens,env}-1$ [%]")
    ax.set_title("Weyl/lensing response inherits the screening")
    ax.grid(alpha=0.25)

    ax = axes[1, 1]
    for name in selected_names:
        model = models[name]
        ax.plot(z_curves[order], model["fsigma8"][z_mask][order], lw=2, color=colors[name], label=name.replace("_", " "))
    ax.plot(z_curves[order], lcdm["fsigma8"][z_mask][order], color="0.45", ls="--", lw=2.2, label="LCDM")
    ax.errorbar(GROWTH_DATA[:, 0], GROWTH_DATA[:, 1], yerr=GROWTH_DATA[:, 2], fmt="o", capsize=3, color="0.15", label="growth data")
    ax.set_xlabel("z")
    ax.set_ylabel(r"$f\sigma_8$")
    ax.set_title("Screened growth histories")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9)

    fig.tight_layout()
    fig.savefig(OUT / "fig1_environmental_screening_summary.png", bbox_inches="tight")
    plt.close(fig)

    fig2, axes2 = plt.subplots(1, 3, figsize=(16.5, 4.8), sharey=True)
    heatmaps = [
        (Sgrid, r"$\mathcal{S}_{env}$", "viridis"),
        (max_mu, r"max $|\mu_{env}-1|$ [%]", "magma"),
        (max_sigma, r"max $|\Sigma_{lens,env}-1|$ [%]", "magma"),
    ]
    for ax, (grid, title, cmap) in zip(axes2, heatmaps):
        im = ax.pcolormesh(Delta_grid, Sigma_grid, grid, shading="auto", cmap=cmap)
        ax.set_xscale("log")
        ax.set_xlabel(r"$\Delta=\rho/\bar{\rho}$")
        ax.set_title(title)
        ax.grid(alpha=0.15)
        for env in ENVIRONMENTS:
            ax.scatter(env.Delta, env.Sigma_env, s=42, edgecolors="white", c="none", linewidth=1.0)
        cbar = fig2.colorbar(im, ax=ax)
        cbar.ax.set_ylabel(title)
    axes2[0].set_ylabel(r"$\Sigma_{env}$")
    fig2.suptitle("Environmental screening phase space")
    fig2.tight_layout()
    fig2.savefig(OUT / "fig2_screening_phase_space.png", bbox_inches="tight")
    plt.close(fig2)

    fig3, ax = plt.subplots(figsize=(11.8, 6.2))
    x = np.arange(len(metrics))
    width = 0.18
    names = [m["name"].replace("_", "\n") for m in metrics]
    bars = [
        ("max_abs_delta_mu_pct", r"$|\mu-1|$"),
        ("max_abs_delta_slip_pct", r"$|\eta_{slip}-1|$"),
        ("max_abs_delta_Sigma_lens_pct", r"$|\Sigma-1|$"),
        ("max_abs_delta_Weyl_proxy_pct", "Weyl proxy"),
    ]
    for i, (key, label) in enumerate(bars):
        ax.bar(x + (i - 1.5) * width, [m[key] for m in metrics], width=width, label=label)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("maximum deviation [%], z <= 3")
    ax.set_title("Environmental hierarchy of MAAT response")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(ncol=4)
    fig3.tight_layout()
    fig3.savefig(OUT / "fig3_environment_response_bars.png", bbox_inches="tight")
    plt.close(fig3)

    fig4, ax = plt.subplots(figsize=(9.2, 6.4))
    ax.loglog([m["Delta"] for m in metrics], [m["S_env"] for m in metrics], "o", ms=8)
    for m in metrics:
        ax.text(m["Delta"] * 1.06, m["S_env"] * 1.02, m["name"].replace("_", " "), fontsize=9)
    ax.plot(Delta_grid, screening_factor(Delta_grid, 0.2), color="black", lw=2, label=r"$\Sigma_{env}=0.2$")
    ax.axhline(0.5, color="tab:red", ls=":", lw=1.2)
    ax.axvline(delta_half, color="tab:red", ls=":", lw=1.2, label=fr"$S_{{env}}\approx0.5$ at $\Delta\approx{delta_half:.2f}$")
    ax.set_xlabel(r"$\Delta=\rho/\bar{\rho}$")
    ax.set_ylabel(r"$\mathcal{S}_{env}$")
    ax.set_title("Screening transition with overdensity")
    ax.grid(alpha=0.25, which="both")
    ax.legend()
    fig4.tight_layout()
    fig4.savefig(OUT / "fig4_screening_transition.png", bbox_inches="tight")
    plt.close(fig4)

    summary = {
        "model": "MAAT Paper 46 Environmental Screening Benchmark",
        "cosmology": asdict(COSMO),
        "projection_params": asdict(PROJ),
        "screening_params": asdict(SCREEN),
        "response_params": {
            "eta_g": ETA_G,
            "beta_s": BETA_S,
        },
        "definitions": {
            "S_env": "[1 + alpha_rho Delta_+^n + alpha_sigma Sigma_env^m]^(-1)",
            "Delta_plus": "max(Delta - 1, 0)",
            "C_env": "C_hat_proj(z) * S_env(Delta, Sigma_env)",
            "mu_env": "1 + eta_g * C_env",
            "eta_slip_env": "1 + beta_s * C_env",
            "Sigma_lens_env": "mu_env * (1 + eta_slip_env) / 2",
        },
        "environment_metrics": metrics,
        "key_results": {
            "stable_environments": int(sum(m["stable"] for m in metrics)),
            "total_environments": int(len(metrics)),
            "S_env_void": float(next(m for m in metrics if m["name"] == "void")["S_env"]),
            "S_env_cluster": float(next(m for m in metrics if m["name"] == "cluster")["S_env"]),
            "S_env_local_dense": float(next(m for m in metrics if m["name"] == "local_dense")["S_env"]),
            "void_max_delta_Sigma_lens_pct": float(next(m for m in metrics if m["name"] == "void")["max_abs_delta_Sigma_lens_pct"]),
            "cluster_max_delta_Sigma_lens_pct": float(next(m for m in metrics if m["name"] == "cluster")["max_abs_delta_Sigma_lens_pct"]),
            "local_dense_max_delta_Sigma_lens_pct": float(next(m for m in metrics if m["name"] == "local_dense")["max_abs_delta_Sigma_lens_pct"]),
            "screening_half_delta_at_sigma_env_0p2": delta_half,
            "screening_half_value": s_half,
        },
        "status": (
            "Effective environmental screening benchmark only; not a microscopic "
            "chameleon, Vainshtein, or symmetron derivation and not a local-gravity test."
        ),
    }
    with open(OUT / "paper46_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\n--- Environment metrics ---")
    for m in metrics:
        print(
            f"{m['name']:>13s}: S_env={m['S_env']:.6f}, "
            f"max|Sigma-1|={m['max_abs_delta_Sigma_lens_pct']:.6f}%, "
            f"max|Weyl-1|={m['max_abs_delta_Weyl_proxy_pct']:.6f}%"
        )

    print("\n--- Screening transition ---")
    print(f"S_env ~= 0.5 at Delta ~= {delta_half:.6f} for Sigma_env=0.2")
    print(f"Outputs written to: {OUT.resolve()}")


if __name__ == "__main__":
    main()
