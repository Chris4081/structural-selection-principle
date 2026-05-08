#!/usr/bin/env python3
"""
MAAT Paper 47 — From Environmental Screening to Observable Void Signatures
==========================================================================

This script turns the effective environmental-screening benchmark of Paper 46
into a more observationally interpretable halo/void-environment prediction.

The main changes relative to Paper 46 are:

1. The environmental suppression variable is written in terms of
   overdensity Delta = rho / rho_bar and a bounded tidal/shear proxy

       Sigma_env = q^2 / (q^2 + q_*^2),

   where q is a dimensionless environment shear amplitude.

2. Instead of only reporting response amplitudes, the script computes the
   lensing-to-growth response ratio

       R_LG(z) = Sigma_lens_env(z) / mu_env(z)
              = [1 + eta_slip_env(z)] / 2,

   and the observable void-cluster contrast

       Delta_LG^(void-cluster) = <R_LG - 1>_void - <R_LG - 1>_cluster.

The goal is not a data fit.  It is a reproducible phenomenological benchmark
that connects MAAT environmental screening to language used in large-scale
structure, void lensing, cluster lensing, and modified-gravity consistency
tests.

Run:
    python3 maat_void_environment_prediction_v01.py
"""

from __future__ import annotations

import csv
import json
import os
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "maat_paper47_matplotlib"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


BASE_DIR = Path(__file__).resolve().parent
OUT = BASE_DIR / "outputs_paper47"
OUT.mkdir(exist_ok=True)

EPS = 1e-12
RNG_SEED = 47


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
    q_star: float = 1.0


@dataclass(frozen=True)
class ResponseParams:
    eta_g: float = 0.02
    beta_s: float = 0.03


@dataclass(frozen=True)
class EnvironmentPopulation:
    name: str
    median_delta: float
    sigma_ln_delta: float
    median_q: float
    sigma_ln_q: float
    description: str


COSMO = Cosmology()
PROJ = ProjectionParams()
SCREEN = ScreeningParams()
RESP = ResponseParams()


ENVIRONMENTS = [
    EnvironmentPopulation(
        "void",
        median_delta=0.20,
        sigma_ln_delta=0.30,
        median_q=0.15,
        sigma_ln_q=0.35,
        description="Underdense low-tidal region; maximally unscreened target.",
    ),
    EnvironmentPopulation(
        "field",
        median_delta=1.00,
        sigma_ln_delta=0.25,
        median_q=0.35,
        sigma_ln_q=0.30,
        description="Mean-density field environment.",
    ),
    EnvironmentPopulation(
        "filament",
        median_delta=5.00,
        sigma_ln_delta=0.35,
        median_q=1.00,
        sigma_ln_q=0.30,
        description="Moderately overdense anisotropic large-scale structure.",
    ),
    EnvironmentPopulation(
        "cluster",
        median_delta=100.00,
        sigma_ln_delta=0.45,
        median_q=2.00,
        sigma_ln_q=0.25,
        description="Dense high-tidal halo/cluster-like environment.",
    ),
    EnvironmentPopulation(
        "local_dense",
        median_delta=1.0e6,
        sigma_ln_delta=0.20,
        median_q=3.00,
        sigma_ln_q=0.20,
        description="Very high-density local-test proxy; near-GR recovery.",
    ),
]


def E_lcdm_z(z: np.ndarray | float, c: Cosmology = COSMO) -> np.ndarray | float:
    z = np.asarray(z)
    return np.sqrt(c.Omega_m0 * (1.0 + z) ** 3 + c.Omega_L0)


def omega_m_z(z: np.ndarray | float, c: Cosmology = COSMO) -> np.ndarray | float:
    z = np.asarray(z)
    return c.Omega_m0 * (1.0 + z) ** 3 / np.maximum(E_lcdm_z(z, c) ** 2, EPS)


def growth_factor_approx(z: np.ndarray | float, c: Cosmology = COSMO) -> np.ndarray | float:
    """Carroll-Press-Turner style normalized LCDM growth approximation."""
    z = np.asarray(z)
    om = omega_m_z(z, c)
    ol = c.Omega_L0 / np.maximum(E_lcdm_z(z, c) ** 2, EPS)
    g = (5.0 * om / 2.0) / (
        om ** (4.0 / 7.0)
        - ol
        + (1.0 + om / 2.0) * (1.0 + ol / 70.0)
    )
    g0 = (5.0 * c.Omega_m0 / 2.0) / (
        c.Omega_m0 ** (4.0 / 7.0)
        - c.Omega_L0
        + (1.0 + c.Omega_m0 / 2.0) * (1.0 + c.Omega_L0 / 70.0)
    )
    return g / (g0 * (1.0 + z))


def latent_depth(z: np.ndarray | float, p: ProjectionParams = PROJ) -> np.ndarray | float:
    z = np.asarray(z)
    return p.latent_floor + (1.0 - p.latent_floor) / (
        1.0 + np.exp(p.latent_sharpness * (z - p.latent_zc))
    )


def projection_raw(z: np.ndarray | float, p: ProjectionParams = PROJ) -> np.ndarray | float:
    z = np.asarray(z)
    b_exp = E_lcdm_z(z) * (1.0 + z)
    d_acc = growth_factor_approx(z) * latent_depth(z, p)
    s_proj = b_exp * np.tanh((b_exp / (p.Bstar_lambda + EPS)) ** p.alpha_lambda)
    return s_proj / (1.0 + p.gamma_lambda * d_acc + EPS)


def projection_kernel(z: np.ndarray, zmax: float = 3.0) -> np.ndarray:
    """Bounded Paper-43/Paper-45 projection kernel on 0 <= z <= zmax."""
    c_raw = projection_raw(z)
    c0 = float(projection_raw(0.0))
    cmax = float(projection_raw(zmax))
    return np.clip((c_raw - c0) / (cmax - c0 + EPS), 0.0, 1.0)


def sigma_env_from_q(q: np.ndarray | float, q_star: float = SCREEN.q_star) -> np.ndarray | float:
    """Bounded tidal/shear proxy used as Paper-47 interpretation of Sigma_env."""
    q = np.asarray(q)
    return q**2 / np.maximum(q**2 + q_star**2, EPS)


def screening_factor(
    Delta: np.ndarray | float,
    Sigma_env: np.ndarray | float,
    s: ScreeningParams = SCREEN,
) -> np.ndarray | float:
    Delta = np.asarray(Delta)
    Sigma_env = np.asarray(Sigma_env)
    Delta_plus = np.maximum(Delta - 1.0, 0.0)
    denom = (
        1.0
        + s.alpha_rho * Delta_plus ** s.n_rho
        + s.alpha_sigma * Sigma_env ** s.m_sigma
    )
    return 1.0 / np.maximum(denom, EPS)


def response_channels(
    z: np.ndarray,
    Delta: np.ndarray | float,
    q: np.ndarray | float,
    resp: ResponseParams = RESP,
) -> dict[str, np.ndarray]:
    Chat = projection_kernel(z)
    Sigma_env = sigma_env_from_q(q)
    Senv = screening_factor(Delta, Sigma_env)
    Cenv = np.asarray(Senv)[..., None] * Chat[None, :]
    mu = 1.0 + resp.eta_g * Cenv
    eta_slip = 1.0 + resp.beta_s * Cenv
    Sigma_lens = mu * (1.0 + eta_slip) / 2.0
    R_LG = Sigma_lens / np.maximum(mu, EPS)
    return {
        "C_hat": Chat,
        "Sigma_env": np.asarray(Sigma_env),
        "S_env": np.asarray(Senv),
        "C_env": Cenv,
        "mu": mu,
        "eta_slip": eta_slip,
        "Sigma_lens": Sigma_lens,
        "R_LG": R_LG,
    }


def draw_environment_samples(
    env: EnvironmentPopulation,
    n: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    delta = rng.lognormal(np.log(env.median_delta), env.sigma_ln_delta, n)
    q = rng.lognormal(np.log(env.median_q), env.sigma_ln_q, n)
    return delta, q


def summarize_population(
    env: EnvironmentPopulation,
    delta: np.ndarray,
    q: np.ndarray,
    z: np.ndarray,
) -> dict[str, float | str]:
    r = response_channels(z, delta, q)
    max_idx = int(np.argmax(r["C_hat"]))
    rlg_at_max = r["R_LG"][:, max_idx] - 1.0
    sigma_at_max = r["Sigma_lens"][:, max_idx] - 1.0
    mu_at_max = r["mu"][:, max_idx] - 1.0
    slip_at_max = r["eta_slip"][:, max_idx] - 1.0

    return {
        "name": env.name,
        "median_Delta": float(np.median(delta)),
        "median_q": float(np.median(q)),
        "median_Sigma_env": float(np.median(r["Sigma_env"])),
        "mean_S_env": float(np.mean(r["S_env"])),
        "median_S_env": float(np.median(r["S_env"])),
        "p16_S_env": float(np.percentile(r["S_env"], 16)),
        "p84_S_env": float(np.percentile(r["S_env"], 84)),
        "mean_max_mu_minus1_pct": float(np.mean(mu_at_max) * 100.0),
        "mean_max_slip_minus1_pct": float(np.mean(slip_at_max) * 100.0),
        "mean_max_Sigma_lens_minus1_pct": float(np.mean(sigma_at_max) * 100.0),
        "mean_max_R_LG_minus1_pct": float(np.mean(rlg_at_max) * 100.0),
        "p16_max_R_LG_minus1_pct": float(np.percentile(rlg_at_max, 16) * 100.0),
        "p84_max_R_LG_minus1_pct": float(np.percentile(rlg_at_max, 84) * 100.0),
        "stable_fraction": float(
            np.mean(
                np.all(r["mu"] > 0.0, axis=1)
                & np.all(r["eta_slip"] > 0.0, axis=1)
                & np.all(r["Sigma_lens"] > 0.0, axis=1)
            )
        ),
        "description": env.description,
    }


def save_dict_rows(path: Path, rows: list[dict]) -> None:
    fields = list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def screening_sensitivity_rows(
    samples: dict[str, tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]],
    z: np.ndarray,
) -> list[dict]:
    """Check whether the void-cluster signature survives ansatz variations."""
    cmax = float(np.max(projection_kernel(z)))
    cases = [
        ("baseline", SCREEN),
        ("weaker_density", ScreeningParams(alpha_rho=0.08, alpha_sigma=1.0, n_rho=0.75, m_sigma=2.0, q_star=1.0)),
        ("stronger_density", ScreeningParams(alpha_rho=0.25, alpha_sigma=1.0, n_rho=0.75, m_sigma=2.0, q_star=1.0)),
        ("soft_density_power", ScreeningParams(alpha_rho=0.15, alpha_sigma=1.0, n_rho=0.50, m_sigma=2.0, q_star=1.0)),
        ("linear_density_power", ScreeningParams(alpha_rho=0.15, alpha_sigma=1.0, n_rho=1.00, m_sigma=2.0, q_star=1.0)),
        ("weaker_tidal", ScreeningParams(alpha_rho=0.15, alpha_sigma=0.4, n_rho=0.75, m_sigma=2.0, q_star=1.0)),
        ("stronger_tidal", ScreeningParams(alpha_rho=0.15, alpha_sigma=1.8, n_rho=0.75, m_sigma=2.0, q_star=1.0)),
        ("linear_tidal_power", ScreeningParams(alpha_rho=0.15, alpha_sigma=1.0, n_rho=0.75, m_sigma=1.0, q_star=1.0)),
        ("sharp_tidal_power", ScreeningParams(alpha_rho=0.15, alpha_sigma=1.0, n_rho=0.75, m_sigma=3.0, q_star=1.0)),
        ("density_only", ScreeningParams(alpha_rho=0.15, alpha_sigma=0.0, n_rho=0.75, m_sigma=2.0, q_star=1.0)),
        ("tidal_only", ScreeningParams(alpha_rho=0.0, alpha_sigma=1.0, n_rho=0.75, m_sigma=2.0, q_star=1.0)),
    ]

    rows: list[dict] = []
    for case_name, params in cases:
        env_means: dict[str, float] = {}
        env_rlg: dict[str, float] = {}
        for env in ENVIRONMENTS:
            delta, q, _ = samples[env.name]
            sigma_env = sigma_env_from_q(q, params.q_star)
            s_env = screening_factor(delta, sigma_env, params)
            env_means[env.name] = float(np.mean(s_env))
            env_rlg[env.name] = float(np.mean(0.5 * RESP.beta_s * cmax * s_env) * 100.0)

        ordered = (
            env_rlg["void"]
            > env_rlg["field"]
            > env_rlg["filament"]
            > env_rlg["cluster"]
            > env_rlg["local_dense"]
        )
        rows.append(
            {
                "case": case_name,
                "alpha_rho": params.alpha_rho,
                "alpha_sigma": params.alpha_sigma,
                "n_rho": params.n_rho,
                "m_sigma": params.m_sigma,
                "q_star": params.q_star,
                "void_mean_S_env": env_means["void"],
                "field_mean_S_env": env_means["field"],
                "filament_mean_S_env": env_means["filament"],
                "cluster_mean_S_env": env_means["cluster"],
                "local_dense_mean_S_env": env_means["local_dense"],
                "void_mean_max_R_LG_minus1_pct": env_rlg["void"],
                "cluster_mean_max_R_LG_minus1_pct": env_rlg["cluster"],
                "local_dense_mean_max_R_LG_minus1_pct": env_rlg["local_dense"],
                "void_cluster_R_LG_contrast_pp": env_rlg["void"] - env_rlg["cluster"],
                "void_local_R_LG_contrast_pp": env_rlg["void"] - env_rlg["local_dense"],
                "positive_void_cluster_contrast": bool(env_rlg["void"] > env_rlg["cluster"]),
                "strict_environment_ordering": bool(ordered),
            }
        )
    return rows


def main() -> None:
    print("=" * 74)
    print("MAAT Paper 47 — Observable Void Signatures from Environmental Screening")
    print("=" * 74)

    rng = np.random.default_rng(RNG_SEED)
    z = np.linspace(0.0, 3.0, 400)
    c_hat = projection_kernel(z)
    z_peak = float(z[int(np.argmax(c_hat))])

    n_per_env = 6000
    sample_rows: list[dict] = []
    summary_rows: list[dict] = []
    curve_rows: list[dict] = []
    samples: dict[str, tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]] = {}

    for env in ENVIRONMENTS:
        delta, q = draw_environment_samples(env, n_per_env, rng)
        r = response_channels(z, delta, q)
        samples[env.name] = (delta, q, r)
        summary_rows.append(summarize_population(env, delta, q, z))

        # Store one compact Monte-Carlo sample table for downstream checks.
        idx_peak = int(np.argmax(c_hat))
        for i in range(n_per_env):
            sample_rows.append(
                {
                    "environment": env.name,
                    "Delta": float(delta[i]),
                    "q_tidal": float(q[i]),
                    "Sigma_env": float(r["Sigma_env"][i]),
                    "S_env": float(r["S_env"][i]),
                    "max_R_LG_minus1_pct": float((r["R_LG"][i, idx_peak] - 1.0) * 100.0),
                    "max_Sigma_lens_minus1_pct": float((r["Sigma_lens"][i, idx_peak] - 1.0) * 100.0),
                }
            )

        for j, zj in enumerate(z):
            curve_rows.append(
                {
                    "environment": env.name,
                    "z": float(zj),
                    "C_hat_proj": float(c_hat[j]),
                    "mean_C_env": float(np.mean(r["C_env"][:, j])),
                    "mean_mu": float(np.mean(r["mu"][:, j])),
                    "mean_eta_slip": float(np.mean(r["eta_slip"][:, j])),
                    "mean_Sigma_lens": float(np.mean(r["Sigma_lens"][:, j])),
                    "mean_R_LG": float(np.mean(r["R_LG"][:, j])),
                    "p16_R_LG": float(np.percentile(r["R_LG"][:, j], 16)),
                    "p84_R_LG": float(np.percentile(r["R_LG"][:, j], 84)),
                }
            )

    save_dict_rows(OUT / "paper47_environment_summary.csv", summary_rows)
    save_dict_rows(OUT / "paper47_monte_carlo_samples.csv", sample_rows)
    save_dict_rows(OUT / "paper47_prediction_curves.csv", curve_rows)

    sensitivity_rows = screening_sensitivity_rows(samples, z)
    save_dict_rows(OUT / "paper47_screening_sensitivity.csv", sensitivity_rows)

    by_name = {row["name"]: row for row in summary_rows}
    void_cluster_contrast = (
        by_name["void"]["mean_max_R_LG_minus1_pct"]
        - by_name["cluster"]["mean_max_R_LG_minus1_pct"]
    )
    void_local_contrast = (
        by_name["void"]["mean_max_R_LG_minus1_pct"]
        - by_name["local_dense"]["mean_max_R_LG_minus1_pct"]
    )
    void_cluster_sigma_contrast = (
        by_name["void"]["mean_max_Sigma_lens_minus1_pct"]
        - by_name["cluster"]["mean_max_Sigma_lens_minus1_pct"]
    )

    # Phase-space grid: converts the Paper-46 ad-hoc axis into a tidal proxy q.
    Delta_grid = np.logspace(np.log10(0.05), 6.0, 180)
    q_grid = np.linspace(0.0, 3.5, 181)
    DD, QQ = np.meshgrid(Delta_grid, q_grid)
    Sigenv_grid = sigma_env_from_q(QQ)
    Sgrid = screening_factor(DD, Sigenv_grid)
    cmax = float(np.max(c_hat))
    max_R_LG = 0.5 * RESP.beta_s * Sgrid * cmax * 100.0
    max_sigma = ((1.0 + RESP.eta_g * Sgrid * cmax) * (1.0 + 0.5 * RESP.beta_s * Sgrid * cmax) - 1.0) * 100.0

    grid_rows = []
    for i in range(DD.shape[0]):
        for j in range(DD.shape[1]):
            grid_rows.append(
                {
                    "Delta": float(DD[i, j]),
                    "q_tidal": float(QQ[i, j]),
                    "Sigma_env": float(Sigenv_grid[i, j]),
                    "S_env": float(Sgrid[i, j]),
                    "max_R_LG_minus1_pct": float(max_R_LG[i, j]),
                    "max_Sigma_lens_minus1_pct": float(max_sigma[i, j]),
                }
            )
    save_dict_rows(OUT / "paper47_tidal_phase_space.csv", grid_rows)

    plt.rcParams.update({"font.size": 11, "figure.dpi": 150, "savefig.dpi": 170})
    palette = {
        "void": "#3fa7ff",
        "field": "#4caf50",
        "filament": "#ff9800",
        "cluster": "#e53935",
        "local_dense": "#111111",
    }

    fig, axes = plt.subplots(2, 2, figsize=(13.4, 9.2))
    ax = axes[0, 0]
    for env in ENVIRONMENTS:
        delta, q, r = samples[env.name]
        ax.scatter(
            delta[::12],
            q[::12],
            s=7,
            alpha=0.24,
            color=palette[env.name],
            label=env.name.replace("_", " "),
        )
    ax.set_xscale("log")
    ax.set_xlabel(r"$\Delta=\rho/\bar{\rho}$")
    ax.set_ylabel(r"tidal proxy $q$")
    ax.set_title("Monte-Carlo halo/void environment populations")
    ax.grid(alpha=0.18)
    ax.legend(markerscale=2.0, fontsize=9)

    ax = axes[0, 1]
    for env in ENVIRONMENTS:
        _, _, r = samples[env.name]
        ax.hist(
            r["S_env"],
            bins=55,
            density=True,
            histtype="step",
            lw=2.0,
            color=palette[env.name],
            label=env.name.replace("_", " "),
        )
    ax.set_xlabel(r"$\mathcal{S}_{env}$")
    ax.set_ylabel("density")
    ax.set_title("Screening distributions")
    ax.grid(alpha=0.22)

    ax = axes[1, 0]
    for env in ENVIRONMENTS:
        curve = [row for row in curve_rows if row["environment"] == env.name]
        zz = np.array([row["z"] for row in curve])
        mean = np.array([row["mean_R_LG"] for row in curve])
        p16 = np.array([row["p16_R_LG"] for row in curve])
        p84 = np.array([row["p84_R_LG"] for row in curve])
        ax.plot(zz, (mean - 1.0) * 100.0, color=palette[env.name], lw=2.2, label=env.name.replace("_", " "))
        ax.fill_between(zz, (p16 - 1.0) * 100.0, (p84 - 1.0) * 100.0, color=palette[env.name], alpha=0.10)
    ax.set_xlabel("z")
    ax.set_ylabel(r"$R_{LG}-1$ [%]")
    ax.set_title("Lensing-to-growth response ratio")
    ax.grid(alpha=0.22)

    ax = axes[1, 1]
    names = [row["name"].replace("_", "\n") for row in summary_rows]
    x = np.arange(len(summary_rows))
    vals_rlg = [row["mean_max_R_LG_minus1_pct"] for row in summary_rows]
    vals_sigma = [row["mean_max_Sigma_lens_minus1_pct"] for row in summary_rows]
    ax.bar(x - 0.18, vals_rlg, width=0.36, label=r"$R_{LG}-1$", color="#00acc1")
    ax.bar(x + 0.18, vals_sigma, width=0.36, label=r"$\Sigma_{lens}-1$", color="#ab47bc")
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("mean maximum response [%]")
    ax.set_title("Observable response hierarchy")
    ax.grid(axis="y", alpha=0.22)
    ax.legend()

    fig.tight_layout()
    fig.savefig(OUT / "fig1_void_lensing_growth_summary.png", bbox_inches="tight")
    plt.close(fig)

    fig2, axes2 = plt.subplots(1, 3, figsize=(16.2, 4.9), sharey=True)
    heatmaps = [
        (Sgrid, r"$\mathcal{S}_{env}$", "viridis"),
        (max_R_LG, r"max $R_{LG}-1$ [%]", "magma"),
        (max_sigma, r"max $\Sigma_{lens}-1$ [%]", "magma"),
    ]
    for ax, (grid, title, cmap) in zip(axes2, heatmaps):
        im = ax.pcolormesh(Delta_grid, q_grid, grid, shading="auto", cmap=cmap)
        ax.set_xscale("log")
        ax.set_xlabel(r"$\Delta=\rho/\bar{\rho}$")
        ax.set_title(title)
        ax.grid(alpha=0.14)
        for env in ENVIRONMENTS:
            ax.scatter(env.median_delta, env.median_q, s=42, edgecolors="white", c="none", linewidth=1.1)
        cbar = fig2.colorbar(im, ax=ax)
        cbar.ax.set_ylabel(title)
    axes2[0].set_ylabel(r"tidal proxy $q$")
    fig2.suptitle("Tidal-environment phase space")
    fig2.tight_layout()
    fig2.savefig(OUT / "fig2_tidal_phase_space.png", bbox_inches="tight")
    plt.close(fig2)

    fig3, ax = plt.subplots(figsize=(9.6, 6.0))
    v_curve = [row for row in curve_rows if row["environment"] == "void"]
    c_curve = [row for row in curve_rows if row["environment"] == "cluster"]
    l_curve = [row for row in curve_rows if row["environment"] == "local_dense"]
    zz = np.array([row["z"] for row in v_curve])
    contrast_vc = (np.array([v["mean_R_LG"] for v in v_curve]) - np.array([c["mean_R_LG"] for c in c_curve])) * 100.0
    contrast_vl = (np.array([v["mean_R_LG"] for v in v_curve]) - np.array([l["mean_R_LG"] for l in l_curve])) * 100.0
    ax.plot(zz, contrast_vc, lw=2.6, color="#1976d2", label="void - cluster")
    ax.plot(zz, contrast_vl, lw=2.6, color="#00695c", label="void - local dense")
    ax.axhline(0.0, color="0.3", lw=1.1)
    ax.set_xlabel("z")
    ax.set_ylabel(r"$\Delta R_{LG}$ [percentage points]")
    ax.set_title("Observable environment contrast")
    ax.grid(alpha=0.24)
    ax.legend()
    fig3.tight_layout()
    fig3.savefig(OUT / "fig3_void_cluster_contrast.png", bbox_inches="tight")
    plt.close(fig3)

    fig4, axes4 = plt.subplots(1, 2, figsize=(13.0, 5.1))
    case_names = [row["case"].replace("_", "\n") for row in sensitivity_rows]
    x = np.arange(len(sensitivity_rows))
    vc = [row["void_cluster_R_LG_contrast_pp"] for row in sensitivity_rows]
    vl = [row["void_local_R_LG_contrast_pp"] for row in sensitivity_rows]

    ax = axes4[0]
    ax.bar(x - 0.18, vc, width=0.36, color="#1976d2", label="void - cluster")
    ax.bar(x + 0.18, vl, width=0.36, color="#00695c", label="void - local dense")
    ax.axhline(0.0, color="0.25", lw=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(case_names, rotation=45, ha="right")
    ax.set_ylabel(r"$\Delta R_{LG}$ [percentage points]")
    ax.set_title("Sensitivity of the observable contrast")
    ax.grid(axis="y", alpha=0.22)
    ax.legend()

    ax = axes4[1]
    ax.plot(x, [row["void_mean_S_env"] for row in sensitivity_rows], "o-", lw=2, label="void")
    ax.plot(x, [row["cluster_mean_S_env"] for row in sensitivity_rows], "o-", lw=2, label="cluster")
    ax.plot(x, [row["local_dense_mean_S_env"] for row in sensitivity_rows], "o-", lw=2, label="local dense")
    ax.set_xticks(x)
    ax.set_xticklabels(case_names, rotation=45, ha="right")
    ax.set_ylabel(r"mean $\mathcal{S}_{env}$")
    ax.set_title("Screening hierarchy under ansatz variations")
    ax.grid(axis="y", alpha=0.22)
    ax.legend()
    fig4.tight_layout()
    fig4.savefig(OUT / "fig4_screening_sensitivity.png", bbox_inches="tight")
    plt.close(fig4)

    summary = {
        "model": "MAAT Paper 47 Void/Environment Observable Signature Benchmark",
        "status": (
            "Phenomenological bridge benchmark; no full halo model, no data fit, "
            "no claim of microscopic screening derivation."
        ),
        "random_seed": RNG_SEED,
        "samples_per_environment": n_per_env,
        "cosmology": asdict(COSMO),
        "projection_params": asdict(PROJ),
        "screening_params": asdict(SCREEN),
        "response_params": asdict(RESP),
        "definitions": {
            "Sigma_env": "q^2 / (q^2 + q_star^2)",
            "S_env": "[1 + alpha_rho Delta_+^n + alpha_sigma Sigma_env^m]^(-1)",
            "C_env": "C_hat_proj(z) * S_env(Delta, Sigma_env)",
            "mu_env": "1 + eta_g * C_env",
            "eta_slip_env": "1 + beta_s * C_env",
            "Sigma_lens_env": "mu_env * (1 + eta_slip_env) / 2",
            "R_LG": "Sigma_lens_env / mu_env = (1 + eta_slip_env) / 2",
        },
        "environment_summary": summary_rows,
        "screening_sensitivity": sensitivity_rows,
        "key_results": {
            "z_peak_projection_kernel": z_peak,
            "void_mean_max_R_LG_minus1_pct": by_name["void"]["mean_max_R_LG_minus1_pct"],
            "cluster_mean_max_R_LG_minus1_pct": by_name["cluster"]["mean_max_R_LG_minus1_pct"],
            "local_dense_mean_max_R_LG_minus1_pct": by_name["local_dense"]["mean_max_R_LG_minus1_pct"],
            "void_cluster_R_LG_contrast_percentage_points": void_cluster_contrast,
            "void_local_dense_R_LG_contrast_percentage_points": void_local_contrast,
            "void_cluster_Sigma_lens_contrast_percentage_points": void_cluster_sigma_contrast,
            "stable_fraction_minimum": float(min(row["stable_fraction"] for row in summary_rows)),
            "sensitivity_cases": int(len(sensitivity_rows)),
            "sensitivity_cases_positive_void_cluster_contrast": int(
                sum(row["positive_void_cluster_contrast"] for row in sensitivity_rows)
            ),
            "sensitivity_void_cluster_contrast_min_pp": float(
                min(row["void_cluster_R_LG_contrast_pp"] for row in sensitivity_rows)
            ),
            "sensitivity_void_cluster_contrast_max_pp": float(
                max(row["void_cluster_R_LG_contrast_pp"] for row in sensitivity_rows)
            ),
        },
    }
    with open(OUT / "paper47_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\n--- Environment summary ---")
    for row in summary_rows:
        print(
            f"{row['name']:>11s}: median Delta={row['median_Delta']:.3g}, "
            f"median q={row['median_q']:.3g}, mean S_env={row['mean_S_env']:.6f}, "
            f"mean max R_LG-1={row['mean_max_R_LG_minus1_pct']:.6f}%"
        )

    print("\n--- Key contrasts ---")
    print(f"z_peak(C_hat_proj) = {z_peak:.4f}")
    print(f"void - cluster R_LG contrast = {void_cluster_contrast:.6f} percentage points")
    print(f"void - local_dense R_LG contrast = {void_local_contrast:.6f} percentage points")
    print(
        "sensitivity positive void-cluster contrast = "
        f"{sum(row['positive_void_cluster_contrast'] for row in sensitivity_rows)}/{len(sensitivity_rows)}"
    )
    print(f"Outputs written to: {OUT.resolve()}")


if __name__ == "__main__":
    main()
