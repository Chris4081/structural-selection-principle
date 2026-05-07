#!/usr/bin/env python3
"""
MAAT Paper 43 — Linear Perturbations and Structure Growth
=========================================================

First perturbation-level benchmark for structural-selection cosmology.

The script solves the linear matter-growth equation

    D_NN + [2 + d ln H / dN] D_N - 3/2 Omega_m(a) mu(a) D = 0,

where N = ln(a) and

    mu(z) = G_eff/G = 1 + eta * C_hat_proj(z).

C_hat_proj is a bounded projection-coupling kernel derived from the
response-derived Paper-42 projection observable. The raw C_proj grows rapidly
with redshift, so the perturbative coupling uses a bounded normalization over
the tested interval 0 <= z <= 3.

Requirements:
    numpy, matplotlib

Run:
    python3 maat_linear_growth_solver_v01.py
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
    str(Path(tempfile.gettempdir()) / "maat_paper43_matplotlib"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


BASE_DIR = Path(__file__).resolve().parent
OUT = BASE_DIR / "outputs_paper43"
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


# Compact f sigma_8 comparison set used in Papers 40 and 42.
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
    return c.Omega_m0 * np.asarray(a) ** (-3) + c.Omega_L0


def omega_m_a(a: np.ndarray | float, c: Cosmology = COSMO) -> np.ndarray | float:
    a = np.asarray(a)
    return c.Omega_m0 * a ** (-3) / np.maximum(E2_of_a(a, c), EPS)


def dlnH_dN(a: float, c: Cosmology = COSMO) -> float:
    # Flat matter+Lambda: d ln H / d ln a = -3 Omega_m(a) / 2.
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
    """
    Bounded perturbative coupling kernel.

    Raw C_proj is a projection-stress diagnostic and grows strongly with
    redshift. For an effective perturbation coupling we use its bounded shape:

        C_hat = (C_raw - C_raw(0)) / (C_raw(zmax) - C_raw(0)).

    This keeps eta interpretable as the maximum fractional modification of
    G_eff/G over the tested interval.
    """
    c_raw = projection_raw(z)
    c0 = float(projection_raw(0.0))
    cmax = float(projection_raw(zmax))
    return np.clip((c_raw - c0) / (cmax - c0 + EPS), 0.0, 1.0)


def mu_of_a(a: float, eta: float) -> float:
    z = 1.0 / a - 1.0
    return float(1.0 + eta * projection_kernel(np.array([z]))[0])


def growth_rhs(N: float, y: np.ndarray, eta: float) -> np.ndarray:
    a = float(np.exp(N))
    delta, delta_N = y
    friction = 2.0 + dlnH_dN(a)
    source = 1.5 * float(omega_m_a(a)) * mu_of_a(a, eta)
    return np.array([delta_N, -friction * delta_N + source * delta])


def integrate_growth(eta: float, a_ini: float = 0.01, n_steps: int = 2600) -> dict:
    N = np.linspace(np.log(a_ini), 0.0, n_steps)
    dN = float(N[1] - N[0])
    y = np.array([a_ini, a_ini], dtype=float)
    delta = np.empty(n_steps)
    delta_N = np.empty(n_steps)
    delta[0], delta_N[0] = y

    for i in range(1, n_steps):
        n0 = float(N[i - 1])
        k1 = growth_rhs(n0, y, eta)
        k2 = growth_rhs(n0 + 0.5 * dN, y + 0.5 * dN * k1, eta)
        k3 = growth_rhs(n0 + 0.5 * dN, y + 0.5 * dN * k2, eta)
        k4 = growth_rhs(n0 + dN, y + dN * k3, eta)
        y = y + (dN / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        delta[i], delta_N[i] = y

    a = np.exp(N)
    z = 1.0 / a - 1.0
    D = delta / np.maximum(delta[-1], EPS)
    f = delta_N / np.maximum(delta, EPS)
    fs8 = f * D * COSMO.sigma8_0
    mu = np.array([mu_of_a(float(ai), eta) for ai in a])
    chat = projection_kernel(z)

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
    print("MAAT Paper 43 — Linear Perturbations and Structure Growth")
    print("=" * 72)

    etas = np.linspace(0.0, 0.08, 41)
    lcdm = integrate_growth(0.0)
    models = {float(eta): integrate_growth(float(eta)) for eta in etas}

    z_mask = lcdm["z"] <= 3.0
    chi2_lcdm = chi2_to_growth_data(lcdm)
    scan_rows = []

    for eta, model in models.items():
        dD = (model["D"] - lcdm["D"]) / np.maximum(lcdm["D"], EPS)
        dfs = (model["fsigma8"] - lcdm["fsigma8"]) / np.maximum(np.abs(lcdm["fsigma8"]), EPS)
        mu = model["mu"]
        stable = bool(
            np.all(mu > 0.0)
            and np.all(model["D"] > 0.0)
            and np.all(model["f"] > 0.0)
            and np.all(np.isfinite(model["fsigma8"]))
        )
        scan_rows.append([
            eta,
            chi2_to_growth_data(model),
            float(np.max(np.abs(dD[z_mask])) * 100),
            float(np.max(np.abs(dfs[z_mask])) * 100),
            float(np.max(np.abs(mu[z_mask] - 1.0)) * 100),
            stable,
        ])

    scan = np.array(scan_rows, dtype=object)
    stable_count = int(np.sum(scan[:, 5].astype(bool)))
    best_idx = int(np.argmin(scan[:, 1].astype(float)))
    best_eta = float(scan[best_idx, 0])
    best = models[best_eta]

    representative_eta = 0.02
    representative = models[representative_eta]
    dD_rep = (representative["D"] - lcdm["D"]) / np.maximum(lcdm["D"], EPS)
    dfs_rep = (
        (representative["fsigma8"] - lcdm["fsigma8"])
        / np.maximum(np.abs(lcdm["fsigma8"]), EPS)
    )

    # CSV outputs
    curve_mask = lcdm["z"] <= 3.0
    curves = np.column_stack([
        lcdm["z"][curve_mask],
        lcdm["a"][curve_mask],
        lcdm["D"][curve_mask],
        representative["D"][curve_mask],
        lcdm["f"][curve_mask],
        representative["f"][curve_mask],
        lcdm["fsigma8"][curve_mask],
        representative["fsigma8"][curve_mask],
        representative["mu"][curve_mask],
        representative["C_hat"][curve_mask],
        dD_rep[curve_mask] * 100,
        dfs_rep[curve_mask] * 100,
    ])
    np.savetxt(
        OUT / "paper43_growth_curves.csv",
        curves,
        delimiter=",",
        header=(
            "z,a,D_LCDM,D_MAAT_eta002,f_LCDM,f_MAAT_eta002,"
            "fsigma8_LCDM,fsigma8_MAAT_eta002,mu_eta002,C_hat,"
            "delta_D_pct,delta_fsigma8_pct"
        ),
        comments="",
    )

    np.savetxt(
        OUT / "paper43_eta_scan.csv",
        scan,
        delimiter=",",
        fmt=["%.8g", "%.10g", "%.10g", "%.10g", "%.10g", "%s"],
        header="eta,chi2,max_abs_delta_D_pct,max_abs_delta_fsigma8_pct,max_abs_delta_mu_pct,stable",
        comments="",
    )

    z_data = GROWTH_DATA[:, 0]
    data_table = np.column_stack([
        z_data,
        GROWTH_DATA[:, 1],
        GROWTH_DATA[:, 2],
        interp_at_z(z_data, lcdm["z"], lcdm["fsigma8"]),
        interp_at_z(z_data, representative["z"], representative["fsigma8"]),
        interp_at_z(z_data, best["z"], best["fsigma8"]),
    ])
    np.savetxt(
        OUT / "paper43_fsigma8_comparison.csv",
        data_table,
        delimiter=",",
        header="z,fsigma8_obs,sigma,fsigma8_LCDM,fsigma8_MAAT_eta002,fsigma8_best_eta",
        comments="",
    )

    # Figures
    plt.rcParams.update({"font.size": 11, "figure.dpi": 150, "savefig.dpi": 160})

    z_plot = lcdm["z"][curve_mask]
    order = np.argsort(z_plot)

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    ax = axes[0, 0]
    ax.plot(z_plot[order], lcdm["D"][curve_mask][order], label="LCDM", lw=2)
    ax.plot(z_plot[order], representative["D"][curve_mask][order], "--", label=r"MAAT $\eta=0.02$", lw=2)
    ax.set_xlabel("z")
    ax.set_ylabel("D(z)")
    ax.set_title("Linear growth factor")
    ax.grid(alpha=0.25)
    ax.legend()

    ax = axes[0, 1]
    ax.errorbar(GROWTH_DATA[:, 0], GROWTH_DATA[:, 1], yerr=GROWTH_DATA[:, 2], fmt="o", capsize=3, label="growth data")
    ax.plot(z_plot[order], lcdm["fsigma8"][curve_mask][order], label="LCDM", lw=2)
    ax.plot(z_plot[order], representative["fsigma8"][curve_mask][order], "--", label=r"MAAT $\eta=0.02$", lw=2)
    ax.set_xlabel("z")
    ax.set_ylabel(r"$f\sigma_8(z)$")
    ax.set_title(r"$f\sigma_8$ comparison")
    ax.grid(alpha=0.25)
    ax.legend()

    ax = axes[1, 0]
    ax.plot(z_plot[order], dD_rep[curve_mask][order] * 100, label=r"$\Delta D/D$", lw=2)
    ax.plot(z_plot[order], dfs_rep[curve_mask][order] * 100, label=r"$\Delta f\sigma_8/f\sigma_8$", lw=2)
    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.set_xlabel("z")
    ax.set_ylabel("relative deviation [%]")
    ax.set_title("Representative perturbative deviation")
    ax.grid(alpha=0.25)
    ax.legend()

    ax = axes[1, 1]
    ax.plot(z_plot[order], representative["C_hat"][curve_mask][order], label=r"$\widehat{C}_{proj}$", lw=2)
    ax.plot(z_plot[order], representative["mu"][curve_mask][order] - 1.0, label=r"$\mu-1$", lw=2)
    ax.set_xlabel("z")
    ax.set_ylabel("bounded coupling")
    ax.set_title("Projection-derived effective coupling")
    ax.grid(alpha=0.25)
    ax.legend()

    fig.tight_layout()
    fig.savefig(OUT / "fig1_growth_perturbation_summary.png", bbox_inches="tight")
    plt.close(fig)

    fig2, axes2 = plt.subplots(1, 2, figsize=(13, 4.8))
    ax = axes2[0]
    ax.plot(scan[:, 0].astype(float), scan[:, 1].astype(float), marker="o")
    ax.axhline(chi2_lcdm, color="k", ls="--", label="LCDM")
    ax.axvline(best_eta, color="tab:red", ls=":", label=fr"best $\eta={best_eta:.3f}$")
    ax.set_xlabel(r"$\eta$")
    ax.set_ylabel(r"$\chi^2$ against compact $f\sigma_8$")
    ax.set_title("Eta scan")
    ax.grid(alpha=0.25)
    ax.legend()

    ax = axes2[1]
    ax.plot(scan[:, 0].astype(float), scan[:, 3].astype(float), marker="o", label=r"max $|\Delta f\sigma_8/f\sigma_8|$")
    ax.plot(scan[:, 0].astype(float), scan[:, 4].astype(float), marker="s", label=r"max $|\mu-1|$")
    ax.set_xlabel(r"$\eta$")
    ax.set_ylabel("maximum deviation [%]")
    ax.set_title("Perturbative response size")
    ax.grid(alpha=0.25)
    ax.legend()
    fig2.tight_layout()
    fig2.savefig(OUT / "fig2_eta_scan.png", bbox_inches="tight")
    plt.close(fig2)

    summary = {
        "model": "MAAT Paper 43 Linear Perturbation Benchmark",
        "cosmology": asdict(COSMO),
        "projection_params": asdict(PROJ),
        "coupling_definition": "mu(z)=1+eta*C_hat_proj(z), C_hat bounded over 0<=z<=3",
        "n_growth_points": int(len(GROWTH_DATA)),
        "eta_scan": {
            "min_eta": float(etas.min()),
            "max_eta": float(etas.max()),
            "n_eta": int(len(etas)),
            "stable_count": stable_count,
            "total_count": int(len(etas)),
            "best_eta": best_eta,
            "chi2_lcdm": chi2_lcdm,
            "chi2_best": float(scan[best_idx, 1]),
        },
        "representative_eta": representative_eta,
        "representative_results": {
            "max_abs_delta_D_pct_z_le_3": float(np.max(np.abs(dD_rep[z_mask])) * 100),
            "max_abs_delta_fsigma8_pct_z_le_3": float(np.max(np.abs(dfs_rep[z_mask])) * 100),
            "max_abs_delta_mu_pct_z_le_3": float(np.max(np.abs(representative["mu"][z_mask] - 1.0)) * 100),
            "chi2": chi2_to_growth_data(representative),
        },
        "stability": {
            "mu_positive_all_scan": bool(stable_count == len(etas)),
            "D_positive_all_scan": bool(all(np.all(m["D"] > 0) for m in models.values())),
            "f_positive_all_scan": bool(all(np.all(m["f"] > 0) for m in models.values())),
        },
        "note": (
            "This is an effective linear-growth benchmark, not a Boltzmann-code "
            "likelihood or a precision cosmological fit."
        ),
    }
    with open(OUT / "paper43_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\n--- Eta scan ---")
    print(f"LCDM chi2: {chi2_lcdm:.6f}")
    print(f"Best eta: {best_eta:.6f}")
    print(f"Best chi2: {float(scan[best_idx, 1]):.6f}")
    print(f"Stable branches: {stable_count}/{len(etas)}")

    print("\n--- Representative eta=0.02 ---")
    print(f"Max |Delta D/D| z<=3: {summary['representative_results']['max_abs_delta_D_pct_z_le_3']:.6f}%")
    print(f"Max |Delta fsigma8/fsigma8| z<=3: {summary['representative_results']['max_abs_delta_fsigma8_pct_z_le_3']:.6f}%")
    print(f"Max |mu-1| z<=3: {summary['representative_results']['max_abs_delta_mu_pct_z_le_3']:.6f}%")
    print(f"chi2 eta=0.02: {summary['representative_results']['chi2']:.6f}")
    print(f"\nOutputs written to: {OUT.resolve()}")


if __name__ == "__main__":
    main()
