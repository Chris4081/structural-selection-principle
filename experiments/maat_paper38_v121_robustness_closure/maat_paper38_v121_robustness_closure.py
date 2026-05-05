#!/usr/bin/env python3
"""
MAAT Paper 38 — v1.2.1 Robustness Closure Reproduction Script
=============================================================

Updates the Paper 35 growth pipeline with the MAAT v1.2.1 emergent
Respect / Robustness closure.

Core v1.2.1 closure:
    R_resp = (H * B * V)^(1/3)
    R_rob  = min(R_resp, (H * B * S * V)^(1/4))
    Stability = R_rob

CCI definitions:
    CCI_min  = S / (H + B + V + eps)
    CCI_diag = S / (H + B + V + R_rob + eps)

R_resp is NOT included separately in the CCI denominator to avoid double counting.

Requirements:
    numpy, matplotlib

Run:
    python3 maat_paper38_v121_robustness_closure.py
"""

from __future__ import annotations

import json
from pathlib import Path
from dataclasses import dataclass, asdict
import os

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/tmp/codex-mpl-cache")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size": 11,
    "figure.dpi": 150,
    "savefig.dpi": 150,
})

OUTDIR = Path("paper38_v121_outputs")
OUTDIR.mkdir(exist_ok=True)

EPS = 1e-10


# ══════════════════════════════════════════════════════════════
# 0. PARAMETERS
# ══════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class LCDMParams:
    H0: float = 67.4
    Omm: float = 0.315
    OmL: float = 0.685
    sigma8_0: float = 0.811
    gamma_growth: float = 0.55


@dataclass(frozen=True)
class MAATParams:
    Omm_MAAT: float = 0.00331
    w_MAAT: float = -0.801


@dataclass(frozen=True)
class DiffusionParams:
    tau_a: float = 1.0
    D_a: float = 0.1
    lam_star: float = 0.5


COSMO = LCDMParams()
MAAT = MAATParams()
DIFF = DiffusionParams()
MAAT_OFF = MAATParams(Omm_MAAT=0.0, w_MAAT=MAAT.w_MAAT)


# ══════════════════════════════════════════════════════════════
# 1. BACKGROUND COSMOLOGY
# ══════════════════════════════════════════════════════════════

def E2(a: np.ndarray, p: LCDMParams = COSMO, m: MAATParams = MAAT) -> np.ndarray:
    """E²(a) = H²(a)/H₀² including MAAT component."""
    return (
        p.Omm * a ** (-3)
        + p.OmL
        + m.Omm_MAAT * a ** (-3 * (1 + m.w_MAAT))
    )


def dE2da(a: np.ndarray, p: LCDMParams = COSMO, m: MAATParams = MAAT) -> np.ndarray:
    """dE²/da."""
    return (
        -3 * p.Omm * a ** (-4)
        - 3 * (1 + m.w_MAAT) * m.Omm_MAAT
        * a ** (-3 * (1 + m.w_MAAT) - 1)
    )


def Omm_eff(a: np.ndarray, p: LCDMParams = COSMO, m: MAATParams = MAAT) -> np.ndarray:
    """Effective matter fraction Ω_m(a)."""
    return p.Omm * a ** (-3) / np.maximum(E2(a, p, m), EPS)


def Omm_maat_eff(a: np.ndarray, p: LCDMParams = COSMO, m: MAATParams = MAAT) -> np.ndarray:
    """Effective MAAT fraction Ω_MAAT(a)."""
    return (
        m.Omm_MAAT * a ** (-3 * (1 + m.w_MAAT))
        / np.maximum(E2(a, p, m), EPS)
    )


# ══════════════════════════════════════════════════════════════
# 2. MODIFIED GROWTH EQUATION
# ══════════════════════════════════════════════════════════════

def growth_rhs(
    lna: float,
    y: list[float],
    include_maat: bool,
    p: LCDMParams = COSMO,
    m: MAATParams = MAAT,
) -> list[float]:
    """
    Growth equation in N = ln(a):

        δ'' + (2 - q) δ' = 3/2 [Ω_m(a) + μ_MAAT(a)] δ

    q = -1 - a/(2E²) dE²/da

    μ_MAAT = Ω_MAAT(a) * (1 - w_MAAT)
    """

    a = float(np.exp(lna))
    m_bg = m if include_maat else MAAT_OFF

    e2 = float(E2(np.array([a]), p, m_bg)[0])
    de = float(dE2da(np.array([a]), p, m_bg)[0])

    q = -1.0 - a * de / (2.0 * e2)

    om = float(Omm_eff(np.array([a]), p, m_bg)[0])
    source = 1.5 * om

    if include_maat:
        om_m = float(Omm_maat_eff(np.array([a]), p, m)[0])
        mu_maat = om_m * (1.0 - m.w_MAAT)
        source += 1.5 * mu_maat

    delta, delta_prime = y
    return [
        delta_prime,
        -(2.0 - q) * delta_prime + source * delta,
    ]


def integrate_growth(
    lna_arr: np.ndarray,
    include_maat: bool = True,
    p: LCDMParams = COSMO,
    m: MAATParams = MAAT,
) -> np.ndarray:
    """RK4 integration of linear growth."""

    y = [1.0, 1.0]
    dlna = float(lna_arr[1] - lna_arr[0])
    result = [y[0]]

    for i in range(1, len(lna_arr)):
        ln = float(lna_arr[i - 1])

        k1 = growth_rhs(ln, y, include_maat, p, m)
        k2 = growth_rhs(
            ln + 0.5 * dlna,
            [y[j] + 0.5 * dlna * k1[j] for j in range(2)],
            include_maat,
            p,
            m,
        )
        k3 = growth_rhs(
            ln + 0.5 * dlna,
            [y[j] + 0.5 * dlna * k2[j] for j in range(2)],
            include_maat,
            p,
            m,
        )
        k4 = growth_rhs(
            ln + dlna,
            [y[j] + dlna * k3[j] for j in range(2)],
            include_maat,
            p,
            m,
        )

        y = [
            y[j] + (dlna / 6.0) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j])
            for j in range(2)
        ]

        result.append(y[0])

    return np.array(result)


def fsigma8_from_D(
    lna_arr: np.ndarray,
    D_arr: np.ndarray,
    sigma8_0: float,
) -> np.ndarray:
    """fσ8 = f(z) σ8,0 D(z)/D(0)."""

    D_norm = D_arr / np.maximum(D_arr[-1], EPS)
    f_arr = np.gradient(np.log(np.maximum(D_arr, EPS)), lna_arr)
    return f_arr * sigma8_0 * D_norm


# ══════════════════════════════════════════════════════════════
# 3. MAAT v1.2.1 RESPECT / ROBUSTNESS CLOSURE
# ══════════════════════════════════════════════════════════════

def support_from_defect(d: np.ndarray) -> np.ndarray:
    """Γ = 1/(1+d), with d ≥ 0."""
    return 1.0 / (1.0 + np.maximum(d, 0.0))


def r_resp(H: np.ndarray, B: np.ndarray, V: np.ndarray) -> np.ndarray:
    """R_resp = (H B V)^(1/3)."""
    return np.power(np.maximum(H * B * V, EPS), 1.0 / 3.0)


def r_rob(H: np.ndarray, B: np.ndarray, S: np.ndarray, V: np.ndarray) -> np.ndarray:
    """R_rob = min(R_resp, (H B S V)^(1/4))."""
    resp = r_resp(H, B, V)
    full_bound = np.power(np.maximum(H * B * S * V, EPS), 1.0 / 4.0)
    return np.minimum(resp, full_bound)


def cci_v121(
    H: np.ndarray,
    B: np.ndarray,
    S: np.ndarray,
    V: np.ndarray,
    robustness_weighted: bool = True,
) -> np.ndarray:
    """
    MAAT v1.2.1 CCI.

    Minimal:
        CCI = S / (H+B+V+eps)

    Diagnostic:
        CCI = S / (H+B+V+R_rob+eps)

    R_resp is not included separately.
    """

    if robustness_weighted:
        rob = r_rob(H, B, S, V)
        denom = H + B + V + rob + EPS
    else:
        denom = H + B + V + EPS

    return S / np.maximum(denom, EPS)


def build_v121_structural_supports(
    a_arr: np.ndarray,
    dD_rel: np.ndarray,
    dfs_rel: np.ndarray,
    p: LCDMParams = COSMO,
    m: MAATParams = MAAT,
) -> dict[str, np.ndarray]:
    """
    Builds diagnostic support fields H,B,S,V from existing cosmological proxies.

    These are proxy supports, not fundamental field measurements.

    H defect: growth-amplitude deviation |ΔD/D|
    B defect: growth-rate deviation |Δfσ8/fσ8|
    S defect: MAAT fraction activity proxy, normalized
    V defect: relative growth stress proxy, normalized
    """

    OmM = Omm_maat_eff(a_arr, p, m)

    dH = np.abs(dD_rel)
    dB = np.abs(dfs_rel)

    dS = OmM / np.maximum(np.max(OmM), EPS)

    abs_dfs = np.abs(dfs_rel)
    dV = abs_dfs / np.maximum(np.max(abs_dfs), EPS)

    H = support_from_defect(dH)
    B = support_from_defect(dB)
    S = support_from_defect(dS)
    V = support_from_defect(dV)

    Rresp = r_resp(H, B, V)
    Rrob = r_rob(H, B, S, V)

    CCI_min = cci_v121(H, B, S, V, robustness_weighted=False)
    CCI_diag = cci_v121(H, B, S, V, robustness_weighted=True)

    return {
        "dH": dH,
        "dB": dB,
        "dS": dS,
        "dV": dV,
        "H": H,
        "B": B,
        "S": S,
        "V": V,
        "R_resp": Rresp,
        "R_rob": Rrob,
        "CCI_min": CCI_min,
        "CCI_diag": CCI_diag,
        "Omega_MAAT": OmM,
    }


# ══════════════════════════════════════════════════════════════
# 4. δλ_a PERTURBATION IN FOURIER SPACE
# ══════════════════════════════════════════════════════════════

def delta_lambda_evolution(
    k_vals: list[float],
    t_arr: np.ndarray,
    d: DiffusionParams = DIFF,
    eps_coupling: float | None = None,
    p: LCDMParams = COSMO,
    m: MAATParams = MAAT,
) -> dict[float, dict[str, np.ndarray | float]]:
    """
    tau_a * d_t δλ_a = -(1 + D_a k²) δλ_a + S_k(t)

    Source:
        S_k = eps_coupling * exp(t)

    eps_coupling = Ω_MAAT / Ω_m
    """

    if eps_coupling is None:
        eps_coupling = m.Omm_MAAT / p.Omm

    dt = float(t_arr[1] - t_arr[0])
    results = {}

    for k in k_vals:
        gamma_k = (1.0 + d.D_a * k ** 2) / d.tau_a
        dlam = np.zeros(len(t_arr))
        dlam[0] = 0.0

        for i in range(1, len(t_arr)):
            source = eps_coupling * np.exp(t_arr[i - 1])
            dlam[i] = dlam[i - 1] + dt * ((-gamma_k * dlam[i - 1] + source) / d.tau_a)

        results[k] = {
            "gamma_k": gamma_k,
            "dlam_arr": dlam,
            "dlam_t5": float(dlam[int(len(t_arr) * 0.5)]),
        }

    return results


# ══════════════════════════════════════════════════════════════
# 5. POSITIVITY STRESS TEST
# ══════════════════════════════════════════════════════════════

def positivity_stress_test(d: DiffusionParams = DIFF) -> dict[str, dict[str, float | bool]]:
    """
    PDE:
        tau_a * d_t lambda = D_a ∇² lambda - lambda + lambda*
    """

    N = 200
    L = 10.0
    dx = L / N
    dt = 0.001
    T = 20.0
    steps = int(T / dt)
    x = np.linspace(0, L, N)

    np.random.seed(42)

    scenarios = {
        "near-zero": 0.01 * np.ones(N),
        "spike-down": 0.5 * np.ones(N) - 0.49 * np.sin(2 * np.pi * x / L) ** 2,
        "random-small": 0.05 * np.abs(np.random.randn(N)),
    }

    results = {}

    for name, lam_init in scenarios.items():
        lam = lam_init.copy()
        min_ever = float(np.min(lam))

        for _ in range(steps):
            laplacian = (np.roll(lam, -1) - 2 * lam + np.roll(lam, 1)) / dx ** 2
            dlam_dt = (d.D_a * laplacian - lam + d.lam_star) / d.tau_a
            lam = lam + dt * dlam_dt

            val = float(np.min(lam))
            if val < min_ever:
                min_ever = val

        results[name] = {
            "init_min": float(np.min(lam_init)),
            "min_ever": min_ever,
            "final_min": float(np.min(lam)),
            "final_mean": float(np.mean(lam)),
            "positive": bool(min_ever >= -1e-9),
        }

    return results


# ══════════════════════════════════════════════════════════════
# 6. FIXED-POINT STABILITY
# ══════════════════════════════════════════════════════════════

def fixedpoint_stability(
    k_vals: list[float],
    t_arr: np.ndarray,
    d: DiffusionParams = DIFF,
) -> dict[float, dict[str, np.ndarray | float]]:
    """
    tau_a * d_t δ_a = D_a ∇² δ_a - δ_a

    δ_k(t) = δ_k(0) exp(-γ_k t)
    γ_k = (1 + D_a k²)/tau_a
    """

    results = {}

    for k in k_vals:
        gamma_k = (1.0 + d.D_a * k ** 2) / d.tau_a
        delta_t = np.exp(-gamma_k * t_arr)

        results[k] = {
            "gamma_k": gamma_k,
            "delta_t": delta_t,
            "decay_time": 1.0 / gamma_k,
        }

    return results


# ══════════════════════════════════════════════════════════════
# 7. ENERGY CONDITIONS
# ══════════════════════════════════════════════════════════════

def energy_conditions(m: MAATParams = MAAT) -> dict[str, float | bool]:
    rho = 1.0
    p = m.w_MAAT * rho

    return {
        "w_MAAT": m.w_MAAT,
        "NEC": rho + p,
        "WEC_rho": rho >= 0,
        "WEC_NEC": (rho + p) >= 0,
        "SEC": rho + 3 * p,
        "DEC_rho": rho >= 0,
        "DEC_p": abs(p) <= rho,
        "NEC_ok": (rho + p) >= 0,
        "WEC_ok": rho >= 0 and (rho + p) >= 0,
        "SEC_ok": (rho + 3 * p) >= 0,
        "DEC_ok": rho >= 0 and abs(p) <= rho,
    }


# ══════════════════════════════════════════════════════════════
# 8. PLOTTING
# ══════════════════════════════════════════════════════════════

def make_growth_figures(
    z_arr: np.ndarray,
    a_arr: np.ndarray,
    D_lcdm_n: np.ndarray,
    D_maat_n: np.ndarray,
    fs8_lcdm: np.ndarray,
    fs8_maat: np.ndarray,
    dD_rel: np.ndarray,
    dfs_rel: np.ndarray,
) -> None:
    fsig8_data = np.array([
        [0.067, 0.423, 0.055],
        [0.17, 0.510, 0.060],
        [0.22, 0.420, 0.070],
        [0.37, 0.460, 0.038],
        [0.57, 0.441, 0.043],
        [0.78, 0.380, 0.040],
        [0.80, 0.470, 0.080],
    ])

    mask = z_arr < 2.5

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    ax = axes[0]
    ax.plot(z_arr[mask], D_lcdm_n[mask], "k-", lw=2, label=r"$\Lambda$CDM")
    ax.plot(z_arr[mask], D_maat_n[mask], "r--", lw=2, label="MAAT")
    ax.set_xlabel("$z$")
    ax.set_ylabel("$D(z)/D(0)$")
    ax.set_title("Linear Growth Factor")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)
    ax.invert_xaxis()

    ax = axes[1]
    ax.errorbar(
        fsig8_data[:, 0],
        fsig8_data[:, 1],
        yerr=fsig8_data[:, 2],
        fmt="o",
        color="steelblue",
        capsize=3,
        label="Data",
        zorder=3,
    )
    ax.plot(z_arr[mask], fs8_lcdm[mask], "k-", lw=2, label=r"$\Lambda$CDM")
    ax.plot(z_arr[mask], fs8_maat[mask], "r--", lw=2, label="MAAT")
    ax.set_xlabel("$z$")
    ax.set_ylabel(r"$f\sigma_8(z)$")
    ax.set_title(r"$f\sigma_8$ — MAAT vs $\Lambda$CDM")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    ax = axes[2]
    ax.plot(z_arr[mask], dfs_rel[mask] * 100, "b-", lw=2, label=r"$\Delta f\sigma_8/f\sigma_8$ [%]")
    ax.plot(z_arr[mask], dD_rel[mask] * 100, "r--", lw=1.8, label=r"$\Delta D/D$ [%]")
    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.set_xlabel("$z$")
    ax.set_ylabel("Relative difference [%]")
    ax.set_title("MAAT modification to growth")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(OUTDIR / "fig_paper38_growth.png", bbox_inches="tight")
    plt.close(fig)


def make_structural_figures(z_arr: np.ndarray, structural: dict[str, np.ndarray]) -> None:
    mask = z_arr < 2.5

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    ax = axes[0]
    for key in ["H", "B", "S", "V"]:
        ax.plot(z_arr[mask], structural[key][mask], lw=1.8, label=key)
    ax.set_xlabel("$z$")
    ax.set_ylabel("support")
    ax.set_title("MAAT v1.2.1 support fields")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    ax = axes[1]
    ax.plot(z_arr[mask], structural["R_resp"][mask], lw=2, label=r"$R_{\rm resp}$")
    ax.plot(z_arr[mask], structural["R_rob"][mask], lw=2, ls="--", label=r"$R_{\rm rob}$")
    ax.set_xlabel("$z$")
    ax.set_ylabel("closure score")
    ax.set_title("Respect / Robustness closure")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    ax = axes[2]
    ax.plot(z_arr[mask], structural["CCI_min"][mask], lw=2, label=r"$CCI_{\rm min}$")
    ax.plot(z_arr[mask], structural["CCI_diag"][mask], lw=2, ls="--", label=r"$CCI_{\rm diag}$")
    ax.set_xlabel("$z$")
    ax.set_ylabel("CCI")
    ax.set_title("CCI v1.2.1")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(OUTDIR / "fig_paper38_v121_closure.png", bbox_inches="tight")
    plt.close(fig)


def make_perturbation_figures(
    k_vals: list[float],
    t_arr: np.ndarray,
    dlam_results: dict,
    fp_results: dict,
) -> None:
    colors = ["b", "r", "g", "orange"]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    for i, k in enumerate(k_vals):
        ax.plot(
            t_arr,
            dlam_results[k]["dlam_arr"],
            color=colors[i],
            lw=1.8,
            label=f"k={k} h/Mpc",
        )
    ax.set_xlabel("t [arb. units]")
    ax.set_ylabel(r"$\delta\lambda_a(k,t)$")
    ax.set_title(r"$\delta\lambda_a$ evolution")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    ax = axes[1]
    for i, k in enumerate(k_vals):
        ax.semilogy(
            t_arr,
            np.abs(fp_results[k]["delta_t"]),
            color=colors[i],
            lw=1.8,
            label=f"k={k}, γ={fp_results[k]['gamma_k']:.3f}",
        )
    ax.set_xlabel("t [arb. units]")
    ax.set_ylabel(r"$|\delta\lambda_a(k,t)|$")
    ax.set_title("Fixed-point stability")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(OUTDIR / "fig_paper38_perturbations.png", bbox_inches="tight")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════
# 9. MAIN
# ══════════════════════════════════════════════════════════════

def main() -> None:
    print("=" * 72)
    print("MAAT Paper 38 — v1.2.1 Robustness Closure")
    print("=" * 72)

    lna = np.linspace(np.log(0.01), 0.0, 500)
    a_arr = np.exp(lna)
    z_arr = 1.0 / a_arr - 1.0

    D_lcdm = integrate_growth(lna, include_maat=False)
    D_maat = integrate_growth(lna, include_maat=True)

    D_lcdm_n = D_lcdm / np.maximum(D_lcdm[-1], EPS)
    D_maat_n = D_maat / np.maximum(D_maat[-1], EPS)

    fs8_lcdm = fsigma8_from_D(lna, D_lcdm, COSMO.sigma8_0)
    fs8_maat = fsigma8_from_D(lna, D_maat, COSMO.sigma8_0)

    dD_rel = (D_maat_n - D_lcdm_n) / np.maximum(D_lcdm_n, EPS)
    dfs_rel = (fs8_maat - fs8_lcdm) / np.maximum(np.abs(fs8_lcdm), EPS)

    mask_z2 = z_arr < 2.0
    max_dD = float(np.max(np.abs(dD_rel[mask_z2])) * 100)
    max_dfs = float(np.max(np.abs(dfs_rel[mask_z2])) * 100)

    print("\n--- Growth: MAAT vs ΛCDM ---")
    print(f"{'z':>6}  {'D_MAAT/D_LCDM':>14}  {'ΔD/D [%]':>10}  {'Δfσ8/fσ8 [%]':>14}")

    for zc in [0.0, 0.3, 0.57, 0.85, 1.0, 2.0]:
        idx = int(np.argmin(np.abs(z_arr - zc)))
        print(
            f"{zc:6.2f}  "
            f"{D_maat_n[idx] / np.maximum(D_lcdm_n[idx], EPS):14.6f}  "
            f"{dD_rel[idx] * 100:10.4f}  "
            f"{dfs_rel[idx] * 100:14.4f}"
        )

    print(f"\nMax |ΔD/D|      = {max_dD:.4f} %")
    print(f"Max |Δfσ8/fσ8| = {max_dfs:.4f} %")

    structural = build_v121_structural_supports(a_arr, dD_rel, dfs_rel)

    print("\n--- MAAT v1.2.1 Respect / Robustness Closure ---")
    print(f"<H>        = {np.mean(structural['H']):.6f}")
    print(f"<B>        = {np.mean(structural['B']):.6f}")
    print(f"<S>        = {np.mean(structural['S']):.6f}")
    print(f"<V>        = {np.mean(structural['V']):.6f}")
    print(f"<R_resp>   = {np.mean(structural['R_resp']):.6f}")
    print(f"<R_rob>    = {np.mean(structural['R_rob']):.6f}")
    print(f"<CCI_min>  = {np.mean(structural['CCI_min']):.6f}")
    print(f"<CCI_diag> = {np.mean(structural['CCI_diag']):.6f}")
    print("Stability = R_rob")

    k_vals = [0.01, 0.10, 0.30, 1.00]
    t_arr = np.linspace(0.0, 10.0, 1000)
    eps_c = MAAT.Omm_MAAT / COSMO.Omm

    dlam_results = delta_lambda_evolution(k_vals, t_arr, eps_coupling=eps_c)

    print(f"\n--- δλ_a Perturbation (ε_coupling = {eps_c:.5f}) ---")
    print(f"{'k [h/Mpc]':>12}  {'γ_k':>8}  {'δλ_a(t=5)':>12}  {'decay time':>12}")

    for k in k_vals:
        r = dlam_results[k]
        print(
            f"{k:12.2f}  "
            f"{r['gamma_k']:8.4f}  "
            f"{r['dlam_t5']:12.6f}  "
            f"{1 / r['gamma_k']:12.4f}"
        )

    pos_results = positivity_stress_test()
    all_positive = all(bool(r["positive"]) for r in pos_results.values())

    print("\n--- Positivity Stress Test ---")
    print(f"{'Scenario':>15}  {'init min':>10}  {'min ever':>12}  {'positive':>10}")

    for name, r in pos_results.items():
        print(
            f"{name:>15}  "
            f"{r['init_min']:10.4f}  "
            f"{r['min_ever']:12.8f}  "
            f"{'✓' if r['positive'] else '✗':>10}"
        )

    fp_results = fixedpoint_stability(k_vals, t_arr)

    print("\n--- Fixed-Point Stability ---")
    print(f"{'k [h/Mpc]':>12}  {'γ_k':>8}  {'decay time':>12}")

    for k in k_vals:
        r = fp_results[k]
        print(f"{k:12.2f}  {r['gamma_k']:8.4f}  {r['decay_time']:12.4f}")

    ec = energy_conditions()

    print(f"\n--- Energy Conditions (w_MAAT = {ec['w_MAAT']}) ---")
    print(f"NEC: {ec['NEC']:.4f}  {'✓' if ec['NEC_ok'] else '✗'}")
    print(f"WEC: {'✓' if ec['WEC_ok'] else '✗'}")
    print(f"SEC: {ec['SEC']:.4f}  {'✓' if ec['SEC_ok'] else '✗ expected for dark-energy-like sector'}")
    print(f"DEC: {'✓' if ec['DEC_ok'] else '✗'}")

    # CSV exports
    growth_data = np.column_stack([
        z_arr,
        a_arr,
        D_lcdm_n,
        D_maat_n,
        fs8_lcdm,
        fs8_maat,
        dD_rel * 100,
        dfs_rel * 100,
    ])

    np.savetxt(
        OUTDIR / "paper38_growth_curves.csv",
        growth_data,
        delimiter=",",
        header="z,a,D_LCDM,D_MAAT,fsigma8_LCDM,fsigma8_MAAT,dD_pct,dfsigma8_pct",
        comments="",
    )

    structural_data = np.column_stack([
        z_arr,
        a_arr,
        structural["dH"],
        structural["dB"],
        structural["dS"],
        structural["dV"],
        structural["H"],
        structural["B"],
        structural["S"],
        structural["V"],
        structural["R_resp"],
        structural["R_rob"],
        structural["CCI_min"],
        structural["CCI_diag"],
        structural["Omega_MAAT"],
    ])

    np.savetxt(
        OUTDIR / "paper38_v121_structural_closure.csv",
        structural_data,
        delimiter=",",
        header=(
            "z,a,dH,dB,dS,dV,H,B,S,V,"
            "R_resp,R_rob,CCI_min,CCI_diag,Omega_MAAT"
        ),
        comments="",
    )

    dlam_rows = []
    for k in k_vals:
        for i, t in enumerate(t_arr):
            dlam_rows.append([k, t, dlam_results[k]["dlam_arr"][i]])

    np.savetxt(
        OUTDIR / "paper38_dlambda_evolution.csv",
        np.array(dlam_rows),
        delimiter=",",
        header="k_hMpc,t,delta_lambda_a",
        comments="",
    )

    make_growth_figures(
        z_arr,
        a_arr,
        D_lcdm_n,
        D_maat_n,
        fs8_lcdm,
        fs8_maat,
        dD_rel,
        dfs_rel,
    )

    make_structural_figures(z_arr, structural)
    make_perturbation_figures(k_vals, t_arr, dlam_results, fp_results)

    summary = {
        "paper": "MAAT Paper 38 — v1.2.1 Robustness Closure",
        "cosmology": asdict(COSMO),
        "maat": asdict(MAAT),
        "diffusion": asdict(DIFF),
        "growth_results": {
            "max_dD_pct_z_lt_2": max_dD,
            "max_dfsigma8_pct_z_lt_2": max_dfs,
        },
        "maat_v121_closure": {
            "definition": {
                "R_resp": "(H*B*V)^(1/3)",
                "R_rob": "min(R_resp, (H*B*S*V)^(1/4))",
                "Stability": "R_rob",
                "CCI_min": "S/(H+B+V+eps)",
                "CCI_diag": "S/(H+B+V+R_rob+eps)",
                "note": "R_resp is not included separately in CCI denominator to avoid double counting.",
            },
            "means": {
                "H": float(np.mean(structural["H"])),
                "B": float(np.mean(structural["B"])),
                "S": float(np.mean(structural["S"])),
                "V": float(np.mean(structural["V"])),
                "R_resp": float(np.mean(structural["R_resp"])),
                "R_rob": float(np.mean(structural["R_rob"])),
                "CCI_min": float(np.mean(structural["CCI_min"])),
                "CCI_diag": float(np.mean(structural["CCI_diag"])),
            },
            "mins": {
                "R_resp": float(np.min(structural["R_resp"])),
                "R_rob": float(np.min(structural["R_rob"])),
                "CCI_min": float(np.min(structural["CCI_min"])),
                "CCI_diag": float(np.min(structural["CCI_diag"])),
            },
            "maxs": {
                "R_resp": float(np.max(structural["R_resp"])),
                "R_rob": float(np.max(structural["R_rob"])),
                "CCI_min": float(np.max(structural["CCI_min"])),
                "CCI_diag": float(np.max(structural["CCI_diag"])),
            },
        },
        "dlambda_results": {
            str(k): {
                "gamma_k": float(v["gamma_k"]),
                "dlam_t5": float(v["dlam_t5"]),
            }
            for k, v in dlam_results.items()
        },
        "positivity_results": pos_results,
        "fixedpoint_stability": {
            str(k): {
                "gamma_k": float(v["gamma_k"]),
                "decay_time": float(v["decay_time"]),
            }
            for k, v in fp_results.items()
        },
        "energy_conditions": {
            k: float(v) if isinstance(v, float) else bool(v)
            for k, v in ec.items()
        },
    }

    with open(OUTDIR / "paper38_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 72)
    print("PAPER 38 KEY RESULTS")
    print("=" * 72)
    print(f"Max |ΔD/D|          = {max_dD:.4f} %")
    print(f"Max |Δfσ8/fσ8|      = {max_dfs:.4f} %")
    print(f"Mean R_resp         = {np.mean(structural['R_resp']):.6f}")
    print(f"Mean R_rob          = {np.mean(structural['R_rob']):.6f}")
    print(f"Mean CCI_min        = {np.mean(structural['CCI_min']):.6f}")
    print(f"Mean CCI_diag       = {np.mean(structural['CCI_diag']):.6f}")
    print(f"λ_a positivity      = {'✓ ALL SCENARIOS' if all_positive else '✗ VIOLATION'}")
    print("Fixed-point stable  = ✓ all γ_k > 0")
    print(f"\nOutputs: {OUTDIR.resolve()}")
    print("Run: python3 maat_paper38_v121_robustness_closure.py")


if __name__ == "__main__":
    main()
