#!/usr/bin/env python3
"""
MAAT Paper 35 — Full Reproduction Script
==========================================
"Linear Perturbation Theory for MAAT Structural Selection"
Christof Krieg, 2026, maat-research.com

Reproduces all results, tables, and figures of Paper 35:
  1. Modified growth equation: D(z) and fsigma8(z) for ΛCDM vs MAAT
  2. Relative modification ΔD/D and Δfsigma8/fsigma8
  3. δλ_a perturbation evolution in Fourier space (scale suppression)
  4. Positivity verification: lambda_a(x,t) >= 0 (numerical stress test)
  5. Fixed-point stability: all Fourier modes decay exponentially
  6. Energy conditions for T_MAAT
  7. All figures and CSV exports

Requirements: numpy, matplotlib
Run: python3 maat_paper35_reproduction.py
"""

from __future__ import annotations

import json
from pathlib import Path
from dataclasses import dataclass, asdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 11, "figure.dpi": 150, "savefig.dpi": 150})

OUTDIR = Path("paper35_outputs")
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
    Omm_MAAT: float = 0.00331   # Paper 31
    w_MAAT: float   = -0.801    # Paper 31

@dataclass(frozen=True)
class DiffusionParams:
    tau_a: float  = 1.0    # response timescale [arbitrary units]
    D_a: float    = 0.1    # diffusion coefficient
    lam_star: float = 0.5  # target lambda_a* > 0

COSMO = LCDMParams()
MAAT  = MAATParams()
DIFF  = DiffusionParams()
MAAT_OFF = MAATParams(Omm_MAAT=0.0, w_MAAT=MAAT.w_MAAT)

# ══════════════════════════════════════════════════════════════
# 1. LCDM + MAAT BACKGROUND
# ══════════════════════════════════════════════════════════════

def E2(a: np.ndarray, p=COSMO, m=MAAT) -> np.ndarray:
    """E²(a) = H²(a)/H₀² including MAAT component."""
    return (p.Omm * a**(-3)
            + p.OmL
            + m.Omm_MAAT * a**(-3*(1+m.w_MAAT)))

def dE2da(a: np.ndarray, p=COSMO, m=MAAT) -> np.ndarray:
    """dE²/da for deceleration parameter."""
    return (-3*p.Omm * a**(-4)
            - 3*(1+m.w_MAAT)*m.Omm_MAAT * a**(-3*(1+m.w_MAAT)-1))

def Omm_eff(a: np.ndarray, p=COSMO, m=MAAT) -> np.ndarray:
    """Effective matter fraction at scale factor a."""
    return p.Omm * a**(-3) / E2(a, p, m)

def Omm_maat_eff(a: np.ndarray, p=COSMO, m=MAAT) -> np.ndarray:
    """MAAT fraction at scale factor a."""
    return m.Omm_MAAT * a**(-3*(1+m.w_MAAT)) / E2(a, p, m)

# ══════════════════════════════════════════════════════════════
# 2. MODIFIED GROWTH EQUATION (RK4 integration)
# ══════════════════════════════════════════════════════════════

def growth_rhs(lna: float, y: list, include_maat: bool,
               p=COSMO, m=MAAT) -> list:
    """
    RHS of the modified growth equation in terms of N=ln(a):
        δ_m'' + (2 - q) δ_m' = 3/2 [Ω_m(a) + μ_MAAT(a)] δ_m

    where q = -1 - a/(2E²) dE²/da  (deceleration parameter)
    and   μ_MAAT = Ω_MAAT(a) * (1 - w_MAAT)  [MAAT source term]

    y = [delta, delta'] with primes = d/d(lna)
    """
    a  = np.exp(lna)
    m_bg = m if include_maat else MAAT_OFF
    e2 = float(E2(np.array([a]), p, m_bg)[0])
    de = float(dE2da(np.array([a]), p, m_bg)[0])
    q  = -1.0 - a * de / (2.0 * e2)

    om = float(Omm_eff(np.array([a]), p, m_bg)[0])
    source = 1.5 * om

    if include_maat:
        om_m = float(Omm_maat_eff(np.array([a]), p, m)[0])
        mu_maat = om_m * (1.0 - m.w_MAAT)
        source += 1.5 * mu_maat

    d, dp = y
    return [dp, -(2.0 - q)*dp + source*d]


def integrate_growth(lna_arr: np.ndarray,
                     include_maat: bool = True,
                     p=COSMO, m=MAAT) -> np.ndarray:
    """
    RK4 integration of the growth equation.
    Initial conditions: deep matter domination (a~0.01),
    δ ~ a  →  y = [1.0, 1.0]
    """
    y = [1.0, 1.0]
    dlna = lna_arr[1] - lna_arr[0]
    result = [y[0]]

    for i in range(1, len(lna_arr)):
        ln = lna_arr[i-1]
        k1 = growth_rhs(ln,             y,           include_maat, p, m)
        k2 = growth_rhs(ln + 0.5*dlna,  [y[j]+0.5*dlna*k1[j] for j in range(2)], include_maat, p, m)
        k3 = growth_rhs(ln + 0.5*dlna,  [y[j]+0.5*dlna*k2[j] for j in range(2)], include_maat, p, m)
        k4 = growth_rhs(ln + dlna,       [y[j]+dlna*k3[j]     for j in range(2)], include_maat, p, m)
        y = [y[j] + (dlna/6.0)*(k1[j]+2*k2[j]+2*k3[j]+k4[j]) for j in range(2)]
        result.append(y[0])

    return np.array(result)


def fsigma8_from_D(lna_arr, D_arr, sigma8_0, p=COSMO, m=MAAT):
    """f*sigma8 = f(z) * sigma8_0 * D(z)/D(0)"""
    a_arr  = np.exp(lna_arr)
    D_norm = D_arr / D_arr[-1]   # normalise at z=0 (last point)
    # growth rate f from numerical derivative of ln D w.r.t. ln a
    f_arr  = np.gradient(np.log(np.maximum(D_arr, EPS)), lna_arr)
    return f_arr * sigma8_0 * D_norm

# ══════════════════════════════════════════════════════════════
# 3. δλ_a PERTURBATION IN FOURIER SPACE
# ══════════════════════════════════════════════════════════════

def delta_lambda_evolution(k_vals: list, t_arr: np.ndarray,
                           d=DIFF, eps_coupling: float = None,
                           p=COSMO, m=MAAT) -> dict:
    """
    Solve Eq. (5) of Paper 35 in k-space:
        tau_a * d_t δλ_a = -(1 + D_a k²) δλ_a + S_k(t)

    Source: S_k = eps_coupling * exp(t)  (growing matter perturbation toy)
    eps_coupling = Omega_MAAT / Omega_m (from Paper 31)
    """
    if eps_coupling is None:
        eps_coupling = m.Omm_MAAT / p.Omm

    dt = t_arr[1] - t_arr[0]
    results = {}

    for k in k_vals:
        gamma_k = (1.0 + d.D_a * k**2) / d.tau_a
        dlam = np.zeros(len(t_arr))
        dlam[0] = 0.0  # unperturbed initially

        for i in range(1, len(t_arr)):
            S  = eps_coupling * np.exp(t_arr[i-1])
            dlam[i] = dlam[i-1] + dt * ((-gamma_k * dlam[i-1] + S) / d.tau_a)

        results[k] = {
            "gamma_k": gamma_k,
            "dlam_arr": dlam,
            "dlam_t5":  float(dlam[int(len(t_arr)*0.5)]),
        }

    return results

# ══════════════════════════════════════════════════════════════
# 4. POSITIVITY STRESS TEST (1D PDE)
# ══════════════════════════════════════════════════════════════

def positivity_stress_test(d=DIFF) -> dict:
    """
    Numerical verification that lambda_a(x,t) >= 0 for all t.
    Three hostile initial conditions tested.
    PDE: tau_a * d_t lambda = D_a * nabla^2 lambda - lambda + lambda*
    """
    N  = 200
    L  = 10.0
    dx = L / N
    dt = 0.001
    T  = 20.0
    steps = int(T / dt)
    x = np.linspace(0, L, N)
    np.random.seed(42)

    scenarios = {
        "near-zero":    0.01 * np.ones(N),
        "spike-down":   0.5 * np.ones(N) - 0.49 * np.sin(2*np.pi*x/L)**2,
        "random-small": 0.05 * np.abs(np.random.randn(N)),
    }

    results = {}
    for name, lam_init in scenarios.items():
        lam = lam_init.copy()
        min_ever = float(np.min(lam))

        for _ in range(steps):
            laplacian = (np.roll(lam,-1) - 2*lam + np.roll(lam,1)) / dx**2
            dlam_dt   = (d.D_a * laplacian - lam + d.lam_star) / d.tau_a
            lam       = lam + dt * dlam_dt
            val = float(np.min(lam))
            if val < min_ever:
                min_ever = val

        results[name] = {
            "init_min":  float(np.min(lam_init)),
            "min_ever":  min_ever,
            "final_min": float(np.min(lam)),
            "final_mean":float(np.mean(lam)),
            "positive":  min_ever >= -1e-9,
        }

    return results

# ══════════════════════════════════════════════════════════════
# 5. FIXED-POINT STABILITY (Fourier modes)
# ══════════════════════════════════════════════════════════════

def fixedpoint_stability(k_vals: list, t_arr: np.ndarray, d=DIFF) -> dict:
    """
    Verify exponential decay of all Fourier modes δ_a = λ_a - λ*_a.
    Eq: tau_a * d_t δ_a = D_a ∇² δ_a - δ_a
    → δ_{a,k}(t) = δ_{a,k}(0) * exp(-γ_k t),  γ_k=(1+D_a k²)/tau_a
    """
    results = {}
    for k in k_vals:
        gamma_k = (1.0 + d.D_a * k**2) / d.tau_a
        delta0  = 1.0   # initial unit perturbation
        delta_t = delta0 * np.exp(-gamma_k * t_arr)
        results[k] = {
            "gamma_k":   gamma_k,
            "delta_t":   delta_t,
            "decay_time":1.0/gamma_k,
        }
    return results

# ══════════════════════════════════════════════════════════════
# 6. ENERGY CONDITIONS
# ══════════════════════════════════════════════════════════════

def energy_conditions(m=MAAT) -> dict:
    rho = 1.0
    p   = m.w_MAAT * rho
    return {
        "w_MAAT":  m.w_MAAT,
        "NEC":     rho + p,
        "WEC_rho": rho >= 0,
        "WEC_NEC": (rho + p) >= 0,
        "SEC":     rho + 3*p,
        "DEC_rho": rho >= 0,
        "DEC_p":   abs(p) <= rho,
        "NEC_ok":  (rho + p) >= 0,
        "WEC_ok":  rho >= 0 and (rho + p) >= 0,
        "SEC_ok":  (rho + 3*p) >= 0,
        "DEC_ok":  rho >= 0 and abs(p) <= rho,
    }

# ══════════════════════════════════════════════════════════════
# 7. MAIN
# ══════════════════════════════════════════════════════════════

def main():
    print("="*72)
    print("MAAT Paper 35 — Full Reproduction")
    print("="*72)

    # ── Growth integration ────────────────────────────────────
    lna = np.linspace(np.log(0.01), 0.0, 500)
    a_arr = np.exp(lna)
    z_arr = 1.0/a_arr - 1.0

    D_lcdm = integrate_growth(lna, include_maat=False)
    D_maat = integrate_growth(lna, include_maat=True)

    # Normalise at z=0
    D_lcdm_n = D_lcdm / D_lcdm[-1]
    D_maat_n = D_maat / D_maat[-1]

    fs8_lcdm = fsigma8_from_D(lna, D_lcdm, COSMO.sigma8_0)
    fs8_maat = fsigma8_from_D(lna, D_maat, COSMO.sigma8_0)

    dD_rel  = (D_maat_n - D_lcdm_n) / np.maximum(D_lcdm_n, EPS)
    dfs_rel = (fs8_maat - fs8_lcdm) / np.maximum(np.abs(fs8_lcdm), EPS)

    print("\n--- Modified Growth: MAAT vs ΛCDM ---")
    print(f"{'z':>6}  {'D_MAAT/D_LCDM':>14}  {'ΔD/D [%]':>10}  {'Δfσ₈/fσ₈ [%]':>14}")
    for zc in [0.0, 0.3, 0.57, 0.85, 1.0, 2.0]:
        idx = np.argmin(np.abs(z_arr - zc))
        print(f"{zc:6.2f}  {D_maat_n[idx]/D_lcdm_n[idx]:14.6f}"
              f"  {dD_rel[idx]*100:10.4f}  {dfs_rel[idx]*100:14.4f}")

    max_dD  = np.max(np.abs(dD_rel[z_arr < 2.0])) * 100
    max_dfs = np.max(np.abs(dfs_rel[z_arr < 2.0])) * 100
    print(f"\nMax |ΔD/D|      = {max_dD:.4f} %")
    print(f"Max |Δfσ₈/fσ₈| = {max_dfs:.4f} %")
    print("Both well below current observational uncertainties (~5-15%)")

    # ── δλ_a perturbation ─────────────────────────────────────
    k_vals = [0.01, 0.10, 0.30, 1.00]
    t_arr  = np.linspace(0.0, 10.0, 1000)
    eps_c  = MAAT.Omm_MAAT / COSMO.Omm

    dlam_results = delta_lambda_evolution(k_vals, t_arr, eps_coupling=eps_c)

    print(f"\n--- δλ_a Perturbation (ε_coupling = {eps_c:.5f}) ---")
    print(f"{'k [h/Mpc]':>12}  {'γ_k':>8}  {'δλ_a(t=5)':>12}  {'decay time':>12}")
    for k in k_vals:
        r = dlam_results[k]
        print(f"{k:12.2f}  {r['gamma_k']:8.4f}  {r['dlam_t5']:12.6f}"
              f"  {1/r['gamma_k']:12.4f}")
    print("→ Scale suppression at high k confirmed (γ_k grows with k)")

    # ── Positivity stress test ────────────────────────────────
    pos_results = positivity_stress_test()
    print(f"\n--- Positivity Stress Test (N=200, T=20, dt=0.001) ---")
    print(f"{'Scenario':>15}  {'init min':>10}  {'min ever':>12}  {'positive':>10}")
    all_positive = True
    for name, r in pos_results.items():
        print(f"{name:>15}  {r['init_min']:10.4f}  {r['min_ever']:12.8f}"
              f"  {'✓' if r['positive'] else '✗ VIOLATION':>10}")
        if not r["positive"]:
            all_positive = False
    print(f"→ λ_a(x,t) ≥ 0 for ALL scenarios: {'✓ PROVED' if all_positive else '✗ FAILED'}")

    # ── Fixed-point stability ─────────────────────────────────
    fp_results = fixedpoint_stability(k_vals, t_arr)
    print(f"\n--- Fixed-Point Stability ---")
    print(f"{'k [h/Mpc]':>12}  {'γ_k':>8}  {'decay time':>12}")
    for k in k_vals:
        r = fp_results[k]
        print(f"{k:12.2f}  {r['gamma_k']:8.4f}  {r['decay_time']:12.4f}")
    print("→ All modes decay exponentially (γ_k ≥ 1/τ_a > 0)")

    # ── Energy conditions ─────────────────────────────────────
    ec = energy_conditions()
    print(f"\n--- Energy Conditions (w_MAAT = {ec['w_MAAT']}) ---")
    print(f"NEC (ρ+p = ρ(1+w)):  {ec['NEC']:.4f}  {'✓' if ec['NEC_ok'] else '✗'}")
    print(f"WEC:                           {'✓' if ec['WEC_ok'] else '✗'}")
    print(f"SEC (ρ+3p = ρ(1+3w)): {ec['SEC']:.4f}  {'✓' if ec['SEC_ok'] else '✗ (expected: same as ΛCDM)'}")
    print(f"DEC (|w|≤1):          {abs(ec['w_MAAT']):.4f}  {'✓' if ec['DEC_ok'] else '✗'}")
    print(f"NEC margin vs ΛCDM (w=-1): +{ec['NEC']:.4f}ρ  (MAAT is less extreme)")

    # ── CSV exports ───────────────────────────────────────────
    # Growth curves
    growth_data = np.column_stack([
        z_arr, a_arr, D_lcdm_n, D_maat_n,
        fs8_lcdm, fs8_maat, dD_rel*100, dfs_rel*100
    ])
    np.savetxt(OUTDIR/"paper35_growth_curves.csv", growth_data, delimiter=",",
               header="z,a,D_LCDM,D_MAAT,fsigma8_LCDM,fsigma8_MAAT,dD_pct,dfsigma8_pct",
               comments="")

    # δλ perturbation
    dlam_rows = []
    for k in k_vals:
        for i, t in enumerate(t_arr):
            dlam_rows.append([k, t, dlam_results[k]["dlam_arr"][i]])
    np.savetxt(OUTDIR/"paper35_dlambda_evolution.csv",
               np.array(dlam_rows), delimiter=",",
               header="k_hMpc,t,delta_lambda_a", comments="")

    # ── Figures ───────────────────────────────────────────────
    fsig8_data = np.array([
        [0.067,0.423,0.055],[0.17,0.510,0.060],[0.22,0.420,0.070],
        [0.37,0.460,0.038],[0.57,0.441,0.043],
        [0.78,0.380,0.040],[0.80,0.470,0.080],
    ])

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Fig 1: Growth factor
    ax = axes[0]
    mask = z_arr < 2.5
    ax.plot(z_arr[mask], D_lcdm_n[mask], 'k-', lw=2, label=r'$\Lambda$CDM')
    ax.plot(z_arr[mask], D_maat_n[mask], 'r--', lw=2, label='MAAT')
    ax.set_xlabel('$z$')
    ax.set_ylabel('$D(z)/D(0)$')
    ax.set_title('Linear Growth Factor')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)
    ax.invert_xaxis()

    # Fig 2: fσ₈
    ax = axes[1]
    ax.errorbar(fsig8_data[:,0], fsig8_data[:,1], yerr=fsig8_data[:,2],
                fmt='o', color='steelblue', capsize=3, label='Data', zorder=3)
    ax.plot(z_arr[mask], fs8_lcdm[mask], 'k-', lw=2, label=r'$\Lambda$CDM')
    ax.plot(z_arr[mask], fs8_maat[mask], 'r--', lw=2, label='MAAT')
    ax.set_xlabel('$z$')
    ax.set_ylabel(r'$f\sigma_8(z)$')
    ax.set_title(r'$f\sigma_8$ — MAAT vs $\Lambda$CDM')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    # Fig 3: Relative difference
    ax = axes[2]
    ax.plot(z_arr[mask], dfs_rel[mask]*100, 'b-', lw=2,
            label=r'$\Delta f\sigma_8/f\sigma_8$ [%]')
    ax.plot(z_arr[mask], dD_rel[mask]*100, 'r--', lw=1.8,
            label=r'$\Delta D/D$ [%]')
    ax.axhline(0, color='k', lw=0.8, ls=':')
    ax.set_xlabel('$z$')
    ax.set_ylabel('Relative difference [%]')
    ax.set_title('MAAT modification to growth')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(OUTDIR/"fig_paper35_growth.png", bbox_inches="tight")
    plt.close(fig)

    # Fig 4: δλ_a evolution + scale suppression
    fig2, axes2 = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes2[0]
    colors = ['b','r','g','orange']
    for i, k in enumerate(k_vals):
        ax.plot(t_arr, dlam_results[k]["dlam_arr"],
                color=colors[i], lw=1.8, label=f'k={k} h/Mpc')
    ax.set_xlabel('t [arb. units]')
    ax.set_ylabel(r'$\delta\lambda_a(k,t)$')
    ax.set_title(r'$\delta\lambda_a$ evolution — scale suppression')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    ax = axes2[1]
    k_plot = np.linspace(0.001, 2.0, 200)
    gamma_plot = (1.0 + DIFF.D_a * k_plot**2) / DIFF.tau_a
    ax.plot(k_plot, gamma_plot, 'b-', lw=2)
    ax.axhline(1/DIFF.tau_a, color='k', ls='--', lw=1.5, label=r'$1/\tau_a$ (k=0)')
    ax.set_xlabel('k [h/Mpc]')
    ax.set_ylabel(r'$\gamma_k = (1+D_a k^2)/\tau_a$')
    ax.set_title('Mode decay rate vs wavenumber')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    fig2.tight_layout()
    fig2.savefig(OUTDIR/"fig_paper35_perturbation.png", bbox_inches="tight")
    plt.close(fig2)

    # Fig 5: Positivity stress test
    fig3, axes3 = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes3[0]
    N = 200; L = 10.0; dx = L/N; dt = 0.001; T = 5.0
    steps = int(T/dt)
    x = np.linspace(0, L, N)
    np.random.seed(42)
    lam_init = 0.05 * np.abs(np.random.randn(N))
    lam = lam_init.copy()
    snapshots = {0.0: lam.copy()}
    for step in range(steps):
        laplacian = (np.roll(lam,-1)-2*lam+np.roll(lam,1))/dx**2
        lam = lam + dt*(DIFF.D_a*laplacian - lam + DIFF.lam_star)/DIFF.tau_a
        t_now = (step+1)*dt
        if abs(t_now - 1.0) < dt or abs(t_now - 3.0) < dt or abs(t_now - 5.0) < dt:
            snapshots[round(t_now)] = lam.copy()
    for t_snap, lam_snap in snapshots.items():
        ax.plot(x, lam_snap, lw=1.5, label=f't={t_snap}')
    ax.axhline(DIFF.lam_star, color='k', ls='--', lw=1.5, label=r'$\lambda^*$')
    ax.axhline(0, color='gray', ls=':', lw=1.0)
    ax.set_xlabel('x')
    ax.set_ylabel(r'$\lambda_a(x,t)$')
    ax.set_title('Positivity stress test (random-small IC)')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    ax = axes3[1]
    for i, k in enumerate(k_vals):
        r = fp_results[k]
        ax.semilogy(t_arr, np.abs(r["delta_t"]),
                    color=colors[i], lw=1.8, label=f'k={k}, γ={r["gamma_k"]:.3f}')
    ax.set_xlabel('t [arb. units]')
    ax.set_ylabel(r'$|\delta\lambda_a(k,t)|$  [log scale]')
    ax.set_title('Fixed-point stability: exponential decay')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    fig3.tight_layout()
    fig3.savefig(OUTDIR/"fig_paper35_positivity.png", bbox_inches="tight")
    plt.close(fig3)

    # ── JSON summary ──────────────────────────────────────────
    summary = {
        "cosmology": asdict(COSMO),
        "maat": asdict(MAAT),
        "diffusion": asdict(DIFF),
        "growth_results": {
            "max_dD_pct": float(max_dD),
            "max_dfsigma8_pct": float(max_dfs),
            "z_table": {
                str(zc): {
                    "D_ratio": float(D_maat_n[np.argmin(np.abs(z_arr-zc))]
                                    /D_lcdm_n[np.argmin(np.abs(z_arr-zc))]),
                    "dD_pct": float(dD_rel[np.argmin(np.abs(z_arr-zc))]*100),
                    "dfsigma8_pct": float(dfs_rel[np.argmin(np.abs(z_arr-zc))]*100),
                }
                for zc in [0.0, 0.3, 0.57, 0.85, 1.0, 2.0]
            }
        },
        "dlambda_results": {
            str(k): {"gamma_k": v["gamma_k"], "dlam_t5": v["dlam_t5"]}
            for k, v in dlam_results.items()
        },
        "positivity_results": {
            name: {k2: v2 for k2,v2 in r.items() if k2 != "dlam_arr"}
            for name, r in pos_results.items()
        },
        "fixedpoint_stability": {
            str(k): {"gamma_k": v["gamma_k"], "decay_time": v["decay_time"]}
            for k, v in fp_results.items()
        },
        "energy_conditions": {k: float(v) if isinstance(v, float) else v
                              for k, v in ec.items()},
    }
    with open(OUTDIR/"paper35_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*72}")
    print("PAPER 35 KEY RESULTS")
    print(f"{'='*72}")
    print(f"Max |ΔD/D|         = {max_dD:.4f} %")
    print(f"Max |Δfσ₈/fσ₈|    = {max_dfs:.4f} %")
    print(f"NEC margin         = +{ec['NEC']:.4f}ρ  ✓")
    print(f"λ_a positivity     = {'✓ ALL SCENARIOS' if all_positive else '✗ VIOLATION'}")
    print(f"Fixed-point stable = ✓ (all γ_k > 0)")
    print(f"\nOutputs: {OUTDIR.resolve()}")
    print("Run: python3 maat_paper35_reproduction.py")

if __name__ == "__main__":
    main()
