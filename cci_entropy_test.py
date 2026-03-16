"""
CCI vs Entropy Production / Structural Information Test
========================================================
Tests Maatis' hypothesis:

    CCI ≈ Ṡ_prod / (I_struct + ε)

For each configuration class (vacuum, kink, localized, chaotic)
we measure:
    1. Ṡ_prod  — coarse-grained Shannon entropy growth rate
    2. I_nn    — nearest-neighbour mutual information
    3. CCI     — Critical Coherence Index
    4. Ratio   — Ṡ_prod / (I_nn + ε)

Then we check: does CCI correlate monotonically with the ratio?

Christof Krieg / MAAT Research / March 2026
"""

from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

SEED = 42
rng = np.random.default_rng(SEED)
OUTDIR = "maat_results_cci_test"
os.makedirs(OUTDIR, exist_ok=True)

# ── Field theory ───────────────────────────────────────────────
N      = 128
DX     = 1.0
DT     = 0.02
STEPS  = 1200
SAVE_EVERY = 10
RUNS_PER_CLASS = 12
CLASSES = ["vacuum", "kink", "localized", "chaotic"]

def V(phi):  return 0.25*(phi**2-1.)**2
def dV(phi): return phi*(phi**2-1.)

def laplacian(phi):
    return (np.roll(phi,-1)+np.roll(phi,1)-2.*phi)/DX**2

def gradient(phi):
    return (np.roll(phi,-1)-np.roll(phi,1))/(2.*DX)

def rk4_step(phi, pi, dt):
    def a(p): return laplacian(p)-dV(p)
    k1p=pi; k1v=a(phi)
    k2p=pi+.5*dt*k1v; k2v=a(phi+.5*dt*k1p)
    k3p=pi+.5*dt*k2v; k3v=a(phi+.5*dt*k2p)
    k4p=pi+dt*k3v;    k4v=a(phi+dt*k3p)
    return (phi+(dt/6.)*(k1p+2*k2p+2*k3p+k4p),
            pi +(dt/6.)*(k1v+2*k2v+2*k3v+k4v))

# ── Initial conditions ─────────────────────────────────────────
def make_initial_condition(init_class):
    x = np.arange(N) - N/2
    if init_class == "vacuum":
        phi = rng.choice([-1.,1.])*np.ones(N)+.02*rng.normal(size=N)
        pi  = .02*rng.normal(size=N)
    elif init_class == "kink":
        w=rng.uniform(4.,10.); c=rng.uniform(-10,10)
        phi=np.tanh((x-c)/w)+.02*rng.normal(size=N)
        if rng.random()<.5: phi=-phi
        pi=.02*rng.normal(size=N)
    elif init_class == "localized":
        amp=rng.uniform(.8,1.8); w=rng.uniform(3.,8.)
        cx=rng.uniform(-15,15)
        phi=amp*np.exp(-(x-cx)**2/(2*w**2))+.05*rng.normal(size=N)
        pi=.05*rng.normal(size=N)
    elif init_class == "chaotic":
        phi=rng.normal(0.,1.,size=N); pi=rng.normal(0.,.7,size=N)
    return phi.astype(float), pi.astype(float)

# ── Entropy production ─────────────────────────────────────────
def coarse_grain_entropy(phi, n_bins=20):
    """Shannon entropy of coarse-grained field value distribution."""
    counts, _ = np.histogram(phi, bins=n_bins, range=(-3., 3.))
    p = counts / counts.sum()
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))

def entropy_production_rate(phi_series, dt_save):
    """
    Estimate Ṡ_prod from entropy growth over trajectory.
    Returns mean |dS/dt| across the trajectory.
    """
    entropies = [coarse_grain_entropy(s) for s in phi_series]
    dS = np.diff(entropies)
    return float(np.mean(np.abs(dS)) / dt_save)

# ── Mutual information ─────────────────────────────────────────
def mutual_information_nn(phi, bins=16):
    """Nearest-neighbour mutual information I(φ_i : φ_{i+1})."""
    x = phi[:-1]; y = phi[1:]
    h, _, _ = np.histogram2d(x, y, bins=bins)
    total = h.sum()
    if total <= 0: return 0.
    p   = h / total
    px  = p.sum(1, keepdims=True)
    py  = p.sum(0, keepdims=True)
    nz  = p > 0
    return float(np.sum(p[nz] * np.log(p[nz] / (px @ py)[nz])))

# ── CCI diagnostics ────────────────────────────────────────────
LAMBDA_STAB = 10.0
LAMBDA_BOUND = 10.0
PHI_MAX = 3.0
ALPHA_GRAD = 0.5
S_STAR = 0.42
S_WIDTH = 0.18
KAPPA = 1.0
EPS = 1e-6

def gamma_coh(phi, pi):
    """Dynamical coherence: inverse residual norm proxy."""
    res = laplacian(phi) - dV(phi) - pi  # simplified residual
    r2  = np.mean(res**2)
    return float(1. - r2 / (LAMBDA_STAB + r2))

def gamma_cons(phi, pi):
    """Conservation: energy variance proxy (low = well-conserved)."""
    e_loc = 0.5*pi**2 + 0.5*gradient(phi)**2 + V(phi)
    var   = np.var(e_loc)
    return float(1. - var / (LAMBDA_STAB + var))

def gamma_corr(phi):
    """Correlation structure via MI."""
    mi = mutual_information_nn(phi)
    return float(mi / (1. + mi))

def gamma_prod(phi, pi):
    """Non-equilibrium activity."""
    raw = np.mean(pi**2) + ALPHA_GRAD*np.mean(gradient(phi)**2)
    return float(raw / (10. + raw))

def gamma_int(phi):
    """Constraint integrity: no boundary violations."""
    bv = np.mean(np.maximum(0., np.abs(phi)-PHI_MAX)**2)
    return float(1. - bv / (LAMBDA_BOUND + bv))

def gamma_inst(phi, pi):
    """Instability proxy: kinetic energy fraction."""
    kin = np.mean(pi**2)
    return float(kin / (LAMBDA_STAB + kin))

def U_struct(phi, pi):
    """Structural imbalance = std dev of five diagnostics."""
    gs = np.array([gamma_coh(phi,pi), gamma_cons(phi,pi),
                   gamma_corr(phi),   gamma_prod(phi,pi),
                   gamma_int(phi)])
    return float(np.std(gs))

def CCI(phi, pi):
    """Critical Coherence Index."""
    num = gamma_inst(phi,pi) * gamma_prod(phi,pi) * (1. + KAPPA*U_struct(phi,pi))
    den = (gamma_coh(phi,pi) + gamma_cons(phi,pi) +
           gamma_corr(phi)   + gamma_int(phi) + EPS)
    return float(num / den)

# ── Run simulation + measure ───────────────────────────────────
def run_one(init_class):
    phi, pi = make_initial_condition(init_class)
    phi_series = []

    for step in range(STEPS):
        phi, pi = rk4_step(phi, pi, DT)
        if step % SAVE_EVERY == 0:
            phi_series.append(phi.copy())

    # Measures at final state
    cci_val  = CCI(phi, pi)
    i_nn     = mutual_information_nn(phi)
    s_dot    = entropy_production_rate(phi_series, dt_save=DT*SAVE_EVERY)
    ratio    = s_dot / (i_nn + EPS)

    # Also record individual diagnostics
    g_coh  = gamma_coh(phi, pi)
    g_cons = gamma_cons(phi, pi)
    g_corr = gamma_corr(phi)
    g_prod = gamma_prod(phi, pi)
    g_int  = gamma_int(phi)
    g_inst = gamma_inst(phi, pi)
    u_str  = U_struct(phi, pi)

    return dict(init_class=init_class,
                cci=cci_val, i_nn=i_nn,
                s_dot=s_dot, ratio=ratio,
                g_coh=g_coh, g_cons=g_cons,
                g_corr=g_corr, g_prod=g_prod,
                g_int=g_int, g_inst=g_inst,
                u_struct=u_str)

# ── Main ───────────────────────────────────────────────────────
def main():
    print("="*65)
    print("CCI vs Entropy Production / Structural Information Test")
    print("Maatis Hypothesis: CCI ≈ Ṡ_prod / (I_struct + ε)")
    print("="*65)
    print(f"Ensemble: {len(CLASSES)} classes × {RUNS_PER_CLASS} runs\n")

    rows = []
    for cls in CLASSES:
        print(f"  Running {cls}...")
        for _ in range(RUNS_PER_CLASS):
            rows.append(run_one(cls))

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUTDIR, "cci_entropy_test.csv"), index=False)

    # ── Correlation analysis ───────────────────────────────────
    cci    = df["cci"].values
    ratio  = df["ratio"].values
    s_dot  = df["s_dot"].values
    i_nn   = df["i_nn"].values

    sp_cci_ratio, p_sp = spearmanr(cci, ratio)
    pe_cci_ratio, p_pe = pearsonr(cci, ratio)
    sp_cci_sdot,  _    = spearmanr(cci, s_dot)
    sp_cci_inn,   _    = spearmanr(cci, i_nn)

    print("\n" + "="*65)
    print("CORRELATION RESULTS")
    print("="*65)
    print(f"CCI vs Ṡ/I ratio  — Spearman r = {sp_cci_ratio:.4f}  p = {p_sp:.4f}")
    print(f"CCI vs Ṡ/I ratio  — Pearson  r = {pe_cci_ratio:.4f}  p = {p_pe:.4f}")
    print(f"CCI vs Ṡ_prod     — Spearman r = {sp_cci_sdot:.4f}")
    print(f"CCI vs I_nn       — Spearman r = {sp_cci_inn:.4f}")

    # ── Class summary ──────────────────────────────────────────
    print("\n" + "="*65)
    print("CLASS SUMMARY")
    print("="*65)
    summary = df.groupby("init_class")[
        ["cci","s_dot","i_nn","ratio","g_prod","g_corr","u_struct"]
    ].mean().round(4)
    print(summary.to_string())

    # ── Interpretation ─────────────────────────────────────────
    print("\n" + "="*65)
    print("HYPOTHESIS CHECK")
    print("="*65)
    if abs(sp_cci_ratio) > 0.7:
        verdict = "STRONG — CCI correlates well with Ṡ/I ratio ✓"
    elif abs(sp_cci_ratio) > 0.4:
        verdict = "MODERATE — some correlation, needs more work"
    else:
        verdict = "WEAK — CCI and Ṡ/I ratio do not correlate well"
    print(f"Spearman |r| = {abs(sp_cci_ratio):.4f} → {verdict}")

    # ── Plots ──────────────────────────────────────────────────
    colors = {"vacuum":"gray","kink":"green",
              "localized":"blue","chaotic":"red"}

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    # Plot 1: CCI vs Ṡ/I ratio
    ax = axes[0]
    for cls in CLASSES:
        sub = df[df["init_class"]==cls]
        ax.scatter(sub["ratio"], sub["cci"],
                   label=cls, color=colors[cls], alpha=0.7, s=40)
    ax.set_xlabel("Ṡ_prod / (I_nn + ε)", fontsize=11)
    ax.set_ylabel("CCI", fontsize=11)
    ax.set_title(f"CCI vs Entropy/Info ratio\nSpearman r={sp_cci_ratio:.3f}")
    ax.legend(fontsize=9); ax.grid(alpha=0.3)

    # Plot 2: CCI vs Ṡ_prod
    ax = axes[1]
    for cls in CLASSES:
        sub = df[df["init_class"]==cls]
        ax.scatter(sub["s_dot"], sub["cci"],
                   label=cls, color=colors[cls], alpha=0.7, s=40)
    ax.set_xlabel("Ṡ_prod (entropy growth rate)", fontsize=11)
    ax.set_ylabel("CCI", fontsize=11)
    ax.set_title(f"CCI vs Entropy production\nSpearman r={sp_cci_sdot:.3f}")
    ax.legend(fontsize=9); ax.grid(alpha=0.3)

    # Plot 3: CCI vs I_nn
    ax = axes[2]
    for cls in CLASSES:
        sub = df[df["init_class"]==cls]
        ax.scatter(sub["i_nn"], sub["cci"],
                   label=cls, color=colors[cls], alpha=0.7, s=40)
    ax.set_xlabel("I_nn (structural mutual information)", fontsize=11)
    ax.set_ylabel("CCI", fontsize=11)
    ax.set_title(f"CCI vs Structural info\nSpearman r={sp_cci_inn:.3f}")
    ax.legend(fontsize=9); ax.grid(alpha=0.3)

    plt.suptitle("Maatis Hypothesis Test: CCI ≈ Ṡ_prod / (I_struct + ε)",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "cci_correlation_test.png"), dpi=200)
    plt.close()
    print(f"\nPlot saved to {OUTDIR}/cci_correlation_test.png")
    print(f"Data saved to {OUTDIR}/cci_entropy_test.csv")

if __name__ == "__main__":
    main()
