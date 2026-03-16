import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

# ============================================================
# PARAMETERS
# ============================================================

SEED = 42
rng = np.random.default_rng(SEED)

# lattice + dynamics
N = 128
DX = 1.0
DT = 0.02
T_MAX = 24.0
N_STEPS = int(T_MAX / DT)

# sampling
SAMPLE_EVERY = 10  # sample every 10 time steps
EPS = 1e-8

# structural parameters
LAMBDA_STAB = 1.0
LAMBDA_CONS = 0.0   # set to zero for now (not implemented here)
LAMBDA_CORR = 1.0
LAMBDA_DYN = 1.0

# CCI parameters
KAPPA = 1.0
EPS_CCI = 1e-8

# diagnostic scales
LAMBDA_PI = 10.0
LAMBDA_BOUND = 10.0
PHI_MAX = 3.0
ALPHA_GRAD = 0.5
S0 = 1.0
I0 = 1.0

# ensemble
N_RUNS = 160

# coarse-graining / MI
N_BINS_ENTROPY = 32
N_BINS_MI = 16


# ============================================================
# BASIC MODEL
# ============================================================

def Vprime(phi: np.ndarray) -> np.ndarray:
    return phi * (phi**2 - 1.0)

def laplacian(phi: np.ndarray) -> np.ndarray:
    return (np.roll(phi, -1) + np.roll(phi, 1) - 2.0 * phi) / DX**2

def rhs(phi: np.ndarray, pi: np.ndarray):
    """Return time derivatives dphi/dt, dpi/dt."""
    dphi = pi
    dpi = laplacian(phi) - Vprime(phi)
    return dphi, dpi

def rk4_step(phi: np.ndarray, pi: np.ndarray, dt: float):
    k1_phi, k1_pi = rhs(phi, pi)

    k2_phi, k2_pi = rhs(phi + 0.5 * dt * k1_phi, pi + 0.5 * dt * k1_pi)
    k3_phi, k3_pi = rhs(phi + 0.5 * dt * k2_phi, pi + 0.5 * dt * k2_pi)
    k4_phi, k4_pi = rhs(phi + dt * k3_phi, pi + dt * k3_pi)

    phi_new = phi + (dt / 6.0) * (k1_phi + 2*k2_phi + 2*k3_phi + k4_phi)
    pi_new  = pi  + (dt / 6.0) * (k1_pi  + 2*k2_pi  + 2*k3_pi  + k4_pi)

    return phi_new, pi_new


# ============================================================
# INITIAL CONDITIONS: CONTINUOUS MIXED ENSEMBLE
# ============================================================

def make_vacuum_component(x):
    sign = rng.choice([-1.0, 1.0])
    return sign * np.ones_like(x)

def make_kink_component(x):
    x0 = rng.uniform(-0.35 * N, 0.35 * N)
    w = rng.uniform(4.0, 12.0)
    sign = rng.choice([-1.0, 1.0])
    return sign * np.tanh((x - x0) / w)

def make_localized_component(x):
    x0 = rng.uniform(-0.35 * N, 0.35 * N)
    sigma = rng.uniform(3.0, 10.0)
    amp = rng.uniform(0.4, 1.8)
    return amp * np.exp(-0.5 * ((x - x0) / sigma)**2)

def make_noise_component(x):
    return rng.normal(0.0, 1.0, size=len(x))

def sample_continuous_initial_condition():
    """
    Mixture of vacuum, kink, localized packet, and noise.
    This gives a continuum of states, not hard classes.
    """
    x = np.arange(N) - N/2

    a_vac   = rng.uniform(0.0, 1.0)
    a_kink  = rng.uniform(0.0, 1.0)
    a_loc   = rng.uniform(0.0, 1.0)
    a_noise = rng.uniform(0.0, 1.0)

    phi = (
        a_vac   * make_vacuum_component(x) +
        a_kink  * make_kink_component(x) +
        a_loc   * make_localized_component(x) +
        a_noise * make_noise_component(x)
    )

    # normalize gently to avoid absurd amplitudes
    phi = phi / (1.0 + np.std(phi))
    phi = np.clip(phi, -2.5, 2.5)

    # momentum field
    pi = rng.normal(0.0, 0.2 + 0.5 * a_noise, size=N)

    metadata = {
        "a_vac": a_vac,
        "a_kink": a_kink,
        "a_loc": a_loc,
        "a_noise": a_noise
    }
    return phi, pi, metadata


# ============================================================
# INFORMATION / ENTROPY MEASURES
# ============================================================

def histogram_entropy(values: np.ndarray, bins: int = N_BINS_ENTROPY) -> float:
    hist, _ = np.histogram(values, bins=bins, density=False)
    p = hist.astype(float) / max(np.sum(hist), 1.0)
    p = p[p > 0]
    return float(-np.sum(p * np.log(p + EPS)))

def mutual_information_nn(phi: np.ndarray, bins: int = N_BINS_MI) -> float:
    x = phi[:-1]
    y = phi[1:]

    # fixed global-ish range based on field clipping
    hist2d, xedges, yedges = np.histogram2d(
        x, y,
        bins=bins,
        range=[[-PHI_MAX, PHI_MAX], [-PHI_MAX, PHI_MAX]]
    )
    pxy = hist2d / max(np.sum(hist2d), 1.0)
    px = np.sum(pxy, axis=1)
    py = np.sum(pxy, axis=0)

    mi = 0.0
    for i in range(len(px)):
        for j in range(len(py)):
            if pxy[i, j] > 0 and px[i] > 0 and py[j] > 0:
                mi += pxy[i, j] * np.log(pxy[i, j] / (px[i] * py[j] + EPS))
    return float(mi)


# ============================================================
# STRUCTURAL DIAGNOSTICS
# ============================================================

def gamma_inst(phi: np.ndarray, pi: np.ndarray) -> float:
    mean_pi2 = np.mean(pi**2)
    return float(mean_pi2 / (LAMBDA_PI + mean_pi2))

def gamma_prod(phi: np.ndarray, pi: np.ndarray) -> float:
    grad = (np.roll(phi, -1) - np.roll(phi, 1)) / (2.0 * DX)
    activity = np.mean(pi**2) + ALPHA_GRAD * np.mean(grad**2)
    return float(activity / (S0 + activity))

def gamma_corr(phi: np.ndarray) -> float:
    mi = mutual_information_nn(phi)
    return float(mi / (I0 + mi))

def gamma_coh(phi: np.ndarray, pi: np.ndarray) -> float:
    """
    Simple coherence proxy:
    lower spatial irregularity => higher coherence.
    """
    grad = (np.roll(phi, -1) - np.roll(phi, 1)) / (2.0 * DX)
    var_grad = np.var(grad)
    return float(1.0 - var_grad / (1.0 + var_grad))

def gamma_cons(phi: np.ndarray, pi: np.ndarray) -> float:
    """
    Placeholder conservation integrity.
    For now, use a neutral value since we do not implement
    a true flux-conservation defect here.
    """
    return 1.0

def gamma_int(phi: np.ndarray) -> float:
    excess = np.maximum(0.0, np.abs(phi) - PHI_MAX)
    bad = np.mean(excess**2)
    return float(1.0 - bad / (LAMBDA_BOUND + bad))

def structural_imbalance(g_coh, g_cons, g_prod, g_corr, g_int) -> float:
    vals = np.array([g_coh, g_cons, g_prod, g_corr, g_int], dtype=float)
    mu = np.mean(vals)
    return float(np.sqrt(np.mean((vals - mu)**2)))

def compute_cci(phi: np.ndarray, pi: np.ndarray) -> float:
    g_inst = gamma_inst(phi, pi)
    g_prod = gamma_prod(phi, pi)
    g_coh  = gamma_coh(phi, pi)
    g_cons = gamma_cons(phi, pi)
    g_corr = gamma_corr(phi)
    g_int  = gamma_int(phi)

    u_struct = structural_imbalance(g_coh, g_cons, g_prod, g_corr, g_int)

    num = g_inst * g_prod * (1.0 + KAPPA * u_struct)
    den = g_coh + g_cons + g_corr + g_int + EPS_CCI
    return float(num / den)


# ============================================================
# STRUCTURAL FREE ENERGY (OPTIONAL, FOR COMPARISON)
# ============================================================

def structural_free_energy(phi: np.ndarray, pi: np.ndarray) -> float:
    g_inst = gamma_inst(phi, pi)
    g_prod = gamma_prod(phi, pi)
    g_corr = gamma_corr(phi)
    g_int  = gamma_int(phi)

    # stability sector
    E_stab = g_inst + (1.0 - g_int)

    # correlation sector (higher corr should lower energy)
    E_corr = -g_corr

    # activity sector:
    # keep simple here; later can use preferred-window version
    E_dyn = g_prod

    # conservation sector omitted for now
    E_cons = 0.0

    return (
        LAMBDA_STAB * E_stab
        + LAMBDA_CONS * E_cons
        + LAMBDA_CORR * E_corr
        + LAMBDA_DYN * E_dyn
    )


# ============================================================
# RUN ONE TRAJECTORY
# ============================================================

def evolve_and_measure(phi0: np.ndarray, pi0: np.ndarray):
    phi = phi0.copy()
    pi = pi0.copy()

    cci_series = []
    entropy_series = []
    mi_series = []
    fstruct_series = []

    for step in range(N_STEPS):
        phi, pi = rk4_step(phi, pi, DT)

        if step % SAMPLE_EVERY == 0:
            cci_series.append(compute_cci(phi, pi))
            entropy_series.append(histogram_entropy(phi, bins=N_BINS_ENTROPY))
            mi_series.append(mutual_information_nn(phi, bins=N_BINS_MI))
            fstruct_series.append(structural_free_energy(phi, pi))

    entropy_series = np.array(entropy_series, dtype=float)
    mi_series = np.array(mi_series, dtype=float)
    cci_series = np.array(cci_series, dtype=float)
    fstruct_series = np.array(fstruct_series, dtype=float)

    dt_sample = DT * SAMPLE_EVERY
    dS = np.diff(entropy_series) / dt_sample
    dS_pos = np.maximum(0.0, dS)

    mean_cci = float(np.mean(cci_series[:-1])) if len(cci_series) > 1 else float(np.mean(cci_series))
    mean_fstruct = float(np.mean(fstruct_series[:-1])) if len(fstruct_series) > 1 else float(np.mean(fstruct_series))
    mean_dS = float(np.mean(dS_pos)) if len(dS_pos) > 0 else 0.0
    mean_mi = float(np.mean(mi_series[:-1])) if len(mi_series) > 1 else float(np.mean(mi_series))
    ratio_SI = float(mean_dS / (mean_mi + EPS))

    return {
        "mean_cci": mean_cci,
        "mean_fstruct": mean_fstruct,
        "mean_dS": mean_dS,
        "mean_mi": mean_mi,
        "ratio_SI": ratio_SI,
        "final_phi_std": float(np.std(phi)),
        "final_pi_std": float(np.std(pi))
    }


# ============================================================
# MAIN ENSEMBLE
# ============================================================

def run_ensemble(n_runs=N_RUNS):
    rows = []

    for run_id in range(n_runs):
        phi0, pi0, meta = sample_continuous_initial_condition()
        res = evolve_and_measure(phi0, pi0)

        row = {"run_id": run_id}
        row.update(meta)
        row.update(res)
        rows.append(row)

        if (run_id + 1) % 20 == 0:
            print(f"Completed {run_id + 1}/{n_runs}")

    df = pd.DataFrame(rows)

    # MaxEnt weights from mean structural free energy
    beta = 4.0
    shifted = df["mean_fstruct"] - df["mean_fstruct"].min()
    w = np.exp(-beta * shifted.to_numpy())
    w /= np.sum(w)
    df["weight"] = w

    return df


# ============================================================
# ANALYSIS
# ============================================================

def print_correlations(df: pd.DataFrame):
    pairs = [
        ("mean_cci", "ratio_SI"),
        ("mean_cci", "mean_dS"),
        ("mean_cci", "mean_mi"),
        ("mean_cci", "mean_fstruct"),
    ]

    print("\nCorrelation summary")
    print("-" * 60)
    for a, b in pairs:
        sp = spearmanr(df[a], df[b])
        pe = pearsonr(df[a], df[b])
        print(f"{a} vs {b}")
        print(f"  Spearman r = {sp.statistic:.4f}, p = {sp.pvalue:.3e}")
        print(f"  Pearson  r = {pe.statistic:.4f}, p = {pe.pvalue:.3e}")
        print("-" * 60)


# ============================================================
# PLOTS
# ============================================================

def make_plots(df: pd.DataFrame, out_prefix="cci_struct_test"):
    # 1. CCI vs ratio
    plt.figure(figsize=(6, 5))
    plt.scatter(df["ratio_SI"], df["mean_cci"], alpha=0.75)
    plt.xlabel(r"$\overline{\dot S}_{cg} / (\overline{I}_{nn} + \varepsilon)$")
    plt.ylabel(r"mean CCI")
    plt.title("CCI vs structural-information ratio")
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_scatter_ratio.png", dpi=180)
    plt.close()

    # 2. CCI vs dS
    plt.figure(figsize=(6, 5))
    plt.scatter(df["mean_dS"], df["mean_cci"], alpha=0.75)
    plt.xlabel(r"$\overline{\dot S}_{cg}$")
    plt.ylabel(r"mean CCI")
    plt.title("CCI vs coarse-grained entropy production")
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_scatter_dS.png", dpi=180)
    plt.close()

    # 3. CCI vs MI
    plt.figure(figsize=(6, 5))
    plt.scatter(df["mean_mi"], df["mean_cci"], alpha=0.75)
    plt.xlabel(r"$\overline{I}_{nn}$")
    plt.ylabel(r"mean CCI")
    plt.title("CCI vs nearest-neighbour mutual information")
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_scatter_mi.png", dpi=180)
    plt.close()

    # 4. Histogram of mean CCI
    plt.figure(figsize=(6, 5))
    plt.hist(df["mean_cci"], bins=25, alpha=0.85)
    plt.xlabel("mean CCI")
    plt.ylabel("count")
    plt.title("Distribution of mean CCI")
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_hist_cci.png", dpi=180)
    plt.close()

    # 5. Heatmap-like scatter
    plt.figure(figsize=(6, 5))
    sc = plt.scatter(df["mean_mi"], df["mean_dS"], c=df["mean_cci"], alpha=0.85)
    plt.xlabel(r"$\overline{I}_{nn}$")
    plt.ylabel(r"$\overline{\dot S}_{cg}$")
    plt.title("CCI over structural-information plane")
    cbar = plt.colorbar(sc)
    cbar.set_label("mean CCI")
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_heat_scatter.png", dpi=180)
    plt.close()


# ============================================================
# SAVE
# ============================================================

def save_outputs(df: pd.DataFrame, filename="cci_struct_test_results.csv"):
    df.to_csv(filename, index=False)
    print(f"Saved results to {filename}")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    df = run_ensemble(n_runs=N_RUNS)
    save_outputs(df, "cci_struct_test_results.csv")
    print_correlations(df)
    make_plots(df, out_prefix="cci_struct_test")

    print("\nTop 10 by weight")
    print(df.sort_values("weight", ascending=False)[[
        "run_id", "weight", "mean_cci", "ratio_SI", "mean_dS", "mean_mi",
        "a_vac", "a_kink", "a_loc", "a_noise"
    ]].head(10).to_string(index=False))