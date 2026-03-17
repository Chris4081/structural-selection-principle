import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr
from sklearn.linear_model import LinearRegression

# ============================================================
# PARAMETERS
# ============================================================

NX = 48
NY = 48
dx = 1.0
dt = 0.02
Tmax = 12.0
NSTEPS = int(Tmax / dt)

rng = np.random.default_rng(42)

PHI_MAX = 3.0
EPS = 1e-8
N_BINS_MI = 16
N_BINS_ENT = 32

RUNS = 60
SAMPLE_EVERY = 10

ALPHAS = np.linspace(0.5, 3.0, 26)

# ============================================================
# 2D phi^4 MODEL
# ============================================================

def dV(phi):
    return phi**3 - phi

def lap2d(phi):
    return (
        np.roll(phi, -1, axis=0)
        + np.roll(phi, 1, axis=0)
        + np.roll(phi, -1, axis=1)
        + np.roll(phi, 1, axis=1)
        - 4.0 * phi
    ) / dx**2

def grad2_sq(phi):
    gx = (np.roll(phi, -1, axis=0) - np.roll(phi, 1, axis=0)) / (2.0 * dx)
    gy = (np.roll(phi, -1, axis=1) - np.roll(phi, 1, axis=1)) / (2.0 * dx)
    return gx**2 + gy**2

def rhs(phi, pi):
    phi_dot = pi
    pi_dot = lap2d(phi) - dV(phi)
    return phi_dot, pi_dot

def rk4_step(phi, pi, dt):
    k1_phi, k1_pi = rhs(phi, pi)

    k2_phi, k2_pi = rhs(phi + 0.5 * dt * k1_phi, pi + 0.5 * dt * k1_pi)
    k3_phi, k3_pi = rhs(phi + 0.5 * dt * k2_phi, pi + 0.5 * dt * k2_pi)
    k4_phi, k4_pi = rhs(phi + dt * k3_phi, pi + dt * k3_pi)

    phi_new = phi + (dt / 6.0) * (k1_phi + 2 * k2_phi + 2 * k3_phi + k4_phi)
    pi_new  = pi  + (dt / 6.0) * (k1_pi  + 2 * k2_pi  + 2 * k3_pi  + k4_pi)

    return phi_new, pi_new

# ============================================================
# INITIAL CONDITIONS
# ============================================================

def random_initial_2d():
    x = np.arange(NX) - NX / 2
    y = np.arange(NY) - NY / 2
    X, Y = np.meshgrid(x, y, indexing="ij")

    # vacuum background
    vac = rng.choice([-1.0, 1.0]) * np.ones((NX, NY))

    # 2D domain wall
    theta = rng.uniform(0, np.pi)
    width = rng.uniform(3.0, 8.0)
    offset = rng.uniform(-0.2 * NX, 0.2 * NX)
    wall_coord = np.cos(theta) * X + np.sin(theta) * Y - offset
    wall = np.tanh(wall_coord / width)

    # localized lump
    x0 = rng.uniform(-0.2 * NX, 0.2 * NX)
    y0 = rng.uniform(-0.2 * NY, 0.2 * NY)
    sigma = rng.uniform(2.5, 6.0)
    amp = rng.uniform(0.7, 1.5)
    lump = amp * np.exp(-((X - x0) ** 2 + (Y - y0) ** 2) / (2 * sigma**2))

    # noise
    noise = rng.normal(0.0, 0.5, size=(NX, NY))

    a = rng.uniform(0, 1, 4)
    phi = a[0] * vac + a[1] * wall + a[2] * lump + a[3] * noise
    phi /= (1.0 + np.std(phi))
    phi = np.clip(phi, -2.5, 2.5)

    pi = rng.normal(0.0, 0.3, size=(NX, NY))
    return phi, pi

# ============================================================
# DIAGNOSTICS
# ============================================================

def residual(phi):
    return lap2d(phi) - dV(phi)

def coarse_grained_entropy(phi):
    hist, _ = np.histogram(phi, bins=N_BINS_ENT, range=[-PHI_MAX, PHI_MAX])
    p = hist.astype(float) / max(np.sum(hist), 1.0)
    p = p[p > 0]
    return float(-np.sum(p * np.log(p + EPS)))

def mutual_information_pair(a, b, bins=N_BINS_MI):
    hist2d, _, _ = np.histogram2d(
        a, b,
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

def mutual_information_2d(phi):
    """
    Average nearest-neighbour MI in x and y directions.
    """
    x1 = phi.ravel()
    x2 = np.roll(phi, -1, axis=0).ravel()
    y2 = np.roll(phi, -1, axis=1).ravel()

    mi_x = mutual_information_pair(x1, x2)
    mi_y = mutual_information_pair(x1, y2)
    return 0.5 * (mi_x + mi_y)

def Fstruct(phi, pi):
    res = np.mean(np.abs(residual(phi)))
    act = np.mean(pi**2)
    corr = mutual_information_2d(phi)
    return float(res + 0.1 * act - 0.2 * corr)

def CCI(phi, pi):
    prod = np.mean(np.abs(pi))
    coh = np.exp(-np.mean(np.abs(residual(phi))))
    corr = mutual_information_2d(phi)
    return float(prod / (coh + corr + 1e-8))

# ============================================================
# SINGLE TRAJECTORY
# ============================================================

def run_traj_2d():
    phi, pi = random_initial_2d()

    cci_series = []
    f_series = []
    mi_series = []
    s_series = []

    phi_snapshots = []
    times = []

    for step in range(NSTEPS + 1):
        if step % SAMPLE_EVERY == 0:
            cci_series.append(CCI(phi, pi))
            f_series.append(Fstruct(phi, pi))
            mi_series.append(mutual_information_2d(phi))
            s_series.append(coarse_grained_entropy(phi))
            phi_snapshots.append(phi.copy())
            times.append(step * dt)

        if step < NSTEPS:
            phi, pi = rk4_step(phi, pi, dt)

    cci_series = np.array(cci_series)
    f_series = np.array(f_series)
    mi_series = np.array(mi_series)
    s_series = np.array(s_series)
    times = np.array(times)

    dt_sample = SAMPLE_EVERY * dt
    dS = np.diff(s_series) / dt_sample
    dS_pos = np.maximum(0.0, dS)

    mean_cci = float(np.mean(cci_series[:-1]))
    mean_fstruct = float(np.mean(f_series[:-1]))
    mean_mi = float(np.mean(mi_series[:-1]))
    mean_dS_pos = float(np.mean(dS_pos)) if len(dS_pos) > 0 else 0.0

    return {
        "mean_cci": mean_cci,
        "mean_fstruct": mean_fstruct,
        "mean_mi": mean_mi,
        "mean_dS_pos": mean_dS_pos,
        "phi_initial": phi_snapshots[0],
        "phi_final": phi_snapshots[-1],
    }

# ============================================================
# ENSEMBLE
# ============================================================

rows = []
snapshots = []

for run in range(RUNS):
    result = run_traj_2d()
    rows.append({
        "run": run,
        "mean_cci": result["mean_cci"],
        "mean_fstruct": result["mean_fstruct"],
        "mean_mi": result["mean_mi"],
        "mean_dS_pos": result["mean_dS_pos"],
    })
    snapshots.append((result["phi_initial"], result["phi_final"]))

    if (run + 1) % 10 == 0:
        print(f"Completed {run + 1}/{RUNS}")

df = pd.DataFrame(rows)

# ============================================================
# ALPHA SCAN
# ============================================================

alpha_rows = []

for alpha in ALPHAS:
    ratio = df["mean_dS_pos"] / ((df["mean_mi"] + EPS) ** alpha)

    sp = spearmanr(df["mean_cci"], ratio)
    pe = pearsonr(df["mean_cci"], ratio)

    alpha_rows.append({
        "alpha": alpha,
        "spearman_r": sp.statistic,
        "spearman_p": sp.pvalue,
        "pearson_r": pe.statistic,
        "pearson_p": pe.pvalue,
    })

alpha_df = pd.DataFrame(alpha_rows)
best_idx = alpha_df["spearman_r"].idxmax()
best_alpha = float(alpha_df.loc[best_idx, "alpha"])

df["ratio_best"] = df["mean_dS_pos"] / ((df["mean_mi"] + EPS) ** best_alpha)

print("\n=== 2D alpha scan summary ===")
print(alpha_df.sort_values("spearman_r", ascending=False).head(10))
print(f"\nBest alpha (2D, by Spearman): {best_alpha:.3f}")

# ============================================================
# LOG-FIT
# ============================================================

mask = (df["mean_cci"] > 0) & (df["mean_dS_pos"] > 0) & (df["mean_mi"] > 0)

X = np.column_stack([
    np.log(df.loc[mask, "mean_dS_pos"].to_numpy()),
    -np.log(df.loc[mask, "mean_mi"].to_numpy())
])
y = np.log(df.loc[mask, "mean_cci"].to_numpy())

reg = LinearRegression()
reg.fit(X, y)

a0 = reg.intercept_
b_dS = reg.coef_[0]
b_I = reg.coef_[1]
r2 = reg.score(X, y)

print("\n=== 2D log-scaling fit ===")
print(f"log(CCI) = {a0:.4f} + {b_dS:.4f} log(dS_pos) + {b_I:.4f} (-log(MI))")
print(f"R^2 = {r2:.4f}")

# ============================================================
# SAVE
# ============================================================

df.to_csv("2d_cci_entropy_information_test.csv", index=False)
alpha_df.to_csv("2d_cci_alpha_scan.csv", index=False)

print("\nSaved:")
print(" - 2d_cci_entropy_information_test.csv")
print(" - 2d_cci_alpha_scan.csv")

# ============================================================
# PLOTS
# ============================================================

# 1. alpha scan
plt.figure(figsize=(7, 5))
plt.plot(alpha_df["alpha"], alpha_df["spearman_r"], marker="o", label="Spearman")
plt.plot(alpha_df["alpha"], alpha_df["pearson_r"], marker="s", label="Pearson")
plt.axvline(best_alpha, color="black", linestyle="--", label=fr"best $\alpha={best_alpha:.2f}$")
plt.xlabel(r"$\alpha$")
plt.ylabel("correlation with CCI")
plt.title(r"2D test of $CCI \sim \dot S / I^\alpha$")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("2d_cci_alpha_scan.png", dpi=180)
plt.show()

# 2. CCI vs best ratio
plt.figure(figsize=(6, 5))
plt.scatter(df["ratio_best"], df["mean_cci"], alpha=0.8)
plt.xlabel(r"$\overline{\dot S}_{cg}^{(+)} / (\overline{I}_{nn}+\varepsilon)^\alpha$")
plt.ylabel("mean CCI")
plt.title(fr"2D: CCI vs best structural ratio ($\alpha={best_alpha:.2f}$)")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("2d_cci_vs_best_ratio.png", dpi=180)
plt.show()

# 3. CCI vs MI
plt.figure(figsize=(6, 5))
plt.scatter(df["mean_mi"], df["mean_cci"], alpha=0.8)
plt.xlabel(r"$\overline{I}_{nn}$")
plt.ylabel("mean CCI")
plt.title("2D: CCI vs structural information")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("2d_cci_vs_mi.png", dpi=180)
plt.show()

# 4. CCI vs dS
plt.figure(figsize=(6, 5))
plt.scatter(df["mean_dS_pos"], df["mean_cci"], alpha=0.8)
plt.xlabel(r"$\overline{\dot S}_{cg}^{(+)}$")
plt.ylabel("mean CCI")
plt.title("2D: CCI vs entropy-production proxy")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("2d_cci_vs_dS.png", dpi=180)
plt.show()

# 5. show one example trajectory: initial/final field
phi0, phif = snapshots[0]

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
im0 = axes[0].imshow(phi0, origin="lower", aspect="auto")
axes[0].set_title("Initial field")
plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

im1 = axes[1].imshow(phif, origin="lower", aspect="auto")
axes[1].set_title("Final field")
plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

plt.tight_layout()
plt.savefig("2d_example_initial_final.png", dpi=180)
plt.show()
