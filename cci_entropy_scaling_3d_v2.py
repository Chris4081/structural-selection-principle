import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr
from sklearn.linear_model import LinearRegression

# ============================================================
# PARAMETERS
# ============================================================

NX = 24
NY = 24
NZ = 24

dx = 1.0
dt = 0.02
Tmax = 8.0
NSTEPS = int(Tmax / dt)

rng = np.random.default_rng(42)

PHI_MAX = 3.0
EPS = 1e-8
N_BINS_MI = 16
N_BINS_ENT = 32

RUNS = 24
SAMPLE_EVERY = 10

# Same scan window as before for direct comparability.
# If best alpha hits the upper boundary again, interpret as lower bound.
ALPHAS = np.linspace(0.5, 3.5, 31)

# ============================================================
# 3D phi^4 MODEL
# ============================================================

def dV(phi):
    return phi**3 - phi

def lap3d(phi):
    return (
        np.roll(phi, -1, axis=0) + np.roll(phi, 1, axis=0) +
        np.roll(phi, -1, axis=1) + np.roll(phi, 1, axis=1) +
        np.roll(phi, -1, axis=2) + np.roll(phi, 1, axis=2) -
        6.0 * phi
    ) / dx**2

def rhs(phi, pi):
    phi_dot = pi
    pi_dot = lap3d(phi) - dV(phi)
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

def random_initial_3d():
    x = np.arange(NX) - NX / 2
    y = np.arange(NY) - NY / 2
    z = np.arange(NZ) - NZ / 2
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

    # vacuum background
    vac = rng.choice([-1.0, 1.0]) * np.ones((NX, NY, NZ))

    # 3D wall-like structure
    theta = rng.uniform(0, np.pi)
    phi_ang = rng.uniform(0, 2 * np.pi)
    width = rng.uniform(2.5, 6.0)
    offset = rng.uniform(-0.2 * NX, 0.2 * NX)

    nx = np.sin(theta) * np.cos(phi_ang)
    ny = np.sin(theta) * np.sin(phi_ang)
    nz = np.cos(theta)

    wall_coord = nx * X + ny * Y + nz * Z - offset
    wall = np.tanh(wall_coord / width)

    # localized lump
    x0 = rng.uniform(-0.2 * NX, 0.2 * NX)
    y0 = rng.uniform(-0.2 * NY, 0.2 * NY)
    z0 = rng.uniform(-0.2 * NZ, 0.2 * NZ)
    sigma = rng.uniform(2.0, 5.0)
    amp = rng.uniform(0.7, 1.5)

    lump = amp * np.exp(
        -((X - x0) ** 2 + (Y - y0) ** 2 + (Z - z0) ** 2) / (2 * sigma**2)
    )

    # noise
    noise = rng.normal(0.0, 0.5, size=(NX, NY, NZ))

    # mixture
    a = rng.uniform(0, 1, 4)
    phi = a[0] * vac + a[1] * wall + a[2] * lump + a[3] * noise
    phi /= (1.0 + np.std(phi))
    phi = np.clip(phi, -2.5, 2.5)

    pi = rng.normal(0.0, 0.3, size=(NX, NY, NZ))
    return phi, pi

# ============================================================
# DIAGNOSTICS
# ============================================================

def residual(phi):
    return lap3d(phi) - dV(phi)

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

def mutual_information_3d(phi):
    """
    Average nearest-neighbour mutual information in x, y, z directions.
    Returns mean and per-direction values.
    """
    x1 = phi.ravel()
    x2 = np.roll(phi, -1, axis=0).ravel()
    y2 = np.roll(phi, -1, axis=1).ravel()
    z2 = np.roll(phi, -1, axis=2).ravel()

    mi_x = mutual_information_pair(x1, x2)
    mi_y = mutual_information_pair(x1, y2)
    mi_z = mutual_information_pair(x1, z2)

    return (mi_x + mi_y + mi_z) / 3.0, mi_x, mi_y, mi_z

def gradient_anisotropy(phi):
    """
    Per-direction mean squared gradient.
    """
    gx = (np.roll(phi, -1, axis=0) - np.roll(phi, 1, axis=0)) / 2.0
    gy = (np.roll(phi, -1, axis=1) - np.roll(phi, 1, axis=1)) / 2.0
    gz = (np.roll(phi, -1, axis=2) - np.roll(phi, 1, axis=2)) / 2.0
    return float(np.mean(gx**2)), float(np.mean(gy**2)), float(np.mean(gz**2))

def Fstruct(phi, pi):
    """
    Structural free-energy proxy.
    """
    res = np.mean(np.abs(residual(phi)))
    act = np.mean(pi**2)
    corr, _, _, _ = mutual_information_3d(phi)
    return float(res + 0.1 * act - 0.2 * corr)

def CCI(phi, pi):
    """
    Reduced numerical CCI proxy.
    """
    prod = np.mean(np.abs(pi))
    coh = np.exp(-np.mean(np.abs(residual(phi))))
    corr, _, _, _ = mutual_information_3d(phi)
    return float(prod / (coh + corr + EPS))

# ============================================================
# SINGLE TRAJECTORY
# ============================================================

def run_traj_3d():
    phi, pi = random_initial_3d()

    cci_series = []
    f_series = []
    mi_series = []
    mi_x_series = []
    mi_y_series = []
    mi_z_series = []
    grad_x2_series = []
    grad_y2_series = []
    grad_z2_series = []
    s_series = []
    phi_snapshots = []

    for step in range(NSTEPS + 1):
        if step % SAMPLE_EVERY == 0:
            cci_series.append(CCI(phi, pi))
            f_series.append(Fstruct(phi, pi))
            mi_mean, mi_x, mi_y, mi_z = mutual_information_3d(phi)
            mi_series.append(mi_mean)
            mi_x_series.append(mi_x)
            mi_y_series.append(mi_y)
            mi_z_series.append(mi_z)
            gx2, gy2, gz2 = gradient_anisotropy(phi)
            grad_x2_series.append(gx2)
            grad_y2_series.append(gy2)
            grad_z2_series.append(gz2)
            s_series.append(coarse_grained_entropy(phi))
            phi_snapshots.append(phi.copy())

        if step < NSTEPS:
            phi, pi = rk4_step(phi, pi, dt)

    cci_series = np.array(cci_series)
    f_series = np.array(f_series)
    mi_series = np.array(mi_series)
    s_series = np.array(s_series)

    dt_sample = SAMPLE_EVERY * dt
    dS = np.diff(s_series) / dt_sample
    dS_pos = np.maximum(0.0, dS)

    mean_cci = float(np.mean(cci_series[:-1]))
    mean_fstruct = float(np.mean(f_series[:-1]))
    mean_mi = float(np.mean(mi_series[:-1]))
    mean_dS_pos = float(np.mean(dS_pos)) if len(dS_pos) > 0 else 0.0

    mi_x_series = np.array(mi_x_series)
    mi_y_series = np.array(mi_y_series)
    mi_z_series = np.array(mi_z_series)
    grad_x2_series = np.array(grad_x2_series)
    grad_y2_series = np.array(grad_y2_series)
    grad_z2_series = np.array(grad_z2_series)

    return {
        "mean_cci": mean_cci,
        "mean_fstruct": mean_fstruct,
        "mean_mi": mean_mi,
        "mean_mi_x": float(np.mean(mi_x_series[:-1])),
        "mean_mi_y": float(np.mean(mi_y_series[:-1])),
        "mean_mi_z": float(np.mean(mi_z_series[:-1])),
        "mean_grad_x2": float(np.mean(grad_x2_series[:-1])),
        "mean_grad_y2": float(np.mean(grad_y2_series[:-1])),
        "mean_grad_z2": float(np.mean(grad_z2_series[:-1])),
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
    result = run_traj_3d()

    rows.append({
        "run": run,
        "mean_cci": result["mean_cci"],
        "mean_fstruct": result["mean_fstruct"],
        "mean_mi": result["mean_mi"],
        "mean_mi_x": result["mean_mi_x"],
        "mean_mi_y": result["mean_mi_y"],
        "mean_mi_z": result["mean_mi_z"],
        "mean_grad_x2": result["mean_grad_x2"],
        "mean_grad_y2": result["mean_grad_y2"],
        "mean_grad_z2": result["mean_grad_z2"],
        "mean_dS_pos": result["mean_dS_pos"],
    })

    snapshots.append((result["phi_initial"], result["phi_final"]))

    if (run + 1) % 4 == 0:
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
df["log_ratio_best"] = np.log10(df["ratio_best"] + EPS)

print("\n=== 3D alpha scan summary ===")
print(alpha_df.sort_values("spearman_r", ascending=False).head(10))
print(f"\nBest alpha (3D, by Spearman): {best_alpha:.3f}")

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

print("\n=== 3D log-scaling fit ===")
print(f"log(CCI) = {a0:.4f} + {b_dS:.4f} log(dS_pos) + {b_I:.4f} (-log(MI))")
print(f"R^2 = {r2:.4f}")

# ============================================================
# SAVE
# ============================================================

df.to_csv("3d_cci_entropy_information_v2.csv", index=False)
alpha_df.to_csv("3d_cci_alpha_scan.csv", index=False)

print("\nSaved:")
print(" - 3d_cci_entropy_information_test.csv")
print(" - 3d_cci_alpha_scan.csv")

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
plt.title(r"3D test of $CCI \sim \dot S / I^\alpha$")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("3d_cci_alpha_scan.png", dpi=180)
plt.show()

# 2. CCI vs best ratio
plt.figure(figsize=(6, 5))
plt.scatter(df["ratio_best"], df["mean_cci"], alpha=0.8)
plt.xlabel(r"$\overline{\dot S}_{cg}^{(+)} / (\overline{I}_{nn}+\varepsilon)^\alpha$")
plt.ylabel("mean CCI")
plt.title(fr"3D: CCI vs best structural ratio ($\alpha={best_alpha:.2f}$)")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("3d_cci_vs_best_ratio.png", dpi=180)
plt.show()

# 3. CCI vs MI
plt.figure(figsize=(6, 5))
plt.scatter(df["mean_mi"], df["mean_cci"], alpha=0.8)
plt.xlabel(r"$\overline{I}_{nn}$")
plt.ylabel("mean CCI")
plt.title("3D: CCI vs structural information")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("3d_cci_vs_mi.png", dpi=180)
plt.show()

# 4. CCI vs dS
plt.figure(figsize=(6, 5))
plt.scatter(df["mean_dS_pos"], df["mean_cci"], alpha=0.8)
plt.xlabel(r"$\overline{\dot S}_{cg}^{(+)}$")
plt.ylabel("mean CCI")
plt.title("3D: CCI vs entropy-production proxy")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("3d_cci_vs_dS.png", dpi=180)
plt.show()

# 5. Example central slice
phi0, phif = snapshots[0]
zmid = NZ // 2

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

im0 = axes[0].imshow(phi0[:, :, zmid], origin="lower", aspect="auto")
axes[0].set_title("Initial field (central z-slice)")
plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

im1 = axes[1].imshow(phif[:, :, zmid], origin="lower", aspect="auto")
axes[1].set_title("Final field (central z-slice)")
plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

plt.tight_layout()
plt.savefig("3d_example_initial_final_slice.png", dpi=180)
plt.show()