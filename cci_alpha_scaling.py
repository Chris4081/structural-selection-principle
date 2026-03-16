import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.stats import spearmanr, pearsonr
from sklearn.linear_model import LinearRegression

# ============================================================
# PARAMETERS
# ============================================================

N = 128
dx = 1.0
dt = 0.02
Tmax = 24.0
rng = np.random.default_rng(42)

PHI_MAX = 3.0
EPS = 1e-8
N_BINS_MI = 16
N_BINS_ENT = 32
RUNS = 120

# from your current calibration
TAU1 = 0.278899
TAU2 = 0.354858

# alpha scan
ALPHAS = np.linspace(0.5, 2.5, 21)

# ============================================================
# FIELD MODEL
# ============================================================

def dV(phi):
    return phi**3 - phi

def lap(phi):
    return (np.roll(phi, -1) + np.roll(phi, 1) - 2 * phi) / dx**2

def grad(phi):
    return (np.roll(phi, -1) - np.roll(phi, 1)) / (2 * dx)

def rhs(t, y):
    phi = y[:N]
    pi = y[N:]
    return np.concatenate([pi, lap(phi) - dV(phi)])

# ============================================================
# INITIAL CONDITIONS
# ============================================================

def random_initial():
    x = np.arange(N) - N / 2

    vac = rng.choice([-1, 1]) * np.ones(N)
    kink = np.tanh((x - rng.uniform(-N / 4, N / 4)) / rng.uniform(4, 8))
    loc = np.exp(-((x - rng.uniform(-N / 4, N / 4)) / rng.uniform(3, 7)) ** 2)
    noise = rng.normal(0, 0.5, N)

    a = rng.uniform(0, 1, 4)
    phi = a[0] * vac + a[1] * kink + a[2] * loc + a[3] * noise
    phi /= 1 + np.std(phi)
    phi = np.clip(phi, -2.5, 2.5)

    pi = rng.normal(0, 0.3, N)
    return phi, pi

# ============================================================
# DIAGNOSTICS
# ============================================================

def mutual_information(phi):
    x = phi[:-1]
    y = phi[1:]

    hist, _, _ = np.histogram2d(
        x, y,
        bins=N_BINS_MI,
        range=[[-PHI_MAX, PHI_MAX], [-PHI_MAX, PHI_MAX]]
    )

    pxy = hist / max(np.sum(hist), 1.0)
    px = np.sum(pxy, axis=1)
    py = np.sum(pxy, axis=0)

    mi = 0.0
    for i in range(len(px)):
        for j in range(len(py)):
            if pxy[i, j] > 0 and px[i] > 0 and py[j] > 0:
                mi += pxy[i, j] * np.log(pxy[i, j] / (px[i] * py[j] + EPS))
    return float(mi)

def coarse_grained_entropy(phi):
    hist, _ = np.histogram(phi, bins=N_BINS_ENT, range=[-PHI_MAX, PHI_MAX])
    p = hist.astype(float) / max(np.sum(hist), 1.0)
    p = p[p > 0]
    return float(-np.sum(p * np.log(p + EPS)))

def residual(phi):
    return lap(phi) - dV(phi)

def Fstruct(phi, pi):
    res = np.mean(np.abs(residual(phi)))
    act = np.mean(pi**2)
    corr = mutual_information(phi)
    return res + 0.1 * act - 0.2 * corr

def CCI(phi, pi):
    prod = np.mean(np.abs(pi))
    coh = np.exp(-np.mean(np.abs(residual(phi))))
    corr = mutual_information(phi)
    return prod / (coh + corr + 1e-8)

def regime_from_cci(cci):
    if cci < TAU1:
        return "ordered"
    elif cci < TAU2:
        return "critical"
    else:
        return "chaotic"

# ============================================================
# SINGLE TRAJECTORY
# ============================================================

def run_traj():
    phi0, pi0 = random_initial()
    y0 = np.concatenate([phi0, pi0])

    t_eval = np.linspace(0, Tmax, int(Tmax / dt) + 1)

    sol = solve_ivp(
        rhs,
        [0, Tmax],
        y0,
        t_eval=t_eval,
        rtol=1e-6,
        atol=1e-8
    )

    phi_t = sol.y[:N]
    pi_t = sol.y[N:]
    times = sol.t

    cci_series = np.array([CCI(phi_t[:, i], pi_t[:, i]) for i in range(len(times))])
    f_series = np.array([Fstruct(phi_t[:, i], pi_t[:, i]) for i in range(len(times))])
    mi_series = np.array([mutual_information(phi_t[:, i]) for i in range(len(times))])
    s_series = np.array([coarse_grained_entropy(phi_t[:, i]) for i in range(len(times))])

    dS = np.diff(s_series) / dt
    dS_pos = np.maximum(0.0, dS)

    mean_cci = float(np.mean(cci_series[:-1]))
    mean_f = float(np.mean(f_series[:-1]))
    mean_mi = float(np.mean(mi_series[:-1]))
    mean_dS_pos = float(np.mean(dS_pos)) if len(dS_pos) > 0 else 0.0

    return {
        "mean_cci": mean_cci,
        "mean_fstruct": mean_f,
        "mean_mi": mean_mi,
        "mean_dS_pos": mean_dS_pos,
        "regime": regime_from_cci(mean_cci),
    }

# ============================================================
# BUILD ENSEMBLE
# ============================================================

rows = [run_traj() for _ in range(RUNS)]
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
        "pearson_p": pe.pvalue
    })

alpha_df = pd.DataFrame(alpha_rows)
best_idx = alpha_df["spearman_r"].idxmax()
best_alpha = float(alpha_df.loc[best_idx, "alpha"])

df["ratio_best"] = df["mean_dS_pos"] / ((df["mean_mi"] + EPS) ** best_alpha)

print("\n=== Alpha scan summary ===")
print(alpha_df.sort_values("spearman_r", ascending=False).head(10))
print(f"\nBest alpha (by Spearman): {best_alpha:.3f}")

# ============================================================
# LOG-FIT
# log(CCI) = a + b log(dS) - alpha log(I)
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

print("\n=== Log-scaling fit ===")
print(f"log(CCI) = {a0:.4f} + {b_dS:.4f} log(dS_pos) + {b_I:.4f} (-log(MI))")
print(f"R^2 = {r2:.4f}")

# ============================================================
# SAVE
# ============================================================

df.to_csv("cci_entropy_information_test.csv", index=False)
alpha_df.to_csv("cci_alpha_scan.csv", index=False)

print("\nSaved:")
print(" - cci_entropy_information_test.csv")
print(" - cci_alpha_scan.csv")

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
plt.title(r"Testing $CCI \sim \dot S / I^\alpha$")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("cci_alpha_scan.png", dpi=180)
plt.show()

# 2. CCI vs best ratio
plt.figure(figsize=(6, 5))
for reg_name in ["ordered", "critical", "chaotic"]:
    sub = df[df["regime"] == reg_name]
    plt.scatter(sub["ratio_best"], sub["mean_cci"], label=reg_name, alpha=0.8)
plt.xlabel(r"$\overline{\dot S}_{cg}^{(+)} / (\overline{I}_{nn}+\varepsilon)^\alpha$")
plt.ylabel("mean CCI")
plt.title(fr"CCI vs best structural ratio ($\alpha={best_alpha:.2f}$)")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("cci_vs_best_ratio.png", dpi=180)
plt.show()

# 3. CCI vs dS only
plt.figure(figsize=(6, 5))
plt.scatter(df["mean_dS_pos"], df["mean_cci"], alpha=0.8)
plt.xlabel(r"$\overline{\dot S}_{cg}^{(+)}$")
plt.ylabel("mean CCI")
plt.title("CCI vs entropy-production proxy")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("cci_vs_dS.png", dpi=180)
plt.show()

# 4. CCI vs MI only
plt.figure(figsize=(6, 5))
plt.scatter(df["mean_mi"], df["mean_cci"], alpha=0.8)
plt.xlabel(r"$\overline{I}_{nn}$")
plt.ylabel("mean CCI")
plt.title("CCI vs structural information")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("cci_vs_mi.png", dpi=180)
plt.show()