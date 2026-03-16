import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

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

# from your previous calibration
TAU1 = 0.278899
TAU2 = 0.354858

# phase diagram grid
noise_vals = np.linspace(0.0, 1.0, 9)
pi_vals = np.linspace(0.0, 1.0, 9)

# number of runs per grid point
N_ENSEMBLE = 8

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
# INITIAL CONDITIONS WITH CONTROL PARAMETERS
# ============================================================

def controlled_initial(noise_amp, pi_amp):
    x = np.arange(N) - N / 2

    # mixture base
    vac = rng.choice([-1, 1]) * np.ones(N)
    kink = np.tanh((x - rng.uniform(-N / 4, N / 4)) / rng.uniform(4, 8))
    loc = np.exp(-((x - rng.uniform(-N / 4, N / 4)) / rng.uniform(3, 7)) ** 2)
    noise = rng.normal(0, noise_amp, N)

    # keep some structural variety
    a = rng.uniform(0, 1, 3)
    phi = a[0] * vac + a[1] * kink + a[2] * loc + noise
    phi /= (1 + np.std(phi))
    phi = np.clip(phi, -2.5, 2.5)

    pi = rng.normal(0, pi_amp, N)
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

def residual(phi):
    return lap(phi) - dV(phi)

def Fstruct(phi, pi):
    res = np.mean(np.abs(residual(phi)))
    act = np.mean(pi**2)
    corr = mutual_information(phi)
    return res + 0.1 * act - 0.2 * corr

def activity(phi, pi):
    return float(np.mean(pi**2) + 0.5 * np.mean(grad(phi)**2))

def CCI(phi, pi):
    prod = np.mean(np.abs(pi))
    coh = np.exp(-np.mean(np.abs(residual(phi))))
    corr = mutual_information(phi)
    return prod / (coh + corr + 1e-8)

def regime_from_cci(cci):
    if cci < TAU1:
        return 0   # ordered
    elif cci < TAU2:
        return 1   # critical
    else:
        return 2   # chaotic

# ============================================================
# RUN SINGLE TRAJECTORY
# ============================================================

def run_traj(noise_amp, pi_amp):
    phi0, pi0 = controlled_initial(noise_amp, pi_amp)
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

    mean_cci = np.mean([CCI(phi_t[:, i], pi_t[:, i]) for i in range(len(sol.t))])
    mean_f = np.mean([Fstruct(phi_t[:, i], pi_t[:, i]) for i in range(len(sol.t))])

    phi_f = phi_t[:, -1]
    pi_f = pi_t[:, -1]

    return {
        "CCI": mean_cci,
        "Fstruct": mean_f,
        "MI_final": mutual_information(phi_f),
        "Activity_final": activity(phi_f, pi_f),
        "regime": regime_from_cci(mean_cci)
    }

# ============================================================
# PHASE DIAGRAM SWEEP
# ============================================================

rows = []

for noise_amp in noise_vals:
    for pi_amp in pi_vals:
        local_runs = [run_traj(noise_amp, pi_amp) for _ in range(N_ENSEMBLE)]

        mean_cci = np.mean([r["CCI"] for r in local_runs])
        mean_f = np.mean([r["Fstruct"] for r in local_runs])
        mean_mi = np.mean([r["MI_final"] for r in local_runs])
        mean_act = np.mean([r["Activity_final"] for r in local_runs])

        regimes = [r["regime"] for r in local_runs]
        dominant_regime = max(set(regimes), key=regimes.count)

        rows.append({
            "noise_amp": noise_amp,
            "pi_amp": pi_amp,
            "mean_cci": mean_cci,
            "mean_fstruct": mean_f,
            "mean_mi": mean_mi,
            "mean_activity": mean_act,
            "dominant_regime": dominant_regime
        })

df = pd.DataFrame(rows)
df.to_csv("cci_phase_diagram.csv", index=False)
print("Saved: cci_phase_diagram.csv")

# ============================================================
# BUILD MATRICES
# ============================================================

shape = (len(noise_vals), len(pi_vals))

cci_grid = np.full(shape, np.nan)
regime_grid = np.full(shape, np.nan)
mi_grid = np.full(shape, np.nan)
f_grid = np.full(shape, np.nan)

for i, nv in enumerate(noise_vals):
    for j, pv in enumerate(pi_vals):
        sub = df[(df["noise_amp"] == nv) & (df["pi_amp"] == pv)].iloc[0]
        cci_grid[i, j] = sub["mean_cci"]
        regime_grid[i, j] = sub["dominant_regime"]
        mi_grid[i, j] = sub["mean_mi"]
        f_grid[i, j] = sub["mean_fstruct"]

# ============================================================
# PLOTS
# ============================================================

extent = [pi_vals.min(), pi_vals.max(), noise_vals.min(), noise_vals.max()]

# 1. Mean CCI phase diagram
plt.figure(figsize=(7, 5))
plt.imshow(cci_grid, origin="lower", aspect="auto", extent=extent)
plt.colorbar(label="mean CCI")
plt.xlabel(r"initial momentum amplitude $\sigma_\pi$")
plt.ylabel(r"noise amplitude $\sigma_{\rm noise}$")
plt.title("CCI phase diagram")
plt.tight_layout()
plt.savefig("phase_diagram_cci.png", dpi=180)
plt.show()

# 2. Dominant regime phase diagram
plt.figure(figsize=(7, 5))
plt.imshow(regime_grid, origin="lower", aspect="auto", extent=extent, vmin=0, vmax=2)
cbar = plt.colorbar()
cbar.set_ticks([0, 1, 2])
cbar.set_ticklabels(["ordered", "critical", "chaotic"])
plt.xlabel(r"initial momentum amplitude $\sigma_\pi$")
plt.ylabel(r"noise amplitude $\sigma_{\rm noise}$")
plt.title("Dominant regime phase diagram")
plt.tight_layout()
plt.savefig("phase_diagram_regime.png", dpi=180)
plt.show()

# 3. MI phase diagram
plt.figure(figsize=(7, 5))
plt.imshow(mi_grid, origin="lower", aspect="auto", extent=extent)
plt.colorbar(label="mean final MI")
plt.xlabel(r"initial momentum amplitude $\sigma_\pi$")
plt.ylabel(r"noise amplitude $\sigma_{\rm noise}$")
plt.title("Structural information phase diagram")
plt.tight_layout()
plt.savefig("phase_diagram_mi.png", dpi=180)
plt.show()

# 4. Fstruct phase diagram
plt.figure(figsize=(7, 5))
plt.imshow(f_grid, origin="lower", aspect="auto", extent=extent)
plt.colorbar(label="mean Fstruct")
plt.xlabel(r"initial momentum amplitude $\sigma_\pi$")
plt.ylabel(r"noise amplitude $\sigma_{\rm noise}$")
plt.title("Structural free energy phase diagram")
plt.tight_layout()
plt.savefig("phase_diagram_fstruct.png", dpi=180)
plt.show()