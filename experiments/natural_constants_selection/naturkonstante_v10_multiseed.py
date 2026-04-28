import numpy as np
from scipy.optimize import differential_evolution

# --- reuse functions from v9 ---
from naturkonstante_v9_rg_maat import (
    total_score,
    qED_running_alpha,
    derive_mu,
    ALPHA_OBS,
    MU_OBS,
    LAMBDA_OBS,
)

N_SEEDS = 20

results = []

bounds = [
    (np.log(1e-4), np.log(0.3)),       # alpha_uv
    (np.log(1.0), np.log(1e8)),        # RG span
    (np.log(1e-8), np.log(1e-2)),      # y_e
    (np.log(1.0), np.log(1e6)),        # QCD proxy
    (np.log(1e-150), np.log(1e-90)),   # lambda
]

print("\n=== MAAT Constants v10 Multi-Seed RG + MAAT ===\n")

for seed in range(1, N_SEEDS + 1):

    res = differential_evolution(
        total_score,
        bounds,
        seed=seed,
        maxiter=800,
        popsize=20,
        polish=True,
        tol=1e-9,
        workers=1,
    )

    x = res.x

    alpha_uv = np.exp(x[0])
    rg_span = np.exp(x[1])
    y_e = np.exp(x[2])
    qcd = np.exp(x[3])
    lam = np.exp(x[4])

    alpha_ir = qED_running_alpha(alpha_uv, rg_span)
    mu_eff = derive_mu(y_e, qcd)

    results.append([alpha_ir, mu_eff, lam, res.fun])

    print(
        f"seed={seed:2d} "
        f"score={res.fun:.6f} "
        f"alpha={alpha_ir:.6g} "
        f"mu={mu_eff:.6g} "
        f"lambda={lam:.3e} "
        f"log10=({np.log10(alpha_ir):.3f}, {np.log10(mu_eff):.3f}, {np.log10(lam):.3f})"
    )


results = np.array(results)

alpha_vals = results[:, 0]
mu_vals = results[:, 1]
lam_vals = results[:, 2]
scores = results[:, 3]

print("\n--- Summary fitted values ---")

def summary(name, arr):
    print(
        f"{name:15s} mean={arr.mean():.6g} "
        f"std={arr.std():.3e} "
        f"min={arr.min():.6g} "
        f"max={arr.max():.6g}"
    )

summary("alpha", alpha_vals)
summary("mu", mu_vals)
summary("lambda", lam_vals)

print("\n--- Summary log10 values ---")
summary("log10 alpha", np.log10(alpha_vals))
summary("log10 mu", np.log10(mu_vals))
summary("log10 lambda", np.log10(lam_vals))

print("\n--- Summary deviation factors ---")
summary("alpha factor", alpha_vals / ALPHA_OBS)
summary("mu factor", mu_vals / MU_OBS)
summary("lambda factor", lam_vals / LAMBDA_OBS)

best_idx = np.argmin(scores)

print("\n--- Best seed ---")
print("seed:", best_idx + 1)
print("score:", scores[best_idx])
print("alpha:", alpha_vals[best_idx])
print("mu:", mu_vals[best_idx])
print("lambda:", lam_vals[best_idx])

print("\nObserved comparison only:")
print("alpha_obs :", ALPHA_OBS)
print("mu_obs    :", MU_OBS)
print("lambda_obs:", LAMBDA_OBS)
