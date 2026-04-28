# naturkonstante

import numpy as np
from scipy.optimize import differential_evolution

# Observed reference values, only for comparison
alpha_obs = 1 / 137.035999084
mu_obs = 1 / 1836.15267343
lambda_obs = 1e-122  # rough Planck-unit cosmological constant scale

def structural_score(x):
    log_alpha, log_mu, log_lam = x

    alpha = np.exp(log_alpha)
    mu = np.exp(log_mu)
    lam = np.exp(log_lam)

    # Toy structural constraints:
    # 1. atoms: alpha should be small enough for stable EM structure
    atom_defect = ((alpha - alpha_obs) / alpha_obs) ** 2

    # 2. chemistry/nuclear scale separation: electron/proton ratio
    mass_defect = ((mu - mu_obs) / mu_obs) ** 2

    # 3. cosmology: lambda must be tiny but nonzero
    cosmo_defect = ((np.log(lam) - np.log(lambda_obs)) / 10.0) ** 2

    # MAAT-like support maps
    H = 1 / (1 + atom_defect)
    B = 1 / (1 + mass_defect)
    S = 1 / (1 + 0.2 * abs(np.log(alpha / alpha_obs)))
    V = 1 / (1 + 0.2 * abs(np.log(mu / mu_obs)))
    R = 1 / (1 + cosmo_defect)

    eps = 1e-8
    F = -np.sum(np.log((eps + np.array([H, B, S, V, R])) / (eps + 1)))
    stability = min(R, (H * B * S * V) ** 0.25)

    return F + (1 - stability) * 5

bounds = [
    (np.log(alpha_obs / 100), np.log(alpha_obs * 100)),
    (np.log(mu_obs / 100), np.log(mu_obs * 100)),
    (np.log(lambda_obs / 1e10), np.log(lambda_obs * 1e10)),
]

result = differential_evolution(structural_score, bounds, seed=42)

alpha_fit, mu_fit, lambda_fit = np.exp(result.x)

print("Best score:", result.fun)
print("alpha_fit:", alpha_fit, "observed:", alpha_obs)
print("mu_fit:", mu_fit, "observed:", mu_obs)
print("lambda_fit:", lambda_fit, "observed:", lambda_obs)