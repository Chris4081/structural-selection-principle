import numpy as np
from scipy.optimize import differential_evolution

# Observed values only for comparison after optimisation
ALPHA_OBS = 1 / 137.035999084
MU_OBS = 1 / 1836.15267343
LAMBDA_OBS = 1e-122

EPS = 1e-30


def softplus(x):
    return np.log1p(np.exp(np.clip(x, -80, 80)))


def outside_band(log_x, log_min, log_max, width=1.0):
    """Penalty only outside a viable logarithmic band."""
    low = softplus((log_min - log_x) / width) ** 2
    high = softplus((log_x - log_max) / width) ** 2
    return low + high


def structural_score(x):
    """
    x = [log_alpha, log_mu, log_lambda]

    No observed constants are used in the score.
    Only broad physics-inspired viability constraints are used.
    """
    log_alpha, log_mu, log_lam = x

    alpha = np.exp(log_alpha)
    mu = np.exp(log_mu)
    lam = np.exp(log_lam)

    # ------------------------------------------------------------
    # 1. Atomic viability
    # ------------------------------------------------------------
    # Nonrelativistic atoms require roughly alpha << 1.
    # If alpha too large: relativistic/unstable atomic structure.
    atom_relativistic = softplus((alpha - 0.12) / 0.02) ** 2

    # If alpha too tiny: binding becomes too weak structurally.
    atom_inert = softplus((1e-5 - alpha) / 1e-5) ** 2

    # ------------------------------------------------------------
    # 2. Atomic/nuclear scale separation
    # ------------------------------------------------------------
    # mu = me/mp must be much less than 1.
    # Too large: electrons not well separated from baryons.
    hierarchy_fail = softplus((mu - 0.02) / 0.005) ** 2

    # Too small: chemistry becomes dynamically weak/inert.
    hierarchy_inert = softplus((1e-7 - mu) / 1e-7) ** 2

    # ------------------------------------------------------------
    # 3. Chemistry activity window
    # ------------------------------------------------------------
    # Atomic binding scale ~ alpha^2 * me.
    # In proton-mass units this is roughly alpha^2 * mu.
    chemistry_strength = alpha**2 * mu

    # Need enough chemistry, but not too violent.
    chemistry_too_weak = softplus((1e-12 - chemistry_strength) / 1e-12) ** 2
    chemistry_too_strong = softplus((chemistry_strength - 1e-4) / 1e-4) ** 2

    # ------------------------------------------------------------
    # 4. Stellar viability proxy
    # ------------------------------------------------------------
    # Very rough anthropic-style structural condition:
    # alpha should not be too large and mu should provide scale separation.
    # Use product alpha^2 / mu as a proxy for electromagnetic vs mass hierarchy stress.
    stellar_proxy = alpha**2 / max(mu, EPS)

    stellar_too_small = softplus((1e-4 - stellar_proxy) / 1e-4) ** 2
    stellar_too_large = softplus((stellar_proxy - 10.0) / 10.0) ** 2

    # ------------------------------------------------------------
    # 5. Cosmological structure formation
    # ------------------------------------------------------------
    # Lambda must be tiny in Planck units for long-lived structure.
    # This is NOT targeted to observed 1e-122; it is a broad viability band.
    lambda_band = outside_band(
        log_lam,
        np.log(1e-140),
        np.log(1e-80),
        width=6.0,
    )

    # Too large Lambda: vacuum domination before structure formation.
    lambda_too_large = softplus((log_lam - np.log(1e-90)) / 8.0) ** 2

    # Exactly zero is not allowed in this scan; lower bound already handles this.

    # ------------------------------------------------------------
    # MAAT sector defects
    # ------------------------------------------------------------
    H_defect = atom_relativistic + atom_inert
    B_defect = hierarchy_fail + hierarchy_inert
    S_defect = chemistry_too_weak + chemistry_too_strong
    V_defect = stellar_too_small + stellar_too_large
    R_defect = lambda_band + lambda_too_large

    defects = np.array([H_defect, B_defect, S_defect, V_defect, R_defect])
    supports = 1.0 / (1.0 + defects)

    eps = 1e-8
    F_maat = -np.sum(np.log((eps + supports) / (eps + 1.0)))

    stability = min(
        supports[4],
        (supports[0] * supports[1] * supports[2] * supports[3]) ** 0.25,
    )

    # Conservative stability penalty
    return F_maat + 5.0 * (1.0 - stability)


def main():
    bounds = [
        (np.log(1e-7), np.log(0.3)),       # alpha
        (np.log(1e-9), np.log(1e-1)),      # mu = me/mp
        (np.log(1e-160), np.log(1e-50)),   # lambda
    ]

    result = differential_evolution(
        structural_score,
        bounds,
        seed=42,
        maxiter=800,
        popsize=25,
        polish=True,
        tol=1e-10,
        workers=1,
    )

    alpha_fit, mu_fit, lambda_fit = np.exp(result.x)

    print("\n=== MAAT Constants Structure Scan v3 ===")
    print("Best score:", result.fun)

    print("\n--- Fitted structural values ---")
    print("alpha_fit :", alpha_fit)
    print("mu_fit    :", mu_fit)
    print("lambda_fit:", lambda_fit)

    print("\n--- Observed values, comparison only ---")
    print("alpha_obs :", ALPHA_OBS)
    print("mu_obs    :", MU_OBS)
    print("lambda_obs:", LAMBDA_OBS)

    print("\n--- Deviation factors fit / observed ---")
    print("alpha factor :", alpha_fit / ALPHA_OBS)
    print("mu factor    :", mu_fit / MU_OBS)
    print("lambda factor:", lambda_fit / LAMBDA_OBS)

    print("\n--- Log10 values ---")
    print("log10 alpha_fit :", np.log10(alpha_fit))
    print("log10 mu_fit    :", np.log10(mu_fit))
    print("log10 lambda_fit:", np.log10(lambda_fit))

    print("\n--- Derived structural proxies ---")
    print("chemistry_strength alpha^2*mu:", alpha_fit**2 * mu_fit)
    print("stellar_proxy alpha^2/mu     :", alpha_fit**2 / mu_fit)


if __name__ == "__main__":
    main()
