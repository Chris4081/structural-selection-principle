import numpy as np
from scipy.optimize import differential_evolution

ALPHA_OBS = 1 / 137.035999084
MU_OBS = 1 / 1836.15267343
LAMBDA_OBS = 1e-122

EPS = 1e-30


def softplus(x):
    return np.log1p(np.exp(np.clip(x, -80, 80)))


def band_penalty_log(log_x, log_min, log_max, width=1.0):
    low = softplus((log_min - log_x) / width) ** 2
    high = softplus((log_x - log_max) / width) ** 2
    return low + high


def structural_score(x):
    log_alpha, log_mu, log_lam = x

    alpha = np.exp(log_alpha)
    mu = np.exp(log_mu)
    lam = np.exp(log_lam)

    # ----------------------------
    # 1. Atomic stability
    # ----------------------------
    # alpha must be nonrelativistic but not electromagnetically inert.
    atom_band = band_penalty_log(
        log_alpha,
        np.log(1e-3),
        np.log(0.08),
        width=0.8,
    )

    # ----------------------------
    # 2. Mass hierarchy
    # ----------------------------
    # electrons must be much lighter than protons, but not infinitely light.
    mu_band = band_penalty_log(
        log_mu,
        np.log(1e-5),
        np.log(5e-3),
        width=0.8,
    )

    # ----------------------------
    # 3. Chemistry activity
    # ----------------------------
    # Binding scale proxy: E_chem / m_p ~ alpha^2 * mu.
    chemistry = alpha**2 * mu

    chemistry_band = band_penalty_log(
        np.log(max(chemistry, EPS)),
        np.log(1e-8),
        np.log(1e-5),
        width=0.9,
    )

    # ----------------------------
    # 4. Stellar fusion proxy
    # ----------------------------
    # Very rough Coulomb / mass-hierarchy proxy.
    # If alpha too small, fusion too easy/short-lived stars.
    # If alpha too large, fusion too suppressed.
    fusion_proxy = alpha**2 / max(mu, EPS)

    fusion_band = band_penalty_log(
        np.log(max(fusion_proxy, EPS)),
        np.log(1e-2),
        np.log(3.0),
        width=0.9,
    )

    # ----------------------------
    # 5. Structure formation / Lambda
    # ----------------------------
    # Broad non-targeted Planck-unit viability:
    # lambda must be tiny enough for long-lived structure.
    lambda_band = band_penalty_log(
        log_lam,
        np.log(1e-140),
        np.log(1e-110),
        width=5.0,
    )

    # Extra pressure against too-large vacuum domination.
    lambda_domination = softplus((log_lam - np.log(1e-118)) / 5.0) ** 2

    # ----------------------------
    # MAAT sectors
    # ----------------------------
    H_defect = atom_band
    B_defect = mu_band
    S_defect = chemistry_band
    V_defect = fusion_band
    R_defect = lambda_band + 0.5 * lambda_domination

    defects = np.array([H_defect, B_defect, S_defect, V_defect, R_defect])
    supports = 1.0 / (1.0 + defects)

    eps = 1e-8
    F_maat = -np.sum(np.log((eps + supports) / (eps + 1.0)))

    stability = min(
        supports[4],
        (supports[0] * supports[1] * supports[2] * supports[3]) ** 0.25,
    )

    return F_maat + 5.0 * (1.0 - stability)


def main():
    bounds = [
        (np.log(1e-5), np.log(0.2)),       # alpha
        (np.log(1e-7), np.log(1e-2)),      # mu
        (np.log(1e-150), np.log(1e-90)),   # lambda
    ]

    result = differential_evolution(
        structural_score,
        bounds,
        seed=44,
        maxiter=1000,
        popsize=25,
        polish=True,
        tol=1e-10,
        workers=1,
    )

    alpha_fit, mu_fit, lambda_fit = np.exp(result.x)

    print("\n=== MAAT Constants Structure Scan v4 ===")
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

    print("\n--- Derived proxies ---")
    print("chemistry alpha^2*mu:", alpha_fit**2 * mu_fit)
    print("fusion alpha^2/mu   :", alpha_fit**2 / mu_fit)


if __name__ == "__main__":
    main()
