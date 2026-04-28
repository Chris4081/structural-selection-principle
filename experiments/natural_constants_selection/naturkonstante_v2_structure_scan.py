import numpy as np
from scipy.optimize import differential_evolution

# Observed values ONLY for comparison after optimisation
ALPHA_OBS = 1 / 137.035999084
MU_OBS = 1 / 1836.15267343
LAMBDA_OBS = 1e-122

EPS = 1e-12


def safe_log(x):
    return np.log(np.maximum(x, EPS))


def window_defect(log_x, log_min, log_max):
    """
    Zero inside [log_min, log_max], quadratic penalty outside.
    """
    if log_x < log_min:
        return ((log_min - log_x) / max(abs(log_max - log_min), EPS)) ** 2
    if log_x > log_max:
        return ((log_x - log_max) / max(abs(log_max - log_min), EPS)) ** 2
    return 0.0


def band_center_defect(log_x, log_center, width):
    """
    Soft preference around a structural band center.
    """
    return ((log_x - log_center) / width) ** 2


def structural_score(x):
    """
    x = [log_alpha, log_mu, log_lambda]

    No observed constants are used here.
    The score is based only on broad structural viability constraints.
    """
    log_alpha, log_mu, log_lam = x

    alpha = np.exp(log_alpha)
    mu = np.exp(log_mu)
    lam = np.exp(log_lam)

    # -----------------------------
    # 1. Atom stability constraint
    # -----------------------------
    # Alpha must be small enough for nonrelativistic atoms,
    # but not so tiny that electromagnetic binding becomes structurally inert.
    atom_window = window_defect(
        log_alpha,
        np.log(1e-4),
        np.log(0.1),
    )

    # Soft structural preference: alpha neither too weak nor too strong.
    atom_center = band_center_defect(
        log_alpha,
        np.log(1e-2),
        1.5,
    )

    # -----------------------------
    # 2. Mass hierarchy / chemistry
    # -----------------------------
    # mu = me/mp should be small: electrons much lighter than baryons.
    # This enables separated atomic/nuclear scales.
    hierarchy_window = window_defect(
        log_mu,
        np.log(1e-6),
        np.log(1e-2),
    )

    hierarchy_center = band_center_defect(
        log_mu,
        np.log(1e-4),
        2.0,
    )

    # -----------------------------
    # 3. Structure formation / Lambda
    # -----------------------------
    # Lambda should be tiny on Planck scale, otherwise large-scale structure
    # is suppressed. We do not target the observed value directly.
    lambda_window = window_defect(
        log_lam,
        np.log(1e-140),
        np.log(1e-90),
    )

    # Softly prefer "very small but nonzero" lambda.
    lambda_center = band_center_defect(
        log_lam,
        np.log(1e-120),
        12.0,
    )

    # -----------------------------
    # 4. Coupling compatibility
    # -----------------------------
    # Too large alpha with too large mu compresses atomic/nuclear separation.
    coupling_stress = max(0.0, np.log(alpha * mu / 1e-5)) ** 2

    # Too tiny alpha and mu together produces almost no chemistry/activity.
    inertness_stress = max(0.0, np.log(1e-10 / (alpha * mu))) ** 2

    # -----------------------------
    # MAAT-like sector supports
    # -----------------------------
    H_defect = atom_window + 0.35 * atom_center
    B_defect = hierarchy_window + 0.25 * hierarchy_center
    S_defect = 0.15 * inertness_stress
    V_defect = 0.10 * coupling_stress
    R_defect = lambda_window + 0.25 * lambda_center

    defects = np.array([H_defect, B_defect, S_defect, V_defect, R_defect])
    supports = 1 / (1 + defects)

    # Normalised log-sum MAAT functional
    eps = 1e-8
    F_maat = -np.sum(np.log((eps + supports) / (eps + 1)))

    stability = min(
        supports[4],
        (supports[0] * supports[1] * supports[2] * supports[3]) ** 0.25,
    )

    # Add hard-ish stability penalty
    return F_maat + 5.0 * (1.0 - stability)


def main():
    bounds = [
        (np.log(1e-6), np.log(0.3)),       # alpha
        (np.log(1e-8), np.log(1e-1)),      # mu = me/mp
        (np.log(1e-150), np.log(1e-60)),   # lambda
    ]

    result = differential_evolution(
        structural_score,
        bounds,
        seed=42,
        maxiter=600,
        popsize=20,
        polish=True,
        tol=1e-10,
        workers=1,
    )

    alpha_fit, mu_fit, lambda_fit = np.exp(result.x)

    print("\n=== MAAT Constants Structure Scan v2 ===")
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


if __name__ == "__main__":
    main()
