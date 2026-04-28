import numpy as np
from scipy.optimize import differential_evolution

ALPHA_OBS = 1 / 137.035999084
MU_OBS = 1 / 1836.15267343
LAMBDA_OBS = 1e-122

EPS = 1e-30


def softplus(x):
    return np.log1p(np.exp(np.clip(x, -80, 80)))


def band_penalty_log(log_x, log_min, log_max, width=1.0):
    return (
        softplus((log_min - log_x) / width) ** 2
        + softplus((log_x - log_max) / width) ** 2
    )


def qED_running_alpha(alpha_uv, log_scale_ratio):
    """
    One-loop toy QED-like running:
    alpha_IR = alpha_UV / (1 + b * alpha_UV * log_scale_ratio)

    This is not a full Standard Model RG calculation.
    It is a controlled toy map from UV coupling to IR coupling.
    """
    b = 2.0 / (3.0 * np.pi)
    denom = 1.0 + b * alpha_uv * log_scale_ratio
    return alpha_uv / max(denom, EPS)


def derive_mu(y_e, qcd_proxy):
    """
    Toy electron/proton mass ratio:
    mu ~ y_e / qcd_proxy

    Here qcd_proxy represents how strongly QCD enhances the proton mass
    relative to the electroweak/Yukawa scale.
    """
    return y_e / max(qcd_proxy, EPS)


def maat_score_from_effective(alpha, mu, lam):
    log_alpha = np.log(alpha)
    log_mu = np.log(mu)
    log_lam = np.log(lam)

    atom_band = band_penalty_log(
        log_alpha,
        np.log(1e-3),
        np.log(0.08),
        0.8,
    )

    mu_band = band_penalty_log(
        log_mu,
        np.log(1e-5),
        np.log(5e-3),
        0.8,
    )

    chemistry = alpha**2 * mu
    chemistry_band = band_penalty_log(
        np.log(max(chemistry, EPS)),
        np.log(1e-8),
        np.log(1e-5),
        0.9,
    )

    fusion = alpha**2 / max(mu, EPS)
    fusion_band = band_penalty_log(
        np.log(max(fusion, EPS)),
        np.log(1e-2),
        np.log(3.0),
        0.9,
    )

    lambda_band = band_penalty_log(
        log_lam,
        np.log(1e-140),
        np.log(1e-110),
        5.0,
    )

    lambda_domination = softplus((log_lam - np.log(1e-118)) / 5.0) ** 2

    defects = np.array([
        atom_band,
        mu_band,
        chemistry_band,
        fusion_band,
        lambda_band + 0.5 * lambda_domination,
    ])

    supports = 1.0 / (1.0 + defects)

    eps = 1e-8
    f_maat = -np.sum(np.log((eps + supports) / (eps + 1.0)))

    stability = min(
        supports[4],
        (supports[0] * supports[1] * supports[2] * supports[3]) ** 0.25,
    )

    return f_maat + 5.0 * (1.0 - stability), supports, defects


def total_score(x):
    """
    Optimisation variables:
    x[0] = log(alpha_uv)
    x[1] = log(log_scale_ratio)
    x[2] = log(y_e)
    x[3] = log(qcd_proxy)
    x[4] = log(lambda)

    Observed constants are NOT used in the score.
    """
    log_alpha_uv, log_rg_span, log_y_e, log_qcd_proxy, log_lam = x

    alpha_uv = np.exp(log_alpha_uv)
    rg_span = np.exp(log_rg_span)
    y_e = np.exp(log_y_e)
    qcd_proxy = np.exp(log_qcd_proxy)
    lam = np.exp(log_lam)

    alpha_ir = qED_running_alpha(alpha_uv, rg_span)
    mu_eff = derive_mu(y_e, qcd_proxy)

    if not (0 < alpha_ir < 1 and 0 < mu_eff < 1 and 0 < lam < 1):
        return 1e9

    f_maat, supports, defects = maat_score_from_effective(alpha_ir, mu_eff, lam)

    # RG regularisation:
    # avoid absurdly huge UV coupling or completely arbitrary RG span.
    rg_naturalness = 0.02 * band_penalty_log(
        np.log(max(alpha_uv, EPS)),
        np.log(1e-3),
        np.log(0.2),
        1.2,
    )

    rg_span_penalty = 0.01 * band_penalty_log(
        np.log(max(rg_span, EPS)),
        np.log(1.0),
        np.log(1e6),
        2.0,
    )

    # Yukawa/QCD proxy regularisation:
    yukawa_naturalness = 0.02 * band_penalty_log(
        np.log(max(y_e, EPS)),
        np.log(1e-7),
        np.log(1e-2),
        1.5,
    )

    qcd_proxy_naturalness = 0.02 * band_penalty_log(
        np.log(max(qcd_proxy, EPS)),
        np.log(1.0),
        np.log(1e5),
        1.8,
    )

    return (
        f_maat
        + rg_naturalness
        + rg_span_penalty
        + yukawa_naturalness
        + qcd_proxy_naturalness
    )


def main():
    bounds = [
        (np.log(1e-4), np.log(0.3)),       # alpha_uv
        (np.log(1.0), np.log(1e8)),        # RG log scale ratio proxy
        (np.log(1e-8), np.log(1e-2)),      # y_e proxy
        (np.log(1.0), np.log(1e6)),        # QCD mass enhancement proxy
        (np.log(1e-150), np.log(1e-90)),   # lambda
    ]

    result = differential_evolution(
        total_score,
        bounds,
        seed=91,
        maxiter=1200,
        popsize=25,
        polish=True,
        tol=1e-10,
        workers=1,
    )

    log_alpha_uv, log_rg_span, log_y_e, log_qcd_proxy, log_lam = result.x

    alpha_uv = np.exp(log_alpha_uv)
    rg_span = np.exp(log_rg_span)
    y_e = np.exp(log_y_e)
    qcd_proxy = np.exp(log_qcd_proxy)
    lam = np.exp(log_lam)

    alpha_ir = qED_running_alpha(alpha_uv, rg_span)
    mu_eff = derive_mu(y_e, qcd_proxy)

    f_maat, supports, defects = maat_score_from_effective(alpha_ir, mu_eff, lam)

    print("\n=== MAAT Constants v9: RG-Flow + MAAT ===")
    print("Best total score:", result.fun)
    print("MAAT score only :", f_maat)

    print("\n--- UV / generative parameters ---")
    print("alpha_uv     :", alpha_uv)
    print("rg_span      :", rg_span)
    print("y_e_proxy    :", y_e)
    print("qcd_proxy    :", qcd_proxy)
    print("lambda       :", lam)

    print("\n--- Effective IR constants ---")
    print("alpha_ir:", alpha_ir)
    print("mu_eff  :", mu_eff)
    print("lambda  :", lam)

    print("\n--- Observed comparison only ---")
    print("alpha_obs :", ALPHA_OBS)
    print("mu_obs    :", MU_OBS)
    print("lambda_obs:", LAMBDA_OBS)

    print("\n--- Deviation factors fit / observed ---")
    print("alpha factor :", alpha_ir / ALPHA_OBS)
    print("mu factor    :", mu_eff / MU_OBS)
    print("lambda factor:", lam / LAMBDA_OBS)

    print("\n--- Log10 comparison ---")
    print("log10 alpha_ir :", np.log10(alpha_ir))
    print("log10 alpha_obs:", np.log10(ALPHA_OBS))
    print("log10 mu_eff   :", np.log10(mu_eff))
    print("log10 mu_obs   :", np.log10(MU_OBS))
    print("log10 lambda   :", np.log10(lam))
    print("log10 lambdaobs:", np.log10(LAMBDA_OBS))

    print("\n--- MAAT supports H,B,S,V,R ---")
    for name, value in zip(["H", "B", "S", "V", "R"], supports):
        print(f"{name}: {value:.6f}")

    print("\n--- Derived structural proxies ---")
    print("chemistry alpha^2*mu:", alpha_ir**2 * mu_eff)
    print("fusion alpha^2/mu   :", alpha_ir**2 / mu_eff)


if __name__ == "__main__":
    main()
