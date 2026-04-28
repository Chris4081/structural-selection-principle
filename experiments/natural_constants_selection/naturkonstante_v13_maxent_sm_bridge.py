import json
import os
import numpy as np
from scipy.optimize import differential_evolution

ALPHA_OBS = 1 / 137.035999084
MU_OBS = 1 / 1836.15267343
LAMBDA_OBS = 1e-122

EPS = 1e-30

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
LAMBDA_JSON = os.path.join(BASE_DIR, "naturkonstante_v12_maxent_lambda_results.json")

DEFAULT_LAMBDAS = {
    "H": 1.61087307,
    "B": 2.42829723,
    "S": 4.62365496,
    "V": 4.64491214,
    "R": 8.83234363,
}


def load_lambdas():
    if os.path.exists(LAMBDA_JSON):
        with open(LAMBDA_JSON, "r") as f:
            data = json.load(f)
        return dict(zip(data["sectors"], data["lambdas"]))
    return DEFAULT_LAMBDAS.copy()


LAMBDA_SECTOR = load_lambdas()
LAMBDA_SUM = sum(LAMBDA_SECTOR.values())


def softplus(x):
    return np.log1p(np.exp(np.clip(x, -80, 80)))


def band_penalty_log(log_x, log_min, log_max, width=1.0):
    return (
        softplus((log_min - log_x) / width) ** 2
        + softplus((log_x - log_max) / width) ** 2
    )


def qed_running_alpha(alpha_uv, log_scale_ratio):
    b = 2.0 / (3.0 * np.pi)
    denom = 1.0 + b * alpha_uv * log_scale_ratio
    return alpha_uv / max(denom, EPS)


def derive_mu(y_e, qcd_proxy):
    return y_e / max(qcd_proxy, EPS)


def maat_score_from_effective(alpha, mu, lam):
    log_alpha = np.log(alpha)
    log_mu = np.log(mu)
    log_lam = np.log(lam)

    d_H = band_penalty_log(
        log_alpha,
        np.log(1e-3),
        np.log(0.08),
        0.8,
    )

    d_B = band_penalty_log(
        log_mu,
        np.log(1e-5),
        np.log(5e-3),
        0.8,
    )

    chemistry = alpha**2 * mu
    d_S = band_penalty_log(
        np.log(max(chemistry, EPS)),
        np.log(1e-8),
        np.log(1e-5),
        0.9,
    )

    fusion = alpha**2 / max(mu, EPS)
    d_V = band_penalty_log(
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
    d_R = lambda_band + 0.5 * lambda_domination

    defects = np.array([d_H, d_B, d_S, d_V, d_R])
    supports = 1.0 / (1.0 + defects)

    H, B, S, V, R = supports

    eps = 1e-8

    f_weighted = -(
        LAMBDA_SECTOR["H"] * np.log((eps + H) / (eps + 1.0))
        + LAMBDA_SECTOR["B"] * np.log((eps + B) / (eps + 1.0))
        + LAMBDA_SECTOR["S"] * np.log((eps + S) / (eps + 1.0))
        + LAMBDA_SECTOR["V"] * np.log((eps + V) / (eps + 1.0))
        + LAMBDA_SECTOR["R"] * np.log((eps + R) / (eps + 1.0))
    )

    f_maat = f_weighted / LAMBDA_SUM

    stability = min(
        R,
        (H * B * S * V) ** 0.25,
    )

    total = f_maat + 5.0 * (1.0 - stability)

    return total, f_maat, stability, supports, defects


def total_score(x):
    log_alpha_uv, log_rg_span, log_y_e, log_qcd_proxy, log_lam = x

    alpha_uv = np.exp(log_alpha_uv)
    rg_span = np.exp(log_rg_span)
    y_e = np.exp(log_y_e)
    qcd_proxy = np.exp(log_qcd_proxy)
    lam = np.exp(log_lam)

    alpha_ir = qed_running_alpha(alpha_uv, rg_span)
    mu_eff = derive_mu(y_e, qcd_proxy)

    if not (0 < alpha_ir < 1 and 0 < mu_eff < 1 and 0 < lam < 1):
        return 1e9

    score, _, _, _, _ = maat_score_from_effective(alpha_ir, mu_eff, lam)

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
        score
        + rg_naturalness
        + rg_span_penalty
        + yukawa_naturalness
        + qcd_proxy_naturalness
    )


def main():
    bounds = [
        (np.log(1e-4), np.log(0.3)),
        (np.log(1.0), np.log(1e8)),
        (np.log(1e-8), np.log(1e-2)),
        (np.log(1.0), np.log(1e6)),
        (np.log(1e-150), np.log(1e-90)),
    ]

    result = differential_evolution(
        total_score,
        bounds,
        seed=1313,
        maxiter=1400,
        popsize=28,
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

    alpha_ir = qed_running_alpha(alpha_uv, rg_span)
    mu_eff = derive_mu(y_e, qcd_proxy)

    total, f_maat, stability, supports, defects = maat_score_from_effective(
        alpha_ir, mu_eff, lam
    )

    print("\n=== MAAT Constants v13: MaxEnt-lambda SM Bridge ===")

    print("\n--- Loaded MaxEnt lambdas ---")
    for k in ["H", "B", "S", "V", "R"]:
        print(f"lambda_{k}: {LAMBDA_SECTOR[k]:.8f}")
    print("lambda_sum:", LAMBDA_SUM)

    print("\nBest total score:", result.fun)
    print("MAAT score only :", f_maat)
    print("Stability       :", stability)

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

    print("\n--- Defects d_H,d_B,d_S,d_V,d_R ---")
    for name, value in zip(["H", "B", "S", "V", "R"], defects):
        print(f"d_{name}: {value:.8f}")

    print("\n--- Derived structural proxies ---")
    print("chemistry alpha^2*mu:", alpha_ir**2 * mu_eff)
    print("fusion alpha^2/mu   :", alpha_ir**2 / mu_eff)

    output = {
        "lambdas": LAMBDA_SECTOR,
        "lambda_sum": LAMBDA_SUM,
        "best_total_score": float(result.fun),
        "maat_score_only": float(f_maat),
        "stability": float(stability),
        "alpha_uv": float(alpha_uv),
        "rg_span": float(rg_span),
        "y_e_proxy": float(y_e),
        "qcd_proxy": float(qcd_proxy),
        "lambda": float(lam),
        "alpha_ir": float(alpha_ir),
        "mu_eff": float(mu_eff),
        "alpha_obs": float(ALPHA_OBS),
        "mu_obs": float(MU_OBS),
        "lambda_obs": float(LAMBDA_OBS),
        "alpha_factor": float(alpha_ir / ALPHA_OBS),
        "mu_factor": float(mu_eff / MU_OBS),
        "lambda_factor": float(lam / LAMBDA_OBS),
        "supports": {
            k: float(v) for k, v in zip(["H", "B", "S", "V", "R"], supports)
        },
        "defects": {
            k: float(v) for k, v in zip(["H", "B", "S", "V", "R"], defects)
        },
        "chemistry": float(alpha_ir**2 * mu_eff),
        "fusion": float(alpha_ir**2 / mu_eff),
    }

    outpath = os.path.join(BASE_DIR, "naturkonstante_v13_maxent_sm_bridge_results.json")
    with open(outpath, "w") as f:
        json.dump(output, f, indent=2)

    print("\nSaved:", outpath)


if __name__ == "__main__":
    main()
