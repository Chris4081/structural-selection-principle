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

    atom_band = band_penalty_log(log_alpha, np.log(1e-3), np.log(0.08), width=0.8)
    mu_band = band_penalty_log(log_mu, np.log(1e-5), np.log(5e-3), width=0.8)

    chemistry = alpha**2 * mu
    chemistry_band = band_penalty_log(
        np.log(max(chemistry, EPS)),
        np.log(1e-8),
        np.log(1e-5),
        width=0.9,
    )

    fusion_proxy = alpha**2 / max(mu, EPS)
    fusion_band = band_penalty_log(
        np.log(max(fusion_proxy, EPS)),
        np.log(1e-2),
        np.log(3.0),
        width=0.9,
    )

    lambda_band = band_penalty_log(
        log_lam,
        np.log(1e-140),
        np.log(1e-110),
        width=5.0,
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

    return f_maat + 5.0 * (1.0 - stability)


def run_one(seed):
    bounds = [
        (np.log(1e-5), np.log(0.2)),
        (np.log(1e-7), np.log(1e-2)),
        (np.log(1e-150), np.log(1e-90)),
    ]

    result = differential_evolution(
        structural_score,
        bounds,
        seed=seed,
        maxiter=700,
        popsize=20,
        polish=True,
        tol=1e-9,
        workers=1,
    )

    alpha, mu, lam = np.exp(result.x)
    return {
        "seed": seed,
        "score": result.fun,
        "alpha": alpha,
        "mu": mu,
        "lambda": lam,
        "log10_alpha": np.log10(alpha),
        "log10_mu": np.log10(mu),
        "log10_lambda": np.log10(lam),
        "alpha_factor": alpha / ALPHA_OBS,
        "mu_factor": mu / MU_OBS,
        "lambda_factor": lam / LAMBDA_OBS,
        "chemistry": alpha**2 * mu,
        "fusion": alpha**2 / mu,
    }


def summarize(values, name):
    arr = np.array(values, dtype=float)
    print(f"{name:16s} mean={arr.mean(): .6g}  std={arr.std(): .6g}  min={arr.min(): .6g}  max={arr.max(): .6g}")


def main():
    seeds = list(range(1, 31))
    results = []

    print("\n=== MAAT Constants Structure Scan v5 Robustness ===")
    print(f"Running {len(seeds)} seeds...\n")

    for seed in seeds:
        r = run_one(seed)
        results.append(r)
        print(
            f"seed={seed:02d} score={r['score']:.6f} "
            f"alpha={r['alpha']:.5g} mu={r['mu']:.5g} lambda={r['lambda']:.3e} "
            f"log10=({r['log10_alpha']:.3f}, {r['log10_mu']:.3f}, {r['log10_lambda']:.3f})"
        )

    print("\n--- Summary fitted values ---")
    summarize([r["alpha"] for r in results], "alpha")
    summarize([r["mu"] for r in results], "mu")
    summarize([r["lambda"] for r in results], "lambda")

    print("\n--- Summary log10 values ---")
    summarize([r["log10_alpha"] for r in results], "log10 alpha")
    summarize([r["log10_mu"] for r in results], "log10 mu")
    summarize([r["log10_lambda"] for r in results], "log10 lambda")

    print("\n--- Summary deviation factors fit / observed ---")
    summarize([r["alpha_factor"] for r in results], "alpha factor")
    summarize([r["mu_factor"] for r in results], "mu factor")
    summarize([r["lambda_factor"] for r in results], "lambda factor")

    print("\n--- Summary structural proxies ---")
    summarize([r["chemistry"] for r in results], "chemistry")
    summarize([r["fusion"] for r in results], "fusion")

    best = min(results, key=lambda r: r["score"])
    print("\n--- Best seed ---")
    for k, v in best.items():
        print(f"{k}: {v}")

    print("\nObserved comparison only:")
    print("alpha_obs :", ALPHA_OBS, "log10:", np.log10(ALPHA_OBS))
    print("mu_obs    :", MU_OBS, "log10:", np.log10(MU_OBS))
    print("lambda_obs:", LAMBDA_OBS, "log10:", np.log10(LAMBDA_OBS))


if __name__ == "__main__":
    main()
