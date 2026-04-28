# naturkonstante_v6_ablation_scan

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


def make_score(active):
    def structural_score(x):
        log_alpha, log_mu, log_lam = x
        alpha = np.exp(log_alpha)
        mu = np.exp(log_mu)

        atom_band = band_penalty_log(log_alpha, np.log(1e-3), np.log(0.08), 0.8)
        mu_band = band_penalty_log(log_mu, np.log(1e-5), np.log(5e-3), 0.8)

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
            atom_band if active.get("atom", True) else 0.0,
            mu_band if active.get("mass", True) else 0.0,
            chemistry_band if active.get("chemistry", True) else 0.0,
            fusion_band if active.get("fusion", True) else 0.0,
            (lambda_band + 0.5 * lambda_domination) if active.get("lam", True) else 0.0,
        ])

        supports = 1.0 / (1.0 + defects)

        eps = 1e-8
        f_maat = -np.sum(np.log((eps + supports) / (eps + 1.0)))

        stability = min(
            supports[4],
            (supports[0] * supports[1] * supports[2] * supports[3]) ** 0.25,
        )

        penalty = 5.0 * (1.0 - stability) if active.get("stability", True) else 0.0
        return f_maat + penalty

    return structural_score


def constrained_flags(active):
    alpha_constrained = active.get("atom", False) or active.get("chemistry", False) or active.get("fusion", False)
    mu_constrained = active.get("mass", False) or active.get("chemistry", False) or active.get("fusion", False)
    lambda_constrained = active.get("lam", False)
    return alpha_constrained, mu_constrained, lambda_constrained


def run_case(name, active, seed=42):
    bounds = [
        (np.log(1e-5), np.log(0.2)),
        (np.log(1e-7), np.log(1e-2)),
        (np.log(1e-150), np.log(1e-90)),
    ]

    result = differential_evolution(
        make_score(active),
        bounds,
        seed=seed,
        maxiter=700,
        popsize=20,
        polish=True,
        tol=1e-9,
        workers=1,
    )

    alpha, mu, lam = np.exp(result.x)
    alpha_c, mu_c, lam_c = constrained_flags(active)

    return {
        "case": name,
        "active": active,
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
        "alpha_constrained": alpha_c,
        "mu_constrained": mu_c,
        "lambda_constrained": lam_c,
    }


def fmt_value(value, constrained, fmt=".4g"):
    if not constrained:
        return "not constrained"
    return format(value, fmt)


def fmt_factor(value, constrained):
    if not constrained:
        return "n/a"
    return f"{value:.2g}"


def print_row(r):
    alpha_s = fmt_value(r["alpha"], r["alpha_constrained"], ".4g")
    mu_s = fmt_value(r["mu"], r["mu_constrained"], ".4g")
    lam_s = fmt_value(r["lambda"], r["lambda_constrained"], ".2e")

    log_alpha_s = fmt_value(r["log10_alpha"], r["alpha_constrained"], ".2f")
    log_mu_s = fmt_value(r["log10_mu"], r["mu_constrained"], ".2f")
    log_lam_s = fmt_value(r["log10_lambda"], r["lambda_constrained"], ".2f")

    fa = fmt_factor(r["alpha_factor"], r["alpha_constrained"])
    fm = fmt_factor(r["mu_factor"], r["mu_constrained"])
    fl = fmt_factor(r["lambda_factor"], r["lambda_constrained"])

    print(f"\n{r['case']}")
    print(f"  score   : {r['score']:.6g}")
    print(f"  alpha   : {alpha_s} | log10={log_alpha_s} | factor={fa}")
    print(f"  mu      : {mu_s} | log10={log_mu_s} | factor={fm}")
    print(f"  lambda  : {lam_s} | log10={log_lam_s} | factor={fl}")

    if r["alpha_constrained"] and r["mu_constrained"]:
        print(f"  chemistry alpha^2*mu: {r['chemistry']:.4e}")
        print(f"  fusion alpha^2/mu   : {r['fusion']:.4e}")


def main():
    cases = [
        ("full_model", dict(atom=True, mass=True, chemistry=True, fusion=True, lam=True, stability=True)),
        ("no_atom", dict(atom=False, mass=True, chemistry=True, fusion=True, lam=True, stability=True)),
        ("no_mass", dict(atom=True, mass=False, chemistry=True, fusion=True, lam=True, stability=True)),
        ("no_chemistry", dict(atom=True, mass=True, chemistry=False, fusion=True, lam=True, stability=True)),
        ("no_fusion", dict(atom=True, mass=True, chemistry=True, fusion=False, lam=True, stability=True)),
        ("no_lambda", dict(atom=True, mass=True, chemistry=True, fusion=True, lam=False, stability=True)),
        ("no_stability_penalty", dict(atom=True, mass=True, chemistry=True, fusion=True, lam=True, stability=False)),
        ("atom_mass_only", dict(atom=True, mass=True, chemistry=False, fusion=False, lam=False, stability=True)),
        ("chem_fusion_only", dict(atom=False, mass=False, chemistry=True, fusion=True, lam=False, stability=True)),
        ("lambda_only", dict(atom=False, mass=False, chemistry=False, fusion=False, lam=True, stability=True)),
    ]

    print("\n=== MAAT Constants v6b Ablation Scan ===")
    print("Note: parameters are marked 'not constrained' when their controlling sector is disabled.\n")

    for name, active in cases:
        r = run_case(name, active, seed=42)
        print_row(r)

    print("\nObserved comparison only:")
    print("alpha_obs :", ALPHA_OBS, "log10:", np.log10(ALPHA_OBS))
    print("mu_obs    :", MU_OBS, "log10:", np.log10(MU_OBS))
    print("lambda_obs:", LAMBDA_OBS, "log10:", np.log10(LAMBDA_OBS))

    print("\nInterpretation guide:")
    print("- Atom, chemistry and fusion constrain alpha.")
    print("- Mass, chemistry and fusion constrain mu.")
    print("- Lambda sector alone constrains lambda.")
    print("- When a parameter is marked 'not constrained', ignore its numerical optimizer value.")


if __name__ == "__main__":
    main()