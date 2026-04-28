import numpy as np
import matplotlib.pyplot as plt

ALPHA_OBS = 1 / 137.035999084
MU_OBS = 1 / 1836.15267343

# From v5/v6 optimum
LAMBDA_FIXED = 2.384e-129

EPS = 1e-30


def softplus(x):
    return np.log1p(np.exp(np.clip(x, -80, 80)))


def band_penalty_log(log_x, log_min, log_max, width=1.0):
    return (
        softplus((log_min - log_x) / width) ** 2
        + softplus((log_x - log_max) / width) ** 2
    )


def structural_score_alpha_mu(alpha, mu, lam=LAMBDA_FIXED):
    log_alpha = np.log(alpha)
    log_mu = np.log(mu)
    log_lam = np.log(lam)

    atom_band = band_penalty_log(
        log_alpha,
        np.log(1e-3),
        np.log(0.08),
        width=0.8,
    )

    mu_band = band_penalty_log(
        log_mu,
        np.log(1e-5),
        np.log(5e-3),
        width=0.8,
    )

    chemistry = alpha**2 * mu
    chemistry_band = band_penalty_log(
        np.log(max(chemistry, EPS)),
        np.log(1e-8),
        np.log(1e-5),
        width=0.9,
    )

    fusion = alpha**2 / max(mu, EPS)
    fusion_band = band_penalty_log(
        np.log(max(fusion, EPS)),
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


def main():
    # Scan area around physically interesting region
    log_alpha_vals = np.linspace(-4.0, -0.8, 320)
    log_mu_vals = np.linspace(-5.5, -2.0, 320)

    Z = np.zeros((len(log_mu_vals), len(log_alpha_vals)))

    best_score = float("inf")
    best_alpha = None
    best_mu = None

    for i, log_mu in enumerate(log_mu_vals):
        mu = 10 ** log_mu
        for j, log_alpha in enumerate(log_alpha_vals):
            alpha = 10 ** log_alpha
            score = structural_score_alpha_mu(alpha, mu)
            Z[i, j] = score

            if score < best_score:
                best_score = score
                best_alpha = alpha
                best_mu = mu

    print("\n=== MAAT Constants v7 Landscape Heatmap ===")
    print("Fixed lambda:", LAMBDA_FIXED)
    print("Best grid score:", best_score)
    print("Best alpha:", best_alpha, "observed:", ALPHA_OBS, "factor:", best_alpha / ALPHA_OBS)
    print("Best mu   :", best_mu, "observed:", MU_OBS, "factor:", best_mu / MU_OBS)
    print("Best log10 alpha:", np.log10(best_alpha))
    print("Best log10 mu   :", np.log10(best_mu))

    # Clip for visibility
    Z_plot = np.clip(Z, 0, 2.0)

    plt.figure(figsize=(12, 7))
    extent = [
        log_alpha_vals.min(),
        log_alpha_vals.max(),
        log_mu_vals.min(),
        log_mu_vals.max(),
    ]

    im = plt.imshow(
        Z_plot,
        origin="lower",
        aspect="auto",
        extent=extent,
        interpolation="nearest",
    )

    plt.colorbar(im, label="MAAT structural score F")

    # Observed values
    plt.scatter(
        [np.log10(ALPHA_OBS)],
        [np.log10(MU_OBS)],
        marker="*",
        s=180,
        label="Observed constants",
        edgecolors="black",
    )

    # Best grid point
    plt.scatter(
        [np.log10(best_alpha)],
        [np.log10(best_mu)],
        marker="o",
        s=100,
        label="MAAT optimum",
        edgecolors="black",
    )

    # v5/v6 optimum reference
    plt.scatter(
        [np.log10(0.0122007)],
        [np.log10(0.0006831)],
        marker="x",
        s=120,
        label="v5/v6 optimum",
    )

    plt.xlabel("log10(alpha)")
    plt.ylabel("log10(mu = m_e / m_p)")
    plt.title("MAAT Constants Landscape: Structural Score over alpha and mu")
    plt.legend()
    plt.tight_layout()

    out = "/Volumes/MAATSSD/TOE/maat_constants_v7_heatmap.png"
    plt.savefig(out, dpi=220)
    print("\nSaved heatmap:", out)


if __name__ == "__main__":
    main()
