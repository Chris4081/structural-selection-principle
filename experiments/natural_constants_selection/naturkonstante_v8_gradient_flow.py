import numpy as np
import matplotlib.pyplot as plt

ALPHA_OBS = 1 / 137.035999084
MU_OBS = 1 / 1836.15267343

LAMBDA_FIXED = 2.384e-129
EPS = 1e-30


def softplus(x):
    return np.log1p(np.exp(np.clip(x, -80, 80)))


def band_penalty_log(log_x, log_min, log_max, width=1.0):
    return (
        softplus((log_min - log_x) / width) ** 2
        + softplus((log_x - log_max) / width) ** 2
    )


def structural_score(log_alpha10, log_mu10):
    alpha = 10 ** log_alpha10
    mu = 10 ** log_mu10
    lam = LAMBDA_FIXED

    log_alpha = np.log(alpha)
    log_mu = np.log(mu)
    log_lam = np.log(lam)

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

    lambda_band = band_penalty_log(log_lam, np.log(1e-140), np.log(1e-110), 5.0)
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


def numerical_gradient(x, y, h=1e-4):
    dfdx = (structural_score(x + h, y) - structural_score(x - h, y)) / (2 * h)
    dfdy = (structural_score(x, y + h) - structural_score(x, y - h)) / (2 * h)
    return np.array([dfdx, dfdy])


def gradient_flow(start, lr=0.05, steps=350):
    path = [np.array(start, dtype=float)]

    for _ in range(steps):
        p = path[-1]
        grad = numerical_gradient(p[0], p[1])

        # gradient clipping for stability
        norm = np.linalg.norm(grad)
        if norm > 5:
            grad = grad / norm * 5

        new_p = p - lr * grad

        # keep inside plotting domain
        new_p[0] = np.clip(new_p[0], -4.0, -0.8)
        new_p[1] = np.clip(new_p[1], -5.5, -2.0)

        path.append(new_p)

        if np.linalg.norm(path[-1] - path[-2]) < 1e-6:
            break

    return np.array(path)


def main():
    log_alpha_vals = np.linspace(-4.0, -0.8, 300)
    log_mu_vals = np.linspace(-5.5, -2.0, 300)

    Z = np.zeros((len(log_mu_vals), len(log_alpha_vals)))

    best_score = float("inf")
    best = None

    for i, y in enumerate(log_mu_vals):
        for j, x in enumerate(log_alpha_vals):
            s = structural_score(x, y)
            Z[i, j] = s
            if s < best_score:
                best_score = s
                best = (x, y)

    starts = [
        (-3.7, -5.1),
        (-3.4, -2.4),
        (-2.8, -4.8),
        (-2.4, -2.2),
        (-1.3, -4.6),
        (-1.1, -2.3),
        (-3.8, -3.0),
        (-1.6, -5.0),
    ]

    paths = [gradient_flow(s, lr=0.06, steps=500) for s in starts]

    print("\n=== MAAT Constants v8 Gradient Flow ===")
    print("Fixed lambda:", LAMBDA_FIXED)
    print("Best grid score:", best_score)
    print("Best log10 alpha:", best[0])
    print("Best log10 mu   :", best[1])
    print("Best alpha:", 10 ** best[0])
    print("Best mu   :", 10 ** best[1])
    print("Observed log10 alpha:", np.log10(ALPHA_OBS))
    print("Observed log10 mu   :", np.log10(MU_OBS))

    for idx, p in enumerate(paths, start=1):
        end = p[-1]
        print(
            f"path {idx}: start=({p[0,0]:.2f},{p[0,1]:.2f}) "
            f"end=({end[0]:.4f},{end[1]:.4f}) "
            f"score_end={structural_score(end[0], end[1]):.6f}"
        )

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

    # Contours
    cs = plt.contour(
        log_alpha_vals,
        log_mu_vals,
        Z_plot,
        levels=[0.1, 0.2, 0.4, 0.8, 1.2, 1.6],
        colors="black",
        linewidths=0.6,
        alpha=0.65,
    )
    plt.clabel(cs, inline=True, fontsize=8)

    # Gradient paths
    for p in paths:
        plt.plot(p[:, 0], p[:, 1], linewidth=2)
        plt.scatter(p[0, 0], p[0, 1], marker="s", s=40)
        plt.scatter(p[-1, 0], p[-1, 1], marker="o", s=55)

    # Observed
    plt.scatter(
        [np.log10(ALPHA_OBS)],
        [np.log10(MU_OBS)],
        marker="*",
        s=220,
        edgecolors="black",
        label="Observed constants",
    )

    # Best grid
    plt.scatter(
        [best[0]],
        [best[1]],
        marker="X",
        s=160,
        edgecolors="black",
        label="MAAT basin minimum",
    )

    plt.xlabel("log10(alpha)")
    plt.ylabel("log10(mu = m_e / m_p)")
    plt.title("MAAT Constants Landscape: Gradient Flow into Structural Basin")
    plt.legend()
    plt.tight_layout()

    out = "/Volumes/MAATSSD/TOE/maat_constants_v8_gradient_flow.png"
    plt.savefig(out, dpi=220)
    print("\nSaved gradient-flow figure:", out)


if __name__ == "__main__":
    main()
