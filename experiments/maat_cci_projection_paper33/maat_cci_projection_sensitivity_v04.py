#!/usr/bin/env python3
# MAAT CCI Projection Sensitivity Test v0.4
# Tests whether the projection-transition redshift is robust under parameter changes.
#
# Builds on v0.3:
# - LCDM expansion E(z)
# - Carroll-Press-Turner growth proxy D(z)
# - stabilized projection observable
# - scans alpha, latent_zc, sharpness, depth_floor
# - outputs CSV, JSON summary, heatmaps and histograms

from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUTDIR = Path("maat_cci_projection_sensitivity_v04")
OUTDIR.mkdir(exist_ok=True)

# Planck-like LCDM reference
H0 = 67.4
OMEGA_M0 = 0.315
OMEGA_L0 = 0.685
EPS = 1e-9


def E_lcdm(z):
    return np.sqrt(OMEGA_M0 * (1.0 + z) ** 3 + OMEGA_L0)


def omega_m_z(z):
    Ez2 = E_lcdm(z) ** 2
    return OMEGA_M0 * (1.0 + z) ** 3 / Ez2


def omega_l_z(z):
    Ez2 = E_lcdm(z) ** 2
    return OMEGA_L0 / Ez2


def growth_suppression_g(z):
    om = omega_m_z(z)
    ol = omega_l_z(z)
    return (5.0 * om / 2.0) / (
        om ** (4.0 / 7.0)
        - ol
        + (1.0 + om / 2.0) * (1.0 + ol / 70.0)
        + EPS
    )


def growth_D(z):
    g = growth_suppression_g(z)
    g0 = growth_suppression_g(0.0)
    return (g / (1.0 + z)) / (g0 + EPS)


def matter_lambda_equality_z():
    return (OMEGA_L0 / OMEGA_M0) ** (1.0 / 3.0) - 1.0


def latent_depth(z, zc=1.1, sharpness=3.5, floor=0.20):
    raw = 1.0 / (1.0 + np.exp(sharpness * (z - zc)))
    return floor + (1.0 - floor) * raw


def compute_projection(
    z,
    alpha=1.0,
    latent_zc=1.1,
    sharpness=3.5,
    depth_floor=0.20,
    gamma=1.0,
):
    E = E_lcdm(z)
    D = growth_D(z)
    L = latent_depth(z, zc=latent_zc, sharpness=sharpness, floor=depth_floor)

    breadth = E * (1.0 + z)
    accessible_depth = D * L

    raw_compression = (breadth / breadth[0]) ** alpha
    compression = np.tanh(raw_compression / (np.nanmax(raw_compression) + EPS))

    projected_surface = breadth * compression
    cci = projected_surface / (1.0 + gamma * accessible_depth)
    cci_norm = cci / max(cci[0], EPS)

    residual = np.abs(projected_surface - accessible_depth) / (
        projected_surface + accessible_depth + EPS
    )

    return {
        "z": z,
        "breadth": breadth,
        "accessible_depth": accessible_depth,
        "cci_norm": cci_norm,
        "residual": residual,
        "growth_D": D,
        "latent_depth": L,
        "compression": compression,
    }


def transition_estimate(z, residual, cci_norm):
    """
    Transition = strongest curvature of bounded residual,
    with edge exclusion to avoid boundary artifacts.
    """
    dy = np.gradient(residual, z)
    ddy = np.gradient(dy, z)

    n = len(z)
    lo = int(0.08 * n)
    hi = int(0.92 * n)

    idx = lo + int(np.nanargmax(np.abs(ddy[lo:hi])))

    return {
        "z_transition": float(z[idx]),
        "residual_at_transition": float(residual[idx]),
        "cci_norm_at_transition": float(cci_norm[idx]),
        "curvature_strength": float(abs(ddy[idx])),
    }


def classify_transition(z_transition):
    if 0.5 <= z_transition <= 1.1:
        return "cosmologically_relevant_mid_z"
    if z_transition < 0.5:
        return "low_z"
    return "high_z"


def main():
    z = np.linspace(0.0, 3.0, 1000)
    z_eq = matter_lambda_equality_z()

    # Sensitivity grid
    alpha_values = np.linspace(0.2, 2.0, 10)
    latent_zc_values = np.linspace(0.6, 1.6, 11)
    sharpness_values = np.linspace(2.0, 6.0, 9)
    depth_floor_values = np.linspace(0.10, 0.40, 7)

    rows = []

    total = (
        len(alpha_values)
        * len(latent_zc_values)
        * len(sharpness_values)
        * len(depth_floor_values)
    )

    count = 0
    for alpha in alpha_values:
        for latent_zc in latent_zc_values:
            for sharpness in sharpness_values:
                for depth_floor in depth_floor_values:
                    count += 1

                    proj = compute_projection(
                        z,
                        alpha=alpha,
                        latent_zc=latent_zc,
                        sharpness=sharpness,
                        depth_floor=depth_floor,
                    )

                    trans = transition_estimate(
                        z,
                        proj["residual"],
                        proj["cci_norm"],
                    )

                    zt = trans["z_transition"]

                    rows.append({
                        "alpha": float(alpha),
                        "latent_zc": float(latent_zc),
                        "sharpness": float(sharpness),
                        "depth_floor": float(depth_floor),
                        **trans,
                        "delta_z_vs_matter_lambda": float(zt - z_eq),
                        "abs_delta_z_vs_matter_lambda": float(abs(zt - z_eq)),
                        "final_cci_norm": float(proj["cci_norm"][-1]),
                        "max_residual": float(np.nanmax(proj["residual"])),
                        "regime": classify_transition(zt),
                    })

    df = pd.DataFrame(rows)
    df.to_csv(OUTDIR / "cci_projection_sensitivity_results_v04.csv", index=False)

    # Summary stats
    mid = df[df["regime"] == "cosmologically_relevant_mid_z"]
    summary = {
        "model": "MAAT CCI Projection Sensitivity Test v0.4",
        "total_parameter_points": int(total),
        "matter_lambda_equality_z": float(z_eq),
        "transition_z_mean": float(df["z_transition"].mean()),
        "transition_z_median": float(df["z_transition"].median()),
        "transition_z_std": float(df["z_transition"].std()),
        "transition_z_min": float(df["z_transition"].min()),
        "transition_z_max": float(df["z_transition"].max()),
        "fraction_mid_z_0p5_to_1p1": float(len(mid) / len(df)),
        "best_match_to_matter_lambda": df.sort_values("abs_delta_z_vs_matter_lambda").head(5).to_dict(orient="records"),
        "note": (
            "Sensitivity scan of toy projection model. "
            "A robust mid-z concentration supports the projection-observable interpretation, "
            "but this is not a data fit."
        ),
    }

    with open(OUTDIR / "cci_projection_sensitivity_summary_v04.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Plot 1: histogram of transition redshifts
    plt.figure(figsize=(8.8, 5.4))
    plt.hist(df["z_transition"], bins=35, alpha=0.85)
    plt.axvline(z_eq, linestyle="--", label="matter–Λ equality")
    plt.axvspan(0.5, 1.1, alpha=0.15, label="mid-z relevant band")
    plt.xlabel("estimated transition redshift")
    plt.ylabel("count")
    plt.title("Distribution of CCI Projection Transition Redshifts")
    plt.legend()
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig1_transition_histogram_v04.png", dpi=260)
    plt.close()

    # Plot 2: alpha vs transition
    plt.figure(figsize=(8.8, 5.4))
    plt.scatter(df["alpha"], df["z_transition"], s=12, alpha=0.35)
    plt.axhline(z_eq, linestyle="--", label="matter–Λ equality")
    plt.axhspan(0.5, 1.1, alpha=0.15, label="mid-z relevant band")
    plt.xlabel("projection strength alpha")
    plt.ylabel("transition redshift")
    plt.title("Sensitivity: Projection Strength vs Transition")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig2_alpha_vs_transition_v04.png", dpi=260)
    plt.close()

    # Plot 3: latent_zc vs transition
    plt.figure(figsize=(8.8, 5.4))
    plt.scatter(df["latent_zc"], df["z_transition"], s=12, alpha=0.35)
    plt.axhline(z_eq, linestyle="--", label="matter–Λ equality")
    plt.axhspan(0.5, 1.1, alpha=0.15, label="mid-z relevant band")
    plt.xlabel("latent depth transition parameter zc")
    plt.ylabel("transition redshift")
    plt.title("Sensitivity: Latent Depth Scale vs Transition")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig3_latent_zc_vs_transition_v04.png", dpi=260)
    plt.close()

    # Plot 4: depth_floor vs transition
    plt.figure(figsize=(8.8, 5.4))
    plt.scatter(df["depth_floor"], df["z_transition"], s=12, alpha=0.35)
    plt.axhline(z_eq, linestyle="--", label="matter–Λ equality")
    plt.axhspan(0.5, 1.1, alpha=0.15, label="mid-z relevant band")
    plt.xlabel("depth floor")
    plt.ylabel("transition redshift")
    plt.title("Sensitivity: Depth Floor vs Transition")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig4_depth_floor_vs_transition_v04.png", dpi=260)
    plt.close()

    # Plot 5: heatmap alpha-latent_zc, averaged over sharpness/floor
    pivot = df.groupby(["latent_zc", "alpha"])["z_transition"].mean().reset_index()
    heat = pivot.pivot(index="latent_zc", columns="alpha", values="z_transition")

    plt.figure(figsize=(9.2, 5.8))
    im = plt.imshow(
        heat.values,
        origin="lower",
        aspect="auto",
        extent=[
            alpha_values.min(),
            alpha_values.max(),
            latent_zc_values.min(),
            latent_zc_values.max(),
        ],
    )
    plt.colorbar(im, label="mean transition redshift")
    plt.xlabel("projection strength alpha")
    plt.ylabel("latent depth zc")
    plt.title("Mean Transition Redshift Across Parameter Grid")
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig5_transition_heatmap_alpha_zc_v04.png", dpi=260)
    plt.close()

    # Plot 6: final CCI vs transition
    plt.figure(figsize=(8.8, 5.4))
    plt.scatter(df["z_transition"], df["final_cci_norm"], s=12, alpha=0.35)
    plt.yscale("log")
    plt.axvline(z_eq, linestyle="--", label="matter–Λ equality")
    plt.axvspan(0.5, 1.1, alpha=0.15, label="mid-z relevant band")
    plt.xlabel("transition redshift")
    plt.ylabel("final normalized CCI, log scale")
    plt.title("Final Projection Strength vs Transition Location")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig6_final_cci_vs_transition_v04.png", dpi=260)
    plt.close()

    print("\n=== MAAT CCI Projection Sensitivity Test v0.4 ===")
    for k, v in summary.items():
        if k == "best_match_to_matter_lambda":
            print(f"{k}:")
            for row in v:
                print("  ", row)
        else:
            print(f"{k}: {v}")

    print(f"\nSaved outputs to: {OUTDIR.resolve()}")
    print("Main files:")
    print(" - cci_projection_sensitivity_results_v04.csv")
    print(" - cci_projection_sensitivity_summary_v04.json")
    print(" - fig1_transition_histogram_v04.png")
    print(" - fig2_alpha_vs_transition_v04.png")
    print(" - fig3_latent_zc_vs_transition_v04.png")
    print(" - fig4_depth_floor_vs_transition_v04.png")
    print(" - fig5_transition_heatmap_alpha_zc_v04.png")
    print(" - fig6_final_cci_vs_transition_v04.png")


if __name__ == "__main__":
    main()