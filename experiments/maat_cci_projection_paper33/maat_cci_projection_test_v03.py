#!/usr/bin/env python3
# MAAT CCI Projection Test v0.3
# Realistic cosmology proxy version:
# CCI as projection observable of latent structural depth under expansion.
#
# Improvements over v0.2:
# - Uses LCDM E(z)
# - Uses Carroll-Press-Turner growth proxy D(z)
# - Compares projection transition to matter-Lambda equality
# - Generates paper-ready plots and summary files

from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUTDIR = Path("maat_cci_projection_test_v03")
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
    """
    Carroll-Press-Turner growth suppression approximation.
    """
    om = omega_m_z(z)
    ol = omega_l_z(z)
    return (5.0 * om / 2.0) / (
        om ** (4.0 / 7.0)
        - ol
        + (1.0 + om / 2.0) * (1.0 + ol / 70.0)
        + EPS
    )


def growth_D(z):
    """
    Normalized linear growth proxy D(z), D(0)=1.
    """
    g = growth_suppression_g(z)
    g0 = growth_suppression_g(0.0)
    return (g / (1.0 + z)) / (g0 + EPS)


def latent_depth(z, zc=1.1, sharpness=3.5, floor=0.20):
    """
    Hidden structural depth availability.
    Smoothly decreases as expansion/redshift increases.
    """
    raw = 1.0 / (1.0 + np.exp(sharpness * (z - zc)))
    return floor + (1.0 - floor) * raw


def saturation(x):
    return np.tanh(x)


def compute_projection_cci(
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

    # Breadth = expansion-driven configuration breadth
    breadth = E * (1.0 + z)

    # Depth = physically accessible coherent structure proxy
    accessible_depth = D * L

    # Compression = expansion-induced projection pressure, bounded
    raw_compression = (breadth / breadth[0]) ** alpha
    compression = saturation(raw_compression / np.nanmax(raw_compression))

    # Stabilized projection observable
    projected_surface = breadth * compression
    cci_proj = projected_surface / (1.0 + gamma * accessible_depth)
    cci_norm = cci_proj / max(cci_proj[0], EPS)

    # Bounded surface imbalance
    residual = np.abs(projected_surface - accessible_depth) / (
        projected_surface + accessible_depth + EPS
    )

    depth_loss = 1.0 - accessible_depth / max(accessible_depth[0], EPS)

    return pd.DataFrame({
        "z": z,
        "E_LCDM": E,
        "Omega_m_z": omega_m_z(z),
        "Omega_L_z": omega_l_z(z),
        "growth_D": D,
        "latent_depth": L,
        "breadth": breadth,
        "accessible_depth": accessible_depth,
        "compression": compression,
        "projected_surface": projected_surface,
        "CCI_projection": cci_proj,
        "CCI_projection_norm": cci_norm,
        "projection_residual": residual,
        "depth_loss": depth_loss,
    })


def matter_lambda_equality_z():
    """
    Omega_m (1+z)^3 = Omega_L
    """
    return (OMEGA_L0 / OMEGA_M0) ** (1.0 / 3.0) - 1.0


def transition_estimate(df):
    """
    Estimate transition from maximum curvature of bounded residual.
    Ignore edges to avoid boundary artifacts.
    """
    z = df["z"].to_numpy()
    y = df["projection_residual"].to_numpy()

    dy = np.gradient(y, z)
    ddy = np.gradient(dy, z)

    n = len(z)
    lo = int(0.08 * n)
    hi = int(0.92 * n)

    idx = lo + int(np.nanargmax(np.abs(ddy[lo:hi])))

    return {
        "z_transition_estimate": float(z[idx]),
        "residual_at_transition": float(y[idx]),
        "CCI_norm_at_transition": float(df["CCI_projection_norm"].iloc[idx]),
        "max_abs_curvature": float(abs(ddy[idx])),
    }


def save_plots(df, z_eq):
    # 1 breadth-depth
    plt.figure(figsize=(8.8, 5.4))
    plt.plot(df["z"], df["breadth"] / df["breadth"].iloc[0], label="Expansion breadth")
    plt.plot(df["z"], df["accessible_depth"] / df["accessible_depth"].iloc[0], label="Accessible structural depth")
    plt.axvline(z_eq, linestyle="--", label="matter–Λ equality")
    plt.xlabel("z")
    plt.ylabel("normalized value")
    plt.title("Breadth–Depth Projection Structure with LCDM Growth Proxy")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig1_breadth_depth_lcdm_growth.png", dpi=260)
    plt.close()

    # 2 CCI
    plt.figure(figsize=(8.8, 5.4))
    plt.plot(df["z"], df["CCI_projection_norm"])
    plt.axvline(z_eq, linestyle="--", label="matter–Λ equality")
    plt.xlabel("z")
    plt.ylabel("CCI projection, normalized at z=0")
    plt.title("CCI as Cosmological Projection Observable")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig2_cci_projection_lcdm_growth.png", dpi=260)
    plt.close()

    # 3 residual
    plt.figure(figsize=(8.8, 5.4))
    plt.plot(df["z"], df["projection_residual"])
    plt.axvline(z_eq, linestyle="--", label="matter–Λ equality")
    plt.xlabel("z")
    plt.ylabel("bounded projection residual")
    plt.title("Surface Imbalance Under Cosmological Expansion")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig3_projection_residual_lcdm_growth.png", dpi=260)
    plt.close()

    # 4 components
    plt.figure(figsize=(8.8, 5.4))
    plt.plot(df["z"], df["growth_D"], label="Growth proxy D(z)")
    plt.plot(df["z"], df["latent_depth"], label="Latent structural depth")
    plt.plot(df["z"], df["compression"], label="Saturating compression")
    plt.plot(df["z"], df["depth_loss"], label="Accessible depth loss")
    plt.axvline(z_eq, linestyle="--", label="matter–Λ equality")
    plt.xlabel("z")
    plt.ylabel("bounded proxy")
    plt.title("Growth, Latent Depth, Compression, and Depth Loss")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig4_components_lcdm_growth.png", dpi=260)
    plt.close()

    # 5 Omega components
    plt.figure(figsize=(8.8, 5.4))
    plt.plot(df["z"], df["Omega_m_z"], label="Ωm(z)")
    plt.plot(df["z"], df["Omega_L_z"], label="ΩΛ(z)")
    plt.axvline(z_eq, linestyle="--", label="matter–Λ equality")
    plt.xlabel("z")
    plt.ylabel("density fraction")
    plt.title("Reference LCDM Background Fractions")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig5_lcdm_density_fractions.png", dpi=260)
    plt.close()


def alpha_sweep():
    z = np.linspace(0.0, 3.0, 900)
    rows = []

    for alpha in np.linspace(0.0, 2.0, 31):
        df = compute_projection_cci(z, alpha=alpha)
        trans = transition_estimate(df)

        rows.append({
            "alpha": float(alpha),
            **trans,
            "final_CCI_norm": float(df["CCI_projection_norm"].iloc[-1]),
            "max_projection_residual": float(df["projection_residual"].max()),
        })

    sweep = pd.DataFrame(rows)
    sweep.to_csv(OUTDIR / "alpha_sweep_results_v03.csv", index=False)

    plt.figure(figsize=(8.8, 5.4))
    plt.plot(sweep["alpha"], sweep["z_transition_estimate"], marker="o")
    plt.axhline(matter_lambda_equality_z(), linestyle="--", label="matter–Λ equality")
    plt.xlabel("projection strength alpha")
    plt.ylabel("estimated transition redshift")
    plt.title("Projection Strength vs Transition Estimate")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig6_alpha_transition_sweep_v03.png", dpi=260)
    plt.close()

    plt.figure(figsize=(8.8, 5.4))
    plt.plot(sweep["alpha"], sweep["final_CCI_norm"], marker="o")
    plt.xlabel("projection strength alpha")
    plt.ylabel("final normalized CCI")
    plt.title("Projection Strength vs Final CCI")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig7_alpha_final_cci_v03.png", dpi=260)
    plt.close()

    return sweep


def main():
    z = np.linspace(0.0, 3.0, 1000)

    df = compute_projection_cci(
        z,
        alpha=1.0,
        latent_zc=1.1,
        sharpness=3.5,
        depth_floor=0.20,
        gamma=1.0,
    )

    z_eq = matter_lambda_equality_z()
    trans = transition_estimate(df)
    sweep = alpha_sweep()

    df.to_csv(OUTDIR / "cci_projection_timeseries_v03.csv", index=False)
    save_plots(df, z_eq)

    summary = {
        "model": "MAAT CCI Projection Test v0.3",
        "interpretation": "CCI as stabilized projection observable using LCDM expansion and growth proxy.",
        "redshift_range": [float(z.min()), float(z.max())],
        "H0": H0,
        "Omega_m0": OMEGA_M0,
        "Omega_L0": OMEGA_L0,
        "matter_lambda_equality_z": float(z_eq),
        "baseline_projection_alpha": 1.0,
        "depth_floor": 0.20,
        "compression": "saturating tanh of expansion breadth",
        "transition_estimate": trans,
        "delta_z_transition_minus_matter_lambda_equality": float(trans["z_transition_estimate"] - z_eq),
        "final_CCI_norm": float(df["CCI_projection_norm"].iloc[-1]),
        "max_projection_residual": float(df["projection_residual"].max()),
        "note": "Uses LCDM expansion and an analytic growth proxy; not a full cosmological parameter fit."
    }

    with open(OUTDIR / "cci_projection_summary_v03.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\n=== MAAT CCI Projection Test v0.3 ===")
    for k, v in summary.items():
        print(f"{k}: {v}")

    print(f"\nSaved outputs to: {OUTDIR.resolve()}")
    print("Main files:")
    print(" - cci_projection_timeseries_v03.csv")
    print(" - alpha_sweep_results_v03.csv")
    print(" - cci_projection_summary_v03.json")
    print(" - fig1_breadth_depth_lcdm_growth.png")
    print(" - fig2_cci_projection_lcdm_growth.png")
    print(" - fig3_projection_residual_lcdm_growth.png")
    print(" - fig4_components_lcdm_growth.png")
    print(" - fig5_lcdm_density_fractions.png")
    print(" - fig6_alpha_transition_sweep_v03.png")
    print(" - fig7_alpha_final_cci_v03.png")


if __name__ == "__main__":
    main()