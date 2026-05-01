#!/usr/bin/env python3
# MAAT CCI Projection Test v0.2
# Stabilized toy model: CCI as projection observable of latent structural depth
# under cosmological expansion.
#
# Key improvement vs v0.1:
# - no division-by-zero blow-up
# - accessible depth has a floor
# - compression saturates
# - transition is estimated from curvature / maximum response, not boundary divergence

from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUTDIR = Path("maat_cci_projection_test_v02")
OUTDIR.mkdir(exist_ok=True)

H0 = 67.4
OMEGA_M = 0.315
OMEGA_L = 0.685
EPS = 1e-9


def E_lcdm(z):
    return np.sqrt(OMEGA_M * (1.0 + z) ** 3 + OMEGA_L)


def growth_depth_proxy(z):
    """
    Accessible structure-growth depth proxy.
    Not a Boltzmann solution.
    """
    return 1.0 / (1.0 + z)


def latent_depth(z, transition_z=1.1, sharpness=4.0, floor=0.25):
    """
    Latent structural depth with a nonzero floor.
    The floor prevents artificial singularity when hidden structure becomes
    weakly accessible.
    """
    raw = 1.0 / (1.0 + np.exp(sharpness * (z - transition_z)))
    return floor + (1.0 - floor) * raw


def compression_factor(z, alpha=1.0):
    """
    Saturating expansion-compression factor.
    Prevents runaway divergence while preserving monotonic projection pressure.
    """
    raw = (1.0 + z) ** alpha
    return np.tanh(raw / raw.max())


def compute_projection_cci(
    z,
    alpha=1.0,
    transition_z=1.1,
    sharpness=4.0,
    depth_floor=0.25,
    gamma=1.0,
):
    E = E_lcdm(z)
    D = growth_depth_proxy(z)
    L = latent_depth(
        z,
        transition_z=transition_z,
        sharpness=sharpness,
        floor=depth_floor,
    )

    breadth = E * (1.0 + z)
    accessible_depth = D * L

    compression = compression_factor(z, alpha=alpha)

    # Stabilized projection observable:
    # expansion breadth is projected through saturating compression and moderated
    # by accessible depth.
    projected_surface = breadth * compression
    CCI_proj = projected_surface / (1.0 + gamma * accessible_depth)

    # Normalize to z=0 for readability.
    CCI_norm = CCI_proj / max(CCI_proj[0], EPS)

    # Bounded imbalance observable in [0,1]
    projection_residual = np.abs(projected_surface - accessible_depth) / (
        projected_surface + accessible_depth + EPS
    )

    # Depth compression contrast
    depth_loss = 1.0 - accessible_depth / max(accessible_depth[0], EPS)

    return pd.DataFrame({
        "z": z,
        "E": E,
        "growth_depth_D": D,
        "latent_depth_L": L,
        "breadth": breadth,
        "accessible_depth": accessible_depth,
        "compression": compression,
        "projected_surface": projected_surface,
        "CCI_projection": CCI_proj,
        "CCI_projection_norm": CCI_norm,
        "projection_residual": projection_residual,
        "depth_loss": depth_loss,
    })


def transition_estimate(df):
    """
    Estimate transition from strongest curvature of bounded residual,
    ignoring edges.
    """
    z = df["z"].to_numpy()
    y = df["projection_residual"].to_numpy()

    dy = np.gradient(y, z)
    ddy = np.gradient(dy, z)

    # ignore outer 8% to avoid edge artifacts
    n = len(z)
    lo = int(0.08 * n)
    hi = int(0.92 * n)

    idx_local = np.argmax(np.abs(ddy[lo:hi]))
    idx = lo + idx_local

    return {
        "z_transition_estimate": float(z[idx]),
        "residual_at_transition": float(y[idx]),
        "max_abs_curvature": float(abs(ddy[idx])),
        "CCI_norm_at_transition": float(df["CCI_projection_norm"].iloc[idx]),
    }


def save_plots(df):
    plt.figure(figsize=(8.5, 5.3))
    plt.plot(df["z"], df["breadth"] / df["breadth"].iloc[0], label="Expansion breadth")
    plt.plot(df["z"], df["accessible_depth"] / df["accessible_depth"].iloc[0], label="Accessible structural depth")
    plt.xlabel("z")
    plt.ylabel("normalized value")
    plt.title("Breadth–Depth Projection Structure")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig1_breadth_depth_stabilized.png", dpi=260)
    plt.close()

    plt.figure(figsize=(8.5, 5.3))
    plt.plot(df["z"], df["CCI_projection_norm"])
    plt.xlabel("z")
    plt.ylabel("CCI projection, normalized at z=0")
    plt.title("CCI as Stabilized Projection Observable")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig2_stabilized_CCI_projection.png", dpi=260)
    plt.close()

    plt.figure(figsize=(8.5, 5.3))
    plt.plot(df["z"], df["projection_residual"])
    plt.xlabel("z")
    plt.ylabel("bounded projection residual")
    plt.title("Surface Imbalance Under Expansion")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig3_bounded_projection_residual.png", dpi=260)
    plt.close()

    plt.figure(figsize=(8.5, 5.3))
    plt.plot(df["z"], df["latent_depth_L"], label="Latent structural depth")
    plt.plot(df["z"], df["compression"], label="Saturating compression")
    plt.plot(df["z"], df["depth_loss"], label="Accessible depth loss")
    plt.xlabel("z")
    plt.ylabel("bounded proxy")
    plt.title("Latent Depth, Compression, and Depth Loss")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig4_latent_depth_compression_loss.png", dpi=260)
    plt.close()


def alpha_sweep():
    rows = []
    z = np.linspace(0.0, 3.0, 700)

    for alpha in np.linspace(0.0, 2.0, 25):
        df = compute_projection_cci(z, alpha=alpha)
        trans = transition_estimate(df)

        rows.append({
            "alpha": float(alpha),
            **trans,
            "final_CCI_norm": float(df["CCI_projection_norm"].iloc[-1]),
            "max_projection_residual": float(df["projection_residual"].max()),
        })

    sweep = pd.DataFrame(rows)
    sweep.to_csv(OUTDIR / "alpha_sweep_results_v02.csv", index=False)

    plt.figure(figsize=(8.5, 5.3))
    plt.plot(sweep["alpha"], sweep["z_transition_estimate"], marker="o")
    plt.xlabel("projection strength alpha")
    plt.ylabel("estimated transition redshift")
    plt.title("Projection Strength vs Transition Estimate")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig5_alpha_transition_sweep_v02.png", dpi=260)
    plt.close()

    plt.figure(figsize=(8.5, 5.3))
    plt.plot(sweep["alpha"], sweep["final_CCI_norm"], marker="o")
    plt.xlabel("projection strength alpha")
    plt.ylabel("final normalized CCI")
    plt.title("Projection Strength vs Final CCI")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTDIR / "fig6_alpha_final_CCI_v02.png", dpi=260)
    plt.close()

    return sweep


def main():
    z = np.linspace(0.0, 3.0, 900)

    df = compute_projection_cci(
        z,
        alpha=1.0,
        transition_z=1.1,
        sharpness=4.0,
        depth_floor=0.25,
        gamma=1.0,
    )

    trans = transition_estimate(df)
    sweep = alpha_sweep()

    df.to_csv(OUTDIR / "cci_projection_timeseries_v02.csv", index=False)
    save_plots(df)

    summary = {
        "model": "MAAT CCI Projection Test v0.2",
        "interpretation": "Stabilized toy model of CCI as projection observable of latent structural depth under expansion.",
        "redshift_range": [float(z.min()), float(z.max())],
        "baseline_projection_alpha": 1.0,
        "depth_floor": 0.25,
        "compression": "saturating tanh",
        "transition_estimate": trans,
        "final_CCI_norm": float(df["CCI_projection_norm"].iloc[-1]),
        "max_CCI_norm": float(df["CCI_projection_norm"].max()),
        "max_projection_residual": float(df["projection_residual"].max()),
        "note": "Toy projection model only; not a cosmological parameter fit."
    }

    with open(OUTDIR / "cci_projection_summary_v02.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\n=== MAAT CCI Projection Test v0.2 ===")
    for k, v in summary.items():
        print(f"{k}: {v}")

    print(f"\nSaved outputs to: {OUTDIR.resolve()}")
    print("Main files:")
    print(" - cci_projection_timeseries_v02.csv")
    print(" - alpha_sweep_results_v02.csv")
    print(" - cci_projection_summary_v02.json")
    print(" - fig1_breadth_depth_stabilized.png")
    print(" - fig2_stabilized_CCI_projection.png")
    print(" - fig3_bounded_projection_residual.png")
    print(" - fig4_latent_depth_compression_loss.png")
    print(" - fig5_alpha_transition_sweep_v02.png")
    print(" - fig6_alpha_final_CCI_v02.png")


if __name__ == "__main__":
    main()