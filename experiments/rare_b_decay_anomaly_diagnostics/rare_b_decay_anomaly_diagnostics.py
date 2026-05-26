#!/usr/bin/env python3
"""Toy structural diagnostics for rare B-decay anomaly patterns.

This script does not use LHCb or HEPData measurements.  It creates synthetic
pull patterns that illustrate a methodological distinction used in the paper:

1. a coherent Wilson-coefficient-like shift across q^2 bins and observables;
2. a localized hadronic/charm-like projection distortion;
3. a mixed case.

The goal is only to make the proposed MAAT v1.6 diagnostic reproducible.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


EPS = 1e-12
RNG_SEED = 4081
OBSERVABLES = ("P5prime", "BR", "FL")
N_SHUFFLES = 200


@dataclass(frozen=True)
class Diagnostic:
    H: float
    B: float
    S_eff: float
    V: float
    R_resp: float
    R_rob: float
    F_diag: float
    globality: float
    locality: float


def make_q2_bins() -> np.ndarray:
    """Representative low- and high-q^2 bin centres outside resonance windows."""
    return np.array([0.5, 1.1, 2.0, 3.2, 4.4, 5.6, 6.8, 8.2, 11.2, 12.4, 15.2, 17.4])


def build_synthetic_patterns(q2: np.ndarray, seed: int = RNG_SEED) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)

    # Observable-dependent response vector, loosely mimicking that angular and
    # rate observables can respond with different signs/magnitudes.
    response = np.array([1.0, 0.55, -0.35])
    noise = rng.normal(0.0, 0.08, size=(q2.size, len(response)))

    coherent = -1.15 * np.outer(np.ones_like(q2), response) + noise

    local_peak = np.exp(-0.5 * ((q2 - 8.4) / 0.85) ** 2)
    second_peak = 0.55 * np.exp(-0.5 * ((q2 - 12.6) / 0.75) ** 2)
    local_shape = local_peak - 0.7 * second_peak
    localized = -1.75 * np.outer(local_shape, response) + rng.normal(
        0.0, 0.08, size=(q2.size, len(response))
    )

    mixed = 0.58 * coherent + 0.50 * localized + rng.normal(
        0.0, 0.06, size=(q2.size, len(response))
    )

    standard_model_like = rng.normal(0.0, 0.18, size=(q2.size, len(response)))

    return {
        "standard_model_like_noise": standard_model_like,
        "coherent_wilson_shift": coherent,
        "localized_charm_like_projection": localized,
        "mixed_coherent_plus_local": mixed,
    }


def fit_constant_mode(pattern: np.ndarray) -> np.ndarray:
    """Best constant q^2-independent mode for each observable."""
    return np.mean(pattern, axis=0, keepdims=True) * np.ones_like(pattern)


def support_diagnostics(pattern: np.ndarray) -> Diagnostic:
    n_bins, n_obs = pattern.shape
    constant_fit = fit_constant_mode(pattern)
    residual_after_global = pattern - constant_fit

    rms = float(np.sqrt(np.mean(pattern**2)))
    global_rms = float(np.sqrt(np.mean(residual_after_global**2)))
    if rms < EPS:
        global_fraction = 0.0
    else:
        global_fraction = max(0.0, 1.0 - global_rms / (rms + EPS))

    # H: support for a coherent effective-Hamiltonian deformation.
    H = 1.0 / (1.0 + global_rms)

    # B: support for smooth q^2 behaviour; second differences penalise local
    # hadronic bumps more strongly than broad coherent shifts.
    second_diff = np.diff(pattern, n=2, axis=0)
    curvature = float(np.mean(np.abs(second_diff)))
    B = 1.0 / (1.0 + curvature)

    # S_eff: bounded anomaly activity. Very tiny and very extreme patterns are
    # both less useful for stable diagnostic interpretation.
    activity = float(np.mean(np.abs(pattern)))
    S_eff = float(np.exp(-0.5 * ((activity - 0.85) / 0.45) ** 2))

    # V: cross-observable connectivity. We measure whether observables share a
    # coherent q^2 mode after centering.
    centered = pattern - np.mean(pattern, axis=0, keepdims=True)
    corr_vals: list[float] = []
    for i in range(n_obs):
        for j in range(i + 1, n_obs):
            xi = centered[:, i]
            xj = centered[:, j]
            denom = float(np.linalg.norm(xi) * np.linalg.norm(xj))
            if denom > EPS:
                corr_vals.append(abs(float(np.dot(xi, xj) / denom)))
    V = float(np.mean(corr_vals)) if corr_vals else 0.0

    R_resp = float((H * B * V) ** (1.0 / 3.0))
    R_rob = float(min(R_resp, (H * B * S_eff * V) ** 0.25))

    F_diag = float(
        -np.log(EPS + H)
        - np.log(EPS + B)
        - np.log(EPS + S_eff)
        - np.log(EPS + V)
        - np.log(EPS + R_rob)
    )

    # Locality responds to large bin-to-bin residual after removing the global
    # mode. It is not a hadronic calculation, only a toy pattern diagnostic.
    tail95 = float(np.quantile(np.abs(residual_after_global), 0.95))
    locality = float(tail95 / (tail95 + rms + EPS))

    return Diagnostic(
        H=H,
        B=B,
        S_eff=S_eff,
        V=V,
        R_resp=R_resp,
        R_rob=R_rob,
        F_diag=F_diag,
        globality=global_fraction,
        locality=locality,
    )


def shuffled_bin_null(
    patterns: dict[str, np.ndarray],
    n_shuffles: int = N_SHUFFLES,
    seed: int = RNG_SEED + 1,
) -> pd.DataFrame:
    """Shuffle q^2 bins independently per observable while preserving amplitudes."""
    rng = np.random.default_rng(seed)
    rows = []
    for scenario, pattern in patterns.items():
        baseline = support_diagnostics(pattern)
        for k in range(n_shuffles):
            shuffled = pattern.copy()
            for col in range(shuffled.shape[1]):
                shuffled[:, col] = shuffled[rng.permutation(shuffled.shape[0]), col]
            diag = support_diagnostics(shuffled)
            rows.append(
                {
                    "scenario": scenario,
                    "shuffle_id": k,
                    "F_diag": diag.F_diag,
                    "globality": diag.globality,
                    "locality": diag.locality,
                    "delta_F_vs_real": diag.F_diag - baseline.F_diag,
                    "delta_globality_vs_real": diag.globality - baseline.globality,
                    "delta_locality_vs_real": diag.locality - baseline.locality,
                }
            )
    return pd.DataFrame(rows)


def save_plots(q2: np.ndarray, patterns: dict[str, np.ndarray], table: pd.DataFrame, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True, sharey=True)
    for ax, (name, pattern) in zip(axes.flat, patterns.items()):
        for idx, obs in enumerate(OBSERVABLES):
            ax.plot(q2, pattern[:, idx], marker="o", label=obs)
        ax.axhline(0, color="black", lw=0.8, alpha=0.5)
        ax.set_title(name.replace("_", " "))
        ax.set_ylabel("synthetic pull")
        ax.grid(alpha=0.25)
    axes[-1, 0].set_xlabel(r"$q^2$ bin centre [GeV$^2$]")
    axes[-1, 1].set_xlabel(r"$q^2$ bin centre [GeV$^2$]")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "fig1_synthetic_residual_patterns.png", dpi=180)
    plt.close(fig)

    support_cols = ["H", "B", "S_eff", "V", "R_rob"]
    fig, ax = plt.subplots(figsize=(10, 4.6))
    x = np.arange(len(table))
    width = 0.15
    for k, col in enumerate(support_cols):
        ax.bar(x + (k - 2) * width, table[col], width=width, label=col)
    ax.set_xticks(x)
    ax.set_xticklabels(table["scenario"], rotation=18, ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("support")
    ax.set_title("MAAT v1.6 toy supports for rare B-decay residual patterns")
    ax.legend(frameon=False, ncol=len(support_cols))
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out / "fig2_support_comparison.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.4, 6.2))
    sc = ax.scatter(
        table["globality"],
        table["locality"],
        c=table["F_diag"],
        s=180,
        cmap="viridis_r",
        edgecolor="black",
    )
    labels = {
        "standard_model_like_noise": "SM-like\nnoise",
        "coherent_wilson_shift": "coherent\nWilson\nshift",
        "localized_charm_like_projection": "localized\ncharm-like\nprojection",
        "mixed_coherent_plus_local": "mixed\ncoherent\n+ local",
    }
    offsets = {
        "standard_model_like_noise": (16, -22),
        "coherent_wilson_shift": (-54, 4),
        "localized_charm_like_projection": (14, 10),
        "mixed_coherent_plus_local": (14, 8),
    }
    for _, row in table.iterrows():
        scenario = str(row["scenario"])
        ax.annotate(
            labels.get(scenario, scenario.replace("_", "\n")),
            (row["globality"], row["locality"]),
            xytext=offsets.get(scenario, (8, 8)),
            textcoords="offset points",
            fontsize=8.5,
            ha="left",
            va="center",
            bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="none", alpha=0.78),
        )
    ax.set_xlim(-0.04, 1.02)
    ax.set_ylim(0.11, 0.74)
    ax.set_xlabel("global coherent-mode support")
    ax.set_ylabel("localized residual pressure")
    ax.set_title("Diagnostic plane: global Wilson-like vs local projection-like")
    ax.grid(alpha=0.25)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label(r"$F_{\rm diag}$ (lower is more supported)")
    fig.tight_layout()
    fig.savefig(out / "fig3_diagnostic_plane.png", dpi=180)
    plt.close(fig)


def save_null_plot(table: pd.DataFrame, null_summary: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.2, 4.6))
    x = np.arange(len(table))
    width = 0.36
    real_lookup = table.set_index("scenario")
    ordered = table["scenario"].tolist()
    real_cost = [real_lookup.loc[s, "F_diag"] for s in ordered]
    null_cost = [
        null_summary.loc[null_summary["scenario"] == s, "F_diag_mean"].iloc[0]
        for s in ordered
    ]
    ax.bar(x - width / 2, real_cost, width=width, label="real toy pattern")
    ax.bar(x + width / 2, null_cost, width=width, label="bin-shuffled null mean")
    ax.set_xticks(x)
    ax.set_xticklabels(ordered, rotation=18, ha="right")
    ax.set_ylabel(r"$F_{\rm diag}$ (lower is more supported)")
    ax.set_title("q²-bin shuffled null: same amplitudes, higher structural stress")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out / "fig4_bin_shuffled_null.png", dpi=180)
    plt.close(fig)


def main() -> None:
    base = Path(__file__).resolve().parent
    out = base / "outputs_phenomenological"
    out.mkdir(parents=True, exist_ok=True)

    q2 = make_q2_bins()
    patterns = build_synthetic_patterns(q2)

    rows = []
    residual_rows = []
    for scenario, pattern in patterns.items():
        diag = support_diagnostics(pattern)
        rows.append({"scenario": scenario, **diag.__dict__})
        for i, q in enumerate(q2):
            row = {"scenario": scenario, "q2": q}
            for j, obs in enumerate(OBSERVABLES):
                row[obs] = pattern[i, j]
            residual_rows.append(row)

    table = pd.DataFrame(rows)
    residuals = pd.DataFrame(residual_rows)
    null_samples = shuffled_bin_null(patterns)
    null_summary = (
        null_samples.groupby("scenario")
        .agg(
            F_diag_mean=("F_diag", "mean"),
            F_diag_std=("F_diag", "std"),
            globality_mean=("globality", "mean"),
            globality_std=("globality", "std"),
            locality_mean=("locality", "mean"),
            locality_std=("locality", "std"),
            delta_F_vs_real_mean=("delta_F_vs_real", "mean"),
            delta_globality_vs_real_mean=("delta_globality_vs_real", "mean"),
            delta_locality_vs_real_mean=("delta_locality_vs_real", "mean"),
        )
        .reset_index()
    )
    table.to_csv(out / "toy_diagnostic_summary.csv", index=False)
    residuals.to_csv(out / "toy_residual_patterns.csv", index=False)
    null_samples.to_csv(out / "toy_bin_shuffled_null_samples.csv", index=False)
    null_summary.to_csv(out / "toy_bin_shuffled_null_summary.csv", index=False)

    summary = {
        "status": "synthetic_toy_only_no_lhcb_data",
        "rng_seed": RNG_SEED,
        "n_shuffles": N_SHUFFLES,
        "n_q2_bins": int(q2.size),
        "observables": list(OBSERVABLES),
        "best_supported_scenario_by_F_diag": str(table.loc[table["F_diag"].idxmin(), "scenario"]),
        "highest_globality_scenario": str(table.loc[table["globality"].idxmax(), "scenario"]),
        "highest_locality_scenario": str(table.loc[table["locality"].idxmax(), "scenario"]),
        "rows": table.to_dict(orient="records"),
        "bin_shuffled_null_summary": null_summary.to_dict(orient="records"),
    }
    (out / "toy_diagnostic_summary.json").write_text(json.dumps(summary, indent=2))
    save_plots(q2, patterns, table, out)
    save_null_plot(table, null_summary, out)

    print(table.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print("\nBin-shuffled null summary:")
    print(null_summary.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print(f"\nOutputs written to: {out}")


if __name__ == "__main__":
    main()
