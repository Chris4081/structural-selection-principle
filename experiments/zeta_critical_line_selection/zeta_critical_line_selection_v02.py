#!/usr/bin/env python3
"""Extra Paper: MAAT diagnostic test on zeta-zero critical-line stability.

This is a diagnostic toy benchmark, not a Riemann-Hypothesis proof.

The benchmark asks whether known Riemann zeta zeros have a lower MAAT-style
structural cost on the critical line than under small off-line shifts.  It
also runs controls and ablations:

* random critical-line points,
* jittered zeros,
* without the explicit balance-to-critical-line term,
* without the direct zeta-null support term,
* parameter sweeps over alpha and robustness displacement.

The key safety principle is that the paper must report when the minimum is
partly built into the diagnostic.
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).parent / ".mplconfig"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import pandas as pd


EPS = 1.0e-12
SEED = 4081
MP_DPS = 50
DEFAULT_ALPHA = 10.0
DEFAULT_DSIGMA = 0.002
DEFAULT_PAIR_ZEROS = 200

DELTAS = np.array(
    [-0.050, -0.020, -0.010, -0.005, -0.002, 0.0, 0.002, 0.005, 0.010, 0.020, 0.050],
    dtype=float,
)


@dataclass(frozen=True)
class ScoreConfig:
    alpha: float
    dsigma: float
    use_H: bool = True
    use_B: bool = True
    use_R: bool = True

    @property
    def label(self) -> str:
        parts = []
        if self.use_H:
            parts.append("H")
        if self.use_B:
            parts.append("B")
        if self.use_R:
            parts.append("R")
        return "".join(parts) or "none"


def normalize01(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    lo = np.nanmin(x)
    hi = np.nanmax(x)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo < EPS:
        return np.zeros_like(x, dtype=float)
    return (x - lo) / (hi - lo)


def zeta_abs(sigma: float, t: float) -> float:
    return float(abs(mp.zeta(mp.mpc(sigma, t))))


def H_support(sigma: float, t: float) -> float:
    return 1.0 / (1.0 + zeta_abs(sigma, t))


def B_balance(sigma: float, *, alpha: float) -> float:
    return 1.0 / (1.0 + alpha * abs(sigma - 0.5))


def R_robustness(sigma: float, t: float, *, dsigma: float) -> float:
    z_plus = mp.zeta(mp.mpc(sigma + dsigma, t))
    z_minus = mp.zeta(mp.mpc(sigma - dsigma, t))
    return 1.0 / (1.0 + float(abs(z_plus - z_minus)))


def F_maat(sigma: float, t: float, config: ScoreConfig) -> tuple[float, dict[str, float]]:
    supports: dict[str, float] = {
        "H": H_support(sigma, t),
        "B": B_balance(sigma, alpha=config.alpha),
        "R": R_robustness(sigma, t, dsigma=config.dsigma),
    }
    selected = []
    if config.use_H:
        selected.append(supports["H"])
    if config.use_B:
        selected.append(supports["B"])
    if config.use_R:
        selected.append(supports["R"])
    if not selected:
        return 0.0, supports
    F = float(-sum(np.log(EPS + np.asarray(selected, dtype=float))))
    return F, supports


def generate_t_values(n_zeros: int, n_random: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    zeros = [float(mp.im(mp.zetazero(k))) for k in range(1, n_zeros + 1)]
    lo, hi = min(zeros), max(zeros)

    rows: list[dict] = []
    for idx, t in enumerate(zeros, start=1):
        rows.append({"set": "known_zero", "source_n": idx, "t": t})
        rows.append(
            {
                "set": "jittered_zero",
                "source_n": idx,
                "t": float(t + rng.normal(0.0, 0.12)),
            }
        )

    random_ts = rng.uniform(lo, hi, size=n_random)
    for idx, t in enumerate(random_ts, start=1):
        rows.append({"set": "random_critical_line", "source_n": idx, "t": float(t)})
    return pd.DataFrame(rows)


def evaluate_grid(
    t_df: pd.DataFrame,
    config: ScoreConfig,
    *,
    deltas: np.ndarray = DELTAS,
) -> pd.DataFrame:
    rows: list[dict] = []
    for _, item in t_df.iterrows():
        for delta in deltas:
            sigma = 0.5 + float(delta)
            F, supports = F_maat(sigma, float(item["t"]), config)
            rows.append(
                {
                    "set": item["set"],
                    "source_n": int(item["source_n"]),
                    "t": float(item["t"]),
                    "delta": float(delta),
                    "sigma": sigma,
                    "config": config.label,
                    "alpha": config.alpha,
                    "dsigma": config.dsigma,
                    "F_MAAT": F,
                    **supports,
                }
            )
    return pd.DataFrame(rows)


def summarize_delta(df: pd.DataFrame) -> pd.DataFrame:
    summary = (
        df.groupby(["set", "config", "alpha", "dsigma", "delta"], as_index=False)
        .agg(
            mean_F=("F_MAAT", "mean"),
            median_F=("F_MAAT", "median"),
            std_F=("F_MAAT", "std"),
            mean_H=("H", "mean"),
            mean_B=("B", "mean"),
            mean_R=("R", "mean"),
        )
        .sort_values(["set", "config", "alpha", "dsigma", "delta"])
    )
    return summary


def summarize_minima(summary: pd.DataFrame) -> pd.DataFrame:
    rows: list[pd.Series] = []
    group_cols = ["set", "config", "alpha", "dsigma"]
    for _, group in summary.groupby(group_cols):
        best = group.loc[group["mean_F"].idxmin()].copy()
        f0 = group.loc[np.isclose(group["delta"], 0.0), "mean_F"]
        best["mean_F_at_delta0"] = float(f0.iloc[0]) if len(f0) else np.nan
        best["delta0_margin_vs_best"] = float(best["mean_F_at_delta0"] - best["mean_F"])
        rows.append(best)
    return pd.DataFrame(rows).reset_index(drop=True)


def run_parameter_sweep(t_df: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    alpha_values = [0.0, 2.0, 5.0, 10.0, 20.0]
    dsigma_values = [0.001, 0.002, 0.005]
    configs: list[ScoreConfig] = []
    for alpha in alpha_values:
        for dsigma in dsigma_values:
            configs.append(ScoreConfig(alpha=alpha, dsigma=dsigma, use_H=True, use_B=True, use_R=True))
            configs.append(ScoreConfig(alpha=alpha, dsigma=dsigma, use_H=True, use_B=False, use_R=True))
            configs.append(ScoreConfig(alpha=alpha, dsigma=dsigma, use_H=False, use_B=True, use_R=True))

    pieces = []
    for config in configs:
        pieces.append(evaluate_grid(t_df, config))
    sweep_df = pd.concat(pieces, ignore_index=True)
    sweep_summary = summarize_delta(sweep_df)
    sweep_minima = summarize_minima(sweep_summary)
    sweep_df.to_csv(output_dir / "zeta_sweep_grid.csv", index=False)
    sweep_summary.to_csv(output_dir / "zeta_sweep_delta_summary.csv", index=False)
    sweep_minima.to_csv(output_dir / "zeta_sweep_minima.csv", index=False)
    return sweep_minima


def gue_wigner_surmise_density(s: np.ndarray) -> np.ndarray:
    """Unit-mean GUE nearest-neighbour spacing proxy.

    This is the Wigner-surmise approximation, not the full Fredholm-determinant
    spacing law.  It is sufficient for a conservative diagnostic benchmark.
    """

    s = np.asarray(s, dtype=float)
    return (32.0 / np.pi**2) * s**2 * np.exp(-4.0 * s**2 / np.pi)


def spacing_probabilities(spacings: np.ndarray, bins: np.ndarray) -> np.ndarray:
    hist, _ = np.histogram(spacings, bins=bins)
    probs = hist.astype(float)
    if probs.sum() <= EPS:
        return np.zeros(len(bins) - 1, dtype=float)
    return probs / probs.sum()


def pair_correlation_structural_test(
    *,
    n_pair_zeros: int,
    seed: int,
    output_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Compare normalized zeta-zero spacings with a GUE spacing proxy."""

    rng = np.random.default_rng(seed + 2026)
    zeros = np.array([float(mp.im(mp.zetazero(k))) for k in range(1, n_pair_zeros + 1)])
    raw_spacing = np.diff(zeros)
    local_density = np.log(zeros[:-1] / (2.0 * np.pi)) / (2.0 * np.pi)
    zeta_spacing = raw_spacing * local_density
    zeta_spacing = zeta_spacing / np.mean(zeta_spacing)

    poisson_spacing = rng.exponential(scale=1.0, size=len(zeta_spacing))

    uniform_points = np.sort(rng.uniform(zeros[0], zeros[-1], size=n_pair_zeros))
    uniform_spacing = np.diff(uniform_points)
    uniform_spacing = uniform_spacing / np.mean(uniform_spacing)

    jittered = np.sort(zeros + rng.normal(0.0, 0.04, size=len(zeros)))
    jittered_spacing = np.diff(jittered)
    jittered_density = np.log(jittered[:-1] / (2.0 * np.pi)) / (2.0 * np.pi)
    jittered_spacing = jittered_spacing * jittered_density
    jittered_spacing = jittered_spacing / np.mean(jittered_spacing)

    bins = np.linspace(0.0, 4.0, 33)
    centers = 0.5 * (bins[:-1] + bins[1:])
    width = bins[1] - bins[0]
    gue_prob = gue_wigner_surmise_density(centers) * width
    gue_prob = gue_prob / gue_prob.sum()

    series = {
        "zeta_normalized_spacings": zeta_spacing,
        "jittered_zeta_spacings": jittered_spacing,
        "poisson_unit_mean_control": poisson_spacing,
        "uniform_points_control": uniform_spacing,
    }

    rows: list[dict] = []
    summary_rows: list[dict] = []
    for name, spacings in series.items():
        probs = spacing_probabilities(spacings, bins)
        abs_diff = np.abs(probs - gue_prob)
        V_bins = 1.0 / (1.0 + abs_diff)
        V_weighted = float(np.sum(gue_prob * V_bins))
        V_mean = float(np.mean(V_bins))
        l1 = float(np.sum(abs_diff))
        l2 = float(np.sqrt(np.sum((probs - gue_prob) ** 2)))
        F_pair = float(-np.log(EPS + V_weighted))
        for center, p_zeta, p_gue, diff, v_bin in zip(centers, probs, gue_prob, abs_diff, V_bins):
            rows.append(
                {
                    "series": name,
                    "spacing_center": float(center),
                    "P_empirical": float(p_zeta),
                    "P_GUE_proxy": float(p_gue),
                    "abs_diff": float(diff),
                    "V_bin": float(v_bin),
                }
            )
        summary_rows.append(
            {
                "series": name,
                "n_spacings": int(len(spacings)),
                "mean_spacing": float(np.mean(spacings)),
                "std_spacing": float(np.std(spacings)),
                "L1_to_GUE": l1,
                "L2_to_GUE": l2,
                "V_pair_weighted": V_weighted,
                "V_pair_mean": V_mean,
                "F_pair": F_pair,
            }
        )

    pair_df = pd.DataFrame(rows)
    pair_summary = pd.DataFrame(summary_rows).sort_values("F_pair")
    pair_df.to_csv(output_dir / "zeta_pair_correlation_bins.csv", index=False)
    pair_summary.to_csv(output_dir / "zeta_pair_correlation_summary.csv", index=False)

    metadata = {
        "n_pair_zeros": n_pair_zeros,
        "n_spacings": int(len(zeta_spacing)),
        "spacing_normalisation": "Delta t_n * log(t_n/(2pi))/(2pi), rescaled to mean 1",
        "gue_reference": "Wigner-surmise unit-mean GUE nearest-neighbour spacing proxy",
    }
    return pair_df, pair_summary, metadata


def make_figures(
    grid_df: pd.DataFrame,
    delta_summary: pd.DataFrame,
    minima: pd.DataFrame,
    sweep_minima: pd.DataFrame,
    pair_df: pd.DataFrame,
    pair_summary: pd.DataFrame,
    output_dir: Path,
) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")

    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    for set_name, group in delta_summary[delta_summary["config"] == "HBR"].groupby("set"):
        ax.plot(group["delta"], group["mean_F"], marker="o", label=set_name.replace("_", " "))
    ax.axvline(0.0, color="black", linewidth=1, linestyle="--")
    ax.set_xlabel(r"off-line shift $\delta=\sigma-1/2$")
    ax.set_ylabel(r"mean $F_{\mathrm{MAAT}}$")
    ax.set_title("Full HBR diagnostic: known zeros versus controls")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "fig1_full_hbr_delta_scan.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 5.0))
    config_order = ["HBR", "HR", "BR"]
    for config in config_order:
        subset = delta_summary[
            (delta_summary["set"] == "known_zero")
            & (delta_summary["config"] == config)
            & np.isclose(delta_summary["alpha"], DEFAULT_ALPHA)
            & np.isclose(delta_summary["dsigma"], DEFAULT_DSIGMA)
        ]
        if len(subset):
            ax.plot(subset["delta"], subset["mean_F"], marker="o", label=config)
    ax.axvline(0.0, color="black", linewidth=1, linestyle="--")
    ax.set_xlabel(r"off-line shift $\delta$")
    ax.set_ylabel(r"mean $F_{\mathrm{MAAT}}$")
    ax.set_title("Ablation scan on known zeta zeros")
    ax.legend(title="supports")
    fig.tight_layout()
    fig.savefig(output_dir / "fig2_known_zero_ablation_scan.png", dpi=180)
    plt.close(fig)

    pivot = minima[
        (np.isclose(minima["alpha"], DEFAULT_ALPHA))
        & (np.isclose(minima["dsigma"], DEFAULT_DSIGMA))
    ].pivot(index="config", columns="set", values="delta")
    fig, ax = plt.subplots(figsize=(8.2, 3.8))
    im = ax.imshow(np.abs(pivot.values), cmap="viridis_r", vmin=0, vmax=float(np.max(np.abs(DELTAS))))
    ax.set_xticks(np.arange(len(pivot.columns)), labels=[c.replace("_", "\n") for c in pivot.columns])
    ax.set_yticks(np.arange(len(pivot.index)), labels=pivot.index)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            ax.text(j, i, f"{pivot.iloc[i, j]:+.3f}", ha="center", va="center", color="white", fontsize=9)
    ax.set_title("Best mean off-line shift by control set and ablation")
    fig.colorbar(im, ax=ax, label=r"$|\delta_{\mathrm{best}}|$")
    fig.tight_layout()
    fig.savefig(output_dir / "fig3_best_delta_heatmap.png", dpi=180)
    plt.close(fig)

    known = sweep_minima[sweep_minima["set"] == "known_zero"].copy()
    known["delta_is_zero"] = np.isclose(known["delta"], 0.0)
    rate = (
        known.groupby(["config", "alpha"], as_index=False)["delta_is_zero"]
        .mean()
        .rename(columns={"delta_is_zero": "fraction_delta0_best"})
    )
    fig, ax = plt.subplots(figsize=(8.2, 4.5))
    for config, group in rate.groupby("config"):
        ax.plot(group["alpha"], group["fraction_delta0_best"], marker="o", label=config)
    ax.set_xlabel(r"balance strength $\alpha$")
    ax.set_ylabel(r"fraction of $\Delta\sigma_R$ sweeps with $\delta=0$ best")
    ax.set_ylim(-0.05, 1.05)
    ax.set_title("Parameter-sweep stability of the critical-line minimum")
    ax.legend(title="supports")
    fig.tight_layout()
    fig.savefig(output_dir / "fig4_parameter_sweep_delta0_fraction.png", dpi=180)
    plt.close(fig)

    supports = ["H", "B", "R"]
    subset = grid_df[
        (grid_df["set"] == "known_zero")
        & (grid_df["config"] == "HBR")
        & np.isclose(grid_df["delta"], 0.0)
    ]
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    values = [subset[s].mean() for s in supports]
    ax.bar(supports, values, color=["#1f77b4", "#ff7f0e", "#2ca02c"], alpha=0.85)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("mean support at known zeros")
    ax.set_title("Support decomposition on the critical line")
    fig.tight_layout()
    fig.savefig(output_dir / "fig5_support_decomposition.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 5.0))
    for name, group in pair_df.groupby("series"):
        if name == "zeta_normalized_spacings":
            lw = 2.8
            alpha = 0.95
        else:
            lw = 1.7
            alpha = 0.72
        ax.plot(
            group["spacing_center"],
            group["P_empirical"],
            marker="o",
            linewidth=lw,
            alpha=alpha,
            label=name.replace("_", " "),
        )
    gue = pair_df[pair_df["series"] == "zeta_normalized_spacings"]
    ax.plot(
        gue["spacing_center"],
        gue["P_GUE_proxy"],
        color="black",
        linewidth=2.2,
        linestyle="--",
        label="GUE proxy",
    )
    ax.set_xlabel("normalised nearest-neighbour spacing")
    ax.set_ylabel("bin probability")
    ax.set_title("Pair-correlation diagnostic: spacing distribution vs GUE proxy")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output_dir / "fig6_pair_correlation_spacing_distribution.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 4.4))
    ordered = pair_summary.sort_values("F_pair")
    ax.bar(ordered["series"].str.replace("_", "\n"), ordered["V_pair_weighted"], color="#1f77b4", alpha=0.85)
    ax.set_ylabel(r"$V_{\mathrm{pair}}$ weighted support")
    ax.set_ylim(0.0, 1.02)
    ax.set_title("Structural support against GUE spacing proxy")
    fig.tight_layout()
    fig.savefig(output_dir / "fig7_pair_correlation_support.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent / "outputs")
    parser.add_argument("--n-zeros", type=int, default=40)
    parser.add_argument("--n-random", type=int, default=40)
    parser.add_argument("--n-pair-zeros", type=int, default=DEFAULT_PAIR_ZEROS)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()

    mp.mp.dps = MP_DPS
    args.output_dir.mkdir(parents=True, exist_ok=True)

    t_df = generate_t_values(args.n_zeros, args.n_random, args.seed)
    t_df.to_csv(args.output_dir / "zeta_t_sets.csv", index=False)

    configs = [
        ScoreConfig(DEFAULT_ALPHA, DEFAULT_DSIGMA, True, True, True),
        ScoreConfig(DEFAULT_ALPHA, DEFAULT_DSIGMA, True, False, True),
        ScoreConfig(DEFAULT_ALPHA, DEFAULT_DSIGMA, False, True, True),
    ]
    grid_df = pd.concat([evaluate_grid(t_df, config) for config in configs], ignore_index=True)
    delta_summary = summarize_delta(grid_df)
    minima = summarize_minima(delta_summary)

    grid_df.to_csv(args.output_dir / "zeta_default_grid.csv", index=False)
    delta_summary.to_csv(args.output_dir / "zeta_default_delta_summary.csv", index=False)
    minima.to_csv(args.output_dir / "zeta_default_minima.csv", index=False)

    sweep_minima = run_parameter_sweep(t_df, args.output_dir)
    pair_df, pair_summary, pair_metadata = pair_correlation_structural_test(
        n_pair_zeros=args.n_pair_zeros,
        seed=args.seed,
        output_dir=args.output_dir,
    )
    make_figures(grid_df, delta_summary, minima, sweep_minima, pair_df, pair_summary, args.output_dir)

    full_known = minima[
        (minima["set"] == "known_zero")
        & (minima["config"] == "HBR")
        & np.isclose(minima["alpha"], DEFAULT_ALPHA)
        & np.isclose(minima["dsigma"], DEFAULT_DSIGMA)
    ].iloc[0]
    no_b_known = minima[
        (minima["set"] == "known_zero")
        & (minima["config"] == "HR")
        & np.isclose(minima["alpha"], DEFAULT_ALPHA)
        & np.isclose(minima["dsigma"], DEFAULT_DSIGMA)
    ].iloc[0]
    no_h_known = minima[
        (minima["set"] == "known_zero")
        & (minima["config"] == "BR")
        & np.isclose(minima["alpha"], DEFAULT_ALPHA)
        & np.isclose(minima["dsigma"], DEFAULT_DSIGMA)
    ].iloc[0]

    known_sweep = sweep_minima[sweep_minima["set"] == "known_zero"].copy()
    sweep_stats = (
        known_sweep.assign(delta_is_zero=np.isclose(known_sweep["delta"], 0.0))
        .groupby("config", as_index=False)
        .agg(
            delta0_best_fraction=("delta_is_zero", "mean"),
            mean_abs_best_delta=("delta", lambda x: float(np.mean(np.abs(x)))),
        )
    )

    summary = {
        "seed": args.seed,
        "n_zeros": args.n_zeros,
        "n_random": args.n_random,
        "n_pair_zeros": args.n_pair_zeros,
        "deltas": DELTAS.tolist(),
        "default_alpha": DEFAULT_ALPHA,
        "default_dsigma": DEFAULT_DSIGMA,
        "headline": {
            "known_zero_full_HBR_best_delta": float(full_known["delta"]),
            "known_zero_full_HBR_mean_F": float(full_known["mean_F"]),
            "known_zero_without_B_HR_best_delta": float(no_b_known["delta"]),
            "known_zero_without_H_BR_best_delta": float(no_h_known["delta"]),
            "delta0_margin_vs_best_full_HBR": float(full_known["delta0_margin_vs_best"]),
        },
        "sweep_stats_known_zero": sweep_stats.to_dict(orient="records"),
        "pair_correlation": {
            **pair_metadata,
            "summary": pair_summary.to_dict(orient="records"),
        },
    }
    with open(args.output_dir / "zeta_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("[zeta-extra] default minima")
    print(minima.to_string(index=False))
    print("\n[zeta-extra] known-zero sweep stats")
    print(sweep_stats.to_string(index=False))
    print("\n[zeta-extra] pair-correlation summary")
    print(pair_summary.to_string(index=False))
    print(f"\n[zeta-extra] outputs written to {args.output_dir}")


if __name__ == "__main__":
    main()
