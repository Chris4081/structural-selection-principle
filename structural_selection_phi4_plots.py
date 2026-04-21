#!/usr/bin/env python3
"""Plot and summarize the static phi^4 structural-selection benchmark results."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path.cwd() / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path.cwd() / ".cache"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate plots for the static phi^4 structural-selection benchmark."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("structural_selection_phi4_results.json"),
        help="Path to the JSON results file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("structural_selection_phi4_plots"),
        help="Directory where plots and summary files are written.",
    )
    return parser.parse_args()


def load_results(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def configure_style() -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update(
        {
            "figure.dpi": 140,
            "savefig.dpi": 180,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "legend.fontsize": 9,
            "font.size": 10,
        }
    )


def save_protocol_a(results: dict, output_dir: Path) -> None:
    protocol_a = results["protocol_a"]
    perturbations = protocol_a["perturbations"]
    eta = np.array([row["eta"] for row in perturbations], dtype=float)
    f_values = np.array([row["F_struct"] for row in perturbations], dtype=float)
    ordering_pass = np.array([row["ordering_pass"] for row in perturbations], dtype=bool)

    vac_f = protocol_a["exact_states"]["vacuum_plus"]["F_struct"]
    saddle_f = protocol_a["exact_states"]["saddle"]["F_struct"]
    eta_max = float(protocol_a["moderate_eta_max"])

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(eta, f_values, color="#1f77b4", lw=2.2, marker="o", label="Perturbed vacuum")
    ax.scatter(
        eta[ordering_pass],
        f_values[ordering_pass],
        color="#2ca02c",
        s=38,
        zorder=3,
        label="Ordering pass",
    )
    ax.scatter(
        eta[~ordering_pass],
        f_values[~ordering_pass],
        color="#d62728",
        s=38,
        zorder=3,
        label="Ordering fail",
    )
    ax.axhline(vac_f, color="#444444", lw=1.6, linestyle="--", label="Exact vacuum")
    ax.axhline(saddle_f, color="#9467bd", lw=1.8, linestyle=":", label="False saddle")
    ax.axvspan(eta.min(), eta_max, color="#f2f2f2", alpha=0.65, label="Moderate range")
    ax.set_title("Protocol A: Vacuum Perturbations")
    ax.set_xlabel(r"Perturbation amplitude $\eta$")
    ax.set_ylabel(r"$F_{\mathrm{struct}}$")
    ax.legend(loc="upper left", frameon=True)
    ax.text(
        0.98,
        0.03,
        (
            f"vacuum_vs_saddle_pass = {protocol_a['vacuum_vs_saddle_pass']}\n"
            f"moderate_perturbation_pass = {protocol_a['moderate_perturbation_pass']}"
        ),
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.9, "edgecolor": "#cccccc"},
    )
    fig.tight_layout()
    fig.savefig(output_dir / "protocol_a_vacuum_discrimination.png", bbox_inches="tight")
    plt.close(fig)


def save_protocol_b(results: dict, output_dir: Path) -> None:
    protocol_b = results["protocol_b"]
    distortions = protocol_b["distortions"]
    kink_f = protocol_b["kink"]["F_struct"]

    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    colors = {1: "#1f77b4", 2: "#ff7f0e", 3: "#2ca02c"}
    for mode in (1, 2, 3):
        subset = [row for row in distortions if row["mode"] == mode]
        epsilon = np.array([row["epsilon"] for row in subset], dtype=float)
        f_values = np.array([row["F_struct"] for row in subset], dtype=float)
        ax.plot(
            epsilon,
            f_values,
            marker="o",
            lw=1.8,
            color=colors[mode],
            label=f"Mode {mode}",
        )
    ax.axhline(kink_f, color="#222222", lw=2.0, linestyle="--", label="Discrete kink")
    ax.set_title("Protocol B: Kink Versus Sector-Matched Distortions")
    ax.set_xlabel(r"Distortion amplitude $\varepsilon$")
    ax.set_ylabel(r"$F_{\mathrm{struct}}$")
    ax.legend(loc="upper left", ncol=2, frameon=True)
    ax.text(
        0.98,
        0.03,
        f"$p_K$ = {protocol_b['p_k']:.3f}\npass = {protocol_b['pass']}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.9, "edgecolor": "#cccccc"},
    )
    fig.tight_layout()
    fig.savefig(output_dir / "protocol_b_kink_discrimination.png", bbox_inches="tight")
    plt.close(fig)


def save_protocol_c(results: dict, output_dir: Path) -> None:
    protocol_c = results["protocol_c"]
    rows = protocol_c["rank_comparison"]
    percentiles = [f"{row['percentile']}%" for row in rows]
    d_struct = np.array([row["mean_distance_structural"] for row in rows], dtype=float)
    d_energy = np.array([row["mean_distance_energy"] for row in rows], dtype=float)
    x = np.arange(len(rows))
    width = 0.34

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    bars1 = ax.bar(x - width / 2, d_struct, width, color="#1f77b4", label="Structural ranking")
    bars2 = ax.bar(x + width / 2, d_energy, width, color="#ff7f0e", label="Energy ranking")

    for row, bar in zip(rows, bars1):
        if row["structural_better"]:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                "better",
                ha="center",
                va="bottom",
                fontsize=8,
                color="#1f77b4",
            )

    ax.set_title("Protocol C: Ranking Quality Near the Kink")
    ax.set_xlabel("Lowest-ranked subset")
    ax.set_ylabel("Mean distance to discrete kink")
    ax.set_xticks(x)
    ax.set_xticklabels(percentiles)
    ax.legend(loc="upper left", frameon=True)
    ax.text(
        0.98,
        0.03,
        f"any_structural_gain = {protocol_c['any_structural_gain']}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.9, "edgecolor": "#cccccc"},
    )
    fig.tight_layout()
    fig.savefig(output_dir / "protocol_c_ranking_comparison.png", bbox_inches="tight")
    plt.close(fig)


def save_distortion_scatter(results: dict, output_dir: Path) -> None:
    distortions = results["protocol_b"]["distortions"]
    f_values = np.array([row["F_struct"] for row in distortions], dtype=float)
    e_values = np.array([row["E_stat"] for row in distortions], dtype=float)
    distances = np.array([row["distance_to_kink"] for row in distortions], dtype=float)
    modes = np.array([row["mode"] for row in distortions], dtype=int)

    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    scatter = ax.scatter(
        e_values,
        f_values,
        c=distances,
        cmap="viridis",
        s=46,
        alpha=0.9,
        edgecolors="none",
    )
    for mode, marker in zip((1, 2, 3), ("o", "s", "^")):
        subset = modes == mode
        ax.scatter(
            e_values[subset],
            f_values[subset],
            facecolors="none",
            edgecolors="#333333",
            marker=marker,
            s=70,
            linewidths=0.8,
            label=f"Mode {mode}",
        )
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("Distance to discrete kink")
    ax.set_title("Distortion Ensemble: Structural Score Versus Static Energy")
    ax.set_xlabel(r"$E_{\mathrm{stat}}$")
    ax.set_ylabel(r"$F_{\mathrm{struct}}$")
    ax.legend(loc="upper left", frameon=True)
    fig.tight_layout()
    fig.savefig(output_dir / "distortion_energy_vs_structural_scatter.png", bbox_inches="tight")
    plt.close(fig)


def write_summary(results: dict, output_dir: Path) -> None:
    protocol_a = results["protocol_a"]
    protocol_b = results["protocol_b"]
    protocol_c = results["protocol_c"]

    summary = f"""# Structural Selection phi4 Plot Summary

- Protocol A
  - vacuum_vs_saddle_pass: {protocol_a['vacuum_vs_saddle_pass']}
  - moderate_perturbation_pass: {protocol_a['moderate_perturbation_pass']}
  - delta_vac: {protocol_a['delta_vac']:.6f}
- Protocol B
  - p_K: {protocol_b['p_k']:.6f}
  - threshold: {protocol_b['success_threshold']:.2f}
  - pass: {protocol_b['pass']}
- Protocol C
  - any_structural_gain: {protocol_c['any_structural_gain']}
"""
    for row in protocol_c["rank_comparison"]:
        summary += (
            f"  - percentile {row['percentile']}%: "
            f"d_F={row['mean_distance_structural']:.6f}, "
            f"d_E={row['mean_distance_energy']:.6f}, "
            f"structural_better={row['structural_better']}\n"
        )

    (output_dir / "summary.md").write_text(summary, encoding="utf-8")


def main() -> None:
    args = parse_args()
    configure_style()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

    results = load_results(args.input)
    save_protocol_a(results, args.output_dir)
    save_protocol_b(results, args.output_dir)
    save_protocol_c(results, args.output_dir)
    save_distortion_scatter(results, args.output_dir)
    write_summary(results, args.output_dir)

    print(f"Wrote plots to {args.output_dir}")


if __name__ == "__main__":
    main()
