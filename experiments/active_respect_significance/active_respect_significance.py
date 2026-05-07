#!/usr/bin/env python3
"""
Active Respect / Structural Significance Toy Simulation
======================================================

This script tests an extension of the MAAT v1.2.1 closure idea:

    H, B, S, V  = primitive support sectors
    R_resp      = passive structural adherence
    R_rob       = activity-sensitive robustness
    R_sig       = active structural significance

The goal is not to introduce R as a fifth primitive field. Instead, the
simulation treats respect/robustness/significance as emergent functions of the
four primitive support sectors and a controlled-activity window.

The toy model is intentionally simple and deterministic once the seed is fixed.
It is meant as a reproducible sanity check for the idea that a dormant but
coherent system can have passive structural adherence, while an actively
meaningful system requires controlled activity near an optimal window.
"""

from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path

import os

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "maat_active_respect_matplotlib"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT = BASE_DIR / "outputs"
OUT.mkdir(exist_ok=True)

SEED = 42
EPS = 1e-12
A_STAR = 1.0
SIGMA_A = 0.28
ALPHA_SIG = 0.45


def clip01(x: np.ndarray | float) -> np.ndarray | float:
    return np.clip(x, 1e-6, 1.0)


def sigmoid(x: np.ndarray | float) -> np.ndarray | float:
    return 1.0 / (1.0 + np.exp(-x))


def controlled_activity_support(
    activity: np.ndarray | float,
    a_star: float = A_STAR,
    sigma_a: float = SIGMA_A,
) -> np.ndarray | float:
    """Activity support: high only inside the controlled creative window."""
    activity = np.asarray(activity)
    return np.exp(-0.5 * ((activity - a_star) / sigma_a) ** 2)


def primitive_supports(activity: np.ndarray | float) -> tuple[np.ndarray, ...]:
    """
    Smooth toy supports for H, B, S_eff, V as a function of raw activity.

    Interpretation:
    - H remains possible in quiet systems but degrades in chaotic overdrive.
    - B is high for controlled systems and falls when activity overshoots.
    - S_eff is the controlled-activity window, not raw activity.
    - V grows with activity/coupling but can be damaged by overdrive.
    """
    activity = np.asarray(activity)
    s_eff = controlled_activity_support(activity)
    overdrive = sigmoid((activity - 1.55) / 0.16)

    h = 0.65 + 0.35 * s_eff - 0.38 * overdrive
    b = 0.75 + 0.25 * s_eff - 0.48 * overdrive
    v = (0.25 + 0.75 * sigmoid((activity - 0.25) / 0.18)) * (1.0 - 0.22 * overdrive)

    return clip01(h), clip01(b), clip01(s_eff), clip01(v)


def closures(
    h: np.ndarray,
    b: np.ndarray,
    s_eff: np.ndarray,
    v: np.ndarray,
    alpha: float = ALPHA_SIG,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Emergent respect hierarchy.

    R_resp is the v1.2.1 passive adherence closure.
    R_rob is the activity-sensitive robustness boundary.
    R_sig is the proposed active-significance field.
    """
    r_resp = np.power(h * b * v, 1.0 / 3.0)
    r_rob = np.minimum(r_resp, np.power(h * b * s_eff * v, 1.0 / 4.0))
    r_sig = np.power(r_resp, 1.0 - alpha) * np.power(s_eff, alpha)
    return clip01(r_resp), clip01(r_rob), clip01(r_sig)


def rankdata(x: np.ndarray) -> np.ndarray:
    """Simple average-rank implementation without scipy."""
    order = np.argsort(x)
    ranks = np.empty_like(order, dtype=float)
    sorted_x = x[order]
    n = len(x)
    i = 0
    while i < n:
        j = i + 1
        while j < n and sorted_x[j] == sorted_x[i]:
            j += 1
        ranks[order[i:j]] = 0.5 * (i + j - 1) + 1.0
        i = j
    return ranks


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx = rankdata(np.asarray(x))
    ry = rankdata(np.asarray(y))
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    return float(np.dot(rx, ry) / (np.linalg.norm(rx) * np.linalg.norm(ry) + EPS))


def archetype_table() -> list[dict]:
    archetypes = {
        "dormant_coherent": 0.15,
        "stable_star": 1.00,
        "chaotic_burst": 1.85,
    }
    rows: list[dict] = []
    for name, activity in archetypes.items():
        h, b, s_eff, v = primitive_supports(np.array([activity]))
        r_resp, r_rob, r_sig = closures(h, b, s_eff, v)
        rows.append(
            {
                "name": name,
                "activity": activity,
                "H": float(h[0]),
                "B": float(b[0]),
                "S_eff": float(s_eff[0]),
                "V": float(v[0]),
                "R_resp": float(r_resp[0]),
                "R_rob": float(r_rob[0]),
                "R_sig": float(r_sig[0]),
            }
        )
    return rows


def main() -> None:
    rng = np.random.default_rng(SEED)

    activity_grid = np.linspace(0.0, 2.2, 700)
    h, b, s_eff, v = primitive_supports(activity_grid)
    r_resp, r_rob, r_sig = closures(h, b, s_eff, v)

    # Ensemble with small support noise to test whether the hierarchy survives
    # beyond a single smooth curve.
    n = 6000
    activity = rng.uniform(0.0, 2.2, n)
    h_e, b_e, s_e, v_e = primitive_supports(activity)
    h_e = clip01(h_e + rng.normal(0.0, 0.035, n))
    b_e = clip01(b_e + rng.normal(0.0, 0.035, n))
    s_e = clip01(s_e + rng.normal(0.0, 0.025, n))
    v_e = clip01(v_e + rng.normal(0.0, 0.035, n))
    r_resp_e, r_rob_e, r_sig_e = closures(h_e, b_e, s_e, v_e)

    regimes = {
        "dormant": activity < 0.35,
        "controlled_active": (activity >= 0.75) & (activity <= 1.25),
        "overdriven": activity > 1.55,
    }

    regime_summary = {}
    for name, mask in regimes.items():
        regime_summary[name] = {
            "count": int(mask.sum()),
            "mean_activity": float(activity[mask].mean()),
            "mean_R_resp": float(r_resp_e[mask].mean()),
            "mean_R_rob": float(r_rob_e[mask].mean()),
            "mean_R_sig": float(r_sig_e[mask].mean()),
            "mean_S_eff": float(s_e[mask].mean()),
        }

    high_resp_low_sig = (r_resp_e > 0.60) & (r_sig_e < 0.25)
    top_sig = r_sig_e >= np.quantile(r_sig_e, 0.90)

    summary = {
        "model": "Active Respect / Structural Significance Toy Simulation",
        "seed": SEED,
        "n_ensemble": n,
        "definitions": {
            "R_resp": "(H B V)^(1/3)",
            "R_rob": "min(R_resp, (H B S_eff V)^(1/4))",
            "R_sig": "R_resp^(1-alpha) S_eff^alpha",
            "S_eff": "exp(-0.5*((A-A_star)/sigma_A)^2)",
        },
        "parameters": {
            "A_star": A_STAR,
            "sigma_A": SIGMA_A,
            "alpha_sig": ALPHA_SIG,
        },
        "curve_optima": {
            "A_at_max_R_resp": float(activity_grid[np.argmax(r_resp)]),
            "max_R_resp": float(np.max(r_resp)),
            "A_at_max_R_rob": float(activity_grid[np.argmax(r_rob)]),
            "max_R_rob": float(np.max(r_rob)),
            "A_at_max_R_sig": float(activity_grid[np.argmax(r_sig)]),
            "max_R_sig": float(np.max(r_sig)),
        },
        "regime_summary": regime_summary,
        "archetypes": archetype_table(),
        "ensemble_tests": {
            "fraction_high_R_resp_low_R_sig": float(high_resp_low_sig.mean()),
            "top_10pct_R_sig_mean_activity": float(activity[top_sig].mean()),
            "top_10pct_R_sig_std_activity": float(activity[top_sig].std()),
            "spearman_R_sig_activity": spearman(r_sig_e, activity),
            "spearman_R_sig_S_eff": spearman(r_sig_e, s_e),
            "spearman_R_sig_R_resp": spearman(r_sig_e, r_resp_e),
        },
        "interpretation": (
            "The simulation supports treating R as an emergent hierarchy: "
            "passive adherence can be nonzero in quiet systems, but active "
            "significance peaks only when controlled activity is present."
        ),
    }

    with open(OUT / "active_respect_significance_results.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    csv = np.column_stack([activity_grid, h, b, s_eff, v, r_resp, r_rob, r_sig])
    np.savetxt(
        OUT / "active_respect_curves.csv",
        csv,
        delimiter=",",
        header="activity,H,B,S_eff,V,R_resp,R_rob,R_sig",
        comments="",
    )

    ensemble_csv = np.column_stack([activity, h_e, b_e, s_e, v_e, r_resp_e, r_rob_e, r_sig_e])
    np.savetxt(
        OUT / "active_respect_ensemble.csv",
        ensemble_csv,
        delimiter=",",
        header="activity,H,B,S_eff,V,R_resp,R_rob,R_sig",
        comments="",
    )

    plt.rcParams.update({"font.size": 11, "figure.dpi": 150, "savefig.dpi": 170})

    fig, ax = plt.subplots(figsize=(11.5, 6.4))
    ax.plot(activity_grid, h, label="H coherence", lw=2)
    ax.plot(activity_grid, b, label="B balance", lw=2)
    ax.plot(activity_grid, s_eff, label=r"$S_{\rm eff}$ controlled activity", lw=2)
    ax.plot(activity_grid, v, label="V coupling", lw=2)
    ax.plot(activity_grid, r_resp, label=r"$R_{\rm resp}$", lw=3, ls="--")
    ax.plot(activity_grid, r_rob, label=r"$R_{\rm rob}$", lw=3, ls="-.")
    ax.plot(activity_grid, r_sig, label=r"$R_{\rm sig}$", lw=3, color="black")
    ax.axvline(A_STAR, color="black", lw=1, ls=":", alpha=0.8)
    ax.text(A_STAR + 0.02, 0.05, r"$A_\star$", fontsize=12)
    ax.set_xlabel("raw activity A")
    ax.set_ylabel("support / emergent field value")
    ax.set_title("Emergent respect hierarchy from primitive H, B, S, V")
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, frameon=True)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_respect_hierarchy_vs_activity.png", bbox_inches="tight")
    plt.close(fig)

    alpha_grid = np.linspace(0.05, 0.95, 180)
    phase = np.array([
        closures(h, b, s_eff, v, alpha=float(alpha))[2] for alpha in alpha_grid
    ])
    fig, ax = plt.subplots(figsize=(11.5, 6.2))
    im = ax.imshow(
        phase,
        origin="lower",
        aspect="auto",
        extent=[activity_grid.min(), activity_grid.max(), alpha_grid.min(), alpha_grid.max()],
        cmap="viridis",
    )
    ax.axvline(A_STAR, color="white", ls=":", lw=1.5)
    ax.set_xlabel("raw activity A")
    ax.set_ylabel(r"activity weight $\alpha$")
    ax.set_title(r"Phase map of active significance $R_{\rm sig}$")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$R_{\rm sig}$")
    fig.tight_layout()
    fig.savefig(OUT / "fig2_rsig_phase_map.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.5, 7.2))
    sc = ax.scatter(r_resp_e, r_sig_e, c=activity, cmap="plasma", s=12, alpha=0.55, edgecolors="none")
    rows = archetype_table()
    for row in rows:
        ax.scatter(row["R_resp"], row["R_sig"], s=160, marker="*", edgecolor="black", linewidth=1.1)
        ax.text(row["R_resp"] + 0.015, row["R_sig"] + 0.015, row["name"].replace("_", " "), fontsize=10)
    ax.set_xlabel(r"passive adherence $R_{\rm resp}$")
    ax.set_ylabel(r"active significance $R_{\rm sig}$")
    ax.set_title("Passive respect is not the same as active significance")
    ax.grid(alpha=0.25)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("raw activity A")
    fig.tight_layout()
    fig.savefig(OUT / "fig3_passive_vs_active_significance.png", bbox_inches="tight")
    plt.close(fig)

    labels = [row["name"].replace("_", "\n") for row in rows]
    metrics = ["H", "B", "S_eff", "V", "R_resp", "R_rob", "R_sig"]
    values = np.array([[row[m] for m in metrics] for row in rows])
    x = np.arange(len(labels))
    width = 0.11
    fig, ax = plt.subplots(figsize=(12.2, 6.4))
    for i, metric in enumerate(metrics):
        ax.bar(x + (i - 3) * width, values[:, i], width=width, label=metric)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("value")
    ax.set_title("Archetype comparison: dormant coherence, stable star, chaotic burst")
    ax.grid(axis="y", alpha=0.22)
    ax.legend(ncol=4)
    fig.tight_layout()
    fig.savefig(OUT / "fig4_archetype_supports.png", bbox_inches="tight")
    plt.close(fig)

    print("=" * 72)
    print("Active Respect / Structural Significance Toy Simulation")
    print("=" * 72)
    print(f"A at max R_sig: {summary['curve_optima']['A_at_max_R_sig']:.4f}")
    print(f"max R_sig: {summary['curve_optima']['max_R_sig']:.4f}")
    print(f"fraction high R_resp but low R_sig: {summary['ensemble_tests']['fraction_high_R_resp_low_R_sig']:.4f}")
    print(f"top 10% R_sig mean activity: {summary['ensemble_tests']['top_10pct_R_sig_mean_activity']:.4f}")
    print(f"outputs written to: {OUT.resolve()}")


if __name__ == "__main__":
    main()
