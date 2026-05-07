#!/usr/bin/env python3
"""
Societal Critical Coherence Index Toy Framework
===============================================

This script builds a synthetic, ethically bounded toy benchmark for applying
the MAAT/CCI logic to social-system states.

It does NOT rank people, parties, countries, or real communities.
All inputs are synthetic archetypes chosen to illustrate regime structure:

    - stagnant order
    - constructive reform
    - polarized mobilization
    - fragmented activism
    - authoritarian stability
    - creative democratic renewal

The purpose is to test whether a societal CCI and active significance index can
distinguish raw activity, controlled transformation, robustness, and connected
social coherence.
"""

from __future__ import annotations

import json
import tempfile
from dataclasses import dataclass, asdict
from pathlib import Path
import os

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "maat_societal_cci_matplotlib"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT = BASE_DIR / "outputs"
OUT.mkdir(exist_ok=True)

EPS = 1e-9
A_STAR = 1.0
SIGMA_A = 0.30
ALPHA_SIG = 0.45


@dataclass(frozen=True)
class Archetype:
    name: str
    H: float
    B: float
    A_raw: float
    V: float
    description: str


ARCHETYPES = [
    Archetype(
        "stagnant_order",
        H=0.82,
        B=0.76,
        A_raw=0.18,
        V=0.70,
        description="Stable but low-transformational system with limited civic activity.",
    ),
    Archetype(
        "constructive_reform",
        H=0.78,
        B=0.74,
        A_raw=0.96,
        V=0.76,
        description="Moderate-to-high activity integrated by coherence, balance, and trust.",
    ),
    Archetype(
        "polarized_mobilization",
        H=0.42,
        B=0.28,
        A_raw=1.42,
        V=0.35,
        description="High mobilization with weak balance and low shared coherence.",
    ),
    Archetype(
        "fragmented_activism",
        H=0.38,
        B=0.46,
        A_raw=1.10,
        V=0.22,
        description="Activity near the constructive window but poor connectedness.",
    ),
    Archetype(
        "authoritarian_stability",
        H=0.70,
        B=0.32,
        A_raw=0.24,
        V=0.30,
        description="Low activity and restricted connectedness with apparent order.",
    ),
    Archetype(
        "creative_democratic_renewal",
        H=0.84,
        B=0.81,
        A_raw=1.02,
        V=0.86,
        description="High coherence, balance, trust, and controlled transformation.",
    ),
]


def controlled_activity(a_raw: np.ndarray | float) -> np.ndarray | float:
    a_raw = np.asarray(a_raw)
    return np.exp(-0.5 * ((a_raw - A_STAR) / SIGMA_A) ** 2)


def compute_metrics(H: float, B: float, A_raw: float, V: float) -> dict:
    S_eff = float(controlled_activity(A_raw))
    R_resp = float((H * B * V) ** (1.0 / 3.0))
    R_rob = float(min(R_resp, (H * B * S_eff * V) ** (1.0 / 4.0)))
    R_sig = float((R_resp ** (1.0 - ALPHA_SIG)) * (S_eff ** ALPHA_SIG))
    CCI_soc = float(A_raw / (H + B + V + R_rob + EPS))
    stress_soc = CCI_soc
    ASI_soc = float(R_sig / (1.0 + CCI_soc))
    stagnation = float((1.0 - S_eff) * (H + B + V) / 3.0)
    fragmentation = float(A_raw * (1.0 - V) * (1.0 - B))
    return {
        "S_eff": S_eff,
        "R_resp_soc": R_resp,
        "R_rob_soc": R_rob,
        "R_sig_soc": R_sig,
        "CCI_soc": CCI_soc,
        "Stress_soc": stress_soc,
        "ASI_soc": ASI_soc,
        "Stagnation_proxy": stagnation,
        "Fragmentation_proxy": fragmentation,
    }


def classify(row: dict) -> str:
    A = row["A_raw"]
    H, B, V = row["H"], row["B"], row["V"]
    R_sig, CCI = row["R_sig_soc"], row["CCI_soc"]
    if R_sig > 0.75 and CCI < 0.45:
        return "constructive transformation"
    if A < 0.40 and min(H, B, V) > 0.60:
        return "stable but stagnant"
    if A > 1.20 and B < 0.45:
        return "polarization risk"
    if A > 0.80 and V < 0.35:
        return "fragmentation risk"
    if A < 0.45 and (B < 0.45 or V < 0.45):
        return "constrained stability"
    return "mixed regime"


def archetype_rows() -> list[dict]:
    rows = []
    for a in ARCHETYPES:
        row = asdict(a)
        row.update(compute_metrics(a.H, a.B, a.A_raw, a.V))
        row["regime"] = classify(row)
        rows.append(row)
    return rows


def random_ensemble(n: int = 5000, seed: int = 44) -> list[dict]:
    rng = np.random.default_rng(seed)
    H = rng.beta(3.0, 2.0, n)
    B = rng.beta(2.4, 2.4, n)
    V = rng.beta(2.6, 2.2, n)
    A_raw = rng.uniform(0.0, 1.8, n)
    rows = []
    for i in range(n):
        row = {
            "H": float(H[i]),
            "B": float(B[i]),
            "A_raw": float(A_raw[i]),
            "V": float(V[i]),
        }
        row.update(compute_metrics(row["H"], row["B"], row["A_raw"], row["V"]))
        row["regime"] = classify(row)
        rows.append(row)
    return rows


def save_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(",".join(fields) + "\n")
        for row in rows:
            values = []
            for field in fields:
                value = row[field]
                if isinstance(value, str):
                    values.append('"' + value.replace('"', '""') + '"')
                else:
                    values.append(f"{float(value):.10g}")
            f.write(",".join(values) + "\n")


def main() -> None:
    archetypes = archetype_rows()
    ensemble = random_ensemble()

    best_asi = max(archetypes, key=lambda r: r["ASI_soc"])
    highest_stress = max(archetypes, key=lambda r: r["CCI_soc"])
    highest_rsig = max(archetypes, key=lambda r: r["R_sig_soc"])

    regime_counts: dict[str, int] = {}
    for row in ensemble:
        regime_counts[row["regime"]] = regime_counts.get(row["regime"], 0) + 1

    summary = {
        "model": "Societal Critical Coherence Index Toy Framework",
        "status": "synthetic toy model; no real people, parties, countries, or communities are evaluated",
        "parameters": {
            "A_star": A_STAR,
            "sigma_A": SIGMA_A,
            "alpha_sig": ALPHA_SIG,
            "epsilon": EPS,
        },
        "definitions": {
            "S_eff": "exp[-0.5 ((A_raw - A_star)/sigma_A)^2]",
            "R_resp_soc": "(H B V)^(1/3)",
            "R_rob_soc": "min(R_resp_soc, (H B S_eff V)^(1/4))",
            "R_sig_soc": "R_resp_soc^(1-alpha) S_eff^alpha",
            "CCI_soc": "A_raw / (H + B + V + R_rob_soc + epsilon)",
            "ASI_soc": "R_sig_soc / (1 + CCI_soc)",
        },
        "archetypes": archetypes,
        "key_results": {
            "best_ASI_archetype": best_asi["name"],
            "best_ASI_value": best_asi["ASI_soc"],
            "highest_CCI_archetype": highest_stress["name"],
            "highest_CCI_value": highest_stress["CCI_soc"],
            "highest_R_sig_archetype": highest_rsig["name"],
            "highest_R_sig_value": highest_rsig["R_sig_soc"],
            "constructive_reform_ASI": next(r for r in archetypes if r["name"] == "constructive_reform")["ASI_soc"],
            "creative_democratic_renewal_ASI": next(r for r in archetypes if r["name"] == "creative_democratic_renewal")["ASI_soc"],
            "polarized_mobilization_CCI": next(r for r in archetypes if r["name"] == "polarized_mobilization")["CCI_soc"],
            "authoritarian_stability_R_sig": next(r for r in archetypes if r["name"] == "authoritarian_stability")["R_sig_soc"],
        },
        "ensemble": {
            "n": len(ensemble),
            "seed": 44,
            "regime_counts": regime_counts,
            "mean_CCI_soc": float(np.mean([r["CCI_soc"] for r in ensemble])),
            "mean_ASI_soc": float(np.mean([r["ASI_soc"] for r in ensemble])),
            "corr_CCI_ASI": float(np.corrcoef([r["CCI_soc"] for r in ensemble], [r["ASI_soc"] for r in ensemble])[0, 1]),
        },
    }

    with open(OUT / "societal_cci_results.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    fields = [
        "name",
        "H",
        "B",
        "A_raw",
        "S_eff",
        "V",
        "R_resp_soc",
        "R_rob_soc",
        "R_sig_soc",
        "CCI_soc",
        "ASI_soc",
        "regime",
        "description",
    ]
    save_csv(OUT / "societal_cci_archetypes.csv", archetypes, fields)

    ensemble_fields = [
        "H",
        "B",
        "A_raw",
        "S_eff",
        "V",
        "R_resp_soc",
        "R_rob_soc",
        "R_sig_soc",
        "CCI_soc",
        "ASI_soc",
        "regime",
    ]
    save_csv(OUT / "societal_cci_ensemble.csv", ensemble, ensemble_fields)

    plt.rcParams.update({"font.size": 11, "figure.dpi": 150, "savefig.dpi": 170})

    names = [r["name"].replace("_", "\n") for r in archetypes]
    x = np.arange(len(names))

    fig, ax = plt.subplots(figsize=(12.5, 6.4))
    width = 0.16
    metrics = ["CCI_soc", "R_rob_soc", "R_sig_soc", "ASI_soc"]
    labels = [r"$CCI_{\rm soc}$", r"$R_{\rm rob,soc}$", r"$R_{\rm sig,soc}$", r"$ASI_{\rm soc}$"]
    for i, (metric, label) in enumerate(zip(metrics, labels)):
        ax.bar(x + (i - 1.5) * width, [r[metric] for r in archetypes], width=width, label=label)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylim(0, 1.25)
    ax.set_ylabel("index value")
    ax.set_title("Societal structural stress and active significance archetypes")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(ncol=4)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_societal_indices_by_archetype.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 7))
    scatter = ax.scatter(
        [r["CCI_soc"] for r in ensemble],
        [r["ASI_soc"] for r in ensemble],
        c=[r["R_sig_soc"] for r in ensemble],
        s=13,
        alpha=0.45,
        cmap="viridis",
        edgecolors="none",
    )
    for r in archetypes:
        ax.scatter(r["CCI_soc"], r["ASI_soc"], s=150, marker="*", edgecolor="black", linewidth=1.0)
        ax.text(r["CCI_soc"] + 0.01, r["ASI_soc"] + 0.01, r["name"].replace("_", " "), fontsize=9)
    ax.set_xlabel(r"$CCI_{\rm soc}$ structural stress")
    ax.set_ylabel(r"$ASI_{\rm soc}$ constructive significance")
    ax.set_title("Stress is not the same as constructive transformation")
    ax.grid(alpha=0.25)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label(r"$R_{\rm sig,soc}$")
    fig.tight_layout()
    fig.savefig(OUT / "fig2_stress_vs_significance_phase_space.png", bbox_inches="tight")
    plt.close(fig)

    a_grid = np.linspace(0, 1.8, 500)
    s_grid = controlled_activity(a_grid)
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    ax.plot(a_grid, s_grid, lw=3, color="black")
    ax.axvline(A_STAR, color="tab:green", ls=":", lw=2, label=r"$A_\star$")
    ax.fill_between(a_grid, 0, s_grid, where=s_grid > 0.70, color="tab:green", alpha=0.18, label="controlled window")
    ax.set_xlabel("raw social activity / mobilization")
    ax.set_ylabel(r"$S_{\rm eff}$")
    ax.set_title("Controlled activity: constructive transformation is not maximal activity")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "fig3_controlled_activity_window.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11.5, 6.2))
    regimes = list(regime_counts)
    counts = [regime_counts[k] for k in regimes]
    ax.bar(regimes, counts, color="steelblue")
    ax.set_ylabel("synthetic ensemble count")
    ax.set_title("Regime classification in synthetic societal ensemble")
    ax.tick_params(axis="x", rotation=25)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT / "fig4_synthetic_regime_counts.png", bbox_inches="tight")
    plt.close(fig)

    print("=" * 72)
    print("Societal Critical Coherence Index Toy Framework")
    print("=" * 72)
    print(f"best ASI archetype: {best_asi['name']} ({best_asi['ASI_soc']:.4f})")
    print(f"highest stress archetype: {highest_stress['name']} ({highest_stress['CCI_soc']:.4f})")
    print(f"highest R_sig archetype: {highest_rsig['name']} ({highest_rsig['R_sig_soc']:.4f})")
    print(f"outputs written to: {OUT.resolve()}")


if __name__ == "__main__":
    main()
