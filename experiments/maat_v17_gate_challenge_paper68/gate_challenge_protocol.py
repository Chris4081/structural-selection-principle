#!/usr/bin/env python3
"""Paper 68: predeclared MAAT v1.7 Gate Challenge protocol.

This script does not run external benchmarks.  It writes the preregistered
gate architecture, domain matrix, baseline matrix, metrics, and failure rules
that future external validation runs must follow.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs_paper68_gate_challenge_protocol"


SPEC = {
    "paper": 68,
    "title": "The Gate Challenge",
    "subtitle": "A Predeclared Falsification Protocol for MAAT v1.7 Gate-Aware Structural Selection",
    "year": 2026,
    "status": "protocol_only_no_external_results",
    "default_assumption": "MAAT v1.7 is false until supported by predeclared external utility tests.",
    "gate_architecture": {
        "name": "polarity_normalized_response_gate",
        "candidate": "X",
        "supports": ["H", "B", "S", "V"],
        "robustness": {
            "R_resp": "(H*B*V)^(1/3)",
            "R_rob": "min(R_resp, (H*B*S*V)^(1/4))",
        },
        "calibration": {
            "split": "calibration fold only",
            "center": "median(R_rob)",
            "scale": "MAD(R_rob) + eps",
            "eps": 1e-9,
            "no_test_set_refit": True,
        },
        "polarity": {
            "+1": "high robustness activates selection/preservation modes",
            "-1": "low robustness activates warning/refinement modes",
            "predeclared_per_domain": True,
        },
        "response_terms": {
            "z_R": "p_D * (R_rob - median_cal(R_rob)) / (MAD_cal(R_rob) + eps)",
            "z_t": "p_D * d_t R_rob standardized on calibration fold",
            "z_x": "p_D * grad R_rob standardized on calibration fold",
            "a_t": 0.50,
            "a_x": 0.25,
        },
        "gate_signal": "g = z_R + a_t*z_t + a_x*z_x",
        "soft_gate": "G = 1 / (1 + exp(-g/tau))",
        "hard_gate": "G_hard = 1[G >= 0.5]",
        "tau": 1.0,
        "score_comparison": {
            "score_only": "F_score(H,B,S,V)",
            "score_with_R": "F_score(H,B,S,V) weighted by continuous R_rob",
            "gate_with_R": "G * F_score(H,B,S,V) or gate-triggered action",
        },
    },
    "success_criteria": {
        "minimum_support": (
            "Gate beats score_with_R in at least 2 of 3 primary external domains "
            "with paired bootstrap 95% CI lower bound > 0, and no primary domain "
            "shows a CI upper bound < 0 against score_with_R."
        ),
        "strong_support": (
            "Gate beats both score_with_R and the strongest classical baseline in "
            "all 3 primary external domains with paired bootstrap 95% CI lower bound > 0."
        ),
        "failure": (
            "Gate fails if it does not beat score_with_R in at least 2 primary domains, "
            "or if improvements vanish under shuffled/null controls, or if gains require "
            "post-hoc threshold tuning."
        ),
        "v2_claim_rule": "No v2.0 claim is allowed unless strong_support is met on external benchmarks.",
    },
    "null_controls": [
        "shuffled_R_rob_within_domain",
        "shuffled_gate_signal_preserving_marginal_distribution",
        "random_polarity_assignment",
        "score_label_shuffle_where_applicable",
    ],
    "ablation_rules": [
        "remove_response_terms_use_z_R_only",
        "remove_R_value_use_response_terms_only",
        "soft_gate_vs_hard_gate",
        "score_only_vs_score_with_R_vs_gate_with_R",
    ],
}


DOMAINS = [
    {
        "domain": "Fluid / JHTDB turbulence",
        "external_status": "primary_external",
        "polarity": -1,
        "decision": "trigger adaptive refinement or early warning",
        "primary_metric": "lead_coverage_utility_at_matched_false_alarm_rate",
        "secondary_metrics": "detection_rate, median_lead_time, false_alarm_rate, bootstrap_CI",
        "classical_baselines": "high-k monitor, palinstrophy, max vorticity, RMS vorticity",
        "MAAT_score": "continuous W_MAAT score",
        "MAAT_gate": "polarity_normalized_response_gate on R_rob",
        "external_data_requirement": "JHTDB public DNS sub-cubes or equivalent public turbulence benchmark",
    },
    {
        "domain": "SAT / CDCL",
        "external_status": "primary_external",
        "polarity": 1,
        "decision": "activate structural branching mode or fall back to classical heuristic",
        "primary_metric": "paired_compute_to_solution_regret",
        "secondary_metrics": "conflicts, decisions, solved_within_budget, PAR2, bootstrap_CI",
        "classical_baselines": "VSIDS, LRB, MOMS, Jeroslow-Wang",
        "MAAT_score": "continuous Mode-A structural score",
        "MAAT_gate": "polarity_normalized_response_gate deciding Mode-A activation",
        "external_data_requirement": "SAT Competition, SATLIB, or other public CNF families",
    },
    {
        "domain": "Open quantum systems",
        "external_status": "primary_external_or_public_recipe",
        "polarity": 1,
        "decision": "select/predict robust pointer-state or entanglement-preserving candidate",
        "primary_metric": "ranking_utility_or_top_k_enrichment_against_final_robustness",
        "secondary_metrics": "Spearman, AUROC, calibration, bootstrap_CI",
        "classical_baselines": "decoherence rate, purity, fidelity, stationarity baseline",
        "MAAT_score": "continuous HBSV or field score",
        "MAAT_gate": "polarity_normalized_response_gate on correlated stationarity robustness",
        "external_data_requirement": "public Lindblad benchmark recipes or independently specified open-system ensemble",
    },
    {
        "domain": "Optional public ML hardness",
        "external_status": "secondary_external",
        "polarity": -1,
        "decision": "flag hard or unstable samples",
        "primary_metric": "sample_hardness_ranking_utility",
        "secondary_metrics": "Spearman, AUROC, top-k enrichment, bootstrap_CI",
        "classical_baselines": "local label disagreement, kNN density, centroid distance, PCA error",
        "MAAT_score": "sample-level structural score",
        "MAAT_gate": "polarity_normalized_response_gate",
        "external_data_requirement": "public sklearn/OpenML datasets",
    },
]


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("No rows to write")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def make_protocol_figure(path: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - fallback for minimal envs
        fallback = path.with_suffix(".txt")
        fallback.write_text(f"matplotlib unavailable: {exc}\n", encoding="utf-8")
        return

    labels = [row["domain"].split(" / ")[0] for row in DOMAINS[:3]]
    polarity = [row["polarity"] for row in DOMAINS[:3]]
    colors = ["#2a9d8f" if p > 0 else "#e76f51" for p in polarity]

    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    bars = ax.bar(labels, [abs(p) for p in polarity], color=colors, edgecolor="#1f2933")
    ax.set_ylim(0, 1.35)
    ax.set_ylabel("predeclared gate active")
    ax.set_title("Paper 68 Gate Challenge: same gate architecture, predeclared task polarity")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["off", "on"])
    for bar, p in zip(bars, polarity):
        text = "high-R selects" if p > 0 else "low-R warns"
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            1.05,
            text,
            ha="center",
            va="bottom",
            fontsize=10,
        )
    ax.text(
        0.5,
        -0.22,
        "Thresholds are calibrated only from calibration folds; no test-set retuning.",
        transform=ax.transAxes,
        ha="center",
        fontsize=9,
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    spec_path = OUT / "gate_challenge_preregistration.json"
    spec_path.write_text(json.dumps(SPEC, indent=2, sort_keys=True), encoding="utf-8")

    write_csv(OUT / "gate_challenge_domain_matrix.csv", DOMAINS)

    baseline_rows = []
    for row in DOMAINS:
        for baseline in [b.strip() for b in row["classical_baselines"].split(",")]:
            baseline_rows.append(
                {
                    "domain": row["domain"],
                    "baseline": baseline,
                    "role": "strong_reasonable_baseline",
                    "comparison": "must be compared against MAAT_score and MAAT_gate",
                }
            )
    write_csv(OUT / "gate_challenge_baseline_matrix.csv", baseline_rows)

    metric_rows = [
        {
            "metric": "delta_utility_gate_vs_score",
            "definition": "U(MAAT_gate) - U(MAAT_score_with_R)",
            "primary": True,
        },
        {
            "metric": "delta_utility_gate_vs_best_classical",
            "definition": "U(MAAT_gate) - U(best_classical_baseline)",
            "primary": True,
        },
        {
            "metric": "paired_bootstrap_CI_95",
            "definition": "paired bootstrap confidence interval for utility deltas",
            "primary": True,
        },
        {
            "metric": "null_survival",
            "definition": "real gate improvement exceeds shuffled/null controls",
            "primary": True,
        },
    ]
    write_csv(OUT / "gate_challenge_metric_registry.csv", metric_rows)

    make_protocol_figure(OUT / "fig1_gate_challenge_protocol.png")

    readme = f"""# Paper 68 Gate Challenge Protocol

This folder contains the preregistered protocol artifacts for Paper 68.
It does not run the external benchmarks.  It freezes the gate equation,
calibration rule, domains, baselines, metrics, null controls, and failure
criteria before future external validation.

Run:

```bash
python3 gate_challenge_protocol.py
```

Outputs are written to:

```text
{OUT.name}/
```

Key files:

- `gate_challenge_preregistration.json`
- `gate_challenge_domain_matrix.csv`
- `gate_challenge_baseline_matrix.csv`
- `gate_challenge_metric_registry.csv`
- `fig1_gate_challenge_protocol.png`

Scientific status: protocol only, no external validation result.
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")

    summary = {
        "output_dir": str(OUT),
        "files": sorted(p.name for p in OUT.iterdir()),
        "domains": [row["domain"] for row in DOMAINS],
        "success_criteria": SPEC["success_criteria"],
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
