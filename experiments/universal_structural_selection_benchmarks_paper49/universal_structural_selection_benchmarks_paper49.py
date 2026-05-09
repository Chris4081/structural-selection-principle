#!/usr/bin/env python3
"""
Paper 49: Universal Structural Selection Benchmarks.

This script is a reproducibility aggregator, not a new physics solver.  It
collects key metrics from existing MAAT/Structural Selection experiments and
turns them into a benchmark-readiness registry:

    complete = an explicit competitor/null test exists
    partial  = a diagnostic or indirect competitor exists
    missing  = the benchmark class is still absent
    n/a      = the competitor class is not natural for the domain

Run from this folder:

    python3 universal_structural_selection_benchmarks_paper49.py
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]
OUTDIR = SCRIPT_PATH.parent / "outputs_paper49"


def read_json(rel_path: str) -> dict[str, Any]:
    path = REPO_ROOT / rel_path
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def status_value(status: str) -> float:
    return {
        "missing": 0.0,
        "partial": 1.0,
        "complete": 2.0,
        "n/a": np.nan,
    }[status]


@dataclass
class DomainEvidence:
    domain: str
    benchmark: str
    evidence_status: str
    main_competitor: str
    primary_metric: str
    primary_value: str
    interpretation: str
    source: str


def fixed_energy_evidence() -> DomainEvidence:
    data = read_json(
        "experiments/fixed_energy_structural_selection/"
        "fixed_energy_structural_selection_results.json"
    )
    eq = data["equal_energy_phi4_1d"]
    results = eq["results"]
    perturbed = [k for k in results if k != "kink_reference"]
    energies = np.array([results[k]["Energy"] for k in perturbed], dtype=float)
    f_scores = {k: float(results[k]["F_MAAT"]) for k in perturbed}
    f_best = min(f_scores.values())
    f_worst = max(f_scores.values())
    spread_pct = 100.0 * (float(np.max(energies)) - float(np.min(energies))) / float(eq["target_energy"])
    ratio = f_worst / f_best
    return DomainEvidence(
        domain="nonlinear_fields",
        benchmark="fixed-energy phi4 / sine-Gordon structural ranking",
        evidence_status="complete",
        main_competitor="energy-only",
        primary_metric="equal-energy F_MAAT separation",
        primary_value=f"F_worst/F_best={ratio:.3f} at energy spread={spread_pct:.4f}%",
        interpretation=(
            "Near-equal-energy configurations are separated strongly by "
            "structural score, so the ranking is not reducible to total energy."
        ),
        source="experiments/fixed_energy_structural_selection/",
    )


def string_path_evidence() -> DomainEvidence:
    data = read_json(
        "experiments/string_landscape_path_selection/outputs_path_selection/"
        "path_selection_summary.json"
    )
    overlap = data["edge_statistics"]["top20_overlap_structural_energy"]
    rho = data["spearman_abs_energy_vs_F_bridge_minima"]
    cost_struct = data["best_energy_to_best_bridge_structural_path"]["cost"]
    cost_energy = data["best_energy_to_best_bridge_energy_action_path"]["cost"]
    return DomainEvidence(
        domain="string_landscape",
        benchmark="reduced KKLT vacuum/path selection",
        evidence_status="complete",
        main_competitor="energy-only path action",
        primary_metric="top-20 edge overlap and path cost",
        primary_value=(
            f"top20 overlap={overlap}/20, "
            f"rho(|E|,F_bridge)={rho:.4f}, "
            f"A_struct={cost_struct:.6f}, A_energy={cost_energy:.6f}"
        ),
        interpretation=(
            "Structural flow and energy flow select disjoint top transition "
            "edges and different path costs in the reduced KKLT graph."
        ),
        source="experiments/string_landscape_path_selection/",
    )


def cosmology_evidence() -> DomainEvidence:
    data = read_json(
        "experiments/maat_paper42_blind_projection_test/outputs_paper42/"
        "paper42_summary.json"
    )
    null = data["main_null_tests"]["CCI_diag_vs_abs_residual_sigma"]
    rho = null["spearman"]
    p_perm = null["p_perm"]
    return DomainEvidence(
        domain="cosmology",
        benchmark="blind projection residual test",
        evidence_status="partial",
        main_competitor="redshift-shuffled null / LCDM residual baseline",
        primary_metric="blind Spearman residual correlation",
        primary_value=f"rho={rho:.4f}, permutation p={p_perm:.4f}",
        interpretation=(
            "The blind projection test provides a reproducible null protocol, "
            "but the current signal is suggestive rather than a detection."
        ),
        source="experiments/maat_paper42_blind_projection_test/",
    )


def sat_evidence() -> DomainEvidence:
    data = read_json(
        "experiments/maat_v121_observables_stability_paper37/sat_validation/"
        "maat_v121_sat_summary.json"
    )
    base = data["base_cv_r2_mean"]
    v121 = data["v121_cv_r2_mean"]
    delta = data["delta_cv_r2"]
    return DomainEvidence(
        domain="sat_constraints",
        benchmark="SAT/runtime validation using v1.2.1 fields",
        evidence_status="partial",
        main_competitor="base field regression",
        primary_metric="cross-validated R2 change",
        primary_value=f"base R2={base:.4f}, v1.2.1 R2={v121:.4f}, delta={delta:.4f}",
        interpretation=(
            "The existing SAT validation is a stress test rather than a success: "
            "the reported CV performance is poor and marks SAT as an open "
            "universality target."
        ),
        source=(
            "experiments/maat_v121_observables_stability_paper37/sat_validation/"
        ),
    )


def safety_evidence() -> DomainEvidence:
    data = read_json(
        "experiments/boundary_aware_lambda_calibration/"
        "closed_maat_lambda_fit_results.json"
    )
    r_share = data["shares"]["R"]
    lam_r = data["lambdas"]["R"]
    loss = data["loss"]
    return DomainEvidence(
        domain="ai_safety_boundary",
        benchmark="boundary-aware response-weight calibration",
        evidence_status="partial",
        main_competitor="unconstrained / non-boundary lambda balance",
        primary_metric="robustness weight share",
        primary_value=f"R share={r_share:.4f}, lambda_R={lam_r:.4f}, loss={loss:.6f}",
        interpretation=(
            "Constraint-boundary data drive robustness dominance, but this is "
            "not yet an external AI-system generalization benchmark."
        ),
        source="experiments/boundary_aware_lambda_calibration/",
    )


def sm_bridge_evidence() -> DomainEvidence:
    data = read_json(
        "experiments/standard_model_bridge/"
        "standard_model_rg_maat_v11_holdout_results.json"
    )
    rows = data["summary_rows"]
    holdouts = [r for r in rows if r["held_out"] != "none"]
    mean_log_error = float(np.mean([r["log10_error_abs"] for r in holdouts]))
    max_factor = float(max(max(r["factor"], 1.0 / r["factor"]) for r in holdouts))
    return DomainEvidence(
        domain="sm_like_constants",
        benchmark="SM-like RG bridge direct-term holdout",
        evidence_status="partial",
        main_competitor="random UV ensemble / direct-term ablation",
        primary_metric="holdout predictive error",
        primary_value=f"mean |log10 error|={mean_log_error:.4f}, worst factor={max_factor:.3f}",
        interpretation=(
            "Holdouts retain basin-level proximity for some observables but "
            "are not a precision derivation of constants."
        ),
        source="experiments/standard_model_bridge/",
    )


def quantum_evidence() -> DomainEvidence:
    return DomainEvidence(
        domain="quantum_measurement",
        benchmark="pointer-state / decoherence branch selection",
        evidence_status="missing",
        main_competitor="decoherence-only / stability-only pointer ranking",
        primary_metric="not implemented",
        primary_value="missing",
        interpretation=(
            "A quantum measurement benchmark is required for universality but "
            "is not yet present in the repository."
        ),
        source="future benchmark",
    )


def build_status_matrix() -> tuple[list[str], list[str], np.ndarray]:
    domains = [
        "fields",
        "string",
        "cosmology",
        "SAT",
        "AI/safety",
        "SM bridge",
        "quantum",
    ]
    tests = [
        "energy-only",
        "entropy/activity",
        "stability-only",
        "random/null",
        "domain baseline",
        "cross-domain",
    ]
    matrix_labels = [
        ["complete", "partial", "partial", "missing", "complete", "partial"],
        ["complete", "missing", "partial", "missing", "complete", "partial"],
        ["n/a", "missing", "partial", "complete", "complete", "partial"],
        ["n/a", "missing", "partial", "missing", "partial", "missing"],
        ["n/a", "missing", "complete", "missing", "partial", "missing"],
        ["n/a", "missing", "partial", "complete", "partial", "partial"],
        ["missing", "missing", "missing", "missing", "missing", "missing"],
    ]
    matrix = np.array([[status_value(v) for v in row] for row in matrix_labels], dtype=float)
    return domains, tests, matrix


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_status_matrix(domains: list[str], tests: list[str], matrix: np.ndarray) -> None:
    masked = np.ma.masked_invalid(matrix)
    cmap = plt.colormaps.get_cmap("viridis").resampled(3).copy()
    cmap.set_bad(color="#d9d9d9")
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    im = ax.imshow(masked, vmin=0, vmax=2, cmap=cmap, aspect="auto")
    ax.set_xticks(range(len(tests)), tests, rotation=35, ha="right")
    ax.set_yticks(range(len(domains)), domains)
    labels = {0: "missing", 1: "partial", 2: "complete"}
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if np.isnan(matrix[i, j]):
                text = "n/a"
                color = "#333333"
            else:
                text = labels[int(matrix[i, j])]
                color = "white" if matrix[i, j] > 0.6 else "black"
            ax.text(j, i, text, ha="center", va="center", fontsize=8.5, color=color)
    cbar = fig.colorbar(im, ax=ax, ticks=[0, 1, 2], fraction=0.04, pad=0.03)
    cbar.ax.set_yticklabels(["missing", "partial", "complete"])
    ax.set_title("Paper 49 Benchmark Readiness Matrix")
    ax.set_xlabel("Competitor / validation class")
    ax.set_ylabel("Domain")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig1_benchmark_readiness_matrix.png", dpi=260)
    plt.close(fig)


def plot_evidence_status(evidence: list[DomainEvidence]) -> None:
    values = {"missing": 0, "partial": 1, "complete": 2}
    domains = [e.domain for e in evidence]
    scores = [values[e.evidence_status] for e in evidence]
    colors = ["#d73027" if s == 0 else "#fee08b" if s == 1 else "#1a9850" for s in scores]
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ax.bar(domains, scores, color=colors, edgecolor="#222222", linewidth=0.8)
    ax.set_ylim(0, 2.25)
    ax.set_yticks([0, 1, 2], ["missing", "partial", "complete"])
    ax.set_ylabel("Current evidence status")
    ax.set_title("Domain-Level Universal Structural Selection Evidence")
    ax.tick_params(axis="x", rotation=35)
    for i, e in enumerate(evidence):
        ax.text(i, scores[i] + 0.06, e.evidence_status, ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig2_domain_evidence_status.png", dpi=260)
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    evidence = [
        fixed_energy_evidence(),
        string_path_evidence(),
        cosmology_evidence(),
        sat_evidence(),
        safety_evidence(),
        sm_bridge_evidence(),
        quantum_evidence(),
    ]
    evidence_rows = [asdict(e) for e in evidence]
    write_csv(OUTDIR / "paper49_domain_evidence_registry.csv", evidence_rows)

    domains, tests, matrix = build_status_matrix()
    matrix_rows: list[dict[str, Any]] = []
    for domain, row in zip(domains, matrix):
        entry = {"domain": domain}
        for test, value in zip(tests, row):
            if np.isnan(value):
                entry[test] = "n/a"
            elif value == 0:
                entry[test] = "missing"
            elif value == 1:
                entry[test] = "partial"
            else:
                entry[test] = "complete"
        matrix_rows.append(entry)
    write_csv(OUTDIR / "paper49_benchmark_readiness_matrix.csv", matrix_rows)

    finite = matrix[np.isfinite(matrix)]
    complete = int(np.sum(finite == 2))
    partial = int(np.sum(finite == 1))
    missing = int(np.sum(finite == 0))
    readiness_index = float(np.sum(finite) / (2.0 * len(finite)))
    summary = {
        "model": "Paper 49 Universal Structural Selection Benchmark Registry",
        "status": (
            "Meta-benchmark aggregator. It audits existing evidence and "
            "identifies missing competitor tests; it is not a new solver."
        ),
        "domains": [e.domain for e in evidence],
        "complete_domain_evidence": sum(e.evidence_status == "complete" for e in evidence),
        "partial_domain_evidence": sum(e.evidence_status == "partial" for e in evidence),
        "missing_domain_evidence": sum(e.evidence_status == "missing" for e in evidence),
        "matrix_complete_cells": complete,
        "matrix_partial_cells": partial,
        "matrix_missing_cells": missing,
        "benchmark_readiness_index": readiness_index,
        "key_findings": [
            "Energy-only baselines are already defeated in field and string benchmarks.",
            "Cosmology has null/permutation protocols but remains at diagnostic evidence level.",
            "SAT and quantum measurement are the main open universality stress tests.",
            "Entropy/activity-only and stability-only competitors are not yet standardized across all domains.",
        ],
        "outputs": [
            "paper49_domain_evidence_registry.csv",
            "paper49_benchmark_readiness_matrix.csv",
            "paper49_summary.json",
            "fig1_benchmark_readiness_matrix.png",
            "fig2_domain_evidence_status.png",
        ],
    }
    with (OUTDIR / "paper49_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    plot_status_matrix(domains, tests, matrix)
    plot_evidence_status(evidence)

    print("Paper 49 universal benchmark registry complete.")
    print(f"Output directory: {OUTDIR}")
    print(f"Complete/partial/missing domain evidence: {summary['complete_domain_evidence']}/"
          f"{summary['partial_domain_evidence']}/{summary['missing_domain_evidence']}")
    print(f"Readiness index: {readiness_index:.4f}")


if __name__ == "__main__":
    main()
