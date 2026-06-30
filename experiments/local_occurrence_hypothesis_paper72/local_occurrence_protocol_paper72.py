#!/usr/bin/env python3
"""Paper 72 preregistration artifact for the Local Occurrence Hypothesis.

This script intentionally does not execute a SAT solver, load CNFs, compute
features, or tune parameters. It writes the frozen protocol that a future Paper
73 execution must follow.
"""

from __future__ import annotations

import json
from pathlib import Path


OUT_DIR = Path("outputs_paper72_local_occurrence_protocol")
OUT_FILE = OUT_DIR / "paper72_protocol.json"


def build_protocol() -> dict:
    return {
        "paper": 72,
        "title": "The Local Occurrence Hypothesis",
        "subtitle": "A Preregistered SAT-Specific Test Protocol after Gate Challenge Cycle I",
        "year": 2026,
        "status": "preregistration_only_no_solver_execution",
        "series_archive": {
            "name": "MAAT Research Series II -- Gate Challenge Cycle I",
            "doi": "10.5281/zenodo.21062386",
            "url": "https://doi.org/10.5281/zenodo.21062386",
        },
        "guardrails": [
            "no_retroactive_change_to_paper68",
            "no_new_maat_version",
            "no_post_hoc_retuning",
            "no_new_gate_claim",
            "paper72_defines_protocol_only",
            "gate_plus_L_solver_run_reserved_for_paper73",
        ],
        "background": {
            "paper69": "Frozen MAAT v1.7 gate was not supported in the SATLIB smoke execution and lost to MOMS.",
            "paper70": "MOMS advantage was diagnosed as decision-local structural information.",
            "paper71": "Most local channels were reconstructible from H,B,S,V,R_rob,G_gate; only weak occurrence/degree residuals remained.",
        },
        "hypothesis": (
            "Decision-local occurrence/degree residuals provide additional "
            "SAT-specific decision information beyond the frozen v1.7 gate."
        ),
        "null_hypothesis": (
            "The residual local channels provide no robust improvement over "
            "frozen gate, score-with-R, or classical SAT heuristics."
        ),
        "frozen_support_basis": ["H", "B", "S", "V", "R_rob", "G_gate"],
        "local_channel_definitions": {
            "variable_occurrence": {
                "symbol": "O_raw",
                "definition": "std_v(o_v)/(mean_v(o_v)+epsilon), where o_v counts signed or unsigned variable occurrences across clauses.",
            },
            "degree_cv": {
                "symbol": "D_raw",
                "definition": "coefficient of variation of the Paper-69 incidence degree vector.",
            },
            "residualization": {
                "calibration_only": True,
                "model": "linear reconstruction from z-scored H,B,S,V,R_rob,G_gate plus intercept",
                "test_application": "freeze calibration coefficients and calibration z-score statistics; apply without refitting",
                "residual_z_standardization": "after residualization, standardize O_res and D_res using calibration-fold means and standard deviations only; apply frozen statistics to test fold",
            },
            "primary_combined_channel": {
                "symbol": "L_star",
                "formula": "0.5 * (z(D_res) - z(O_res))",
                "z_definition": "post-residualization calibration-fold standardization",
                "orientation_source": "Paper 71 hypothesis-generation result; fixed before Paper 73 execution",
            },
        },
        "gate_plus_L_combination": {
            "base_gate": "frozen Paper-69 gate_v17",
            "primary_formula": "G_72 = sigmoid(logit(G_gate) + L_star)",
            "lambda_L": 1.0,
            "lambda_L_rationale": "fixed at 1.0 to avoid turning Paper 73 into a hyperparameter search over local-channel strength",
            "clipping": "clip G_gate to [1e-6, 1-1e-6] before logit",
            "policy_scope": "replace only the root-level gate activation used by the Paper-69 gate policy; all other solver rules remain frozen",
        },
        "dataset_rules": {
            "sources": ["public SATLIB/DIMACS CNF archives"],
            "raw_data_policy": "raw CNF files excluded from git; manifest and SHA256 required before execution",
            "primary_family_balanced_instances_per_family": 10,
            "calibration_fraction_per_family": 0.20,
            "test_fraction_per_family": 0.80,
            "split_seed": 72072,
            "fallback_if_family_has_fewer_than_10": "use all available validated CNFs and preserve 20/80 split where possible",
        },
        "budgets": {
            "primary_conflict_budget": 5000,
            "runtime_budget_seconds_per_policy_instance": 30,
        },
        "baselines": {
            "classical": ["MOMS", "Jeroslow-Wang", "VSIDS"],
            "maat": ["score_with_R", "frozen_gate_v17"],
            "nulls": ["shuffled_L_within_family"],
        },
        "metrics": {
            "primary": "paired utility delta; positive means lower compute cost for Gate+L than comparator",
            "cost_components": ["conflicts", "decisions", "runtime", "compute_cost"],
            "confidence_interval": {
                "method": "paired bootstrap on instance-level deltas",
                "resamples": 10000,
                "seed": 72072,
                "level": 0.95,
            },
        },
        "success_criteria": {
            "minimum_support": (
                "Gate+L beats frozen gate and score-with-R with 95% CI lower "
                "bound greater than zero."
            ),
            "strong_sat_support": (
                "Gate+L also beats MOMS and Jeroslow-Wang with 95% CI lower "
                "bound greater than zero in the family-balanced external benchmark."
            ),
        },
        "failure_criteria": {
            "no_minimum_support": [
                "Gate+L fails to beat frozen gate with 95% CI lower bound greater than zero.",
                "Gate+L fails to beat score-with-R with 95% CI lower bound greater than zero.",
                "The apparent MAAT-baseline gain disappears under shuffled-L within-family null.",
            ],
            "no_strong_sat_support": [
                "Gate+L fails to beat MOMS with 95% CI lower bound greater than zero.",
                "Gate+L fails to beat Jeroslow-Wang with 95% CI lower bound greater than zero.",
            ],
            "minimum_only_result": (
                "If Gate+L beats frozen gate and score-with-R but does not beat "
                "MOMS or Jeroslow-Wang, report minimum support only, not strong SAT support. "
                "This category must not be promoted rhetorically to strong SAT support."
            ),
        },
        "non_goals": [
            "This protocol does not rescue Paper 69.",
            "This protocol does not modify Paper 68.",
            "This protocol does not define MAAT v1.8.",
            "This protocol does not claim that occurrence residuals are a universal MAAT support.",
        ],
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    protocol = build_protocol()
    OUT_FILE.write_text(json.dumps(protocol, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "PROTOCOL WRITTEN",
        "output": str(OUT_FILE),
        "solver_executed": False,
    }, indent=2))


if __name__ == "__main__":
    main()
