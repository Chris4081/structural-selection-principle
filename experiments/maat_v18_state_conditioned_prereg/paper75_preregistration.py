#!/usr/bin/env python3
"""Write the Paper-75 preregistration. This script cannot execute Paper 76."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "dataset_manifest_paper75.csv"
OUTDIR = ROOT / "outputs_paper75_preregistration"
OUTFILE = OUTDIR / "paper75_preregistration.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def manifest_summary() -> dict:
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {
        "path": MANIFEST.name,
        "sha256": sha256_file(MANIFEST),
        "rows": len(rows),
        "family_counts": dict(sorted(Counter(r["family"] for r in rows).items())),
        "split_counts": dict(sorted(Counter(r["paper75_split"] for r in rows).items())),
        "raw_cnf_in_release": False,
    }


def build_protocol() -> dict:  # noqa: C901 - protocol mirrors the paper.
    return {
        "paper": 75,
        "year": 2026,
        "title": "State-Conditioned Gate Arbitration in SAT/CDCL",
        "subtitle": "A Preregistration for MAAT v1.8 Gate Challenge Cycle II",
        "status": "preregistration_only_no_external_execution",
        "execution_reserved_for": "Paper 76 after Paper 75 has a public Zenodo DOI",
        "provenance_archives": {
            "gate_challenge_cycle_I_papers_69_71": "10.5281/zenodo.21062386",
            "paper72_preregistration": "10.5281/zenodo.21063444",
            "paper73_execution": "10.5281/zenodo.21123988",
        },
        "guardrails": {
            "external_state_conditioned_run_before_doi": False,
            "permitted_pre_release_testing": "synthetic or explicitly excluded development fixtures only",
            "retuning_after_release": False,
            "uncertainty_role": "secondary_report_only_no_primary_gate_retuning",
            "raw_external_data_in_release": False,
        },
        "state_X_t": {
            "evaluation_boundary": "after unit propagation reaches a conflict-free fixpoint and immediately before each branch decision, including the root decision",
            "tuple": ["original_clauses", "learned_clauses", "partial_assignment", "decision_levels", "trail", "reason_clauses", "VSIDS_activity", "saved_phase"],
            "residual_formula": "drop satisfied clauses and remove falsified literals from all remaining original and persistent learned clauses",
            "terminal_states": "conflicting or fully satisfied states are never gate-evaluation states",
            "clause_deletion": "disabled for the primary transparent backend",
            "time_index": "t increments exactly once per branch-decision boundary",
        },
        "supports": {
            "population": "unassigned variables appearing in at least one nonempty residual clause",
            "H": "1/(1+mean_v(abs(pos_v-neg_v)/(pos_v+neg_v+1e-9)))",
            "B": "1/(1+std_v(deg_v)/(mean_v(deg_v)+1e-9))",
            "S": "mean_v(short_v)/(max_v(short_v)+1e-9), short_v=sum_{C contains v, |C|<=3}1/|C|; S=0 if max short_v=0",
            "V": "clip(0.65*largest_component_fraction+0.35*mean_local_clustering,0,1) on the residual variable co-occurrence graph",
            "R_resp": "(H*B*V)^(1/3)",
            "R_struct": "(H*B*S*V)^(1/4)",
            "R_rob": "min(R_resp,R_struct)",
            "raw_support_range": "clip each H,B,S,V to [1e-9,1] before closure",
        },
        "normalization_and_gate": {
            "polarity": 1,
            "epsilon": 1e-9,
            "z_R": "(R_rob-median_cal(R_rob))/(MAD_cal(R_rob)+1e-9)",
            "delta_R": "R_rob[X_t]-R_rob[X_previous_evaluated_state_on_same_run] with 0 at the root",
            "z_t": "(delta_R-median_cal(delta_R))/(MAD_cal(delta_R)+1e-9)",
            "z_x": 0.0,
            "z_x_rationale": "no canonical SAT spatial derivative was frozen; Paper-68 missing-response rule sets it to zero",
            "gate_signal": "g=z_R+0.50*z_t+0.25*z_x",
            "soft_gate": "G=1/(1+exp(-g))",
            "activation_budget_q": 0.25,
            "threshold": "weighted calibration-state 0.75 quantile of G using lower interpolation; frozen after calibration",
            "tie_fraction": "f_tie=clip((0.25-W_cal[G>theta_q])/W_cal[G=theta_q],0,1); set f_tie=0 if denominator is zero",
            "primary_activation": "G>theta_q activates; G<theta_q does not; G=theta_q activates iff HashToUnit(75075,gate_tie,dataset_id,decision_index,state_fingerprint)<f_tie",
            "hard_threshold_0_5_role": "reported Paper-68 ablation only, not the primary budget-matched arbitration",
        },
        "policy": {
            "fallback": "MOMS",
            "active_policy": "Paper-63 Mode A on the current residual frontier",
            "mode_A": "softmax_beta8(progress_v-0.38*h_MAAT_v)",
            "progress": "0.52*VSIDS_norm+0.26*JW_norm+0.14*unit_gain_norm+0.08*degree_norm",
            "h_MAAT": "min-max normalized -sum(log(1e-9+q)) over candidate-local H,B,S,V,R_rob",
            "polarity": "larger sum of occurrence, JW and binary-clause unit-gain scores; existing saved phase is reused with probability 0.25 exactly as in Paper 69",
            "arbitration": "Mode A if G>=theta_q, otherwise MOMS",
            "common_random_numbers": "u_i_t=HashToUnit(paper75|75075|primary|dataset_id|decision_index|state_fingerprint); every directly compared stochastic primary policy uses the same u at the same state; policy identity is excluded",
            "null_random_namespaces": ["gate_tie", "random_null", "shuffled_gate_null", "debug_only"],
        },
        "causal_update_and_rollback": {
            "allowed_information": "current and past solver state only; no trial assignment, rollout, oracle label or test-set aggregate",
            "enqueue": "apply assignment, remove newly satisfied clauses from active counters, remove newly falsified literals, and log inverse deltas",
            "learn": "insert each learned clause after conflict analysis; learned clauses persist across backtracking",
            "backtrack": "undo assignment-induced counter deltas in reverse trail order above target level; do not remove learned clauses",
            "verification": "at every synthetic test event and every 100th Paper-76 decision, incremental supports must equal full recomputation within 1e-12; mismatch aborts execution",
        },
        "evaluation_and_calibration": {
            "frequency": "every branch-decision boundary after propagation and before candidate scoring",
            "trajectory_policy": "MOMS-only",
            "trajectory_rule": "all calibration decision states come exclusively from frozen MOMS-only traces under Paper-75 budgets and primary common-random seed; arbitration, test policies and test instances contribute no calibration state",
            "calibration_states_per_instance": 128,
            "sampling": "retain the 128 smallest SHA256 priorities of paper75|dataset_id|decision_index|state_fingerprint; retain all if fewer",
            "state_weight": "1/(n_families*n_instances_in_family*n_retained_states_for_instance)",
            "instance_weighting": "each family has total weight 1/n_families and each instance has equal weight within family",
            "test_threshold_refit": False,
        },
        "frozen_constants": {
            "activation_budget_q": 0.25,
            "no_harm_delta_relative": 0.05,
            "distinctness_epsilon": 0.01,
            "reconstruction_rho_star": 0.90,
            "conflict_budget": 5000,
            "runtime_budget_seconds_per_policy_instance": 30,
            "bootstrap_resamples": 10000,
            "bootstrap_seed": 75075,
            "minimum_test_activation": 0.125,
        },
        "distinctness": {
            "metric": "absolute weighted Spearman correlation on retained calibration states",
            "degenerate_if": "1-abs(rho)<0.01 for any pair among H,B,S,V or weighted standard deviation<1e-8",
            "outcome": "not_testable_as_constructed",
        },
        "classical_reconstruction": {
            "candidate_score_vectors": {
                "MOMS": "m_v=p_v+n_v+2*p_v*n_v on shortest residual clauses",
                "Jeroslow_Wang": "j_v=jw_pos_v+jw_neg_v with clause weight 2^(-|C|)",
                "VSIDS": "current raw VSIDS activity a_v for every unassigned active variable",
            },
            "aggregations_per_vector": {
                "max": "max(s_v)",
                "mean": "mean(s_v)",
                "std": "population std(s_v)",
                "cv": "std(s_v)/(abs(mean(s_v))+1e-9)",
                "top2_gap": "(s_(1)-s_(2))/(abs(s_(1))+1e-9), descending order; 0 if fewer than 2 candidates",
                "normalized_entropy": "-sum(p_v*log(p_v))/log(n), p_v=max(s_v,0)/sum max(s_u,0); 0 if n<=1 or sum=0",
            },
            "predictor_count": 18,
            "targets": ["H", "B", "S", "V", "R_rob", "G"],
            "model": "ordinary least squares with intercept",
            "validation": "five-fold group cross-validation grouped by calibration instance; deterministic fold hash seed 75075",
            "metric": "family-and-instance-weighted out-of-fold R_squared",
            "threshold": 0.90,
            "primary_gate_coordinates": ["R_rob", "G"],
            "outcome_if_any_primary_gate_coordinate_reaches_threshold": "convergent_rediscovery",
        },
        "dataset": manifest_summary(),
        "dataset_rules": {
            "source": "public SATLIB/DIMACS archives listed in the frozen manifest",
            "split": "within each family sort SHA256(paper75|75075|dataset_id|instance_sha256); first max(1,floor(0.20*n_family)) calibration, remainder test",
            "gatekeeper": ["manifest columns", "file existence", "exact SHA256", "DIMACS parseability", "no unexpected CNFs"],
            "execution_if_gatekeeper_fails": False,
            "raw_data_policy": "raw CNFs and archives excluded; users obtain them from SATLIB and verify hashes",
        },
        "overhead_and_cost": {
            "search_cost": "conflicts+0.30*decisions+0.004*propagations+timeout*5000",
            "timers": "perf_counter_ns around support update, support evaluation, gate calculation and arbitration separately",
            "kappa_cal": "sum calibration MOMS search_cost/sum calibration MOMS wall_seconds",
            "primary_total_cost": "search_cost+kappa_cal*structural_overhead_seconds",
            "relative_regret": "(C_policy-C_MOMS)/max(C_MOMS,1e-9)",
            "no_harm_certified": "paired 95% CI upper endpoint of mean relative regret versus MOMS <=0.05",
            "harmful": "paired 95% CI lower endpoint of mean relative regret versus MOMS >0.05",
            "no_harm_not_established": "CI lower endpoint<=0.05 and CI upper endpoint>0.05",
            "runtime_secondary": "wall runtime includes all structural overhead and is reported separately",
        },
        "baselines": {
            "classical": ["MOMS", "Jeroslow-Wang", "VSIDS"],
            "maat": ["score_with_R_Paper69", "static_gate_v17_Paper69", "state_v18_primary"],
            "nulls": ["state_v18_gate_signal_shuffled_within_family_and_split", "random_activation_matched_to_q"],
            "shuffled_gate_exact_rule": "record G sequences on the already-required MOMS baseline trace; within each family and split, map every instance to the next dataset_id in SHA256(seed75075|dataset_id)-sorted cyclic order; at decision t use donor_G[t mod len(donor_G)] and the same theta_q; no target or utility enters the mapping",
            "random_activation_exact_rule": "at each decision activate iff SHA256(seed75075|dataset_id|decision_index|random_null) interpreted as a uniform 256-bit integer is below 0.25*(2^256-1)",
        },
        "ablations": [
            "MOMS_default_only", "state_Mode_A_always", "state_gate_value_only_z_R",
            "state_gate_response_only_0.5z_t", "state_gate_full_z_R_plus_0.5z_t",
            "Paper68_hard_threshold_0.5", "closure_HBV", "closure_HSV",
            "closure_BSV", "closure_HBS",
        ],
        "metrics": {
            "primary_delta": "C_comparator-C_state_v18; positive favors state_v18",
            "bootstrap": "paired resampling of test instances within family, then equal-weight family aggregation; 10000 resamples",
            "reported": ["mean delta", "95% CI", "wins", "losses", "ties", "family-wise deltas", "activation rate", "overhead fraction", "timeouts"],
            "uncertainty_secondary": ["bootstrap distribution of theta_q", "support confidence bands", "gate-margin sensitivity at frozen theta_q"],
            "uncertainty_cannot_change": ["theta_q", "f_tie", "q", "gate equation", "status axes"],
        },
        "status_axes": {
            "execution_status": {
                "dataset_gatekeeper_fail": "any manifest, hash, parseability or unexpected-file check fails; no execution",
                "executed": "dataset gatekeeper passes and primary execution completes",
            },
            "construct_status": {
                "distinct": "distinctness passes and both primary gate coordinates have weighted OOF reconstruction R_squared<0.90",
                "not_testable_as_constructed": "distinctness precondition fails; primary utility execution stops",
                "convergent_rediscovery": "R_rob or G has weighted OOF reconstruction R_squared>=0.90; utility remains independently reportable",
            },
            "activation_status": {
                "adequate_activation": "family-instance-weighted test decision activation rate>=0.125",
                "activation_collapse": "family-instance-weighted test decision activation rate<0.125",
            },
            "safety_status": {
                "no_harm_certified": "relative-regret CI upper<=0.05",
                "no_harm_not_established": "relative-regret CI lower<=0.05 and CI upper>0.05",
                "harmful": "relative-regret CI lower>0.05",
            },
            "utility_status": {
                "negative_result": "CI upper<=0 against static v1.7 or score-with-R",
                "inconclusive": "minimum-support conditions are not met and negative_result does not apply",
                "minimum_support": "CI lower>0 against static v1.7 and score-with-R and both null advantages have CI lower>0",
                "routing_support": "minimum_support, adequate_activation and no_harm_certified, plus CI lower>0 versus MOMS among test instances with activation rate>=0.125",
                "strong_sat_support": "minimum_support, adequate_activation and no_harm_certified, plus overall CI lower>0 versus both MOMS and Jeroslow-Wang",
            },
            "reporting_rule": "report every axis; construct status never suppresses safety or utility status; activation_collapse forbids routing_support and strong_sat_support",
        },
        "release_condition": "Paper 76 is forbidden until this exact artifact has a public Zenodo DOI.",
    }


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    protocol = build_protocol()
    OUTFILE.write_text(json.dumps(protocol, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PREREGISTRATION WRITTEN", "solver_executed": False, "output": str(OUTFILE.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
