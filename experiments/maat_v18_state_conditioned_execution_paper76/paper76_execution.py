#!/usr/bin/env python3
"""Paper 76: frozen execution of the Paper-75 state-conditioned protocol.

No parameter in this file is fitted on test utility.  Calibration uses only
the manifest-labelled calibration split and MOMS-only trajectories.  The
external run is refused unless the published Paper-75 DOI, protocol hash,
manifest hash, and every raw CNF hash pass their gatekeepers.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from paper76_core import (
    EPS,
    IncrementalResidualTracker,
    SupportSnapshot,
    classify_safety,
    closure,
    hash_to_unit,
    max_snapshot_difference,
    random_activation,
    sigmoid,
    softmax_choice,
    state_priority,
    weighted_corr,
    weighted_mad,
    weighted_mean,
    weighted_quantile_lower,
    weighted_spearman,
    weighted_std,
)


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "outputs_paper76_state_conditioned_execution"
MANIFEST = ROOT / "dataset_manifest_paper75.csv"
P75_PROTOCOL = ROOT / "paper75_preregistration_frozen.json"
P75_RELEASE = ROOT / "paper75_release_verification.json"

PAPER75_DOI = "10.5281/zenodo.21767033"
PAPER75_PROTOCOL_SHA256 = "9e86b1e1f8ddcd3228c90f62cc9465d3dccbd9049d5815c5d7dfda0ffd2f3b9e"
PAPER75_MANIFEST_SHA256 = "941fad18c5ea4839ee3aac4b5d7f1e3d5dd0cf8ed8c8f37f3ec9a4a9a74df3ad"

POLICIES = [
    "moms",
    "jeroslow_wang",
    "vsids",
    "score_with_R_Paper69",
    "static_gate_v17_Paper69",
    "state_v18_primary",
    "state_v18_gate_signal_shuffled",
    "random_activation_matched_q",
]

ABLATIONS = [
    "state_Mode_A_always",
    "state_gate_value_only_z_R",
    "state_gate_response_only_0_5z_t",
    "Paper68_hard_threshold_0_5",
    "closure_HBV",
    "closure_HSV",
    "closure_BSV",
    "closure_HBS",
]

PRIMARY = "state_v18_primary"
COMPARATORS = [
    "static_gate_v17_Paper69",
    "score_with_R_Paper69",
    "state_v18_gate_signal_shuffled",
    "random_activation_matched_q",
    "moms",
    "jeroslow_wang",
    "vsids",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def find_paper69_dir(explicit: Path | None = None) -> Path:
    candidates = []
    if explicit:
        candidates.append(explicit)
    if os.environ.get("MAAT_SSL_REPO"):
        candidates.append(Path(os.environ["MAAT_SSL_REPO"]) / "experiments" / "gate_challenge_sat_paper69")
    candidates.extend(
        [
            ROOT.parent / "gate_challenge_sat_paper69",
            ROOT.parents[1] / "experiments" / "gate_challenge_sat_paper69" if len(ROOT.parents) > 1 else ROOT,
        ]
    )
    for candidate in candidates:
        if (candidate / "gate_challenge_sat_paper69.py").exists():
            return candidate.resolve()
    raise RuntimeError("Cannot locate the frozen Paper-69 transparent CDCL backend")


P69: Any = None


def load_p69(path: Path) -> Any:
    global P69
    if P69 is None:
        paper63_candidates = []
        if os.environ.get("MAAT_SSL_REPO"):
            paper63_candidates.append(Path(os.environ["MAAT_SSL_REPO"]) / "experiments" / "sat_cdcl_structure_gated_mode_a_paper63")
        paper63_candidates.extend(
            [
                path.parent / "sat_cdcl_structure_gated_mode_a_paper63",
            ]
        )
        for candidate in paper63_candidates:
            if (candidate / "sat_cdcl_structure_gated_mode_a_paper63.py").exists():
                sys.path.insert(0, str(candidate))
                break
        sys.path.insert(0, str(path))
        import gate_challenge_sat_paper69 as module  # type: ignore

        P69 = module
    return P69


@dataclass
class Instance:
    dataset_id: str
    family: str
    local_path: str
    sha256: str
    split: str
    instance_weight: float
    n_vars: int
    n_clauses: int
    clauses: list[tuple[int, ...]]


@dataclass
class TraceState:
    dataset_id: str
    family: str
    decision_index: int
    fingerprint: str
    priority: str
    H: float
    B: float
    S: float
    V: float
    predictors: tuple[float, ...]
    R_variants: dict[str, float]
    delta_variants: dict[str, float]
    G_variants: dict[str, float] = field(default_factory=dict)
    weight: float = 0.0


@dataclass
class SolveRecord:
    dataset_id: str
    family: str
    split: str
    policy: str
    n_vars: int
    n_clauses: int
    solved: int
    sat: int
    timeout: int
    decisions: int
    propagations: int
    conflicts: int
    learned_clauses: int
    max_level: int
    runtime_s: float
    search_cost: float
    structural_update_s: float
    support_evaluation_s: float
    gate_calculation_s: float
    arbitration_s: float
    structural_overhead_s: float
    total_cost: float
    activations: int
    activation_rate: float
    audit_max_abs_error: float


def resolve_raw_path(data_dir: Path, manifest_local_path: str) -> Path:
    rel = Path(manifest_local_path)
    if rel.parts and rel.parts[0] == "data_external":
        rel = Path(*rel.parts[1:])
    return data_dir / rel


def read_manifest() -> list[dict[str, str]]:
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def validate_release() -> list[str]:
    failures: list[str] = []
    if sha256_file(P75_PROTOCOL) != PAPER75_PROTOCOL_SHA256:
        failures.append("Paper-75 protocol SHA256 mismatch")
    if sha256_file(MANIFEST) != PAPER75_MANIFEST_SHA256:
        failures.append("Paper-75 manifest SHA256 mismatch")
    if not P75_RELEASE.exists():
        failures.append("missing Paper-75 Zenodo release verification")
    else:
        release = json.loads(P75_RELEASE.read_text(encoding="utf-8"))
        if release.get("doi") != PAPER75_DOI:
            failures.append("Paper-75 DOI mismatch")
        if release.get("status") != "published":
            failures.append("Paper-75 Zenodo record is not published")
        if release.get("publication_date") != "2026-08-03":
            failures.append("unexpected Paper-75 publication date")
    protocol = json.loads(P75_PROTOCOL.read_text(encoding="utf-8"))
    if protocol.get("paper") != 75 or protocol.get("status") != "preregistration_only_no_external_execution":
        failures.append("invalid frozen Paper-75 protocol identity")
    return failures


def validate_dataset(data_dir: Path, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    release_failures = validate_release()
    rows = read_manifest()
    required = {
        "dataset_id", "family", "local_path", "sha256", "paper75_split",
        "paper75_instance_weight", "source_url", "license_or_terms",
    }
    missing_columns = sorted(required - set(rows[0] if rows else {}))
    missing_files: list[str] = []
    hash_mismatches: list[str] = []
    parse_failures: list[str] = []
    expected_paths: set[Path] = set()
    p69 = load_p69(find_paper69_dir())
    for row in rows:
        path = resolve_raw_path(data_dir, row["local_path"])
        expected_paths.add(path.resolve())
        if not path.exists():
            missing_files.append(row["local_path"])
            continue
        if sha256_file(path).lower() != row["sha256"].lower():
            hash_mismatches.append(row["local_path"])
            continue
        try:
            n_vars, clauses = p69.parse_dimacs(path)
            if n_vars <= 0 or not clauses:
                raise ValueError("empty DIMACS")
        except Exception as exc:  # pragma: no cover - external corruption path
            parse_failures.append(f"{row['local_path']}: {exc}")
    actual_paths = {path.resolve() for path in data_dir.rglob("*.cnf")}
    unexpected = sorted(str(path) for path in actual_paths - expected_paths)
    failures = {
        "release": release_failures,
        "missing_columns": missing_columns,
        "missing_files": missing_files,
        "sha256_mismatches": hash_mismatches,
        "parse_failures": parse_failures,
        "unexpected_cnf_files": unexpected,
    }
    passed = not any(failures.values()) and len(rows) == 98
    result = {
        "paper": 76,
        "paper75_preregistration_doi": PAPER75_DOI,
        "validation_type": "release_protocol_manifest_sha256_dimacs_gatekeeper",
        "paths": {"manifest": "./dataset_manifest_paper75.csv", "data": "./data_external"},
        "checks": {
            "Paper75 release": "PASS" if not release_failures else "FAIL",
            "Protocol SHA256": "PASS" if not any("protocol" in x.lower() for x in release_failures) else "FAIL",
            "Manifest": "PASS" if not missing_columns and len(rows) == 98 else "FAIL",
            "CNFs": "PASS" if not missing_files else "FAIL",
            "SHA256": "PASS" if not hash_mismatches else "FAIL",
            "Parseability": "PASS" if not parse_failures else "FAIL",
            "Unexpected Files": "PASS" if not unexpected else "FAIL",
            "Dataset Ready": "YES" if passed else "NO",
        },
        "counts": {
            "manifest_rows": len(rows), "actual_cnf_files": len(actual_paths),
            "verified_files": len(rows) - len(missing_files) - len(hash_mismatches),
        },
        "failures": failures,
        "dataset_ready": passed,
        "status": "VALIDATION PASS" if passed else "VALIDATION FAIL",
    }
    (output_dir / "paper76_dataset_validation.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def load_instances(data_dir: Path) -> list[Instance]:
    p69 = load_p69(find_paper69_dir())
    instances: list[Instance] = []
    for row in read_manifest():
        path = resolve_raw_path(data_dir, row["local_path"])
        n_vars, clauses = p69.parse_dimacs(path)
        instances.append(
            Instance(
                dataset_id=row["dataset_id"], family=row["family"], local_path=row["local_path"],
                sha256=row["sha256"], split=row["paper75_split"],
                instance_weight=float(row["paper75_instance_weight"]), n_vars=n_vars,
                n_clauses=len(clauses), clauses=clauses,
            )
        )
    return instances


def gate_active(value: float, threshold: float, tie_fraction: float, dataset_id: str, decision: int, fingerprint: str) -> bool:
    if value > threshold:
        return True
    if value < threshold:
        return False
    return hash_to_unit(dataset_id, decision, fingerprint, namespace="gate_tie") < tie_fraction


class Paper76Solver:
    """Adapter around the frozen transparent CDCL backend."""

    def __init__(
        self,
        instance: Instance,
        policy: str,
        calibration: dict[str, Any] | None,
        static_gate: float = 0.0,
        donor_g: list[float] | None = None,
        collect_trace: bool = False,
        conflict_budget: int = 5000,
        runtime_budget: float = 30.0,
    ) -> None:
        p69 = load_p69(find_paper69_dir())
        self.backend = p69.Paper69Solver(
            instance.clauses, instance.n_vars, policy="moms", seed=75075,
            structure_gate=0.0, conflict_budget=conflict_budget, beta=8.0, tau=0.38,
        )
        self.instance = instance
        self.policy = policy
        self.calibration = calibration
        self.static_gate = static_gate
        self.donor_g = donor_g or []
        self.collect_trace = collect_trace
        self.runtime_budget = runtime_budget
        self.tracker = IncrementalResidualTracker(instance.clauses, instance.n_vars)
        self.trace: list[TraceState] = []
        self.prev_r = {name: None for name in ["HBSV", "HBV", "HSV", "BSV", "HBS"]}
        self.decision_index = 0
        self.support_ns = 0
        self.gate_ns = 0
        self.arbitration_ns = 0
        self.activations = 0
        self.audit_max_error = 0.0
        self.start_time = 0.0
        self.timed_out = False
        self._patch_backend_hooks()

    def _patch_backend_hooks(self) -> None:
        original_enqueue = self.backend.enqueue
        original_backtrack = self.backend.backtrack

        def enqueue(lit: int, reason: tuple[int, ...] | None) -> bool:
            existed = abs(lit) in self.backend.assignment
            ok = original_enqueue(lit, reason)
            if ok and not existed:
                self.tracker.assign(abs(lit), lit > 0)
            return ok

        def backtrack(level: int) -> None:
            removed = [abs(lit) for lit in self.backend.trail if self.backend.level_of.get(abs(lit), 0) > level]
            for var in reversed(removed):
                self.tracker.unassign(var)
            original_backtrack(level)

        self.backend.enqueue = enqueue
        self.backend.backtrack = backtrack

    def _runtime_exceeded(self) -> bool:
        return time.perf_counter() - self.start_time >= self.runtime_budget

    def _variant_values(self, snapshot: SupportSnapshot) -> tuple[dict[str, float], dict[str, float]]:
        supports = snapshot.as_dict()
        r_values: dict[str, float] = {}
        deltas: dict[str, float] = {}
        for variant in self.prev_r:
            r = snapshot.R_rob if variant == "HBSV" else closure(supports, variant)[2]
            previous = self.prev_r[variant]
            delta = 0.0 if previous is None else r - previous
            r_values[variant] = r
            deltas[variant] = delta
            self.prev_r[variant] = r
        return r_values, deltas

    def _gate_values(self, r_values: dict[str, float], deltas: dict[str, float]) -> dict[str, float]:
        if self.calibration is None:
            return {}
        values: dict[str, float] = {}
        for variant in r_values:
            norm = self.calibration["variants"][variant]
            z_r = (r_values[variant] - norm["median_R"]) / (norm["mad_R"] + EPS)
            z_t = (deltas[variant] - norm["median_delta_R"]) / (norm["mad_delta_R"] + EPS)
            values[variant] = sigmoid(z_r + 0.5 * z_t)
            if variant == "HBSV":
                values["value_only"] = sigmoid(z_r)
                values["response_only"] = sigmoid(0.5 * z_t)
        return values

    def _evaluate(self) -> tuple[SupportSnapshot, str, dict[str, float], dict[str, float], dict[str, float]]:
        start = time.perf_counter_ns()
        snapshot = self.tracker.snapshot(self.backend.vsids)
        self.support_ns += time.perf_counter_ns() - start
        fingerprint = self.tracker.fingerprint(self.backend.decision_level, self.backend.trail)
        r_values, deltas = self._variant_values(snapshot)
        start = time.perf_counter_ns()
        gates = self._gate_values(r_values, deltas)
        self.gate_ns += time.perf_counter_ns() - start
        if self.decision_index % 100 == 0:
            full = self.tracker.full_snapshot(self.backend.vsids)
            error = max_snapshot_difference(snapshot, full)
            self.audit_max_error = max(self.audit_max_error, error)
            if error > 1.0e-12:
                raise RuntimeError(f"incremental/full support audit failed: {error:.3e}")
        if self.collect_trace:
            self.trace.append(
                TraceState(
                    dataset_id=self.instance.dataset_id, family=self.instance.family,
                    decision_index=self.decision_index, fingerprint=fingerprint,
                    priority=state_priority(self.instance.dataset_id, self.decision_index, fingerprint),
                    H=snapshot.H, B=snapshot.B, S=snapshot.S, V=snapshot.V,
                    predictors=snapshot.predictors, R_variants=r_values.copy(),
                    delta_variants=deltas.copy(), G_variants=gates.copy(),
                )
            )
        return snapshot, fingerprint, r_values, deltas, gates

    def _moms_choice(self, scores: dict[int, dict[str, float]], literal: dict[int, dict[str, float]]) -> tuple[int, bool]:
        unassigned = sorted(scores)
        def key(var: int) -> tuple[float, float, float, int]:
            p, n = literal[var]["moms_pos"], literal[var]["moms_neg"]
            return p + n + 2.0 * p * n, p + n, self.backend.vsids[var], -var
        var = max(unassigned, key=key)
        return var, literal[var]["moms_pos"] >= literal[var]["moms_neg"]

    def _mode_a_choice(self, scores: dict[int, dict[str, float]], fingerprint: str) -> tuple[int, bool]:
        uniform = hash_to_unit(self.instance.dataset_id, self.decision_index, fingerprint)
        var = softmax_choice({v: scores[v]["mode_a"] for v in scores}, 8.0, uniform)
        polarity = scores[var]["polarity_true_score"] >= scores[var]["polarity_false_score"]
        if var in self.backend.phase and uniform < 0.25:
            polarity = self.backend.phase[var]
        return var, polarity

    def choose_decision(self) -> int | None:
        unassigned = [v for v in range(1, self.backend.n_vars + 1) if v not in self.backend.assignment]
        if not unassigned:
            return None
        structural = (
            self.collect_trace
            or self.policy.startswith("state_")
            or self.policy.startswith("closure_")
            or self.policy in {"random_activation_matched_q", "Paper68_hard_threshold_0_5"}
        )
        snapshot = None
        fingerprint = self.tracker.fingerprint(self.backend.decision_level, self.backend.trail)
        gates: dict[str, float] = {}
        if structural or self.collect_trace:
            snapshot, fingerprint, _, _, gates = self._evaluate()
        scores = self.backend.score_candidates()
        literal = self.backend.literal_scores_jw_moms()
        start = time.perf_counter_ns()
        active = False
        if self.policy == "moms":
            var, polarity = self._moms_choice(scores, literal)
        elif self.policy == "jeroslow_wang":
            var = max(unassigned, key=lambda v: (literal[v]["jw_pos"] + literal[v]["jw_neg"], self.backend.vsids[v], -v))
            polarity = literal[var]["jw_pos"] >= literal[var]["jw_neg"]
        elif self.policy == "vsids":
            var = max(unassigned, key=lambda v: (self.backend.vsids[v], scores[v]["degree"], -v))
            polarity = self.backend.phase.get(var, True)
        elif self.policy == "score_with_R_Paper69":
            var, polarity = self._mode_a_choice(scores, fingerprint)
            active = True
        elif self.policy == "static_gate_v17_Paper69":
            uniform = hash_to_unit(self.instance.dataset_id, self.decision_index, fingerprint)
            var = softmax_choice(
                {v: scores[v]["progress"] - self.static_gate * 0.38 * scores[v]["h_maat"] for v in scores},
                8.0, uniform,
            )
            polarity = scores[var]["polarity_true_score"] >= scores[var]["polarity_false_score"]
            if var in self.backend.phase and uniform < 0.25:
                polarity = self.backend.phase[var]
            active = self.static_gate > 0
        else:
            if self.calibration is None:
                raise RuntimeError("state policy requires frozen calibration")
            variant = "HBSV"
            gate_key = "HBSV"
            if self.policy.startswith("closure_"):
                variant = self.policy.removeprefix("closure_")
                gate_key = variant
            elif self.policy == "state_gate_value_only_z_R":
                gate_key = "value_only"
            elif self.policy == "state_gate_response_only_0_5z_t":
                gate_key = "response_only"
            gate_value = gates.get(gate_key, 0.0)
            threshold_info = self.calibration["gates"].get(gate_key, self.calibration["gates"]["HBSV"])
            if self.policy == "state_Mode_A_always":
                active = True
            elif self.policy == "Paper68_hard_threshold_0_5":
                active = gate_value >= 0.5
            elif self.policy == "random_activation_matched_q":
                active = random_activation(self.instance.dataset_id, self.decision_index)
            elif self.policy == "state_v18_gate_signal_shuffled":
                gate_value = self.donor_g[self.decision_index % len(self.donor_g)] if self.donor_g else 0.0
                active = gate_active(gate_value, threshold_info["theta"], threshold_info["tie_fraction"], self.instance.dataset_id, self.decision_index, fingerprint)
            else:
                active = gate_active(gate_value, threshold_info["theta"], threshold_info["tie_fraction"], self.instance.dataset_id, self.decision_index, fingerprint)
            if active:
                var, polarity = self._mode_a_choice(scores, fingerprint)
            else:
                var, polarity = self._moms_choice(scores, literal)
        self.arbitration_ns += time.perf_counter_ns() - start
        self.activations += int(active)
        self.decision_index += 1
        return var if polarity else -var

    def _learn(self, learned: tuple[int, ...]) -> None:
        self.backend.clauses.append(learned)
        self.tracker.add_clause(learned)
        self.backend.learned_clauses += 1
        self.backend.bump_vsids(learned)

    def solve(self) -> tuple[bool, bool]:
        self.start_time = time.perf_counter()
        while True:
            if self._runtime_exceeded():
                self.timed_out = True
                return False, False
            conflict = self.backend.propagate()
            if conflict is None:
                break
            self.backend.conflicts += 1
            if self.backend.decision_level == 0:
                return False, True
            learned, backjump = self.backend.analyze_conflict(conflict)
            if not learned:
                return False, True
            self._learn(learned)
            self.backend.backtrack(backjump)
        while True:
            if self._runtime_exceeded():
                self.timed_out = True
                return False, False
            if self.backend.is_satisfied():
                return True, True
            if self.backend.conflicts >= self.backend.conflict_budget:
                return False, False
            decision_lit = self.choose_decision()
            if decision_lit is None:
                return self.backend.is_satisfied(), True
            self.backend.decision_level += 1
            self.backend.max_level = max(self.backend.max_level, self.backend.decision_level)
            self.backend.decisions += 1
            if not self.backend.enqueue(decision_lit, None):
                return False, True
            while True:
                if self._runtime_exceeded():
                    self.timed_out = True
                    return False, False
                conflict = self.backend.propagate()
                if conflict is None:
                    break
                self.backend.conflicts += 1
                if self.backend.decision_level == 0:
                    return False, True
                learned, backjump = self.backend.analyze_conflict(conflict)
                if not learned:
                    return False, True
                self._learn(learned)
                self.backend.backtrack(backjump)
                if self.backend.conflicts >= self.backend.conflict_budget:
                    return False, False

    def record(self, sat: bool, solved: bool, runtime: float, kappa: float = 0.0) -> SolveRecord:
        timeout = int(not solved)
        search_cost = (
            self.backend.conflicts + 0.30 * self.backend.decisions
            + 0.004 * self.backend.propagations + timeout * 5000
        )
        update_s = self.tracker.update_ns / 1.0e9
        support_s = self.support_ns / 1.0e9
        gate_s = self.gate_ns / 1.0e9
        arbitration_s = self.arbitration_ns / 1.0e9
        overhead = update_s + support_s + gate_s + arbitration_s
        charged_policies = {
            "state_v18_primary",
            "state_v18_gate_signal_shuffled",
            "random_activation_matched_q",
            "state_Mode_A_always",
            "state_gate_value_only_z_R",
            "state_gate_response_only_0_5z_t",
            "Paper68_hard_threshold_0_5",
            "closure_HBV",
            "closure_HSV",
            "closure_BSV",
            "closure_HBS",
        }
        total_cost = search_cost + (kappa * overhead if self.policy in charged_policies else 0.0)
        return SolveRecord(
            dataset_id=self.instance.dataset_id, family=self.instance.family, split=self.instance.split,
            policy=self.policy, n_vars=self.instance.n_vars, n_clauses=self.instance.n_clauses,
            solved=int(solved), sat=int(bool(sat)), timeout=timeout,
            decisions=self.backend.decisions, propagations=self.backend.propagations,
            conflicts=self.backend.conflicts, learned_clauses=self.backend.learned_clauses,
            max_level=self.backend.max_level, runtime_s=runtime, search_cost=search_cost,
            structural_update_s=update_s, support_evaluation_s=support_s,
            gate_calculation_s=gate_s, arbitration_s=arbitration_s,
            structural_overhead_s=overhead, total_cost=total_cost,
            activations=self.activations,
            activation_rate=self.activations / max(1, self.backend.decisions),
            audit_max_abs_error=self.audit_max_error,
        )


def run_solver(
    instance: Instance,
    policy: str,
    calibration: dict[str, Any] | None,
    static_gate: float = 0.0,
    donor_g: list[float] | None = None,
    collect_trace: bool = False,
    kappa: float = 0.0,
) -> tuple[SolveRecord, list[TraceState]]:
    solver = Paper76Solver(instance, policy, calibration, static_gate, donor_g, collect_trace)
    start = time.perf_counter()
    sat, solved = solver.solve()
    runtime = time.perf_counter() - start
    return solver.record(sat, solved, runtime, kappa), solver.trace


def retain_calibration_states(traces: dict[str, list[TraceState]], instances: list[Instance]) -> list[TraceState]:
    by_id = {instance.dataset_id: instance for instance in instances}
    family_counts: dict[str, int] = {}
    for instance in instances:
        if instance.split == "calibration":
            family_counts[instance.family] = family_counts.get(instance.family, 0) + 1
    retained: list[TraceState] = []
    n_families = len(family_counts)
    for dataset_id, states in traces.items():
        chosen = sorted(states, key=lambda state: state.priority)[:128]
        instance = by_id[dataset_id]
        weight = 1.0 / (n_families * family_counts[instance.family] * max(1, len(chosen)))
        for state in chosen:
            state.weight = weight
            retained.append(state)
    return retained


def calibrate(states: list[TraceState]) -> dict[str, Any]:
    weights = [state.weight for state in states]
    result: dict[str, Any] = {"variants": {}, "gates": {}}
    for variant in ["HBSV", "HBV", "HSV", "BSV", "HBS"]:
        r_values = [state.R_variants[variant] for state in states]
        deltas = [state.delta_variants[variant] for state in states]
        mad_r = weighted_mad(r_values, weights)
        # Paper 75 freezes weighted medians, not means.
        from paper76_core import weighted_median
        med_r = weighted_median(r_values, weights)
        med_d = weighted_median(deltas, weights)
        mad_d = weighted_mad(deltas, weights)
        result["variants"][variant] = {
            "median_R": med_r, "mad_R": mad_r,
            "median_delta_R": med_d, "mad_delta_R": mad_d,
        }
        gates = [sigmoid((r - med_r) / (mad_r + EPS) + 0.5 * (d - med_d) / (mad_d + EPS)) for r, d in zip(r_values, deltas)]
        for state, gate in zip(states, gates):
            state.G_variants[variant] = gate
        theta = weighted_quantile_lower(gates, weights, 0.75)
        w_above = sum(w for gate, w in zip(gates, weights) if gate > theta)
        w_equal = sum(w for gate, w in zip(gates, weights) if gate == theta)
        tie = max(0.0, min(1.0, (0.25 - w_above) / w_equal)) if w_equal > 0 else 0.0
        result["gates"][variant] = {"theta": theta, "tie_fraction": tie}
    base = result["variants"]["HBSV"]
    for key, formula in [
        ("value_only", lambda state: sigmoid((state.R_variants["HBSV"] - base["median_R"]) / (base["mad_R"] + EPS))),
        ("response_only", lambda state: sigmoid(0.5 * (state.delta_variants["HBSV"] - base["median_delta_R"]) / (base["mad_delta_R"] + EPS))),
    ]:
        gates = [formula(state) for state in states]
        for state, gate in zip(states, gates):
            state.G_variants[key] = gate
        theta = weighted_quantile_lower(gates, weights, 0.75)
        w_above = sum(w for gate, w in zip(gates, weights) if gate > theta)
        w_equal = sum(w for gate, w in zip(gates, weights) if gate == theta)
        tie = max(0.0, min(1.0, (0.25 - w_above) / w_equal)) if w_equal > 0 else 0.0
        result["gates"][key] = {"theta": theta, "tie_fraction": tie}
    return result


def group_fold(dataset_id: str) -> int:
    digest = hashlib.sha256(f"75075|{dataset_id}".encode()).hexdigest()
    return int(digest, 16) % 5


def preconditions(states: list[TraceState], calibration: dict[str, Any]) -> dict[str, Any]:
    weights = [state.weight for state in states]
    supports = {key: [getattr(state, key) for state in states] for key in ["H", "B", "S", "V"]}
    pairs: list[dict[str, Any]] = []
    distinct = True
    for i, left in enumerate(supports):
        if weighted_std(supports[left], weights) < 1.0e-8:
            distinct = False
        for right in list(supports)[i + 1:]:
            rho = weighted_spearman(supports[left], supports[right], weights)
            pairs.append({"left": left, "right": right, "rho": rho})
            if math.isfinite(rho) and 1.0 - abs(rho) < 0.01:
                distinct = False

    X = np.array([state.predictors for state in states], dtype=float)
    targets = {
        "H": np.array(supports["H"]), "B": np.array(supports["B"]),
        "S": np.array(supports["S"]), "V": np.array(supports["V"]),
        "R_rob": np.array([state.R_variants["HBSV"] for state in states]),
        "G": np.array([state.G_variants["HBSV"] for state in states]),
    }
    folds = np.array([group_fold(state.dataset_id) for state in states])
    weight_arr = np.array(weights)
    reconstruction: dict[str, float] = {}
    for name, y in targets.items():
        pred = np.full(len(y), np.nan)
        for fold in range(5):
            train, test = folds != fold, folds == fold
            if not np.any(test) or np.sum(train) <= X.shape[1]:
                continue
            design = np.column_stack([np.ones(np.sum(train)), X[train]])
            coef, *_ = np.linalg.lstsq(design, y[train], rcond=None)
            pred[test] = np.column_stack([np.ones(np.sum(test)), X[test]]) @ coef
        valid = np.isfinite(pred)
        if not np.any(valid):
            reconstruction[name] = float("nan")
            continue
        mean_y = np.average(y[valid], weights=weight_arr[valid])
        ss_res = np.sum(weight_arr[valid] * (y[valid] - pred[valid]) ** 2)
        ss_tot = np.sum(weight_arr[valid] * (y[valid] - mean_y) ** 2)
        reconstruction[name] = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
    if not distinct:
        construct_status = "not_testable_as_constructed"
    elif max(reconstruction.get("R_rob", -math.inf), reconstruction.get("G", -math.inf)) >= 0.90:
        construct_status = "convergent_rediscovery"
    else:
        construct_status = "distinct"
    return {
        "construct_status": construct_status,
        "support_pair_correlations": pairs,
        "support_weighted_std": {key: weighted_std(values, weights) for key, values in supports.items()},
        "reconstruction_oof_weighted_R2": reconstruction,
    }


def static_features(instances: list[Instance]) -> dict[str, float]:
    p69 = load_p69(find_paper69_dir())
    adapted = [
        p69.Instance(
            dataset_id=i.dataset_id, family=i.family, source_name="SATLIB", source_url="",
            license_or_terms="source terms", download_tls_mode="recorded", archive_sha256="",
            local_path=i.local_path, sha256=i.sha256, n_vars=i.n_vars, n_clauses=i.n_clauses,
            clauses=i.clauses,
        ) for i in instances
    ]
    return {feature.dataset_id: feature.G_gate for feature in p69.build_features(adapted)}


def donor_map(instances: list[Instance], traces: dict[str, list[TraceState]]) -> dict[str, list[float]]:
    donors: dict[str, list[float]] = {}
    groups: dict[tuple[str, str], list[Instance]] = {}
    for instance in instances:
        groups.setdefault((instance.family, instance.split), []).append(instance)
    for group in groups.values():
        ordered = sorted(group, key=lambda i: hashlib.sha256(f"75075|{i.dataset_id}".encode()).hexdigest())
        for index, target in enumerate(ordered):
            donor = ordered[(index + 1) % len(ordered)]
            donors[target.dataset_id] = [state.G_variants.get("HBSV", 0.0) for state in traces.get(donor.dataset_id, [])]
    return donors


def family_bootstrap(delta: pd.Series, families: pd.Series, n_boot: int = 10000) -> dict[str, float]:
    frame = pd.DataFrame({"delta": delta.to_numpy(float), "family": families.to_numpy(str)})
    family_names = sorted(frame["family"].unique())
    observed = float(np.mean([frame.loc[frame.family == family, "delta"].mean() for family in family_names]))
    rng = np.random.default_rng(75075)
    samples = np.empty(n_boot)
    for iteration in range(n_boot):
        family_means = []
        for family in family_names:
            values = frame.loc[frame.family == family, "delta"].to_numpy(float)
            family_means.append(float(np.mean(rng.choice(values, size=len(values), replace=True))))
        samples[iteration] = np.mean(family_means)
    low, high = np.quantile(samples, [0.025, 0.975])
    return {"mean": observed, "ci95_low": float(low), "ci95_high": float(high)}


def comparisons(records: pd.DataFrame) -> pd.DataFrame:
    test = records[records.split == "test"]
    pivot = test.pivot(index=["dataset_id", "family"], columns="policy", values="total_cost")
    rows: list[dict[str, Any]] = []
    for comparator in COMPARATORS:
        if PRIMARY not in pivot or comparator not in pivot:
            continue
        pair = pivot[[PRIMARY, comparator]].dropna().reset_index()
        delta = pair[comparator] - pair[PRIMARY]
        boot = family_bootstrap(delta, pair["family"])
        rows.append({
            "comparison": f"{PRIMARY}_vs_{comparator}", **boot,
            "wins_for_state_v18": int((delta > 0).sum()),
            "losses_for_state_v18": int((delta < 0).sum()),
            "ties": int((delta == 0).sum()), "n_pairs": len(delta),
        })
    return pd.DataFrame(rows)


def status_summary(records: pd.DataFrame, comp: pd.DataFrame, construct: str) -> dict[str, Any]:
    primary_rows = records[(records.split == "test") & (records.policy == PRIMARY)]
    family_activation = primary_rows.groupby("family").activation_rate.mean()
    activation_rate = float(family_activation.mean())
    activation_status = "adequate_activation" if activation_rate >= 0.125 else "activation_collapse"
    lookup = comp.set_index("comparison").to_dict("index")
    def row(comparator: str) -> dict[str, float]:
        return lookup.get(f"{PRIMARY}_vs_{comparator}", {"ci95_low": math.nan, "ci95_high": math.nan})
    moms = records[(records.split == "test") & (records.policy == "moms")][["dataset_id", "family", "total_cost"]]
    primary = primary_rows[["dataset_id", "family", "total_cost"]]
    relative = primary.merge(moms, on=["dataset_id", "family"], suffixes=("_primary", "_moms"))
    regret = (relative.total_cost_primary - relative.total_cost_moms) / relative.total_cost_moms.clip(lower=EPS)
    safety_boot = family_bootstrap(regret, relative.family)
    safety = classify_safety(safety_boot["ci95_low"], safety_boot["ci95_high"])
    min_support = all(row(name)["ci95_low"] > 0 for name in [
        "static_gate_v17_Paper69", "score_with_R_Paper69",
        "state_v18_gate_signal_shuffled", "random_activation_matched_q",
    ])
    negative = any(row(name)["ci95_high"] <= 0 for name in ["static_gate_v17_Paper69", "score_with_R_Paper69"])
    strong = min_support and activation_status == "adequate_activation" and safety == "no_harm_certified" and all(
        row(name)["ci95_low"] > 0 for name in ["moms", "jeroslow_wang"]
    )
    eligible_ids = set(primary_rows.loc[primary_rows.activation_rate >= 0.125, "dataset_id"])
    routed = False
    if eligible_ids:
        eligible = records[(records.dataset_id.isin(eligible_ids)) & (records.policy.isin([PRIMARY, "moms"]))]
        ep = eligible.pivot(index=["dataset_id", "family"], columns="policy", values="total_cost").dropna().reset_index()
        if not ep.empty:
            routed = family_bootstrap(ep.moms - ep[PRIMARY], ep.family)["ci95_low"] > 0
    if strong:
        utility = "strong_sat_support"
    elif min_support and activation_status == "adequate_activation" and safety == "no_harm_certified" and routed:
        utility = "routing_support"
    elif min_support:
        utility = "minimum_support"
    elif negative:
        utility = "negative_result"
    else:
        utility = "inconclusive"
    return {
        "execution_status": "executed",
        "construct_status": construct,
        "activation_status": activation_status,
        "safety_status": safety,
        "utility_status": utility,
        "family_instance_weighted_activation_rate": activation_rate,
        "safety_bootstrap_relative_regret_vs_MOMS": safety_boot,
    }


def run_execution(data_dir: Path, output_dir: Path) -> dict[str, Any]:
    validation = validate_dataset(data_dir, output_dir)
    if not validation["dataset_ready"]:
        raise RuntimeError("Paper-76 dataset gatekeeper failed; execution refused")
    instances = load_instances(data_dir)
    calibration_instances = [instance for instance in instances if instance.split == "calibration"]
    test_instances = [instance for instance in instances if instance.split == "test"]

    calibration_records: list[SolveRecord] = []
    calibration_traces: dict[str, list[TraceState]] = {}
    for index, instance in enumerate(calibration_instances, 1):
        record, trace = run_solver(instance, "moms", None, collect_trace=True)
        calibration_records.append(record)
        calibration_traces[instance.dataset_id] = trace
        print(f"[cal {index}/{len(calibration_instances)}] {instance.dataset_id}", flush=True)
    retained = retain_calibration_states(calibration_traces, instances)
    calibration = calibrate(retained)
    precondition = preconditions(retained, calibration)
    calibration["preconditions"] = precondition
    calibration["retained_states"] = len(retained)
    calibration["paper75_doi"] = PAPER75_DOI
    (output_dir / "paper76_calibration.json").write_text(json.dumps(calibration, indent=2), encoding="utf-8")
    pd.DataFrame([
        {**asdict(state), "predictors": json.dumps(state.predictors), "R_variants": json.dumps(state.R_variants),
         "delta_variants": json.dumps(state.delta_variants), "G_variants": json.dumps(state.G_variants)}
        for state in retained
    ]).to_csv(output_dir / "paper76_calibration_states.csv", index=False)
    if precondition["construct_status"] == "not_testable_as_constructed":
        summary = {
            "paper": 76, "paper75_preregistration_doi": PAPER75_DOI,
            "execution_status": "executed", **precondition,
            "conclusion": "Primary utility execution stopped by the frozen distinctness precondition.",
        }
        (output_dir / "paper76_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        return summary

    kappa = sum(record.search_cost for record in calibration_records) / max(EPS, sum(record.runtime_s for record in calibration_records))
    static = static_features(instances)

    # Required MOMS traces are executed first because the shuffled null is a
    # frozen cyclic donor mapping of these already-required trajectories.
    records: list[SolveRecord] = calibration_records.copy()
    moms_test_traces: dict[str, list[TraceState]] = {}
    for index, instance in enumerate(test_instances, 1):
        record, trace = run_solver(instance, "moms", calibration, collect_trace=True, kappa=kappa)
        records.append(record)
        moms_test_traces[instance.dataset_id] = trace
        print(f"[test MOMS {index}/{len(test_instances)}] {instance.dataset_id}", flush=True)
    donors = donor_map(test_instances, moms_test_traces)

    policies = [policy for policy in POLICIES if policy != "moms"] + ABLATIONS
    total = len(test_instances) * len(policies)
    counter = 0
    for instance in test_instances:
        for policy in policies:
            counter += 1
            record, _ = run_solver(
                instance, policy, calibration, static_gate=static[instance.dataset_id],
                donor_g=donors.get(instance.dataset_id), kappa=kappa,
            )
            records.append(record)
            print(f"[policy {counter}/{total}] {instance.dataset_id} :: {policy}", flush=True)

    frame = pd.DataFrame([asdict(record) for record in records])
    frame.to_csv(output_dir / "paper76_solve_records.csv", index=False)
    policy_summary = frame.groupby(["split", "policy"], as_index=False).agg(
        n_instances=("dataset_id", "count"), solve_rate=("solved", "mean"),
        mean_total_cost=("total_cost", "mean"), median_total_cost=("total_cost", "median"),
        mean_runtime_s=("runtime_s", "mean"), mean_activation_rate=("activation_rate", "mean"),
        mean_overhead_fraction=("structural_overhead_s", lambda s: float(np.mean(s / frame.loc[s.index, "runtime_s"].clip(lower=EPS)))),
        timeouts=("timeout", "sum"),
    )
    policy_summary.to_csv(output_dir / "paper76_policy_summary.csv", index=False)
    comp = comparisons(frame)
    comp.to_csv(output_dir / "paper76_primary_comparisons.csv", index=False)
    family = frame[frame.split == "test"].groupby(["family", "policy"], as_index=False).agg(
        mean_total_cost=("total_cost", "mean"), activation_rate=("activation_rate", "mean"),
        solve_rate=("solved", "mean"), timeouts=("timeout", "sum"),
    )
    family.to_csv(output_dir / "paper76_family_results.csv", index=False)
    statuses = status_summary(frame, comp, precondition["construct_status"])
    summary = {
        "paper": 76,
        "title": "State-Conditioned Gate Arbitration Execution",
        "paper75_preregistration_doi": PAPER75_DOI,
        "paper75_protocol_sha256": PAPER75_PROTOCOL_SHA256,
        "manifest_sha256": PAPER75_MANIFEST_SHA256,
        "n_instances": len(instances), "n_calibration_instances": len(calibration_instances),
        "n_test_instances": len(test_instances), "kappa_cal": kappa,
        "preconditions": precondition, "status_axes": statuses,
        "primary_comparisons": comp.to_dict(orient="records"),
        "interpretation_rule": "Only the independent Paper-75 status axes are permitted; no post-hoc retuning or outcome relabelling.",
    }
    (output_dir / "paper76_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (output_dir / "paper76_run_log.json").write_text(json.dumps({
        "status": "executed", "paper75_doi": PAPER75_DOI,
        "protocol_hash": PAPER75_PROTOCOL_SHA256, "manifest_hash": PAPER75_MANIFEST_SHA256,
        "policies": policies, "primary_policy": PRIMARY,
    }, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data_external")
    parser.add_argument("--output-dir", type=Path, default=OUTDIR)
    parser.add_argument("--paper69-dir", type=Path)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    load_p69(find_paper69_dir(args.paper69_dir))
    validation = validate_dataset(args.data_dir, args.output_dir)
    if args.validate_only:
        print(json.dumps(validation, indent=2))
        raise SystemExit(0 if validation["dataset_ready"] else 1)
    summary = run_execution(args.data_dir, args.output_dir)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
