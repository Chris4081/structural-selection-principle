#!/usr/bin/env python3
"""
Paper 63: Structure-gated Mode-A CDCL branching.

Purpose
-------
Paper 62 showed that global h=0 Mode A is not universally compute-positive:
it is useful on some structured families and harmful on others. This benchmark
tests the next v1.6-T hypothesis: gate the structural term by an instance-level
structure diagnostic before inserting it into the CDCL branching policy.

This is not an industrial SAT solver and not a claim against modern CDCL
implementations. The benchmark is deliberately small, reproducible, and
auditable: all policies share the same propagation, learning, restart-free
CDCL loop and differ only in the decision policy.

Tested Paper-63 policy:
    A_g(c) = prog(c) - g(F) * tau * h_MAAT(c)

where g(F) is computed from root-level clause-variable graph diagnostics,
without using family labels.
"""

from __future__ import annotations

import json
import math
import os
import random
import time
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable

_CACHE_ROOT = Path(os.environ.get("TMPDIR", "/tmp")) / "maat_t1_cdcl_cache"
(_CACHE_ROOT / "matplotlib").mkdir(parents=True, exist_ok=True)
(_CACHE_ROOT / "xdg").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE_ROOT / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE_ROOT / "xdg"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SEED = 63063
OUTDIR = Path(__file__).resolve().parent / "outputs_paper63_structure_gated_mode_a"
EPS = 1.0e-9

# Paper-63 predeclared policy parameters.
CONFLICT_BUDGET = 3500
BETA = 8.0
TAU = 0.38
HORIZON_H = 0

POLICIES = ["vsids", "progress_only", "maat_only", "mode_a", "mode_a_gated", "random"]


@dataclass
class SolveRecord:
    family: str
    instance_id: int
    policy: str
    n_vars: int
    n_clauses: int
    alpha: float
    solved: int
    sat: int
    timeout: int
    decisions: int
    propagations: int
    conflicts: int
    learned_clauses: int
    max_level: int
    structure_gate: float
    runtime_s: float
    compute_cost: float


@dataclass
class GateRecord:
    family: str
    instance_id: int
    n_vars: int
    n_clauses: int
    alpha: float
    structure_gate: float
    clustering_mean: float
    degree_cv: float
    v_std: float
    tail_degree_pressure: float
    symmetry_penalty: float


def canonical_clause(clause: Iterable[int]) -> tuple[int, ...] | None:
    lits = tuple(sorted(set(int(x) for x in clause), key=lambda x: (abs(x), x)))
    seen: dict[int, int] = {}
    for lit in lits:
        v = abs(lit)
        if v in seen and seen[v] != lit:
            return None
        seen[v] = lit
    return lits if lits else tuple()


def generate_random_3sat(rng: np.random.Generator, n_vars: int, n_clauses: int) -> list[tuple[int, ...]]:
    clauses: set[tuple[int, ...]] = set()
    attempts = 0
    while len(clauses) < n_clauses and attempts < 150 * n_clauses:
        attempts += 1
        vars_ = rng.choice(np.arange(1, n_vars + 1), size=3, replace=False)
        lits = [int((-1 if rng.random() < 0.5 else 1) * v) for v in vars_]
        c = canonical_clause(lits)
        if c is not None and len(c) == 3:
            clauses.add(c)
    return list(clauses)


def generate_planted_3sat(rng: np.random.Generator, n_vars: int, n_clauses: int) -> list[tuple[int, ...]]:
    assignment = rng.random(n_vars) < 0.5
    clauses: set[tuple[int, ...]] = set()
    attempts = 0
    while len(clauses) < n_clauses and attempts < 180 * n_clauses:
        attempts += 1
        vars_ = rng.choice(np.arange(1, n_vars + 1), size=3, replace=False)
        lits = [int((-1 if rng.random() < 0.5 else 1) * v) for v in vars_]
        if not any((lit > 0 and assignment[abs(lit) - 1]) or (lit < 0 and not assignment[abs(lit) - 1]) for lit in lits):
            j = int(rng.integers(0, 3))
            lits[j] = -lits[j]
        c = canonical_clause(lits)
        if c is not None and len(c) == 3:
            clauses.add(c)
    return list(clauses)


def generate_modular_3sat(rng: np.random.Generator, n_vars: int, n_clauses: int) -> list[tuple[int, ...]]:
    n_blocks = max(3, min(6, n_vars // 6))
    blocks = np.array_split(np.arange(1, n_vars + 1), n_blocks)
    clauses: set[tuple[int, ...]] = set()
    attempts = 0
    while len(clauses) < n_clauses and attempts < 180 * n_clauses:
        attempts += 1
        if rng.random() < 0.80:
            block = blocks[int(rng.integers(0, len(blocks)))]
            pool = block if len(block) >= 3 else np.arange(1, n_vars + 1)
            vars_ = rng.choice(pool, size=3, replace=False)
        else:
            vars_ = rng.choice(np.arange(1, n_vars + 1), size=3, replace=False)
        lits = [int((-1 if rng.random() < 0.5 else 1) * v) for v in vars_]
        c = canonical_clause(lits)
        if c is not None and len(c) == 3:
            clauses.add(c)
    return list(clauses)


def xor3_clauses(a: int, b: int, c: int, rhs: int) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []
    for va in (0, 1):
        for vb in (0, 1):
            for vc in (0, 1):
                if (va ^ vb ^ vc) != rhs:
                    out.append(tuple([(-a if va else a), (-b if vb else b), (-c if vc else c)]))
    return out


def generate_xor3_cnf(rng: np.random.Generator, n_vars: int, n_xors: int) -> list[tuple[int, ...]]:
    clauses: list[tuple[int, ...]] = []
    used: set[tuple[int, int, int, int]] = set()
    attempts = 0
    while len(used) < n_xors and attempts < 150 * n_xors:
        attempts += 1
        a, b, c = sorted(rng.choice(np.arange(1, n_vars + 1), size=3, replace=False))
        rhs = int(rng.integers(0, 2))
        key = (int(a), int(b), int(c), rhs)
        if key in used:
            continue
        used.add(key)
        clauses.extend(xor3_clauses(int(a), int(b), int(c), rhs))
    return clauses


def generate_coloring_cnf(rng: np.random.Generator, n_vertices: int, edge_prob: float, n_colors: int) -> tuple[int, list[tuple[int, ...]]]:
    def var(v: int, color: int) -> int:
        return v * n_colors + color + 1

    clauses: list[tuple[int, ...]] = []
    for v in range(n_vertices):
        clauses.append(tuple(var(v, c) for c in range(n_colors)))
        for c1, c2 in combinations(range(n_colors), 2):
            clauses.append((-var(v, c1), -var(v, c2)))
    for u, v in combinations(range(n_vertices), 2):
        if rng.random() < edge_prob:
            for c in range(n_colors):
                clauses.append((-var(u, c), -var(v, c)))
    return n_vertices * n_colors, clauses


def generate_pigeonhole_cnf(n_holes: int) -> tuple[int, list[tuple[int, ...]]]:
    n_pigeons = n_holes + 1

    def var(p: int, h: int) -> int:
        return p * n_holes + h + 1

    clauses: list[tuple[int, ...]] = []
    for p in range(n_pigeons):
        clauses.append(tuple(var(p, h) for h in range(n_holes)))
    for h in range(n_holes):
        for p1, p2 in combinations(range(n_pigeons), 2):
            clauses.append((-var(p1, h), -var(p2, h)))
    return n_pigeons * n_holes, clauses


def build_instances(rng: np.random.Generator) -> list[tuple[str, int, list[tuple[int, ...]]]]:
    instances: list[tuple[str, int, list[tuple[int, ...]]]] = []
    for n in [18, 22, 26]:
        for alpha in [3.9, 4.25, 4.6]:
            m = int(round(alpha * n))
            for _ in range(3):
                instances.append(("random_3sat", n, generate_random_3sat(rng, n, m)))
                instances.append(("planted_3sat", n, generate_planted_3sat(rng, n, m)))
    for n in [18, 24]:
        for alpha in [4.0, 4.35, 4.7]:
            m = int(round(alpha * n))
            for _ in range(3):
                instances.append(("modular_3sat", n, generate_modular_3sat(rng, n, m)))
    for n in [15, 18, 21]:
        for ratio in [0.85, 1.0, 1.15]:
            for _ in range(2):
                instances.append(("xor3_cnf", n, generate_xor3_cnf(rng, n, int(round(ratio * n)))))
    for vertices in [8, 10, 12]:
        for prob in [0.30, 0.45]:
            n_vars, clauses = generate_coloring_cnf(rng, vertices, prob, 3)
            instances.append(("graph_coloring", n_vars, clauses))
    for holes in [4, 5, 6]:
        n_vars, clauses = generate_pigeonhole_cnf(holes)
        instances.append(("pigeonhole", n_vars, clauses))
    return instances


def formula_structure_gate(clauses: list[tuple[int, ...]], n_vars: int) -> dict[str, float]:
    """Root-level, label-free gate for structural steering.

    The gate is intentionally simple and predeclared:

    - community signal: mean clustering in the variable co-occurrence graph,
    - V-mode spread: variance of local connectivity support,
    - tail pressure: high-degree tail relative to mean degree,
    - symmetry penalty: suppresses globally uniform instances.

    The output g(F) is in [0,1]. Low values reduce Mode A to progress-only;
    high values activate the structural penalty.
    """

    deg = np.zeros(n_vars, dtype=float)
    neighbors: list[set[int]] = [set() for _ in range(n_vars)]
    for clause in clauses:
        vars_ = sorted({abs(lit) - 1 for lit in clause if 1 <= abs(lit) <= n_vars})
        for v in vars_:
            deg[v] += 1.0
        for a, b in combinations(vars_, 2):
            neighbors[a].add(b)
            neighbors[b].add(a)

    graph_deg = np.array([len(n) for n in neighbors], dtype=float)
    max_graph_deg = float(np.max(graph_deg)) + EPS
    v_local = graph_deg / max_graph_deg
    v_std = float(np.std(v_local))
    degree_cv = float(np.std(deg) / (np.mean(deg) + EPS))

    clustering_vals = []
    for v in range(n_vars):
        nb = list(neighbors[v])
        if len(nb) < 2:
            clustering_vals.append(0.0)
            continue
        possible = len(nb) * (len(nb) - 1) / 2.0
        actual = sum(1 for a, b in combinations(nb, 2) if b in neighbors[a])
        clustering_vals.append(actual / (possible + EPS))
    clustering_mean = float(np.mean(clustering_vals)) if clustering_vals else 0.0

    q90 = float(np.quantile(deg, 0.90)) if len(deg) else 0.0
    tail_degree_pressure = float(min(1.0, max(0.0, (q90 / (np.mean(deg) + EPS) - 1.0) / 1.65)))

    community_signal = clustering_mean / (clustering_mean + 0.18 + EPS)
    v_spread_signal = min(1.0, v_std / 0.22)
    symmetry_penalty = math.exp(-degree_cv / 0.20) * math.exp(-v_std / 0.12)

    raw_gate = (
        0.52 * community_signal
        + 0.30 * v_spread_signal
        + 0.18 * tail_degree_pressure
        - 0.55 * symmetry_penalty
    )
    gate = float(np.clip(raw_gate, 0.0, 1.0))
    return {
        "structure_gate": gate,
        "clustering_mean": clustering_mean,
        "degree_cv": degree_cv,
        "v_std": v_std,
        "tail_degree_pressure": tail_degree_pressure,
        "symmetry_penalty": float(symmetry_penalty),
    }


class CDCLSolver:
    def __init__(
        self,
        clauses: list[tuple[int, ...]],
        n_vars: int,
        policy: str,
        seed: int,
        structure_gate: float,
        conflict_budget: int = CONFLICT_BUDGET,
        beta: float = BETA,
        tau: float = TAU,
    ) -> None:
        self.clauses = list(clauses)
        self.n_original = len(clauses)
        self.n_vars = n_vars
        self.policy = policy
        self.rng = random.Random(seed)
        self.structure_gate = float(np.clip(structure_gate, 0.0, 1.0))
        self.conflict_budget = conflict_budget
        self.beta = beta
        self.tau = tau

        self.assignment: dict[int, bool] = {}
        self.level_of: dict[int, int] = {}
        self.reason: dict[int, tuple[int, ...] | None] = {}
        self.trail: list[int] = []
        self.decision_level = 0
        self.phase: dict[int, bool] = {}

        self.vsids = {v: 1.0 for v in range(1, n_vars + 1)}
        for clause in self.clauses:
            for lit in clause:
                self.vsids[abs(lit)] += 0.02
        self.vsids_decay = 0.95

        self.decisions = 0
        self.propagations = 0
        self.conflicts = 0
        self.learned_clauses = 0
        self.max_level = 0

    def lit_value(self, lit: int) -> bool | None:
        val = self.assignment.get(abs(lit))
        if val is None:
            return None
        return val if lit > 0 else not val

    def enqueue(self, lit: int, reason: tuple[int, ...] | None) -> bool:
        var = abs(lit)
        val = lit > 0
        old = self.assignment.get(var)
        if old is not None:
            return old == val
        self.assignment[var] = val
        self.level_of[var] = self.decision_level
        self.reason[var] = reason
        self.phase[var] = val
        self.trail.append(lit)
        if reason is not None:
            self.propagations += 1
        return True

    def clause_state(self, clause: tuple[int, ...]) -> tuple[str, list[int]]:
        unassigned: list[int] = []
        for lit in clause:
            val = self.lit_value(lit)
            if val is True:
                return "sat", []
            if val is None:
                unassigned.append(lit)
        if not unassigned:
            return "conflict", []
        return "open", unassigned

    def propagate(self) -> tuple[int, ...] | None:
        changed = True
        while changed:
            changed = False
            for clause in self.clauses:
                state, unassigned = self.clause_state(clause)
                if state == "sat":
                    continue
                if state == "conflict":
                    return clause
                if len(unassigned) == 1:
                    lit = unassigned[0]
                    if not self.enqueue(lit, clause):
                        return clause
                    changed = True
        return None

    def is_satisfied(self) -> bool:
        return all(self.clause_state(c)[0] == "sat" for c in self.clauses)

    def active_clauses(self) -> list[list[int]]:
        active: list[list[int]] = []
        for clause in self.clauses:
            state, unassigned = self.clause_state(clause)
            if state == "sat":
                continue
            active.append(unassigned)
        return active

    def normalize_map(self, values: dict[int, float], default: float = 0.0) -> dict[int, float]:
        if not values:
            return {}
        lo = min(values.values())
        hi = max(values.values())
        if hi - lo <= EPS:
            return {k: default for k in values}
        return {k: (v - lo) / (hi - lo + EPS) for k, v in values.items()}

    def score_candidates(self) -> dict[int, dict[str, float]]:
        unassigned = [v for v in range(1, self.n_vars + 1) if v not in self.assignment]
        active = [c for c in self.active_clauses() if c]
        pos = {v: 0.0 for v in unassigned}
        neg = {v: 0.0 for v in unassigned}
        jw_pos = {v: 0.0 for v in unassigned}
        jw_neg = {v: 0.0 for v in unassigned}
        short = {v: 0.0 for v in unassigned}
        unit_gain_pos = {v: 0.0 for v in unassigned}
        unit_gain_neg = {v: 0.0 for v in unassigned}
        neigh: dict[int, set[int]] = {v: set() for v in unassigned}

        for clause in active:
            weight = 2.0 ** (-len(clause))
            vars_in = [abs(lit) for lit in clause]
            for lit in clause:
                v = abs(lit)
                if v not in pos:
                    continue
                if lit > 0:
                    pos[v] += 1.0
                    jw_pos[v] += weight
                    if len(clause) == 2:
                        unit_gain_pos[v] += 1.0
                else:
                    neg[v] += 1.0
                    jw_neg[v] += weight
                    if len(clause) == 2:
                        unit_gain_neg[v] += 1.0
                if len(clause) <= 3:
                    short[v] += 1.0 / len(clause)
                neigh[v].update(u for u in vars_in if u != v)

        max_degree = max((pos[v] + neg[v] for v in unassigned), default=1.0) + EPS
        max_jw = max((jw_pos[v] + jw_neg[v] for v in unassigned), default=1.0) + EPS
        max_short = max((short[v] for v in unassigned), default=1.0) + EPS
        max_neigh = max((len(neigh[v]) for v in unassigned), default=1) + EPS
        vs_norm = self.normalize_map({v: self.vsids[v] for v in unassigned}, default=0.5)

        raw: dict[int, dict[str, float]] = {}
        for v in unassigned:
            degree = pos[v] + neg[v]
            sign_total = degree + EPS
            imbalance = abs(pos[v] - neg[v]) / sign_total
            H = 1.0 / (1.0 + imbalance)
            B = 1.0 - imbalance
            S = short[v] / max_short
            V = len(neigh[v]) / max_neigh
            R_resp = (max(EPS, H) * max(EPS, B) * max(EPS, V)) ** (1.0 / 3.0)
            R_rob = min(R_resp, (max(EPS, H) * max(EPS, B) * max(EPS, S) * max(EPS, V)) ** 0.25)
            support_cost = -sum(math.log(EPS + q) for q in (H, B, S, V, R_rob))
            jw = (jw_pos[v] + jw_neg[v]) / max_jw
            degree_norm = degree / max_degree
            unit_gain = max(unit_gain_pos[v], unit_gain_neg[v]) / (max_degree + EPS)
            progress = 0.52 * vs_norm[v] + 0.26 * jw + 0.14 * unit_gain + 0.08 * degree_norm
            raw[v] = {
                "degree": degree,
                "progress": progress,
                "support_cost": support_cost,
                "mode_a": progress - self.tau * support_cost,
                "polarity_true_score": pos[v] + jw_pos[v] + unit_gain_pos[v],
                "polarity_false_score": neg[v] + jw_neg[v] + unit_gain_neg[v],
            }

        # Normalize structural cost across the current decision frontier for the
        # combined policy. This keeps tau interpretable across families.
        if raw:
            costs = {v: raw[v]["support_cost"] for v in raw}
            cost_norm = self.normalize_map(costs, default=0.5)
            for v in raw:
                raw[v]["h_maat"] = cost_norm[v]
                raw[v]["mode_a"] = raw[v]["progress"] - self.tau * cost_norm[v]
                raw[v]["mode_a_gated"] = raw[v]["progress"] - (self.structure_gate * self.tau) * cost_norm[v]
        return raw

    def softmax_choice(self, scores: dict[int, float]) -> int:
        if not scores:
            raise RuntimeError("cannot choose from empty score set")
        keys = list(scores.keys())
        vals = np.array([scores[k] for k in keys], dtype=float)
        vals = vals - float(np.max(vals))
        weights = np.exp(self.beta * vals)
        total = float(np.sum(weights))
        if not np.isfinite(total) or total <= EPS:
            return max(keys, key=lambda k: (scores[k], -k))
        probs = weights / total
        pick = self.rng.random()
        acc = 0.0
        for key, prob in zip(keys, probs):
            acc += float(prob)
            if pick <= acc:
                return int(key)
        return int(keys[-1])

    def choose_decision(self) -> int | None:
        unassigned = [v for v in range(1, self.n_vars + 1) if v not in self.assignment]
        if not unassigned:
            return None
        scores = self.score_candidates()
        if self.policy == "random":
            var = self.rng.choice(unassigned)
        elif self.policy == "vsids":
            var = max(unassigned, key=lambda v: (self.vsids[v], scores.get(v, {}).get("degree", 0.0), -v))
        elif self.policy == "progress_only":
            var = max(unassigned, key=lambda v: (scores[v]["progress"], self.vsids[v], -v))
        elif self.policy == "maat_only":
            # Low structural risk is preferred. This ablation is intentionally
            # goal-blind and is expected to fail when coherence and progress decouple.
            var = self.softmax_choice({v: -scores[v]["h_maat"] for v in unassigned})
        elif self.policy == "mode_a":
            var = self.softmax_choice({v: scores[v]["mode_a"] for v in unassigned})
        elif self.policy == "mode_a_gated":
            var = self.softmax_choice({v: scores[v]["mode_a_gated"] for v in unassigned})
        else:
            raise ValueError(f"unknown policy: {self.policy}")

        sc = scores[var]
        polarity = sc["polarity_true_score"] >= sc["polarity_false_score"]
        if var in self.phase and self.rng.random() < 0.25:
            polarity = self.phase[var]
        return var if polarity else -var

    def resolve(self, c1: tuple[int, ...], c2: tuple[int, ...], var: int) -> tuple[int, ...]:
        lits = [lit for lit in c1 if abs(lit) != var] + [lit for lit in c2 if abs(lit) != var]
        c = canonical_clause(lits)
        if c is None:
            return c1
        return c

    def current_level_count(self, clause: tuple[int, ...]) -> int:
        return sum(1 for lit in clause if self.level_of.get(abs(lit), 0) == self.decision_level)

    def analyze_conflict(self, conflict: tuple[int, ...]) -> tuple[tuple[int, ...], int]:
        learned = conflict
        guard = 0
        while self.current_level_count(learned) > 1 and guard < 3 * max(1, len(self.trail)):
            guard += 1
            pivot_var = None
            for lit in reversed(self.trail):
                v = abs(lit)
                if self.level_of.get(v, 0) == self.decision_level and any(abs(x) == v for x in learned):
                    if self.reason.get(v) is not None:
                        pivot_var = v
                        break
            if pivot_var is None:
                break
            learned = self.resolve(learned, self.reason[pivot_var] or tuple(), pivot_var)

        current_lits = [lit for lit in learned if self.level_of.get(abs(lit), 0) == self.decision_level]
        backjump = 0
        for lit in learned:
            lvl = self.level_of.get(abs(lit), 0)
            if lit not in current_lits:
                backjump = max(backjump, lvl)
        if not current_lits:
            backjump = max(0, self.decision_level - 1)
        return learned, backjump

    def backtrack(self, level: int) -> None:
        keep: list[int] = []
        for lit in self.trail:
            v = abs(lit)
            if self.level_of.get(v, 0) <= level:
                keep.append(lit)
            else:
                self.assignment.pop(v, None)
                self.level_of.pop(v, None)
                self.reason.pop(v, None)
        self.trail = keep
        self.decision_level = level

    def bump_vsids(self, clause: tuple[int, ...]) -> None:
        for lit in clause:
            self.vsids[abs(lit)] = self.vsids.get(abs(lit), 1.0) + 1.0
        if self.conflicts % 64 == 0:
            for v in self.vsids:
                self.vsids[v] *= self.vsids_decay

    def solve(self) -> tuple[bool, bool]:
        while True:
            conflict = self.propagate()
            if conflict is None:
                break
            self.conflicts += 1
            if self.decision_level == 0:
                return False, True
            learned, backjump = self.analyze_conflict(conflict)
            if len(learned) == 0:
                return False, True
            self.clauses.append(learned)
            self.learned_clauses += 1
            self.bump_vsids(learned)
            self.backtrack(backjump)

        while True:
            if self.is_satisfied():
                return True, True
            if self.conflicts >= self.conflict_budget:
                return False, False
            decision_lit = self.choose_decision()
            if decision_lit is None:
                return self.is_satisfied(), True

            self.decision_level += 1
            self.max_level = max(self.max_level, self.decision_level)
            self.decisions += 1
            if not self.enqueue(decision_lit, None):
                return False, True

            while True:
                conflict = self.propagate()
                if conflict is None:
                    break
                self.conflicts += 1
                if self.decision_level == 0:
                    return False, True
                learned, backjump = self.analyze_conflict(conflict)
                if len(learned) == 0:
                    return False, True
                self.clauses.append(learned)
                self.learned_clauses += 1
                self.bump_vsids(learned)
                self.backtrack(backjump)
                if self.conflicts >= self.conflict_budget:
                    return False, False


def solve_instance(
    family: str,
    instance_id: int,
    n_vars: int,
    clauses: list[tuple[int, ...]],
    policy: str,
    structure_gate: float,
) -> SolveRecord:
    start = time.perf_counter()
    solver = CDCLSolver(
        clauses,
        n_vars,
        policy,
        seed=SEED + 1009 * instance_id + 37 * POLICIES.index(policy),
        structure_gate=structure_gate,
    )
    sat, solved = solver.solve()
    runtime_s = time.perf_counter() - start
    timeout = int(not solved)
    # Predeclared compute cost: conflicts dominate, decisions are next, propagation
    # is cheap but not free, and unresolved budgets receive a fixed penalty.
    compute_cost = (
        solver.conflicts
        + 0.30 * solver.decisions
        + 0.004 * solver.propagations
        + timeout * CONFLICT_BUDGET
    )
    return SolveRecord(
        family=family,
        instance_id=instance_id,
        policy=policy,
        n_vars=n_vars,
        n_clauses=len(clauses),
        alpha=len(clauses) / max(1, n_vars),
        solved=int(solved),
        sat=int(bool(sat)),
        timeout=timeout,
        decisions=int(solver.decisions),
        propagations=int(solver.propagations),
        conflicts=int(solver.conflicts),
        learned_clauses=int(solver.learned_clauses),
        max_level=int(solver.max_level),
        structure_gate=float(structure_gate),
        runtime_s=float(runtime_s),
        compute_cost=float(compute_cost),
    )


def paired_regret(df: pd.DataFrame) -> pd.DataFrame:
    pivot = df.pivot_table(index=["family", "instance_id"], columns="policy", values="compute_cost")
    solved = df.pivot_table(index=["family", "instance_id"], columns="policy", values="solved")
    rows: list[dict[str, float | str | int]] = []
    for policy in POLICIES:
        if policy == "vsids" or policy not in pivot.columns or "vsids" not in pivot.columns:
            continue
        comp = (pivot[policy] - pivot["vsids"]).dropna()
        ratio = (pivot[policy] / (pivot["vsids"] + EPS)).replace([np.inf, -np.inf], np.nan).dropna()
        wins = int((comp < 0).sum())
        losses = int((comp > 0).sum())
        ties = int((comp == 0).sum())
        rows.append(
            {
                "policy": policy,
                "n_pairs": int(len(comp)),
                "mean_regret_vs_vsids": float(comp.mean()),
                "median_regret_vs_vsids": float(comp.median()),
                "mean_cost_ratio_vs_vsids": float(ratio.mean()),
                "median_cost_ratio_vs_vsids": float(ratio.median()),
                "wins_vs_vsids": wins,
                "losses_vs_vsids": losses,
                "ties_vs_vsids": ties,
                "solve_rate_delta_vs_vsids": float(solved[policy].mean() - solved["vsids"].mean()),
            }
        )
    return pd.DataFrame(rows)


def paired_regret_against(df: pd.DataFrame, baseline: str) -> pd.DataFrame:
    pivot = df.pivot_table(index=["family", "instance_id"], columns="policy", values="compute_cost")
    rows: list[dict[str, float | str | int]] = []
    for policy in POLICIES:
        if policy == baseline or policy not in pivot.columns or baseline not in pivot.columns:
            continue
        comp = (pivot[policy] - pivot[baseline]).dropna()
        ratio = (pivot[policy] / (pivot[baseline] + EPS)).replace([np.inf, -np.inf], np.nan).dropna()
        rows.append(
            {
                "baseline": baseline,
                "policy": policy,
                "n_pairs": int(len(comp)),
                "mean_regret": float(comp.mean()),
                "median_regret": float(comp.median()),
                "mean_cost_ratio": float(ratio.mean()),
                "median_cost_ratio": float(ratio.median()),
                "wins": int((comp < 0).sum()),
                "losses": int((comp > 0).sum()),
                "ties": int((comp == 0).sum()),
            }
        )
    return pd.DataFrame(rows)


def family_regret(df: pd.DataFrame) -> pd.DataFrame:
    pivot = df.pivot_table(index=["family", "instance_id"], columns="policy", values="compute_cost").reset_index()
    rows = []
    for family, sub in pivot.groupby("family"):
        for policy in POLICIES:
            if policy in {"vsids"} or policy not in sub.columns:
                continue
            comp = (sub[policy] - sub["vsids"]).dropna()
            rows.append(
                {
                    "family": family,
                    "policy": policy,
                    "n_pairs": int(len(comp)),
                    "median_regret_vs_vsids": float(comp.median()),
                    "mean_regret_vs_vsids": float(comp.mean()),
                    "win_rate_vs_vsids": float((comp < 0).mean()),
                }
            )
    return pd.DataFrame(rows)


def run_benchmark() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(SEED)
    instances = build_instances(rng)
    rows: list[SolveRecord] = []
    gate_rows: list[GateRecord] = []
    for instance_id, (family, n_vars, clauses) in enumerate(instances):
        gate_info = formula_structure_gate(clauses, n_vars)
        gate_rows.append(
            GateRecord(
                family=family,
                instance_id=instance_id,
                n_vars=n_vars,
                n_clauses=len(clauses),
                alpha=len(clauses) / max(1, n_vars),
                **gate_info,
            )
        )
        for policy in POLICIES:
            rows.append(
                solve_instance(
                    family,
                    instance_id,
                    n_vars,
                    clauses,
                    policy,
                    structure_gate=gate_info["structure_gate"],
                )
            )
    df = pd.DataFrame([asdict(r) for r in rows])
    gates = pd.DataFrame([asdict(r) for r in gate_rows])
    summary = (
        df.groupby("policy", as_index=False)
        .agg(
            solved_rate=("solved", "mean"),
            timeout_rate=("timeout", "mean"),
            median_conflicts=("conflicts", "median"),
            median_decisions=("decisions", "median"),
            median_compute_cost=("compute_cost", "median"),
            mean_compute_cost=("compute_cost", "mean"),
            median_runtime_s=("runtime_s", "median"),
            n_runs=("solved", "size"),
        )
        .sort_values("median_compute_cost")
    )
    return df, gates, summary, paired_regret(df), family_regret(df), paired_regret_against(df, "mode_a")


def plot_outputs(
    df: pd.DataFrame,
    gates: pd.DataFrame,
    summary: pd.DataFrame,
    regret: pd.DataFrame,
    fam_regret: pd.DataFrame,
    regret_vs_mode_a: pd.DataFrame,
) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")
    colors = {
        "vsids": "#37474f",
        "progress_only": "#1976d2",
        "maat_only": "#8e24aa",
        "mode_a": "#00897b",
        "mode_a_gated": "#f57c00",
        "random": "#9e9e9e",
    }

    fig, ax = plt.subplots(figsize=(8.2, 4.7))
    order = [p for p in POLICIES if p in summary["policy"].tolist()]
    data = summary.set_index("policy").loc[order]
    ax.bar(data.index, data["median_compute_cost"], color=[colors[p] for p in data.index])
    ax.set_ylabel("median compute cost")
    ax.set_title("T1 CDCL policy benchmark")
    ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig1_policy_compute_cost.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 4.7))
    reg = regret.set_index("policy")
    order = [p for p in ["progress_only", "maat_only", "mode_a", "mode_a_gated", "random"] if p in reg.index]
    vals = reg.loc[order, "median_regret_vs_vsids"]
    ax.axhline(0.0, color="#333333", lw=1)
    ax.bar(order, vals, color=[colors[p] for p in order])
    ax.set_ylabel("median regret vs VSIDS")
    ax.set_title("Compute-to-solution regret against VSIDS")
    ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig2_regret_vs_vsids.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.6, 5.2))
    pivot = fam_regret.pivot(index="family", columns="policy", values="median_regret_vs_vsids")
    pivot = pivot[[p for p in ["progress_only", "maat_only", "mode_a", "mode_a_gated", "random"] if p in pivot.columns]]
    pivot.plot(kind="bar", ax=ax, color=[colors[p] for p in pivot.columns])
    ax.axhline(0.0, color="#333333", lw=1)
    ax.set_ylabel("median regret vs VSIDS")
    ax.set_title("Family-level regret reveals where Mode A helps or hurts")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig3_family_regret.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    gate_summary = gates.groupby("family", as_index=False)["structure_gate"].median().sort_values("structure_gate")
    ax.bar(gate_summary["family"], gate_summary["structure_gate"], color="#f57c00")
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("median structure gate g(F)")
    ax.set_title("Instance-level gate by SAT family")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig4_structure_gate_by_family.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.8, 4.6))
    rg = regret_vs_mode_a.set_index("policy")
    if "mode_a_gated" in rg.index:
        vals = rg.loc[["mode_a_gated"], ["median_regret", "mean_regret"]].iloc[0]
        ax.bar(["median", "mean"], [vals["median_regret"], vals["mean_regret"]], color=["#ffb74d", "#f57c00"])
        ax.axhline(0.0, color="#333333", lw=1)
        ax.set_ylabel("regret vs global Mode A")
        ax.set_title("Does gating reduce Mode-A compute cost?")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig5_gated_vs_global_mode_a.png", dpi=180)
    plt.close(fig)

    mode = df[df["policy"] == "mode_a"][["family", "instance_id", "compute_cost"]].rename(
        columns={"compute_cost": "mode_a_cost"}
    )
    gated = df[df["policy"] == "mode_a_gated"][
        ["family", "instance_id", "compute_cost", "structure_gate"]
    ].rename(columns={"compute_cost": "gated_cost"})
    pair = gated.merge(mode, on=["family", "instance_id"], how="inner")
    pair["regret_vs_global_mode_a"] = pair["gated_cost"] - pair["mode_a_cost"]
    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    families = sorted(pair["family"].unique())
    cmap = plt.get_cmap("tab10")
    for idx, family in enumerate(families):
        sub = pair[pair["family"] == family]
        ax.scatter(
            sub["structure_gate"],
            sub["regret_vs_global_mode_a"],
            s=34,
            alpha=0.82,
            label=family,
            color=cmap(idx % 10),
            edgecolor="white",
            linewidth=0.4,
        )
    ax.axhline(0.0, color="#333333", lw=1)
    ax.set_xlim(-0.03, 1.03)
    ax.set_xlabel("structure gate g(F)")
    ax.set_ylabel("gated regret vs global Mode A")
    ax.set_title("Gate strength versus repair of global Mode A")
    ax.legend(fontsize=8, ncol=2, frameon=True)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig6_gate_vs_regret.png", dpi=180)
    plt.close(fig)


def write_readme(summary: pd.DataFrame, regret: pd.DataFrame, regret_vs_mode_a: pd.DataFrame, gates: pd.DataFrame) -> None:
    gated = regret[regret["policy"] == "mode_a_gated"]
    mode_line = "not available"
    if not gated.empty:
        r = gated.iloc[0]
        mode_line = (
            f"median regret {r['median_regret_vs_vsids']:.4g}, "
            f"mean cost ratio {r['mean_cost_ratio_vs_vsids']:.4g}, "
            f"wins/losses/ties {int(r['wins_vs_vsids'])}/"
            f"{int(r['losses_vs_vsids'])}/{int(r['ties_vs_vsids'])}"
        )
    gate_line = gates.groupby("family")["structure_gate"].median().round(4).to_dict()
    vs_global = regret_vs_mode_a[regret_vs_mode_a["policy"] == "mode_a_gated"]
    global_line = "not available"
    if not vs_global.empty:
        rr = vs_global.iloc[0]
        global_line = f"median regret {rr['median_regret']:.4g}, mean regret {rr['mean_regret']:.4g}"
    text = f"""# Paper 63: Structure-Gated Mode-A CDCL Branching

This experiment extends Paper 62 by gating the MAAT v1.6-T Mode-A structural
penalty with a root-level instance diagnostic.

Global Mode A is:

```text
A(c) = prog(c) - tau * h_MAAT(c)
```

Structure-gated Mode A is:

```text
A_g(c) = prog(c) - g(F) * tau * h_MAAT(c)
```

Predeclared parameters:

- seed: `{SEED}`
- conflict budget: `{CONFLICT_BUDGET}`
- beta: `{BETA}`
- tau: `{TAU}`
- rollout horizon h: `{HORIZON_H}`
- policies: `{', '.join(POLICIES)}`

Important scope note: this is not an industrial SAT solver and not a comparison
against CaDiCaL, Kissat, Glucose, or MiniSat. It is a controlled internal policy
test in which VSIDS, global Mode A, and gated Mode A share the same propagation
and conflict-learning machinery.

Headline gated Mode-A result vs VSIDS: {mode_line}.

Gated Mode A vs global Mode A: {global_line}.

Median gate by family: `{gate_line}`.

Outputs:

- `paper63_cdcl_runs.csv`: all per-instance runs
- `paper63_structure_gates.csv`: root-level structure-gate diagnostics
- `paper63_summary_by_policy.csv`: aggregate policy statistics
- `paper63_regret_vs_vsids.csv`: paired regret against VSIDS
- `paper63_regret_vs_mode_a.csv`: paired regret against global Mode A
- `paper63_family_regret_vs_vsids.csv`: family-level paired regret
- `paper63_summary.json`: metadata and compact result table
- `fig1_policy_compute_cost.png`
- `fig2_regret_vs_vsids.png`
- `fig3_family_regret.png`
- `fig4_structure_gate_by_family.png`
- `fig5_gated_vs_global_mode_a.png`
- `fig6_gate_vs_regret.png`

Interpretation rule:

- Negative regret means the policy used less compute than VSIDS on the paired
  instance.
- Positive regret means VSIDS was better.
- A useful gate should reduce global Mode-A tail regret without merely
  collapsing to `progress_only` everywhere.
"""
    (OUTDIR / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df, gates, summary, regret, fam_regret, regret_vs_mode_a = run_benchmark()
    df.to_csv(OUTDIR / "paper63_cdcl_runs.csv", index=False)
    gates.to_csv(OUTDIR / "paper63_structure_gates.csv", index=False)
    summary.to_csv(OUTDIR / "paper63_summary_by_policy.csv", index=False)
    regret.to_csv(OUTDIR / "paper63_regret_vs_vsids.csv", index=False)
    regret_vs_mode_a.to_csv(OUTDIR / "paper63_regret_vs_mode_a.csv", index=False)
    fam_regret.to_csv(OUTDIR / "paper63_family_regret_vs_vsids.csv", index=False)
    plot_outputs(df, gates, summary, regret, fam_regret, regret_vs_mode_a)
    write_readme(summary, regret, regret_vs_mode_a, gates)

    payload = {
        "experiment": "Paper 63 Structure-gated Mode-A CDCL branching",
        "seed": SEED,
        "conflict_budget": CONFLICT_BUDGET,
        "beta": BETA,
        "tau": TAU,
        "horizon_h": HORIZON_H,
        "n_instances": int(df[["family", "instance_id"]].drop_duplicates().shape[0]),
        "policies": POLICIES,
        "summary_by_policy": summary.to_dict(orient="records"),
        "paired_regret_vs_vsids": regret.to_dict(orient="records"),
        "paired_regret_vs_mode_a": regret_vs_mode_a.to_dict(orient="records"),
        "gate_by_family_median": gates.groupby("family")["structure_gate"].median().to_dict(),
    }
    (OUTDIR / "paper63_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("Paper 63 structure-gated Mode-A benchmark complete.")
    print(summary.to_string(index=False))
    print()
    print(regret.to_string(index=False))
    print()
    print(regret_vs_mode_a.to_string(index=False))
    print()
    print(gates.groupby("family")["structure_gate"].median().to_string())
    print(f"\nOutputs written to: {OUTDIR}")


if __name__ == "__main__":
    main()
