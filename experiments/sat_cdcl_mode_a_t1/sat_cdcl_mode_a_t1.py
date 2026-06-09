#!/usr/bin/env python3
"""
MAAT v1.6-T / T1: Mode-A CDCL branching pilot.

Purpose
-------
This script is the first policy-level benchmark for the MAAT v1.6 search
geometry layer. Earlier SAT experiments measured hardness prediction or
DPLL-style branching. Here the structural term is inserted directly into a
small transparent CDCL solver and measured by compute-to-solution regret
against a VSIDS baseline.

This is not an industrial SAT solver and not a claim against modern CDCL
implementations. The benchmark is deliberately small, reproducible, and
auditable: all policies share the same propagation, learning, restart-free
CDCL loop and differ only in the decision policy.

Predeclared T1 policy:
    A(c) = prog(c) - tau * h_MAAT(c)

where prog(c) is a normalized CDCL-local progress score and h_MAAT(c) is a
bounded structural risk score computed from local SAT support coordinates.
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


SEED = 62062
OUTDIR = Path(__file__).resolve().parent / "outputs_t1_cdcl_mode_a"
EPS = 1.0e-9

# T1 predeclared policy parameters.
CONFLICT_BUDGET = 3500
BETA = 8.0
TAU = 0.38
HORIZON_H = 0

POLICIES = ["vsids", "progress_only", "maat_only", "mode_a", "random"]


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
    runtime_s: float
    compute_cost: float


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


class CDCLSolver:
    def __init__(
        self,
        clauses: list[tuple[int, ...]],
        n_vars: int,
        policy: str,
        seed: int,
        conflict_budget: int = CONFLICT_BUDGET,
        beta: float = BETA,
        tau: float = TAU,
    ) -> None:
        self.clauses = list(clauses)
        self.n_original = len(clauses)
        self.n_vars = n_vars
        self.policy = policy
        self.rng = random.Random(seed)
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
) -> SolveRecord:
    start = time.perf_counter()
    solver = CDCLSolver(clauses, n_vars, policy, seed=SEED + 1009 * instance_id + 37 * POLICIES.index(policy))
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


def run_benchmark() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(SEED)
    instances = build_instances(rng)
    rows: list[SolveRecord] = []
    for instance_id, (family, n_vars, clauses) in enumerate(instances):
        for policy in POLICIES:
            rows.append(solve_instance(family, instance_id, n_vars, clauses, policy))
    df = pd.DataFrame([asdict(r) for r in rows])
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
    return df, summary, paired_regret(df), family_regret(df)


def plot_outputs(df: pd.DataFrame, summary: pd.DataFrame, regret: pd.DataFrame, fam_regret: pd.DataFrame) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")
    colors = {
        "vsids": "#37474f",
        "progress_only": "#1976d2",
        "maat_only": "#8e24aa",
        "mode_a": "#00897b",
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
    order = [p for p in ["progress_only", "maat_only", "mode_a", "random"] if p in reg.index]
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
    pivot = pivot[[p for p in ["progress_only", "maat_only", "mode_a", "random"] if p in pivot.columns]]
    pivot.plot(kind="bar", ax=ax, color=[colors[p] for p in pivot.columns])
    ax.axhline(0.0, color="#333333", lw=1)
    ax.set_ylabel("median regret vs VSIDS")
    ax.set_title("Family-level regret reveals where Mode A helps or hurts")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig3_family_regret.png", dpi=180)
    plt.close(fig)


def write_readme(summary: pd.DataFrame, regret: pd.DataFrame) -> None:
    mode_a = regret[regret["policy"] == "mode_a"]
    mode_line = "not available"
    if not mode_a.empty:
        r = mode_a.iloc[0]
        mode_line = (
            f"median regret {r['median_regret_vs_vsids']:.4g}, "
            f"mean cost ratio {r['mean_cost_ratio_vs_vsids']:.4g}, "
            f"wins/losses/ties {int(r['wins_vs_vsids'])}/"
            f"{int(r['losses_vs_vsids'])}/{int(r['ties_vs_vsids'])}"
        )
    text = f"""# MAAT v1.6-T T1: Mode-A CDCL Branching Pilot

This experiment inserts the MAAT v1.6-T Mode-A acquisition rule into a small
transparent CDCL solver and measures compute-to-solution regret against VSIDS.

Mode A is:

```text
A(c) = prog(c) - tau * h_MAAT(c)
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
test in which VSIDS and Mode A share the same propagation and conflict-learning
machinery.

Headline Mode-A result vs VSIDS: {mode_line}.

Outputs:

- `t1_cdcl_runs.csv`: all per-instance runs
- `t1_cdcl_summary_by_policy.csv`: aggregate policy statistics
- `t1_cdcl_regret_vs_vsids.csv`: paired regret against VSIDS
- `t1_cdcl_family_regret_vs_vsids.csv`: family-level paired regret
- `t1_cdcl_summary.json`: metadata and compact result table
- `fig1_policy_compute_cost.png`
- `fig2_regret_vs_vsids.png`
- `fig3_family_regret.png`

Interpretation rule:

- Negative regret means the policy used less compute than VSIDS on the paired
  instance.
- Positive regret means VSIDS was better.
- A useful Mode-A signal must beat `progress_only`; otherwise the structural
  penalty did not earn its compute.
"""
    (OUTDIR / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df, summary, regret, fam_regret = run_benchmark()
    df.to_csv(OUTDIR / "t1_cdcl_runs.csv", index=False)
    summary.to_csv(OUTDIR / "t1_cdcl_summary_by_policy.csv", index=False)
    regret.to_csv(OUTDIR / "t1_cdcl_regret_vs_vsids.csv", index=False)
    fam_regret.to_csv(OUTDIR / "t1_cdcl_family_regret_vs_vsids.csv", index=False)
    plot_outputs(df, summary, regret, fam_regret)
    write_readme(summary, regret)

    payload = {
        "experiment": "MAAT v1.6-T T1 Mode-A CDCL branching pilot",
        "seed": SEED,
        "conflict_budget": CONFLICT_BUDGET,
        "beta": BETA,
        "tau": TAU,
        "horizon_h": HORIZON_H,
        "n_instances": int(df[["family", "instance_id"]].drop_duplicates().shape[0]),
        "policies": POLICIES,
        "summary_by_policy": summary.to_dict(orient="records"),
        "paired_regret_vs_vsids": regret.to_dict(orient="records"),
    }
    (OUTDIR / "t1_cdcl_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("T1 CDCL Mode-A benchmark complete.")
    print(summary.to_string(index=False))
    print()
    print(regret.to_string(index=False))
    print(f"\nOutputs written to: {OUTDIR}")


if __name__ == "__main__":
    main()
