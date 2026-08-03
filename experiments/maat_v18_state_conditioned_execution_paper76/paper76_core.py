#!/usr/bin/env python3
"""Frozen state-conditioned SAT primitives for Paper 76.

This module implements the code-level state semantics preregistered in Paper
75.  It deliberately contains no fitting, threshold search, or result
interpretation.  External execution is orchestrated by ``paper76_execution``.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations
from math import exp, log
from pathlib import Path
from statistics import fmean, pstdev
from time import perf_counter_ns
from typing import Iterable


EPS = 1.0e-9
SEED = 75075


def clip01(value: float) -> float:
    return max(0.0, min(1.0, float(value)))


def sigmoid(value: float) -> float:
    if value >= 0:
        z = exp(-value)
        return 1.0 / (1.0 + z)
    z = exp(value)
    return z / (1.0 + z)


def hash_to_unit(
    dataset_id: str,
    decision_index: int,
    state_fingerprint: str,
    namespace: str = "primary",
) -> float:
    payload = (
        f"paper75|{SEED}|{namespace}|{dataset_id}|{decision_index}|"
        f"{state_fingerprint}"
    )
    return int(sha256(payload.encode("utf-8")).hexdigest(), 16) / float(2**256)


def random_activation(dataset_id: str, decision_index: int) -> bool:
    payload = f"{SEED}|{dataset_id}|{decision_index}|random_null"
    integer = int(sha256(payload.encode("utf-8")).hexdigest(), 16)
    return integer < int(0.25 * (2**256 - 1))


def state_priority(dataset_id: str, decision_index: int, fingerprint: str) -> str:
    payload = f"paper75|{dataset_id}|{decision_index}|{fingerprint}"
    return sha256(payload.encode("utf-8")).hexdigest()


def weighted_quantile_lower(values: list[float], weights: list[float], q: float) -> float:
    if not values or len(values) != len(weights):
        raise ValueError("weighted quantile requires aligned nonempty inputs")
    ordered = sorted(zip(values, weights), key=lambda item: item[0])
    target = q * sum(weight for _, weight in ordered)
    cumulative = 0.0
    for value, weight in ordered:
        cumulative += weight
        if cumulative >= target:
            return float(value)
    return float(ordered[-1][0])


def weighted_median(values: list[float], weights: list[float]) -> float:
    return weighted_quantile_lower(values, weights, 0.5)


def weighted_mad(values: list[float], weights: list[float]) -> float:
    center = weighted_median(values, weights)
    return weighted_median([abs(value - center) for value in values], weights)


def weighted_mean(values: Iterable[float], weights: Iterable[float]) -> float:
    pairs = list(zip(values, weights))
    total = sum(weight for _, weight in pairs)
    if total <= 0:
        return float("nan")
    return sum(value * weight for value, weight in pairs) / total


def weighted_std(values: list[float], weights: list[float]) -> float:
    mean = weighted_mean(values, weights)
    total = sum(weights)
    if total <= 0:
        return float("nan")
    return (sum(w * (v - mean) ** 2 for v, w in zip(values, weights)) / total) ** 0.5


def rank_average(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = 0.5 * (start + 1 + end)
        for pos in range(start, end):
            ranks[order[pos]] = rank
        start = end
    return ranks


def weighted_corr(x: list[float], y: list[float], weights: list[float]) -> float:
    mx = weighted_mean(x, weights)
    my = weighted_mean(y, weights)
    cov = sum(w * (a - mx) * (b - my) for a, b, w in zip(x, y, weights))
    vx = sum(w * (a - mx) ** 2 for a, w in zip(x, weights))
    vy = sum(w * (b - my) ** 2 for b, w in zip(y, weights))
    if vx <= 0 or vy <= 0:
        return float("nan")
    return cov / (vx * vy) ** 0.5


def weighted_spearman(x: list[float], y: list[float], weights: list[float]) -> float:
    return weighted_corr(rank_average(x), rank_average(y), weights)


def safe_entropy(values: list[float]) -> float:
    if len(values) <= 1:
        return 0.0
    nonnegative = [max(0.0, value) for value in values]
    total = sum(nonnegative)
    if total <= 0:
        return 0.0
    probs = [value / total for value in nonnegative if value > 0]
    return -sum(p * log(p) for p in probs) / log(len(values))


def vector_aggregates(values: list[float]) -> list[float]:
    if not values:
        return [0.0] * 6
    ordered = sorted(values, reverse=True)
    mean = fmean(values)
    std = pstdev(values)
    gap = 0.0 if len(ordered) < 2 else (ordered[0] - ordered[1]) / (abs(ordered[0]) + EPS)
    return [
        max(values),
        mean,
        std,
        std / (abs(mean) + EPS),
        gap,
        safe_entropy(values),
    ]


def closure(supports: dict[str, float], variant: str = "HBSV") -> tuple[float, float, float]:
    h, b, s, v = (max(EPS, supports[key]) for key in ("H", "B", "S", "V"))
    r_resp = (h * b * v) ** (1.0 / 3.0)
    if variant == "HBSV":
        r_struct = (h * b * s * v) ** 0.25
    else:
        selected = [max(EPS, supports[key]) for key in variant]
        r_struct = exp(sum(log(value) for value in selected) / len(selected))
    return r_resp, r_struct, min(r_resp, r_struct)


@dataclass
class SupportSnapshot:
    H: float
    B: float
    S: float
    V: float
    R_resp: float
    R_struct: float
    R_rob: float
    predictors: tuple[float, ...]

    def as_dict(self) -> dict[str, float]:
        return {
            "H": self.H,
            "B": self.B,
            "S": self.S,
            "V": self.V,
            "R_resp": self.R_resp,
            "R_struct": self.R_struct,
            "R_rob": self.R_rob,
        }


class IncrementalResidualTracker:
    """Incremental residual-CNF counters with persistent learned clauses."""

    def __init__(self, clauses: list[tuple[int, ...]], n_vars: int) -> None:
        self.n_vars = n_vars
        self.clauses: list[tuple[int, ...]] = []
        self.var_to_clauses: dict[int, set[int]] = {v: set() for v in range(1, n_vars + 1)}
        self.assignment: dict[int, bool] = {}
        self.residual: dict[int, tuple[int, ...] | None] = {}
        self.pos: Counter[int] = Counter()
        self.neg: Counter[int] = Counter()
        self.short: Counter[int] = Counter()
        self.edges: Counter[tuple[int, int]] = Counter()
        self.update_ns = 0
        for clause in clauses:
            self.add_clause(clause, count_time=False)

    def _residual_clause(self, clause: tuple[int, ...]) -> tuple[int, ...] | None:
        open_lits: list[int] = []
        for lit in clause:
            value = self.assignment.get(abs(lit))
            if value is None:
                open_lits.append(lit)
            elif value == (lit > 0):
                return None
        return tuple(open_lits)

    def _apply_contribution(self, residual: tuple[int, ...] | None, sign: int) -> None:
        if residual is None or not residual:
            return
        variables = sorted({abs(lit) for lit in residual})
        for lit in residual:
            target = self.pos if lit > 0 else self.neg
            target[abs(lit)] += sign
            if len(residual) <= 3:
                self.short[abs(lit)] += sign * (1.0 / len(residual))
        for a, b in combinations(variables, 2):
            self.edges[(a, b)] += sign

    def _refresh_clause(self, cid: int) -> None:
        old = self.residual.get(cid)
        self._apply_contribution(old, -1)
        new = self._residual_clause(self.clauses[cid])
        self.residual[cid] = new
        self._apply_contribution(new, +1)

    def add_clause(self, clause: tuple[int, ...], count_time: bool = True) -> None:
        start = perf_counter_ns()
        cid = len(self.clauses)
        canonical = tuple(clause)
        self.clauses.append(canonical)
        for var in {abs(lit) for lit in canonical}:
            self.var_to_clauses.setdefault(var, set()).add(cid)
        residual = self._residual_clause(canonical)
        self.residual[cid] = residual
        self._apply_contribution(residual, +1)
        if count_time:
            self.update_ns += perf_counter_ns() - start

    def assign(self, var: int, value: bool) -> None:
        start = perf_counter_ns()
        if var in self.assignment:
            if self.assignment[var] != value:
                raise ValueError("contradictory tracker assignment")
            return
        affected = tuple(self.var_to_clauses.get(var, ()))
        self.assignment[var] = value
        for cid in affected:
            self._refresh_clause(cid)
        self.update_ns += perf_counter_ns() - start

    def unassign(self, var: int) -> None:
        start = perf_counter_ns()
        if var not in self.assignment:
            return
        affected = tuple(self.var_to_clauses.get(var, ()))
        self.assignment.pop(var)
        for cid in affected:
            self._refresh_clause(cid)
        self.update_ns += perf_counter_ns() - start

    def residual_clauses(self) -> list[tuple[int, ...]]:
        return [res for cid in range(len(self.clauses)) if (res := self.residual[cid]) is not None]

    def has_conflict(self) -> bool:
        return any(res == () for res in self.residual.values())

    def _snapshot_from_counters(
        self,
        pos: Counter[int],
        neg: Counter[int],
        short: Counter[int],
        edges: Counter[tuple[int, int]],
        vsids: dict[int, float],
    ) -> SupportSnapshot:
        variables = sorted(v for v in range(1, self.n_vars + 1) if pos[v] + neg[v] > 0)
        if not variables:
            raise ValueError("terminal state is not an evaluation boundary")
        degrees = [float(pos[v] + neg[v]) for v in variables]
        imbalance = [abs(pos[v] - neg[v]) / (pos[v] + neg[v] + EPS) for v in variables]
        h = 1.0 / (1.0 + fmean(imbalance))
        b = 1.0 / (1.0 + pstdev(degrees) / (fmean(degrees) + EPS))
        short_values = [float(short[v]) for v in variables]
        max_short = max(short_values, default=0.0)
        s = fmean(short_values) / (max_short + EPS) if max_short > 0 else 0.0

        adjacency = {v: set() for v in variables}
        for (a, b_var), multiplicity in edges.items():
            if multiplicity > 0 and a in adjacency and b_var in adjacency:
                adjacency[a].add(b_var)
                adjacency[b_var].add(a)
        seen: set[int] = set()
        component_sizes: list[int] = []
        for start in variables:
            if start in seen:
                continue
            stack = [start]
            seen.add(start)
            size = 0
            while stack:
                node = stack.pop()
                size += 1
                for neighbor in adjacency[node]:
                    if neighbor not in seen:
                        seen.add(neighbor)
                        stack.append(neighbor)
            component_sizes.append(size)
        lcc = max(component_sizes, default=0) / len(variables)
        clustering: list[float] = []
        for var in variables:
            neighbors = sorted(adjacency[var])
            if len(neighbors) < 2:
                clustering.append(0.0)
                continue
            possible = len(neighbors) * (len(neighbors) - 1) / 2.0
            actual = sum(1 for a, b_var in combinations(neighbors, 2) if b_var in adjacency[a])
            clustering.append(actual / possible)
        v_support = clip01(0.65 * lcc + 0.35 * fmean(clustering))
        supports = {"H": clip01(h), "B": clip01(b), "S": clip01(s), "V": v_support}
        supports = {key: max(EPS, value) for key, value in supports.items()}
        r_resp, r_struct, r_rob = closure(supports)

        active_residual = [res for res in self.residual_clauses() if res]
        min_len = min((len(clause) for clause in active_residual), default=0)
        moms: dict[int, float] = {var: 0.0 for var in variables}
        jw: dict[int, float] = {var: 0.0 for var in variables}
        moms_pos: Counter[int] = Counter()
        moms_neg: Counter[int] = Counter()
        for clause in active_residual:
            for lit in clause:
                var = abs(lit)
                jw[var] += 2.0 ** (-len(clause))
                if len(clause) == min_len:
                    if lit > 0:
                        moms_pos[var] += 1
                    else:
                        moms_neg[var] += 1
        for var in variables:
            p, n = moms_pos[var], moms_neg[var]
            moms[var] = p + n + 2.0 * p * n
        predictors = tuple(
            vector_aggregates([moms[var] for var in variables])
            + vector_aggregates([jw[var] for var in variables])
            + vector_aggregates([float(vsids.get(var, 0.0)) for var in variables])
        )
        return SupportSnapshot(
            H=supports["H"], B=supports["B"], S=supports["S"], V=supports["V"],
            R_resp=r_resp, R_struct=r_struct, R_rob=r_rob, predictors=predictors,
        )

    def snapshot(self, vsids: dict[int, float]) -> SupportSnapshot:
        if self.has_conflict():
            raise ValueError("conflicting state is not an evaluation boundary")
        return self._snapshot_from_counters(self.pos, self.neg, self.short, self.edges, vsids)

    def full_snapshot(self, vsids: dict[int, float]) -> SupportSnapshot:
        pos: Counter[int] = Counter()
        neg: Counter[int] = Counter()
        short: Counter[int] = Counter()
        edges: Counter[tuple[int, int]] = Counter()
        for residual in self.residual_clauses():
            if not residual:
                if residual == ():
                    raise ValueError("conflicting state is not an evaluation boundary")
                continue
            variables = sorted({abs(lit) for lit in residual})
            for lit in residual:
                (pos if lit > 0 else neg)[abs(lit)] += 1
                if len(residual) <= 3:
                    short[abs(lit)] += 1.0 / len(residual)
            for a, b in combinations(variables, 2):
                edges[(a, b)] += 1
        return self._snapshot_from_counters(pos, neg, short, edges, vsids)

    def fingerprint(self, decision_level: int, trail: list[int]) -> str:
        residual = sorted(tuple(sorted(clause)) for clause in self.residual_clauses())
        payload = repr((decision_level, tuple(trail), residual))
        return sha256(payload.encode("utf-8")).hexdigest()


def max_snapshot_difference(left: SupportSnapshot, right: SupportSnapshot) -> float:
    keys = ("H", "B", "S", "V", "R_resp", "R_struct", "R_rob")
    return max(abs(getattr(left, key) - getattr(right, key)) for key in keys)


def softmax_choice(scores: dict[int, float], beta: float, uniform: float) -> int:
    if not scores:
        raise ValueError("cannot choose from empty score set")
    keys = sorted(scores)
    peak = max(scores.values())
    weights = [exp(beta * (scores[key] - peak)) for key in keys]
    total = sum(weights)
    if total <= 0:
        return max(keys, key=lambda key: (scores[key], -key))
    target = uniform * total
    cumulative = 0.0
    for key, weight in zip(keys, weights):
        cumulative += weight
        if target <= cumulative:
            return key
    return keys[-1]


def classify_safety(ci_low: float, ci_high: float, delta: float = 0.05) -> str:
    if ci_low > delta:
        return "harmful"
    if ci_high <= delta:
        return "no_harm_certified"
    return "no_harm_not_established"

