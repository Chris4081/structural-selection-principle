"""Pure Paper-75 state semantics used only by preregistration unit tests.

This module is not a Paper-76 solver and does not load external CNFs.  It
defines the frozen residual-CNF state, support equations, deterministic split,
and rollback semantics on small explicitly synthetic fixtures.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
from itertools import combinations
from math import exp
from statistics import fmean, pstdev


EPS = 1e-9
SPLIT_SEED = 75075


def deterministic_split(dataset_id: str, instance_sha256: str) -> str:
    """Return the frozen family-local ordering key, not a global split.

    The manifest freezer sorts this key within each family and assigns the
    first max(1, floor(0.2*n_family)) rows to calibration.
    """

    payload = f"paper75|{SPLIT_SEED}|{dataset_id}|{instance_sha256}"
    return sha256(payload.encode("utf-8")).hexdigest()


def hash_to_unit(
    dataset_id: str,
    decision_index: int,
    state_fingerprint: str,
    namespace: str = "primary",
) -> float:
    """Frozen state-bound random number shared by primary stochastic policies."""

    payload = (
        f"paper75|75075|{namespace}|{dataset_id}|{decision_index}|"
        f"{state_fingerprint}"
    )
    integer = int(sha256(payload.encode("utf-8")).hexdigest(), 16)
    return integer / float(2**256)


def classify_safety(ci_low: float, ci_high: float, delta: float = 0.05) -> str:
    """Classify the independent no-harm axis without conflating uncertainty."""

    if ci_low > delta:
        return "harmful"
    if ci_high <= delta:
        return "no_harm_certified"
    return "no_harm_not_established"


def sigmoid(value: float) -> float:
    if value >= 0:
        z = exp(-value)
        return 1.0 / (1.0 + z)
    z = exp(value)
    return z / (1.0 + z)


@dataclass
class ResidualState:
    """Frozen Paper-75 decision state for semantics and rollback tests.

    Original and learned clauses are persistent. Assignments at levels above a
    backtrack target are removed; learned clauses are not rolled back.
    """

    n_vars: int
    original_clauses: list[tuple[int, ...]]
    learned_clauses: list[tuple[int, ...]] = field(default_factory=list)
    assignment: dict[int, bool] = field(default_factory=dict)
    level_of: dict[int, int] = field(default_factory=dict)
    reason: dict[int, tuple[int, ...] | None] = field(default_factory=dict)
    trail: list[int] = field(default_factory=list)
    decision_level: int = 0
    vsids: dict[int, float] = field(default_factory=dict)
    phase: dict[int, bool] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.vsids:
            self.vsids = {v: 1.0 for v in range(1, self.n_vars + 1)}

    @property
    def clauses(self) -> list[tuple[int, ...]]:
        return self.original_clauses + self.learned_clauses

    def enqueue(self, lit: int, level: int, reason: tuple[int, ...] | None = None) -> None:
        var = abs(lit)
        value = lit > 0
        old = self.assignment.get(var)
        if old is not None and old != value:
            raise ValueError("contradictory synthetic assignment")
        if old is not None:
            return
        self.assignment[var] = value
        self.level_of[var] = level
        self.reason[var] = reason
        self.phase[var] = value
        self.trail.append(lit)
        self.decision_level = max(self.decision_level, level)

    def add_learned_clause(self, clause: tuple[int, ...]) -> None:
        if not clause:
            raise ValueError("empty learned clause is terminal, not a decision state")
        self.learned_clauses.append(tuple(clause))

    def backtrack(self, target_level: int) -> None:
        keep: list[int] = []
        for lit in self.trail:
            var = abs(lit)
            if self.level_of[var] <= target_level:
                keep.append(lit)
            else:
                self.assignment.pop(var, None)
                self.level_of.pop(var, None)
                self.reason.pop(var, None)
        self.trail = keep
        self.decision_level = target_level

    def residual_clauses(self) -> list[tuple[int, ...]]:
        residual: list[tuple[int, ...]] = []
        for clause in self.clauses:
            open_lits: list[int] = []
            satisfied = False
            for lit in clause:
                value = self.assignment.get(abs(lit))
                if value is None:
                    open_lits.append(lit)
                elif value == (lit > 0):
                    satisfied = True
                    break
            if not satisfied:
                residual.append(tuple(open_lits))
        return residual

    def fingerprint(self) -> str:
        residual = sorted(tuple(sorted(c)) for c in self.residual_clauses())
        payload = repr((self.decision_level, tuple(self.trail), residual))
        return sha256(payload.encode("utf-8")).hexdigest()


def state_supports(state: ResidualState) -> dict[str, float]:
    """Compute the exact Paper-75 H, B, S, V and robustness closure.

    Evaluation is valid only at a nonterminal decision boundary: no residual
    empty clause and at least one active variable.
    """

    clauses = state.residual_clauses()
    if any(len(c) == 0 for c in clauses):
        raise ValueError("conflicting state is not an evaluation boundary")
    active = [c for c in clauses if c]
    variables = sorted({abs(lit) for clause in active for lit in clause})
    if not variables:
        raise ValueError("terminal satisfied state is not an evaluation boundary")

    pos = {v: 0.0 for v in variables}
    neg = {v: 0.0 for v in variables}
    short = {v: 0.0 for v in variables}
    neighbors = {v: set() for v in variables}
    for clause in active:
        vars_in = sorted({abs(lit) for lit in clause})
        for lit in clause:
            v = abs(lit)
            if lit > 0:
                pos[v] += 1.0
            else:
                neg[v] += 1.0
            if len(clause) <= 3:
                short[v] += 1.0 / len(clause)
        for a, b in combinations(vars_in, 2):
            neighbors[a].add(b)
            neighbors[b].add(a)

    degree = [pos[v] + neg[v] for v in variables]
    imbalance = fmean(abs(pos[v] - neg[v]) / (pos[v] + neg[v] + EPS) for v in variables)
    degree_cv = pstdev(degree) / (fmean(degree) + EPS)
    h_support = 1.0 / (1.0 + imbalance)
    b_support = 1.0 / (1.0 + degree_cv)

    max_short = max(short.values(), default=0.0)
    s_support = fmean(short.values()) / (max_short + EPS) if max_short > 0 else 0.0

    seen: set[int] = set()
    component_sizes: list[int] = []
    clustering: list[float] = []
    for start in variables:
        if start not in seen:
            stack = [start]
            seen.add(start)
            size = 0
            while stack:
                v = stack.pop()
                size += 1
                for u in neighbors[v]:
                    if u not in seen:
                        seen.add(u)
                        stack.append(u)
            component_sizes.append(size)
        nb = sorted(neighbors[start])
        if len(nb) < 2:
            clustering.append(0.0)
        else:
            possible = len(nb) * (len(nb) - 1) / 2.0
            actual = sum(1 for a, b in combinations(nb, 2) if b in neighbors[a])
            clustering.append(actual / possible)
    lcc = max(component_sizes) / len(variables)
    v_support = 0.65 * lcc + 0.35 * fmean(clustering)

    values = [max(EPS, min(1.0, x)) for x in (h_support, b_support, s_support, v_support)]
    h_support, b_support, s_support, v_support = values
    r_resp = (h_support * b_support * v_support) ** (1.0 / 3.0)
    r_struct = (h_support * b_support * s_support * v_support) ** 0.25
    r_rob = min(r_resp, r_struct)
    return {
        "H": h_support,
        "B": b_support,
        "S": s_support,
        "V": v_support,
        "R_resp": r_resp,
        "R_struct": r_struct,
        "R_rob": r_rob,
        "literal_imbalance": imbalance,
        "degree_cv": degree_cv,
        "n_active_variables": float(len(variables)),
        "n_active_clauses": float(len(active)),
    }
