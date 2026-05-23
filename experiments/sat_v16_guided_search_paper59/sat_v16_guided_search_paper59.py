#!/usr/bin/env python3
"""
Paper 59: MAAT v1.6 Guided SAT Search.

This benchmark tests whether structural-search coordinates can be used as an
active branching heuristic in a transparent DPLL solver with unit propagation.
It is not a proof of P vs NP, not a CDCL competition, and not a replacement for
modern solvers. The purpose is to move from passive SAT-hardness prediction to
active search steering under controlled, reproducible conditions.
"""

from __future__ import annotations

import json
import math
import os
import random
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable

_CACHE_ROOT = Path(os.environ.get("TMPDIR", "/tmp")) / "maat_paper59_cache"
(_CACHE_ROOT / "matplotlib").mkdir(parents=True, exist_ok=True)
(_CACHE_ROOT / "xdg").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE_ROOT / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE_ROOT / "xdg"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SEED = 60060
OUTDIR = Path(__file__).resolve().parent / "outputs_paper59"
DECISION_BUDGET = 5500
EPS = 1.0e-9


@dataclass
class SolveStats:
    family: str
    instance_id: int
    heuristic: str
    n_vars: int
    n_clauses: int
    alpha: float
    solved: int
    sat: int
    timeout: int
    decisions: int
    propagations: int
    conflicts: int
    backtracks: int
    max_depth: int
    search_cost: float


def canonical_clause(clause: Iterable[int]) -> tuple[int, ...] | None:
    lits = tuple(sorted(set(int(x) for x in clause), key=lambda x: (abs(x), x)))
    seen: dict[int, int] = {}
    for lit in lits:
        v = abs(lit)
        if v in seen and seen[v] != lit:
            return None
        seen[v] = lit
    return lits if lits else None


def generate_random_3sat(rng: np.random.Generator, n_vars: int, n_clauses: int) -> list[tuple[int, ...]]:
    clauses: set[tuple[int, ...]] = set()
    attempts = 0
    while len(clauses) < n_clauses and attempts < 120 * n_clauses:
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
    while len(clauses) < n_clauses and attempts < 160 * n_clauses:
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
    while len(clauses) < n_clauses and attempts < 160 * n_clauses:
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
    out = []
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
    while len(used) < n_xors and attempts < 120 * n_xors:
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
        for alpha in [3.8, 4.25, 4.6]:
            m = int(round(alpha * n))
            for _ in range(4):
                instances.append(("random_3sat", n, generate_random_3sat(rng, n, m)))
                instances.append(("planted_3sat", n, generate_planted_3sat(rng, n, m)))
    for n in [18, 24]:
        for alpha in [3.9, 4.25, 4.7]:
            m = int(round(alpha * n))
            for _ in range(4):
                instances.append(("modular_3sat", n, generate_modular_3sat(rng, n, m)))
    for n in [15, 20]:
        for ratio in [0.8, 1.0, 1.2]:
            for _ in range(3):
                instances.append(("xor3_cnf", n, generate_xor3_cnf(rng, n, int(round(ratio * n)))))
    for vertices in [8, 10, 12]:
        for prob in [0.30, 0.45]:
            for colors in [3, 4]:
                for _ in range(2):
                    n_vars, clauses = generate_coloring_cnf(rng, vertices, prob, colors)
                    instances.append(("graph_coloring", n_vars, clauses))
    for holes in [4, 5, 6]:
        n_vars, clauses = generate_pigeonhole_cnf(holes)
        instances.append(("pigeonhole", n_vars, clauses))
    return instances


def lit_value(lit: int, assignment: dict[int, bool]) -> bool | None:
    val = assignment.get(abs(lit))
    if val is None:
        return None
    return val if lit > 0 else not val


def clause_state(clause: tuple[int, ...], assignment: dict[int, bool]) -> tuple[str, list[int]]:
    unassigned: list[int] = []
    for lit in clause:
        val = lit_value(lit, assignment)
        if val is True:
            return "sat", []
        if val is None:
            unassigned.append(lit)
    if not unassigned:
        return "conflict", []
    return "open", unassigned


def unit_propagate(
    clauses: list[tuple[int, ...]],
    assignment: dict[int, bool],
    stats: dict[str, int],
) -> bool:
    changed = True
    while changed:
        changed = False
        for clause in clauses:
            state, unassigned = clause_state(clause, assignment)
            if state == "sat":
                continue
            if state == "conflict":
                stats["conflicts"] += 1
                return False
            if len(unassigned) == 1:
                lit = unassigned[0]
                v = abs(lit)
                val = lit > 0
                existing = assignment.get(v)
                if existing is not None and existing != val:
                    stats["conflicts"] += 1
                    return False
                if existing is None:
                    assignment[v] = val
                    stats["propagations"] += 1
                    changed = True
    return True


def active_clauses(clauses: list[tuple[int, ...]], assignment: dict[int, bool]) -> list[list[int]]:
    out: list[list[int]] = []
    for clause in clauses:
        state, unassigned = clause_state(clause, assignment)
        if state == "sat":
            continue
        if state == "conflict":
            out.append([])
        else:
            out.append(unassigned)
    return out


def all_satisfied(clauses: list[tuple[int, ...]], assignment: dict[int, bool]) -> bool:
    return all(clause_state(c, assignment)[0] == "sat" for c in clauses)


def score_candidates(
    active: list[list[int]],
    unassigned_vars: list[int],
) -> dict[int, dict[str, float]]:
    pos = {v: 0.0 for v in unassigned_vars}
    neg = {v: 0.0 for v in unassigned_vars}
    jw_pos = {v: 0.0 for v in unassigned_vars}
    jw_neg = {v: 0.0 for v in unassigned_vars}
    moms_pos = {v: 0.0 for v in unassigned_vars}
    moms_neg = {v: 0.0 for v in unassigned_vars}
    short = {v: 0.0 for v in unassigned_vars}
    neigh: dict[int, set[int]] = {v: set() for v in unassigned_vars}

    nonempty = [c for c in active if c]
    min_len = min((len(c) for c in nonempty), default=1)
    for clause in nonempty:
        weight = 2.0 ** (-len(clause))
        vars_in = [abs(lit) for lit in clause]
        for lit in clause:
            v = abs(lit)
            if v not in pos:
                continue
            if lit > 0:
                pos[v] += 1.0
                jw_pos[v] += weight
                if len(clause) == min_len:
                    moms_pos[v] += 1.0
            else:
                neg[v] += 1.0
                jw_neg[v] += weight
                if len(clause) == min_len:
                    moms_neg[v] += 1.0
            if len(clause) <= 3:
                short[v] += 1.0 / len(clause)
            neigh[v].update(u for u in vars_in if u != v)

    max_degree = max((pos[v] + neg[v] for v in unassigned_vars), default=1.0) + EPS
    max_jw = max((jw_pos[v] + jw_neg[v] for v in unassigned_vars), default=1.0) + EPS
    max_short = max((short[v] for v in unassigned_vars), default=1.0) + EPS
    max_neigh = max((len(neigh[v]) for v in unassigned_vars), default=1) + EPS

    scores: dict[int, dict[str, float]] = {}
    for v in unassigned_vars:
        degree = pos[v] + neg[v]
        sign_total = degree + EPS
        imbalance = abs(pos[v] - neg[v]) / sign_total
        H = 1.0 / (1.0 + imbalance)
        B = 1.0 - imbalance
        S = short[v] / max_short
        V = len(neigh[v]) / max_neigh
        R_resp = (max(EPS, H) * max(EPS, B) * max(EPS, V)) ** (1.0 / 3.0)
        R_rob = min(R_resp, (max(EPS, H) * max(EPS, B) * max(EPS, S) * max(EPS, V)) ** 0.25)
        verification_gain = max(pos[v], neg[v]) / max_degree
        energy_pressure = (jw_pos[v] + jw_neg[v]) / max_jw
        moms_score = (moms_pos[v] + moms_neg[v]) + (moms_pos[v] * moms_neg[v])
        support_cost = -sum(math.log(EPS + q) for q in (H, B, S, V, R_rob))
        maat = (
            0.44 * verification_gain
            + 0.24 * energy_pressure
            + 0.14 * S
            + 0.13 * V
            + 0.08 * R_rob
            - 0.035 * support_cost
        )
        maat_heavy = (
            0.16 * verification_gain
            + 0.12 * energy_pressure
            + 0.24 * S
            + 0.28 * V
            + 0.24 * R_rob
            - 0.050 * support_cost
        )
        scores[v] = {
            "degree": degree,
            "jw": jw_pos[v] + jw_neg[v],
            "moms": moms_score,
            "H": H,
            "B": B,
            "S": S,
            "V": V,
            "R_rob": R_rob,
            "maat": maat,
            "maat_heavy": maat_heavy,
            "polarity_true_score": pos[v] + jw_pos[v],
            "polarity_false_score": neg[v] + jw_neg[v],
        }
    return scores


def choose_variable(
    clauses: list[tuple[int, ...]],
    n_vars: int,
    assignment: dict[int, bool],
    heuristic: str,
    rng: random.Random,
) -> tuple[int, bool]:
    unassigned = [v for v in range(1, n_vars + 1) if v not in assignment]
    active = active_clauses(clauses, assignment)
    scores = score_candidates(active, unassigned)
    if heuristic == "random":
        v = rng.choice(unassigned)
    elif heuristic == "first":
        v = min(unassigned)
    elif heuristic == "degree":
        v = max(unassigned, key=lambda x: (scores[x]["degree"], -x))
    elif heuristic == "jeroslow_wang":
        v = max(unassigned, key=lambda x: (scores[x]["jw"], scores[x]["degree"], -x))
    elif heuristic == "moms":
        v = max(unassigned, key=lambda x: (scores[x]["moms"], scores[x]["degree"], -x))
    elif heuristic == "maat_v16_branch":
        v = max(unassigned, key=lambda x: (scores[x]["maat"], scores[x]["jw"], -x))
    elif heuristic == "maat_v16_heavy":
        v = max(unassigned, key=lambda x: (scores[x]["maat_heavy"], scores[x]["jw"], -x))
    else:
        raise ValueError(f"unknown heuristic: {heuristic}")
    polarity = scores[v]["polarity_true_score"] >= scores[v]["polarity_false_score"]
    return v, polarity


def dpll_solve(
    clauses: list[tuple[int, ...]],
    n_vars: int,
    heuristic: str,
    seed: int,
    decision_budget: int = DECISION_BUDGET,
) -> tuple[bool, bool, dict[str, int]]:
    rng = random.Random(seed)
    stats = {
        "decisions": 0,
        "propagations": 0,
        "conflicts": 0,
        "backtracks": 0,
        "max_depth": 0,
        "timeout": 0,
    }

    def rec(assignment: dict[int, bool], depth: int) -> bool:
        stats["max_depth"] = max(stats["max_depth"], depth)
        if stats["decisions"] >= decision_budget:
            stats["timeout"] = 1
            return False
        local = dict(assignment)
        if not unit_propagate(clauses, local, stats):
            return False
        if all_satisfied(clauses, local):
            assignment.clear()
            assignment.update(local)
            return True
        if len(local) == n_vars:
            return False
        var, polarity = choose_variable(clauses, n_vars, local, heuristic, rng)
        for val in (polarity, not polarity):
            if stats["decisions"] >= decision_budget:
                stats["timeout"] = 1
                return False
            stats["decisions"] += 1
            nxt = dict(local)
            nxt[var] = val
            if rec(nxt, depth + 1):
                assignment.clear()
                assignment.update(nxt)
                return True
            stats["backtracks"] += 1
        return False

    assignment: dict[int, bool] = {}
    sat = rec(assignment, 0)
    solved = stats["timeout"] == 0
    return sat, solved, stats


def run_benchmark() -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(SEED)
    instances = build_instances(rng)
    heuristics = [
        "random",
        "first",
        "degree",
        "jeroslow_wang",
        "moms",
        "maat_v16_branch",
        "maat_v16_heavy",
    ]
    rows: list[SolveStats] = []

    for instance_id, (family, n_vars, clauses) in enumerate(instances):
        n_clauses = len(clauses)
        alpha = n_clauses / max(1, n_vars)
        for heuristic in heuristics:
            sat, solved, stats = dpll_solve(clauses, n_vars, heuristic, SEED + 997 * instance_id)
            timeout = int(stats["timeout"])
            cost = math.log1p(
                stats["decisions"]
                + 0.35 * stats["conflicts"]
                + 0.10 * stats["propagations"]
                + 0.20 * stats["backtracks"]
            )
            rows.append(
                SolveStats(
                    family=family,
                    instance_id=instance_id,
                    heuristic=heuristic,
                    n_vars=n_vars,
                    n_clauses=n_clauses,
                    alpha=alpha,
                    solved=int(solved),
                    sat=int(sat),
                    timeout=timeout,
                    decisions=int(stats["decisions"]),
                    propagations=int(stats["propagations"]),
                    conflicts=int(stats["conflicts"]),
                    backtracks=int(stats["backtracks"]),
                    max_depth=int(stats["max_depth"]),
                    search_cost=float(cost),
                )
            )

    df = pd.DataFrame([asdict(r) for r in rows])
    summary = (
        df.groupby(["family", "heuristic"], as_index=False)
        .agg(
            solved_rate=("solved", "mean"),
            timeout_rate=("timeout", "mean"),
            median_decisions=("decisions", "median"),
            median_conflicts=("conflicts", "median"),
            median_cost=("search_cost", "median"),
            mean_cost=("search_cost", "mean"),
            n_runs=("solved", "size"),
        )
        .sort_values(["family", "median_cost"])
    )
    return df, summary


def plot_outputs(df: pd.DataFrame, summary: pd.DataFrame) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")
    order = [
        "first",
        "random",
        "degree",
        "jeroslow_wang",
        "moms",
        "maat_v16_branch",
        "maat_v16_heavy",
    ]
    colors = {
        "first": "#9e9e9e",
        "random": "#616161",
        "degree": "#2b6cb0",
        "jeroslow_wang": "#dd6b20",
        "moms": "#805ad5",
        "maat_v16_branch": "#00897b",
        "maat_v16_heavy": "#004d40",
    }

    fig, ax = plt.subplots(figsize=(11, 5.8))
    pivot = summary.pivot(index="family", columns="heuristic", values="median_cost").reindex(columns=order)
    pivot.plot(kind="bar", ax=ax, width=0.84, color=[colors[c] for c in pivot.columns])
    ax.set_ylabel("median log search cost")
    ax.set_xlabel("")
    ax.set_title("Paper 59: active SAT search cost by family")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig1_median_search_cost.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 5.8))
    pivot = summary.pivot(index="family", columns="heuristic", values="timeout_rate").reindex(columns=order)
    pivot.plot(kind="bar", ax=ax, width=0.84, color=[colors[c] for c in pivot.columns])
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("timeout rate")
    ax.set_xlabel("")
    ax.set_title("Timeouts under fixed decision budget")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig2_timeout_rate.png", dpi=180)
    plt.close(fig)

    agg = (
        df.groupby("heuristic", as_index=False)
        .agg(median_cost=("search_cost", "median"), timeout_rate=("timeout", "mean"))
        .set_index("heuristic")
        .reindex(order)
        .reset_index()
    )
    fig, ax1 = plt.subplots(figsize=(9, 5.2))
    ax1.bar(agg["heuristic"], agg["median_cost"], color=[colors[h] for h in agg["heuristic"]])
    ax1.set_ylabel("aggregate median log cost")
    ax1.set_title("Aggregate heuristic comparison")
    ax1.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig3_aggregate_cost.png", dpi=180)
    plt.close(fig)

    pair = df.pivot_table(index=["family", "instance_id"], columns="heuristic", values="search_cost")
    if {"maat_v16_branch", "jeroslow_wang"}.issubset(pair.columns):
        diff = pair["maat_v16_branch"] - pair["jeroslow_wang"]
        fig, ax = plt.subplots(figsize=(8.5, 5.2))
        ax.axvline(0, color="black", lw=1)
        ax.hist(diff.dropna(), bins=28, color="#00897b", alpha=0.8)
        ax.set_xlabel("MAAT v1.6 cost - Jeroslow-Wang cost")
        ax.set_ylabel("instances")
        ax.set_title("Paired cost difference against a strong classical heuristic")
        fig.tight_layout()
        fig.savefig(OUTDIR / "fig4_pairwise_difference.png", dpi=180)
        plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df, summary = run_benchmark()
    plot_outputs(df, summary)
    df.to_csv(OUTDIR / "paper59_solver_runs.csv", index=False)
    summary.to_csv(OUTDIR / "paper59_summary_by_family.csv", index=False)
    aggregate = (
        df.groupby("heuristic", as_index=False)
        .agg(
            solved_rate=("solved", "mean"),
            timeout_rate=("timeout", "mean"),
            median_decisions=("decisions", "median"),
            median_conflicts=("conflicts", "median"),
            median_cost=("search_cost", "median"),
            mean_cost=("search_cost", "mean"),
            n_runs=("solved", "size"),
        )
        .sort_values("median_cost")
    )
    aggregate.to_csv(OUTDIR / "paper59_summary_aggregate.csv", index=False)

    pair = df.pivot_table(index=["family", "instance_id"], columns="heuristic", values="search_cost")
    paired = {}
    for base in ["degree", "jeroslow_wang", "moms"]:
        if {"maat_v16_branch", base}.issubset(pair.columns):
            d = (pair["maat_v16_branch"] - pair[base]).dropna()
            paired[f"maat_minus_{base}_median"] = float(np.median(d))
            paired[f"maat_better_than_{base}_fraction"] = float(np.mean(d < 0))

    payload = {
        "paper": "Paper 59 - MAAT v1.6 Guided SAT Search",
        "status": (
            "Transparent DPLL + unit propagation benchmark. Not a P vs NP claim, "
            "not a CDCL competition, and not a modern SAT-solver replacement."
        ),
        "decision_budget": DECISION_BUDGET,
        "n_instances": int(df[["family", "instance_id"]].drop_duplicates().shape[0]),
        "heuristics": sorted(df["heuristic"].unique().tolist()),
        "families": sorted(df["family"].unique().tolist()),
        "aggregate": aggregate.to_dict(orient="records"),
        "paired_cost_differences": paired,
    }
    with open(OUTDIR / "paper59_summary.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
