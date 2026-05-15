#!/usr/bin/env python3
"""
Paper 55: External SAT/CDCL Hardness Validation.

This benchmark addresses a limitation of Papers 50/50b: the target is no
longer a hand-written DPLL node count, but CDCL solver statistics obtained
from PySAT solvers on canonical SAT benchmark families.

Run:
    python3 experiments/sat_cdcl_external_paper55/sat_cdcl_external_paper55.py
"""

from __future__ import annotations

import json
import math
import os
import time
from collections import deque
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable

SEED = 55055
OUTDIR = Path("outputs_paper55")
EPS = 1.0e-9
CONFLICT_BUDGET = 20000

OUTDIR.mkdir(parents=True, exist_ok=True)
(OUTDIR / ".mplconfig").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUTDIR / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pysat.solvers import Glucose3, Minisat22
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold, LeaveOneGroupOut, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


@dataclass
class InstanceRecord:
    instance_id: int
    family: str
    n_vars: int
    n_clauses: int
    alpha: float
    sat_glucose: int
    sat_minisat: int
    timeout_glucose: int
    timeout_minisat: int
    conflicts_glucose: int
    decisions_glucose: int
    propagations_glucose: int
    conflicts_minisat: int
    decisions_minisat: int
    propagations_minisat: int
    cdcl_hardness: float
    solve_time_s: float
    degree_cv: float
    max_degree_norm: float
    literal_imbalance: float
    mean_var_graph_degree: float
    clustering_mean: float
    component_frac: float
    largest_component_frac: float
    expansion_defect: float
    backdoor_density: float
    propagation_depth_mean: float
    contradiction_rate_mean: float
    covariance_condition_log: float
    H: float
    B: float
    S: float
    V: float
    R_resp: float
    R_rob: float
    F_struct_scalar: float
    density_peak: float
    frust_mean: float
    frust_std: float
    frust_q75: float
    frust_q90: float
    frust_q95: float
    frust_q99: float
    frust_max: float
    frust_tail95_mean: float
    frust_tail_mass: float
    frust_gini: float
    hotspot_fraction: float
    hotspot_cluster_frac: float
    hotspot_edge_density: float
    R_tail: float
    R_var: float
    R_cluster: float
    R_sat: float
    F_mean: float
    F_tail: float
    F_cluster: float
    F_scale: float
    F_struct_multi: float


def canonical_clause(clause: Iterable[int]) -> tuple[int, ...] | None:
    lits = tuple(sorted(set(int(x) for x in clause), key=lambda x: (abs(x), x)))
    vars_seen = {}
    for lit in lits:
        v = abs(lit)
        if v in vars_seen and vars_seen[v] != lit:
            return None
        vars_seen[v] = lit
    return lits if lits else None


def generate_random_3sat(rng: np.random.Generator, n_vars: int, n_clauses: int) -> list[tuple[int, ...]]:
    clauses: set[tuple[int, ...]] = set()
    attempts = 0
    while len(clauses) < n_clauses and attempts < 100 * n_clauses:
        attempts += 1
        vars_ = rng.choice(np.arange(1, n_vars + 1), size=3, replace=False)
        lits = [int((-1 if rng.random() < 0.5 else 1) * v) for v in vars_]
        clause = canonical_clause(lits)
        if clause is not None and len(clause) == 3:
            clauses.add(clause)
    return list(clauses)


def generate_planted_3sat(rng: np.random.Generator, n_vars: int, n_clauses: int) -> list[tuple[int, ...]]:
    assignment = rng.random(n_vars) < 0.5
    clauses: set[tuple[int, ...]] = set()
    attempts = 0
    while len(clauses) < n_clauses and attempts < 140 * n_clauses:
        attempts += 1
        vars_ = rng.choice(np.arange(1, n_vars + 1), size=3, replace=False)
        lits = [int((-1 if rng.random() < 0.5 else 1) * v) for v in vars_]
        if not any((lit > 0 and assignment[abs(lit) - 1]) or (lit < 0 and not assignment[abs(lit) - 1]) for lit in lits):
            j = int(rng.integers(0, 3))
            lits[j] = -lits[j]
        clause = canonical_clause(lits)
        if clause is not None and len(clause) == 3:
            clauses.add(clause)
    return list(clauses)


def generate_modular_3sat(rng: np.random.Generator, n_vars: int, n_clauses: int) -> list[tuple[int, ...]]:
    n_blocks = max(3, min(6, n_vars // 7))
    blocks = np.array_split(np.arange(1, n_vars + 1), n_blocks)
    clauses: set[tuple[int, ...]] = set()
    attempts = 0
    while len(clauses) < n_clauses and attempts < 140 * n_clauses:
        attempts += 1
        if rng.random() < 0.78:
            block = blocks[int(rng.integers(0, len(blocks)))]
            pool = block if len(block) >= 3 else np.arange(1, n_vars + 1)
            vars_ = rng.choice(pool, size=3, replace=False)
        else:
            vars_ = rng.choice(np.arange(1, n_vars + 1), size=3, replace=False)
        lits = [int((-1 if rng.random() < 0.5 else 1) * v) for v in vars_]
        clause = canonical_clause(lits)
        if clause is not None and len(clause) == 3:
            clauses.add(clause)
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
    while len(used) < n_xors and attempts < 100 * n_xors:
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
    for n in [24, 32, 40]:
        for alpha in [3.6, 4.0, 4.25, 4.5]:
            m = int(round(alpha * n))
            for _ in range(8):
                instances.append(("random_3sat", n, generate_random_3sat(rng, n, m)))
                instances.append(("planted_3sat", n, generate_planted_3sat(rng, n, m)))
    for n in [24, 32, 40]:
        for alpha in [3.8, 4.25, 4.7]:
            m = int(round(alpha * n))
            for _ in range(8):
                instances.append(("modular_3sat", n, generate_modular_3sat(rng, n, m)))
    for n in [18, 24, 30]:
        for ratio in [0.8, 1.0, 1.2]:
            for _ in range(6):
                instances.append(("xor3_cnf", n, generate_xor3_cnf(rng, n, int(round(ratio * n)))))
    for vertices in [10, 12, 14, 16]:
        for prob in [0.25, 0.35, 0.45]:
            for colors in [3, 4]:
                for _ in range(3):
                    n_vars, clauses = generate_coloring_cnf(rng, vertices, prob, colors)
                    instances.append(("graph_coloring", n_vars, clauses))
    for holes in [4, 5, 6, 7, 8]:
        n_vars, clauses = generate_pigeonhole_cnf(holes)
        instances.append(("pigeonhole", n_vars, clauses))
    return instances


def solve_cdcl(clauses: list[tuple[int, ...]], solver_cls) -> tuple[int, int, int, int, int]:
    with solver_cls() as solver:
        solver.append_formula([list(c) for c in clauses])
        solver.conf_budget(CONFLICT_BUDGET)
        result = solver.solve_limited(expect_interrupt=True)
        stats = solver.accum_stats()
    timeout = int(result is None)
    sat = -1 if result is None else int(bool(result))
    return (
        sat,
        timeout,
        int(stats.get("conflicts", 0)),
        int(stats.get("decisions", 0)),
        int(stats.get("propagations", 0)),
    )


def unit_propagate(clauses: list[tuple[int, ...]], assignment: dict[int, bool]) -> tuple[bool | None, list[tuple[int, ...]], dict[int, bool], int]:
    assignment = dict(assignment)
    forced = 0
    while True:
        simplified: list[tuple[int, ...]] = []
        for clause in clauses:
            satisfied = False
            remaining = []
            for lit in clause:
                var = abs(lit)
                if var in assignment:
                    val = assignment[var]
                    if (lit > 0 and val) or (lit < 0 and not val):
                        satisfied = True
                        break
                else:
                    remaining.append(lit)
            if satisfied:
                continue
            if not remaining:
                return False, [], assignment, forced
            simplified.append(tuple(remaining))
        if not simplified:
            return True, [], assignment, forced
        units = [c[0] for c in simplified if len(c) == 1]
        if not units:
            return None, simplified, assignment, forced
        changed = False
        for lit in units:
            var = abs(lit)
            val = lit > 0
            if var in assignment:
                if assignment[var] != val:
                    return False, [], assignment, forced
            else:
                assignment[var] = val
                forced += 1
                changed = True
        if not changed:
            return None, simplified, assignment, forced


def graph_arrays(clauses: list[tuple[int, ...]], n_vars: int) -> dict:
    deg = np.zeros(n_vars, dtype=float)
    pos = np.zeros(n_vars, dtype=float)
    neg = np.zeros(n_vars, dtype=float)
    var_to_clauses: list[set[int]] = [set() for _ in range(n_vars)]
    neighbors: list[set[int]] = [set() for _ in range(n_vars)]
    for ci, clause in enumerate(clauses):
        vars_ = sorted({abs(l) - 1 for l in clause})
        for lit in clause:
            idx = abs(lit) - 1
            if 0 <= idx < n_vars:
                deg[idx] += 1
                pos[idx] += lit > 0
                neg[idx] += lit < 0
                var_to_clauses[idx].add(ci)
        for a, b in combinations(vars_, 2):
            neighbors[a].add(b)
            neighbors[b].add(a)
    clustering = np.zeros(n_vars, dtype=float)
    for v in range(n_vars):
        nb = list(neighbors[v])
        if len(nb) < 2:
            continue
        possible = len(nb) * (len(nb) - 1) / 2
        actual = sum(1 for a, b in combinations(nb, 2) if b in neighbors[a])
        clustering[v] = actual / (possible + EPS)
    return {"deg": deg, "pos": pos, "neg": neg, "var_to_clauses": var_to_clauses, "neighbors": neighbors, "clustering": clustering}


def component_stats(neighbors: list[set[int]]) -> tuple[float, float]:
    n = len(neighbors)
    seen = np.zeros(n, dtype=bool)
    sizes = []
    for start in range(n):
        if seen[start]:
            continue
        q = deque([start])
        seen[start] = True
        size = 0
        while q:
            v = q.popleft()
            size += 1
            for u in neighbors[v]:
                if not seen[u]:
                    seen[u] = True
                    q.append(u)
        sizes.append(size)
    return len(sizes) / max(1, n), max(sizes) / max(1, n)


def largest_hotspot_cluster(neighbors: list[set[int]], hotspots: np.ndarray) -> tuple[float, float]:
    idx = set(np.where(hotspots)[0].tolist())
    if not idx:
        return 0.0, 0.0
    seen = set()
    best = 0
    edge_count = 0
    for v in idx:
        edge_count += sum(1 for u in neighbors[v] if u in idx)
        if v in seen:
            continue
        q = deque([v])
        seen.add(v)
        size = 0
        while q:
            cur = q.popleft()
            size += 1
            for u in neighbors[cur]:
                if u in idx and u not in seen:
                    seen.add(u)
                    q.append(u)
        best = max(best, size)
    edge_density = (edge_count / 2.0) / (len(idx) * (len(idx) - 1) / 2.0 + EPS)
    return best / len(neighbors), float(edge_density)


def gini(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    if len(x) == 0 or np.mean(x) <= EPS:
        return 0.0
    diffs = np.abs(x[:, None] - x[None, :])
    return float(np.mean(diffs) / (2.0 * np.mean(x) + EPS))


def extract_features(clauses: list[tuple[int, ...]], n_vars: int, rng: np.random.Generator) -> dict[str, float]:
    arrays = graph_arrays(clauses, n_vars)
    deg = arrays["deg"]
    pos = arrays["pos"]
    neg = arrays["neg"]
    neighbors = arrays["neighbors"]
    clustering = arrays["clustering"]
    var_to_clauses = arrays["var_to_clauses"]
    n_clauses = len(clauses)
    alpha = n_clauses / max(1, n_vars)

    degree_cv = float(np.std(deg) / (np.mean(deg) + EPS))
    max_degree_norm = float(np.max(deg) / (np.mean(deg) + EPS))
    literal_imbalance_vec = np.abs(pos - neg) / (deg + EPS)
    literal_imbalance = float(np.mean(literal_imbalance_vec))
    mean_var_graph_degree = float(np.mean([len(n) for n in neighbors]))
    clustering_mean = float(np.mean(clustering))
    component_frac, largest_component_frac = component_stats(neighbors)

    active = np.where(deg > 0)[0]
    expansion_vals = []
    for size in [3, 5, 8]:
        if len(active) < size:
            continue
        for _ in range(18):
            subset = rng.choice(active, size=size, replace=False)
            touched = set()
            for v in subset:
                touched.update(var_to_clauses[int(v)])
            expansion_vals.append(len(touched) / max(1.0, float(size * np.mean(np.maximum(deg[subset], 1)))))
    expansion_defect = float(max(0.0, 1.0 - (np.mean(expansion_vals) if expansion_vals else 1.0)))

    uncovered = set(range(n_clauses))
    chosen = 0
    while uncovered and (n_clauses - len(uncovered)) / max(1, n_clauses) < 0.80:
        best_v = max(range(n_vars), key=lambda v: len(var_to_clauses[v] & uncovered))
        gain = len(var_to_clauses[best_v] & uncovered)
        if gain == 0:
            break
        uncovered -= var_to_clauses[best_v]
        chosen += 1
    backdoor_density = float(chosen / max(1, n_vars))

    probe_vars = rng.choice(np.arange(1, n_vars + 1), size=min(n_vars, 24), replace=False)
    forced_vals = []
    contradiction_vals = []
    local_forced = np.zeros(n_vars, dtype=float)
    local_contra = np.zeros(n_vars, dtype=float)
    for var in probe_vars:
        vals_for_var = []
        cons_for_var = []
        for val in (False, True):
            status, _, _, forced = unit_propagate(clauses, {int(var): val})
            vals_for_var.append(forced / max(1, n_vars))
            cons_for_var.append(1.0 if status is False else 0.0)
            forced_vals.append(forced / max(1, n_vars))
            contradiction_vals.append(1.0 if status is False else 0.0)
        local_forced[int(var) - 1] = float(np.mean(vals_for_var))
        local_contra[int(var) - 1] = float(np.mean(cons_for_var))
    propagation_depth_mean = float(np.mean(forced_vals)) if forced_vals else 0.0
    contradiction_rate_mean = float(np.mean(contradiction_vals)) if contradiction_vals else 0.0

    local = np.column_stack(
        [
            deg,
            literal_imbalance_vec,
            np.array([len(n) for n in neighbors], dtype=float),
            clustering,
            local_forced,
            local_contra,
        ]
    )
    if local.shape[0] >= 4:
        cov = np.cov(local.T)
        cov += np.eye(cov.shape[0]) * 1e-6
        covariance_condition_log = float(np.log10(np.linalg.cond(cov)))
    else:
        covariance_condition_log = 0.0

    norm_deg = deg / (np.mean(deg) + EPS)
    frust = (
        0.30 * norm_deg
        + 0.20 * literal_imbalance_vec
        + 0.15 * clustering
        + 0.20 * local_forced
        + 0.30 * local_contra
        + 0.10 * np.array([len(n) for n in neighbors], dtype=float) / (np.mean([len(n) for n in neighbors]) + EPS)
    )
    frust_mean = float(np.mean(frust))
    frust_std = float(np.std(frust))
    frust_q75 = float(np.quantile(frust, 0.75))
    frust_q90 = float(np.quantile(frust, 0.90))
    frust_q95 = float(np.quantile(frust, 0.95))
    frust_q99 = float(np.quantile(frust, 0.99))
    frust_max = float(np.max(frust))
    tail = frust[frust >= frust_q95]
    frust_tail95_mean = float(np.mean(tail)) if len(tail) else frust_q95
    frust_tail_mass = float(np.sum(tail) / (np.sum(frust) + EPS))
    frust_gini = gini(frust)
    hotspots = frust >= frust_q90
    hotspot_fraction = float(np.mean(hotspots))
    hotspot_cluster_frac, hotspot_edge_density = largest_hotspot_cluster(neighbors, hotspots)

    H = 1.0 / (1.0 + 2.0 * contradiction_rate_mean + 0.35 * propagation_depth_mean + 0.06 * covariance_condition_log)
    B = 1.0 / (1.0 + degree_cv + literal_imbalance)
    S_phase = math.exp(-0.5 * ((alpha - 4.25) / 0.65) ** 2)
    S_local = math.exp(-0.5 * ((frust_q95 - frust_mean - 0.80) / 0.75) ** 2)
    S = float(math.sqrt(max(EPS, S_phase * S_local)))
    V = 1.0 / (1.0 + component_frac + max(0.0, 0.22 - clustering_mean) + 0.25 * (1.0 - largest_component_frac))
    R_resp = float((H * B * V) ** (1.0 / 3.0))
    R_rob = float(min(R_resp, (H * B * S * V) ** 0.25))
    F_struct_scalar = float(-sum(math.log(EPS + x) for x in [H, B, S, V, R_rob]))
    density_peak = float(math.exp(-0.5 * ((alpha - 4.25) / 0.55) ** 2))
    R_tail = 1.0 / (1.0 + frust_tail95_mean)
    R_var = 1.0 / (1.0 + frust_std + frust_gini)
    R_cluster = 1.0 / (1.0 + hotspot_cluster_frac + hotspot_edge_density)
    R_sat = float(min(R_rob, R_tail, R_var, R_cluster))
    F_mean = float(frust_mean)
    F_tail = float(frust_q95 + frust_tail_mass)
    F_cluster = float(hotspot_cluster_frac + hotspot_edge_density)
    F_scale = float(covariance_condition_log + expansion_defect + propagation_depth_mean)
    F_struct_multi = float(F_mean + F_tail + F_cluster + 0.25 * F_scale)

    return {
        "degree_cv": degree_cv,
        "max_degree_norm": max_degree_norm,
        "literal_imbalance": literal_imbalance,
        "mean_var_graph_degree": mean_var_graph_degree,
        "clustering_mean": clustering_mean,
        "component_frac": component_frac,
        "largest_component_frac": largest_component_frac,
        "expansion_defect": expansion_defect,
        "backdoor_density": backdoor_density,
        "propagation_depth_mean": propagation_depth_mean,
        "contradiction_rate_mean": contradiction_rate_mean,
        "covariance_condition_log": covariance_condition_log,
        "H": H,
        "B": B,
        "S": S,
        "V": V,
        "R_resp": R_resp,
        "R_rob": R_rob,
        "F_struct_scalar": F_struct_scalar,
        "density_peak": density_peak,
        "frust_mean": frust_mean,
        "frust_std": frust_std,
        "frust_q75": frust_q75,
        "frust_q90": frust_q90,
        "frust_q95": frust_q95,
        "frust_q99": frust_q99,
        "frust_max": frust_max,
        "frust_tail95_mean": frust_tail95_mean,
        "frust_tail_mass": frust_tail_mass,
        "frust_gini": frust_gini,
        "hotspot_fraction": hotspot_fraction,
        "hotspot_cluster_frac": hotspot_cluster_frac,
        "hotspot_edge_density": hotspot_edge_density,
        "R_tail": R_tail,
        "R_var": R_var,
        "R_cluster": R_cluster,
        "R_sat": R_sat,
        "F_mean": F_mean,
        "F_tail": F_tail,
        "F_cluster": F_cluster,
        "F_scale": F_scale,
        "F_struct_multi": F_struct_multi,
    }


def make_record(instance_id: int, family: str, n_vars: int, clauses: list[tuple[int, ...]], rng: np.random.Generator) -> InstanceRecord:
    t0 = time.perf_counter()
    g_sat, g_timeout, g_conf, g_dec, g_prop = solve_cdcl(clauses, Glucose3)
    m_sat, m_timeout, m_conf, m_dec, m_prop = solve_cdcl(clauses, Minisat22)
    elapsed = time.perf_counter() - t0
    mean_conf = 0.5 * (g_conf + m_conf)
    mean_dec = 0.5 * (g_dec + m_dec)
    mean_prop = 0.5 * (g_prop + m_prop)
    hardness = float(math.log1p(mean_conf + 0.10 * mean_dec + 0.01 * mean_prop))
    feats = extract_features(clauses, n_vars, rng)
    return InstanceRecord(
        instance_id=instance_id,
        family=family,
        n_vars=n_vars,
        n_clauses=len(clauses),
        alpha=len(clauses) / max(1, n_vars),
        sat_glucose=g_sat,
        sat_minisat=m_sat,
        timeout_glucose=g_timeout,
        timeout_minisat=m_timeout,
        conflicts_glucose=g_conf,
        decisions_glucose=g_dec,
        propagations_glucose=g_prop,
        conflicts_minisat=m_conf,
        decisions_minisat=m_dec,
        propagations_minisat=m_prop,
        cdcl_hardness=hardness,
        solve_time_s=elapsed,
        **feats,
    )


FEATURE_SETS = {
    "density_only": ["n_vars", "n_clauses", "alpha", "density_peak"],
    "standard_graph": [
        "n_vars",
        "n_clauses",
        "alpha",
        "degree_cv",
        "max_degree_norm",
        "literal_imbalance",
        "mean_var_graph_degree",
        "clustering_mean",
        "component_frac",
        "largest_component_frac",
        "expansion_defect",
        "backdoor_density",
    ],
    "scalar_structural": ["F_struct_scalar", "R_rob", "H", "B", "S", "V", "density_peak"],
    "graph_local_structural": [
        "H",
        "B",
        "S",
        "V",
        "R_resp",
        "R_rob",
        "F_struct_scalar",
        "propagation_depth_mean",
        "contradiction_rate_mean",
        "covariance_condition_log",
        "expansion_defect",
        "backdoor_density",
        "density_peak",
    ],
    "v14_defect_field": [
        "H",
        "B",
        "S",
        "V",
        "R_resp",
        "R_rob",
        "R_sat",
        "F_struct_scalar",
        "frust_mean",
        "frust_std",
        "frust_q75",
        "frust_q90",
        "frust_q95",
        "frust_q99",
        "frust_max",
        "frust_tail95_mean",
        "frust_tail_mass",
        "frust_gini",
        "hotspot_fraction",
        "hotspot_cluster_frac",
        "hotspot_edge_density",
        "R_tail",
        "R_var",
        "R_cluster",
        "F_mean",
        "F_tail",
        "F_cluster",
        "F_scale",
        "F_struct_multi",
        "density_peak",
    ],
}


def score_predictions(y: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    rho = spearmanr(y, pred).correlation
    return {
        "r2": float(r2_score(y, pred)),
        "rmse": float(math.sqrt(mean_squared_error(y, pred))),
        "spearman": float(0.0 if np.isnan(rho) else rho),
    }


def evaluate_models(df: pd.DataFrame, rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    y = df["cdcl_hardness"].to_numpy()
    kfold = KFold(n_splits=5, shuffle=True, random_state=SEED)
    rows = []
    predictions = {}
    for name, cols in FEATURE_SETS.items():
        X = df[cols].to_numpy()
        models = {
            "ridge": make_pipeline(StandardScaler(), Ridge(alpha=2.0)),
            "rf": RandomForestRegressor(n_estimators=260, min_samples_leaf=4, random_state=SEED, n_jobs=1),
        }
        for model_name, model in models.items():
            pred = cross_val_predict(model, X, y, cv=kfold, n_jobs=None)
            metrics = score_predictions(y, pred)
            rows.append({"feature_set": name, "model": model_name, "n_features": len(cols), **metrics})
            predictions[(name, model_name)] = pred

    cols = FEATURE_SETS["v14_defect_field"]
    shuffle_scores = []
    for rep in range(16):
        X = df[cols].to_numpy().copy()
        for j in range(X.shape[1]):
            rng.shuffle(X[:, j])
        model = RandomForestRegressor(n_estimators=220, min_samples_leaf=4, random_state=SEED + rep, n_jobs=1)
        pred = cross_val_predict(model, X, y, cv=kfold)
        metrics = score_predictions(y, pred)
        shuffle_scores.append(metrics)
    rows.append(
        {
            "feature_set": "shuffled_v14_null",
            "model": "rf",
            "n_features": len(cols),
            "r2": float(np.mean([m["r2"] for m in shuffle_scores])),
            "rmse": float(np.mean([m["rmse"] for m in shuffle_scores])),
            "spearman": float(np.mean([m["spearman"] for m in shuffle_scores])),
        }
    )
    results = pd.DataFrame(rows).sort_values(["model", "r2"], ascending=[True, False])

    # Leave-family-out is intentionally stricter: it asks whether the feature
    # language extrapolates to an unseen canonical SAT family.
    logo = LeaveOneGroupOut()
    lfo_rows = []
    groups = df["family"].to_numpy()
    for name, cols in FEATURE_SETS.items():
        X = df[cols].to_numpy()
        pred = np.full_like(y, fill_value=np.nan, dtype=float)
        for train, test in logo.split(X, y, groups):
            model = RandomForestRegressor(n_estimators=260, min_samples_leaf=4, random_state=SEED, n_jobs=1)
            model.fit(X[train], y[train])
            pred[test] = model.predict(X[test])
            metrics = score_predictions(y[test], pred[test])
            lfo_rows.append({"feature_set": name, "held_out_family": groups[test][0], "n_test": len(test), **metrics})
        overall = score_predictions(y, pred)
        lfo_rows.append({"feature_set": name, "held_out_family": "ALL_LFO", "n_test": len(y), **overall})
    lfo = pd.DataFrame(lfo_rows)

    best_pred = predictions.get(("v14_defect_field", "rf"))
    pred_df = df[["instance_id", "family", "cdcl_hardness"]].copy()
    pred_df["pred_v14_rf"] = best_pred
    return results, lfo, pred_df


def make_plots(df: pd.DataFrame, results: pd.DataFrame, lfo: pd.DataFrame, pred_df: pd.DataFrame) -> pd.DataFrame:
    plt.style.use("seaborn-v0_8-whitegrid")

    fig, ax = plt.subplots(figsize=(10, 5))
    order = sorted(df["family"].unique())
    ax.boxplot([df.loc[df["family"] == fam, "cdcl_hardness"] for fam in order], tick_labels=order, vert=True, patch_artist=True)
    ax.set_ylabel("CDCL hardness proxy")
    ax.set_title("Paper 55: CDCL hardness by canonical SAT family")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig1_hardness_by_family.png", dpi=180)
    plt.close(fig)

    rf_results = results[results["model"] == "rf"].copy()
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.bar(rf_results["feature_set"], rf_results["r2"], color="#2f6f9f")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("5-fold CV R2")
    ax.set_title("Random-forest model comparison on CDCL hardness")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig2_model_comparison.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5.5))
    for fam, sub in pred_df.groupby("family"):
        ax.scatter(sub["cdcl_hardness"], sub["pred_v14_rf"], s=18, alpha=0.7, label=fam)
    lo = min(pred_df["cdcl_hardness"].min(), pred_df["pred_v14_rf"].min())
    hi = max(pred_df["cdcl_hardness"].max(), pred_df["pred_v14_rf"].max())
    ax.plot([lo, hi], [lo, hi], color="black", linewidth=1.0, linestyle="--")
    ax.set_xlabel("Observed CDCL hardness")
    ax.set_ylabel("Predicted hardness")
    ax.set_title("v1.4 defect-field prediction scatter")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig3_v14_prediction_scatter.png", dpi=180)
    plt.close(fig)

    cols = FEATURE_SETS["v14_defect_field"]
    X = df[cols].to_numpy()
    y = df["cdcl_hardness"].to_numpy()
    model = RandomForestRegressor(n_estimators=360, min_samples_leaf=4, random_state=SEED, n_jobs=1)
    model.fit(X, y)
    perm = permutation_importance(model, X, y, n_repeats=12, random_state=SEED, n_jobs=1)
    importance = pd.DataFrame({"feature": cols, "importance_mean": perm.importances_mean, "importance_std": perm.importances_std})
    importance = importance.sort_values("importance_mean", ascending=False)
    top = importance.head(12).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.barh(top["feature"], top["importance_mean"], xerr=top["importance_std"], color="#7a5195")
    ax.set_xlabel("Permutation importance")
    ax.set_title("Top v1.4 defect-field features")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig4_permutation_importance.png", dpi=180)
    plt.close(fig)

    lfo_all = lfo[lfo["held_out_family"] == "ALL_LFO"].copy()
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.bar(lfo_all["feature_set"], lfo_all["spearman"], color="#ef6f6c")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Leave-family-out Spearman")
    ax.set_title("Strict family-transfer stress test")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig5_leave_family_out.png", dpi=180)
    plt.close(fig)
    return importance


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    raw_instances = build_instances(rng)
    records = []
    for idx, (family, n_vars, clauses) in enumerate(raw_instances):
        records.append(make_record(idx, family, n_vars, clauses, rng))
        if (idx + 1) % 50 == 0:
            print(f"solved {idx + 1}/{len(raw_instances)} instances")

    df = pd.DataFrame([asdict(r) for r in records])
    df.to_csv(OUTDIR / "paper55_sat_cdcl_instances.csv", index=False)

    results, lfo, pred_df = evaluate_models(df, rng)
    results.to_csv(OUTDIR / "paper55_model_results.csv", index=False)
    lfo.to_csv(OUTDIR / "paper55_leave_family_out.csv", index=False)
    pred_df.to_csv(OUTDIR / "paper55_predictions.csv", index=False)
    importance = make_plots(df, results, lfo, pred_df)
    importance.to_csv(OUTDIR / "paper55_permutation_importance.csv", index=False)

    rf = results[results["model"] == "rf"].set_index("feature_set")
    best_feature_set = str(rf["r2"].idxmax())
    summary = {
        "seed": SEED,
        "n_instances": int(len(df)),
        "families": sorted(df["family"].unique().tolist()),
        "conflict_budget": CONFLICT_BUDGET,
        "target": "log1p(mean_conflicts + 0.10*mean_decisions + 0.01*mean_propagations)",
        "best_rf_feature_set": best_feature_set,
        "rf_results": rf[["r2", "rmse", "spearman"]].to_dict(orient="index"),
        "lfo_overall": lfo[lfo["held_out_family"] == "ALL_LFO"].set_index("feature_set")[["r2", "rmse", "spearman"]].to_dict(orient="index"),
        "timeouts_glucose": int(df["timeout_glucose"].sum()),
        "timeouts_minisat": int(df["timeout_minisat"].sum()),
        "sat_agreement_fraction": float(np.mean(df["sat_glucose"] == df["sat_minisat"])),
        "top_importance": importance.head(8).to_dict(orient="records"),
    }
    with (OUTDIR / "paper55_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
