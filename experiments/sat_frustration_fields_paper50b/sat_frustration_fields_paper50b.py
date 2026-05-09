#!/usr/bin/env python3
"""
Paper 50b: SAT Structural Hardness IIb.

This experiment tests the hypothesis suggested by Paper 50:
global MAAT scalar compression is too lossy for SAT hardness, while local
frustration fields, rare-event tails, and hotspot clusters retain signal.

The target is still a deterministic DPLL node-count proxy.  The structural
features are computed from the formula graph and unit-propagation probes, not
from the recursive solver trace.

Run:
    python3 sat_frustration_fields_paper50b.py
"""

from __future__ import annotations

import csv
import json
import math
import time
from collections import deque
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


SEED = 51051
OUTDIR = Path("outputs_paper50b")
EPS = 1.0e-9


@dataclass
class InstanceRecord:
    instance_id: int
    formula_type: str
    n_vars: int
    n_clauses: int
    alpha: float
    satisfiable: bool
    dpll_nodes: int
    dpll_timed_out: bool
    solve_time_s: float
    log_nodes: float
    degree_cv: float
    max_degree_norm: float
    literal_imbalance: float
    interaction_cycle_density: float
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
    F_MAAT: float
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
    F_frust: float
    F_mean: float
    F_tail: float
    F_cluster: float
    F_scale: float
    F_MAAT_multi: float


def generate_random_3sat(rng: np.random.Generator, n_vars: int, n_clauses: int) -> list[tuple[int, int, int]]:
    clauses: set[tuple[int, int, int]] = set()
    attempts = 0
    while len(clauses) < n_clauses and attempts < 80 * n_clauses:
        attempts += 1
        vars_ = rng.choice(np.arange(1, n_vars + 1), size=3, replace=False)
        lits = [int((-1 if rng.random() < 0.5 else 1) * v) for v in vars_]
        key = tuple(sorted(lits, key=lambda x: abs(x)))
        if len({abs(x) for x in key}) == 3:
            clauses.add(key)
    return list(clauses)


def generate_planted_3sat(rng: np.random.Generator, n_vars: int, n_clauses: int) -> list[tuple[int, int, int]]:
    assignment = rng.random(n_vars) < 0.5
    clauses: set[tuple[int, int, int]] = set()
    attempts = 0
    while len(clauses) < n_clauses and attempts < 120 * n_clauses:
        attempts += 1
        vars_ = rng.choice(np.arange(1, n_vars + 1), size=3, replace=False)
        lits = []
        is_satisfied = False
        for v in vars_:
            sign = -1 if rng.random() < 0.5 else 1
            lit = int(sign * v)
            val = bool(assignment[v - 1])
            if (lit > 0 and val) or (lit < 0 and not val):
                is_satisfied = True
            lits.append(lit)
        if not is_satisfied:
            # Flip one literal so the hidden assignment satisfies the clause.
            j = int(rng.integers(0, 3))
            lits[j] = -lits[j]
        key = tuple(sorted(lits, key=lambda x: abs(x)))
        if len({abs(x) for x in key}) == 3:
            clauses.add(key)
    return list(clauses)


def generate_modular_3sat(rng: np.random.Generator, n_vars: int, n_clauses: int, p_internal: float = 0.78) -> list[tuple[int, int, int]]:
    n_blocks = max(3, min(6, n_vars // 7))
    blocks = np.array_split(np.arange(1, n_vars + 1), n_blocks)
    clauses: set[tuple[int, int, int]] = set()
    attempts = 0
    while len(clauses) < n_clauses and attempts < 120 * n_clauses:
        attempts += 1
        if rng.random() < p_internal:
            block = blocks[int(rng.integers(0, len(blocks)))]
            if len(block) < 3:
                pool = np.arange(1, n_vars + 1)
            else:
                pool = block
            vars_ = rng.choice(pool, size=3, replace=False)
        else:
            vars_ = rng.choice(np.arange(1, n_vars + 1), size=3, replace=False)
        lits = [int((-1 if rng.random() < 0.5 else 1) * v) for v in vars_]
        key = tuple(sorted(lits, key=lambda x: abs(x)))
        if len({abs(x) for x in key}) == 3:
            clauses.add(key)
    return list(clauses)


def simplify_clauses(clauses: list[tuple[int, ...]], assignment: dict[int, bool]) -> tuple[bool | None, list[tuple[int, ...]]]:
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
            return False, []
        simplified.append(tuple(remaining))
    if not simplified:
        return True, []
    return None, simplified


def unit_propagate(clauses: list[tuple[int, ...]], assignment: dict[int, bool]) -> tuple[bool | None, list[tuple[int, ...]], dict[int, bool], int]:
    assignment = dict(assignment)
    forced = 0
    while True:
        status, clauses = simplify_clauses(clauses, assignment)
        if status is not None:
            return status, clauses, assignment, forced
        units = [c[0] for c in clauses if len(c) == 1]
        if not units:
            return None, clauses, assignment, forced
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
            return None, clauses, assignment, forced


def choose_branch_var(clauses: list[tuple[int, ...]], assignment: dict[int, bool]) -> int:
    counts: dict[int, int] = {}
    polarity: dict[int, int] = {}
    for clause in clauses:
        for lit in clause:
            var = abs(lit)
            if var in assignment:
                continue
            counts[var] = counts.get(var, 0) + 1
            polarity[var] = polarity.get(var, 0) + (1 if lit > 0 else -1)
    if not counts:
        return 0
    return max(counts, key=lambda v: (counts[v], abs(polarity.get(v, 0)), -v))


def dpll_solve(clauses: list[tuple[int, ...]], node_limit: int = 90000) -> tuple[bool, int, bool]:
    nodes = 0

    def rec(cur_clauses: list[tuple[int, ...]], assignment: dict[int, bool]) -> bool | None:
        nonlocal nodes
        nodes += 1
        if nodes > node_limit:
            return None
        status, simp, assign2, _ = unit_propagate(cur_clauses, assignment)
        if status is not None:
            return status
        var = choose_branch_var(simp, assign2)
        if var == 0:
            return True
        pos = sum(1 for c in simp for lit in c if lit == var)
        neg = sum(1 for c in simp for lit in c if lit == -var)
        first = pos >= neg
        for val in (first, not first):
            assign3 = dict(assign2)
            assign3[var] = val
            res = rec(simp, assign3)
            if res is True:
                return True
            if res is None:
                return None
        return False

    result = rec([tuple(c) for c in clauses], {})
    return bool(result), nodes, result is None


def graph_arrays(clauses: list[tuple[int, int, int]], n_vars: int) -> dict:
    deg = np.zeros(n_vars, dtype=float)
    pos = np.zeros(n_vars, dtype=float)
    neg = np.zeros(n_vars, dtype=float)
    var_to_clauses: list[set[int]] = [set() for _ in range(n_vars)]
    neighbors: list[set[int]] = [set() for _ in range(n_vars)]

    for ci, clause in enumerate(clauses):
        vars_ = [abs(l) - 1 for l in clause]
        for lit in clause:
            idx = abs(lit) - 1
            deg[idx] += 1
            if lit > 0:
                pos[idx] += 1
            else:
                neg[idx] += 1
            var_to_clauses[idx].add(ci)
        for a, b in combinations(vars_, 2):
            neighbors[a].add(b)
            neighbors[b].add(a)

    local_clustering = np.zeros(n_vars, dtype=float)
    tri_total = 0.0
    wedge_total = 0.0
    for v in range(n_vars):
        nb = list(neighbors[v])
        d = len(nb)
        wedges = d * (d - 1) / 2.0
        wedge_total += wedges
        if d >= 2:
            tri = 0
            for a, b in combinations(nb, 2):
                if b in neighbors[a]:
                    tri += 1
            local_clustering[v] = tri / (wedges + EPS)
            tri_total += tri

    degree_cv = float(np.std(deg) / (np.mean(deg) + EPS))
    max_degree_norm = float(np.max(deg) / (np.mean(deg) + EPS))
    literal_imbalance = float(np.mean(np.abs(pos - neg) / (deg + EPS)))
    interaction_cycle_density = float((tri_total / 3.0) / (wedge_total + EPS))

    local = []
    for v in range(n_vars):
        nb = list(neighbors[v])
        nb_degs = deg[nb] if nb else np.array([0.0])
        local.append(
            [
                deg[v],
                abs(pos[v] - neg[v]) / (deg[v] + EPS),
                len(nb),
                float(np.mean(nb_degs)),
                float(np.std(nb_degs)),
                local_clustering[v],
            ]
        )
    local_arr = np.asarray(local, dtype=float)
    cov = np.cov(local_arr.T) + 1.0e-6 * np.eye(local_arr.shape[1])
    eig = np.linalg.eigvalsh(cov)
    covariance_condition_log = float(np.log10((float(np.max(eig)) + EPS) / (float(np.min(eig)) + EPS)))

    return {
        "deg": deg,
        "pos": pos,
        "neg": neg,
        "var_to_clauses": var_to_clauses,
        "neighbors": neighbors,
        "local_clustering": local_clustering,
        "degree_cv": degree_cv,
        "max_degree_norm": max_degree_norm,
        "literal_imbalance": literal_imbalance,
        "interaction_cycle_density": interaction_cycle_density,
        "covariance_condition_log": covariance_condition_log,
    }


def greedy_backdoor_marks(var_to_clauses: list[set[int]], n_clauses: int, target_fraction: float = 0.80) -> tuple[np.ndarray, float]:
    n_vars = len(var_to_clauses)
    uncovered = set(range(n_clauses))
    selected: list[int] = []
    while uncovered and (n_clauses - len(uncovered)) / max(1, n_clauses) < target_fraction:
        best_v = max(range(n_vars), key=lambda v: len(var_to_clauses[v] & uncovered))
        gain = len(var_to_clauses[best_v] & uncovered)
        if gain == 0:
            break
        selected.append(best_v)
        uncovered -= var_to_clauses[best_v]
    marks = np.zeros(n_vars, dtype=float)
    marks[selected] = 1.0
    return marks, float(len(selected) / max(1, n_vars))


def propagation_fields(clauses: list[tuple[int, int, int]], n_vars: int) -> tuple[np.ndarray, np.ndarray]:
    forced = np.zeros(n_vars, dtype=float)
    contradictions = np.zeros(n_vars, dtype=float)
    for var in range(1, n_vars + 1):
        fvals = []
        cvals = []
        for val in (True, False):
            status, _, _, count = unit_propagate([tuple(c) for c in clauses], {var: val})
            fvals.append(count / max(1, n_vars))
            cvals.append(1.0 if status is False else 0.0)
        forced[var - 1] = float(np.mean(fvals))
        contradictions[var - 1] = float(np.mean(cvals))
    return forced, contradictions


def gini(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    if len(x) == 0:
        return 0.0
    x = np.clip(x, 0.0, None)
    if np.sum(x) < EPS:
        return 0.0
    sx = np.sort(x)
    n = len(sx)
    return float((2.0 * np.sum((np.arange(1, n + 1)) * sx) / (n * np.sum(sx))) - (n + 1) / n)


def largest_hotspot_cluster(hotspot: np.ndarray, neighbors: list[set[int]]) -> tuple[float, float]:
    idxs = set(np.where(hotspot)[0].tolist())
    if not idxs:
        return 0.0, 0.0
    seen: set[int] = set()
    largest = 0
    internal_edges = 0
    possible_edges = len(idxs) * (len(idxs) - 1) / 2.0
    for v in idxs:
        internal_edges += sum(1 for u in neighbors[v] if u in idxs)
    internal_edges /= 2.0
    for start in list(idxs):
        if start in seen:
            continue
        q = deque([start])
        seen.add(start)
        size = 0
        while q:
            v = q.popleft()
            size += 1
            for u in neighbors[v]:
                if u in idxs and u not in seen:
                    seen.add(u)
                    q.append(u)
        largest = max(largest, size)
    return float(largest / len(hotspot)), float(internal_edges / (possible_edges + EPS))


def compute_features(clauses: list[tuple[int, int, int]], n_vars: int, alpha: float) -> dict[str, float]:
    g = graph_arrays(clauses, n_vars)
    backdoor_marks, backdoor_density = greedy_backdoor_marks(g["var_to_clauses"], len(clauses))
    prop_depth, contradiction = propagation_fields(clauses, n_vars)

    deg = g["deg"]
    deg_z = np.abs(deg - np.mean(deg)) / (np.std(deg) + EPS)
    degree_anomaly = np.tanh(deg_z / 2.0)
    polarity = np.abs(g["pos"] - g["neg"]) / (deg + EPS)
    clustering = g["local_clustering"]

    local_frust = (
        0.22 * degree_anomaly
        + 0.16 * polarity
        + 0.18 * clustering
        + 2.10 * prop_depth
        + 0.58 * contradiction
        + 0.22 * backdoor_marks
    )
    local_frust = np.clip(local_frust, 0.0, None)

    q75, q90, q95, q99 = np.percentile(local_frust, [75, 90, 95, 99])
    tail = local_frust[local_frust >= q95]
    hotspot = local_frust >= max(q90, float(np.mean(local_frust) + np.std(local_frust)))
    hotspot_cluster_frac, hotspot_edge_density = largest_hotspot_cluster(hotspot, g["neighbors"])
    tail_mass = float(np.sum(local_frust[hotspot]) / (np.sum(local_frust) + EPS))

    frust_std = float(np.std(local_frust))
    R_tail = 1.0 / (1.0 + float(q95))
    R_var = math.exp(-frust_std)
    R_cluster = math.exp(-2.0 * hotspot_cluster_frac)
    R_sat = (R_tail * R_var * R_cluster) ** (1.0 / 3.0)
    F_frust = -sum(math.log(EPS + r) for r in (R_tail, R_var, R_cluster, R_sat))
    # Paper-50c multi-scale defect-field functional components.
    F_mean = float(np.mean(local_frust))
    F_tail = float(q95 + 0.50 * q99 + 0.25 * tail_mass)
    F_cluster = float(hotspot_cluster_frac + 0.50 * hotspot_edge_density)
    F_scale = float(0.45 * frust_std + 0.35 * gini(local_frust) + 0.20 * abs(float(q95) - float(q75)))
    F_MAAT_multi = F_mean + F_tail + F_cluster + F_scale

    # Paper-50-style global support compression, retained for ablation.
    d_H = 2.0 * float(np.mean(contradiction)) + 0.25 * g["covariance_condition_log"]
    d_B = 1.3 * g["literal_imbalance"] + 0.75 * g["degree_cv"]
    s_eff = math.exp(-0.5 * ((float(np.mean(prop_depth)) - 0.18) / 0.10) ** 2)
    d_S = (1.0 - s_eff) + 0.20 * g["max_degree_norm"] / 6.0
    expansion_defect = 1.0 / (1.0 + float(np.mean([len(n) for n in g["neighbors"]])) / (n_vars / 3.0 + EPS))
    d_V = 1.5 * expansion_defect + 2.0 * g["interaction_cycle_density"] + 0.8 * backdoor_density
    H = 1.0 / (1.0 + d_H)
    B = 1.0 / (1.0 + d_B)
    S = 1.0 / (1.0 + d_S)
    V = 1.0 / (1.0 + d_V)
    R_resp = (H * B * V) ** (1.0 / 3.0)
    R_rob = min(R_resp, (H * B * S * V) ** 0.25)
    F_MAAT = -sum(math.log(EPS + x) for x in (H, B, S, V, R_rob))
    density_peak = math.exp(-0.5 * ((alpha - 4.26) / 0.38) ** 2)

    return {
        "degree_cv": g["degree_cv"],
        "max_degree_norm": g["max_degree_norm"],
        "literal_imbalance": g["literal_imbalance"],
        "interaction_cycle_density": g["interaction_cycle_density"],
        "backdoor_density": backdoor_density,
        "propagation_depth_mean": float(np.mean(prop_depth)),
        "contradiction_rate_mean": float(np.mean(contradiction)),
        "covariance_condition_log": g["covariance_condition_log"],
        "H": H,
        "B": B,
        "S": S,
        "V": V,
        "R_resp": R_resp,
        "R_rob": R_rob,
        "F_MAAT": F_MAAT,
        "density_peak": density_peak,
        "frust_mean": float(np.mean(local_frust)),
        "frust_std": frust_std,
        "frust_q75": float(q75),
        "frust_q90": float(q90),
        "frust_q95": float(q95),
        "frust_q99": float(q99),
        "frust_max": float(np.max(local_frust)),
        "frust_tail95_mean": float(np.mean(tail)) if len(tail) else 0.0,
        "frust_tail_mass": tail_mass,
        "frust_gini": gini(local_frust),
        "hotspot_fraction": float(np.mean(hotspot)),
        "hotspot_cluster_frac": hotspot_cluster_frac,
        "hotspot_edge_density": hotspot_edge_density,
        "R_tail": R_tail,
        "R_var": R_var,
        "R_cluster": R_cluster,
        "R_sat": R_sat,
        "F_frust": F_frust,
        "F_mean": F_mean,
        "F_tail": F_tail,
        "F_cluster": F_cluster,
        "F_scale": F_scale,
        "F_MAAT_multi": F_MAAT_multi,
    }


def generate_dataset() -> list[InstanceRecord]:
    rng = np.random.default_rng(SEED)
    records: list[InstanceRecord] = []
    plan = [("random", 330), ("planted", 110), ("modular", 110)]
    alpha_grid = np.linspace(3.35, 5.10, 10)

    idx = 0
    for formula_type, count in plan:
        for _ in range(count):
            n_vars = int(rng.integers(24, 39))
            alpha = float(alpha_grid[idx % len(alpha_grid)] + rng.normal(0.0, 0.04))
            n_clauses = int(round(alpha * n_vars))
            if formula_type == "random":
                clauses = generate_random_3sat(rng, n_vars, n_clauses)
            elif formula_type == "planted":
                clauses = generate_planted_3sat(rng, n_vars, n_clauses)
            elif formula_type == "modular":
                clauses = generate_modular_3sat(rng, n_vars, n_clauses)
            else:
                raise ValueError(formula_type)

            t0 = time.perf_counter()
            sat, nodes, timed_out = dpll_solve(clauses)
            dt = time.perf_counter() - t0
            features = compute_features(clauses, n_vars, alpha)
            records.append(
                InstanceRecord(
                    instance_id=idx,
                    formula_type=formula_type,
                    n_vars=n_vars,
                    n_clauses=n_clauses,
                    alpha=alpha,
                    satisfiable=sat,
                    dpll_nodes=nodes,
                    dpll_timed_out=timed_out,
                    solve_time_s=dt,
                    log_nodes=math.log1p(nodes),
                    **features,
                )
            )
            idx += 1
    return records


def one_hot_types(type_values: list[str]) -> tuple[np.ndarray, list[str]]:
    enc = OneHotEncoder(sparse_output=False, drop=None)
    arr = enc.fit_transform(np.asarray(type_values).reshape(-1, 1))
    names = [f"type_{cat}" for cat in enc.categories_[0]]
    return arr, names


def feature_matrix(rows: list[dict], feature_names: list[str]) -> tuple[np.ndarray, list[str]]:
    numeric = [f for f in feature_names if f != "formula_type"]
    X_num = np.asarray([[row[f] for f in numeric] for row in rows], dtype=float) if numeric else np.empty((len(rows), 0))
    names = list(numeric)
    if "formula_type" in feature_names:
        X_type, type_names = one_hot_types([row["formula_type"] for row in rows])
        X_num = np.column_stack([X_num, X_type])
        names += type_names
    return X_num, names


def corr(x: np.ndarray, y: np.ndarray) -> float:
    if np.std(x) < EPS or np.std(y) < EPS:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def rankdata(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x)
    ranks = np.empty_like(order, dtype=float)
    i = 0
    while i < len(x):
        j = i + 1
        while j < len(x) and x[order[j]] == x[order[i]]:
            j += 1
        ranks[order[i:j]] = 0.5 * (i + j - 1) + 1.0
        i = j
    return ranks


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return corr(rankdata(x), rankdata(y))


def evaluate_model(rows: list[dict], model_name: str, feature_names: list[str], estimator: str) -> dict[str, float | str | int]:
    X, expanded_names = feature_matrix(rows, feature_names)
    y = np.asarray([row["log_nodes"] for row in rows], dtype=float)
    preds = np.zeros_like(y)
    kf = KFold(n_splits=5, shuffle=True, random_state=SEED + 7)
    for train, test in kf.split(X):
        if estimator == "ridge":
            model = make_pipeline(StandardScaler(), Ridge(alpha=1.0))
        elif estimator == "random_forest":
            model = RandomForestRegressor(
                n_estimators=220,
                max_depth=9,
                min_samples_leaf=5,
                random_state=SEED + 11,
                n_jobs=-1,
            )
        else:
            raise ValueError(estimator)
        model.fit(X[train], y[train])
        preds[test] = model.predict(X[test])
    return {
        "model": model_name,
        "estimator": estimator,
        "n_features": len(expanded_names),
        "features": "+".join(feature_names),
        "cv_r2": float(r2_score(y, preds)),
        "pearson": corr(preds, y),
        "spearman": spearman(preds, y),
        "rmse": float(math.sqrt(mean_squared_error(y, preds))),
    }


def feature_importance(rows: list[dict], feature_names: list[str]) -> list[dict[str, float | str]]:
    X, expanded_names = feature_matrix(rows, feature_names)
    y = np.asarray([row["log_nodes"] for row in rows], dtype=float)
    model = RandomForestRegressor(
        n_estimators=280,
        max_depth=9,
        min_samples_leaf=5,
        random_state=SEED + 19,
        n_jobs=-1,
    )
    model.fit(X, y)
    perm = permutation_importance(model, X, y, n_repeats=12, random_state=SEED + 23, n_jobs=-1)
    rows_out = []
    for name, imp, std, builtin in sorted(
        zip(expanded_names, perm.importances_mean, perm.importances_std, model.feature_importances_),
        key=lambda t: t[1],
        reverse=True,
    ):
        rows_out.append(
            {
                "feature": name,
                "permutation_importance_mean": float(imp),
                "permutation_importance_std": float(std),
                "rf_builtin_importance": float(builtin),
            }
        )
    return rows_out


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_results(rows: list[dict], model_rows: list[dict], importances: list[dict]) -> None:
    y = np.asarray([r["log_nodes"] for r in rows], dtype=float)
    q95 = np.asarray([r["frust_q95"] for r in rows], dtype=float)
    cluster = np.asarray([r["hotspot_cluster_frac"] for r in rows], dtype=float)
    alpha = np.asarray([r["alpha"] for r in rows], dtype=float)
    type_color = {"random": "#0ea5e9", "planted": "#22c55e", "modular": "#f97316"}
    colors = [type_color[r["formula_type"]] for r in rows]

    rf_rows = [r for r in model_rows if r["estimator"] == "random_forest"]
    fig, ax = plt.subplots(figsize=(9.8, 5.2))
    names = [r["model"] for r in rf_rows]
    vals = [r["cv_r2"] for r in rf_rows]
    ax.bar(names, vals, color="#334155", edgecolor="#0f172a")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel(r"5-fold CV $R^2$")
    ax.set_title("Paper 50b compression ladder: random forest")
    ax.tick_params(axis="x", rotation=32)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig1_compression_ladder_rf.png", dpi=260)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    ax.scatter(q95, y, c=colors, s=26, alpha=0.78, edgecolor="none")
    ax.set_xlabel(r"local frustration $q_{95}$")
    ax.set_ylabel(r"$\log(1+\mathrm{DPLL\ nodes})$")
    ax.set_title("Rare-event frustration tails vs SAT hardness")
    for name, color in type_color.items():
        ax.scatter([], [], c=color, label=name)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig2_tail_q95_vs_nodes.png", dpi=260)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    sc = ax.scatter(cluster, y, c=alpha, s=28, cmap="magma", alpha=0.80, edgecolor="none")
    ax.set_xlabel("largest hotspot cluster fraction")
    ax.set_ylabel(r"$\log(1+\mathrm{DPLL\ nodes})$")
    ax.set_title("Frustration-hotspot clusters vs SAT hardness")
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(r"clause density $\alpha$")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig3_hotspot_cluster_vs_nodes.png", dpi=260)
    plt.close(fig)

    top = importances[:14]
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    ax.barh([r["feature"] for r in reversed(top)], [r["permutation_importance_mean"] for r in reversed(top)], color="#14b8a6")
    ax.set_xlabel("permutation importance")
    ax.set_title("Top random-forest features for local frustration model")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig4_feature_importance.png", dpi=260)
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    records = generate_dataset()
    rows = [asdict(r) for r in records]
    write_csv(OUTDIR / "paper50b_sat_instances.csv", rows)

    models = {
        "density_only": ["alpha", "density_peak"],
        "density_plus_type": ["alpha", "density_peak", "formula_type"],
        "standard_graph": [
            "alpha",
            "degree_cv",
            "max_degree_norm",
            "literal_imbalance",
            "interaction_cycle_density",
            "formula_type",
        ],
        "global_MAAT_scalar": ["F_MAAT", "R_rob", "H", "B", "S", "V"],
        "frustration_tail_only": ["frust_q90", "frust_q95", "frust_q99", "frust_tail_mass", "hotspot_cluster_frac"],
        "frustration_robustness": ["R_tail", "R_var", "R_cluster", "R_sat", "F_frust"],
        "multi_scalar_50c": ["F_MAAT_multi"],
        "multi_components_50c": ["F_mean", "F_tail", "F_cluster", "F_scale"],
        "local_frustration_fields": [
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
            "R_sat",
            "F_frust",
        ],
        "local_frustration_plus_density": [
            "alpha",
            "density_peak",
            "formula_type",
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
            "R_sat",
            "F_frust",
        ],
        "multi_plus_density_50c": [
            "alpha",
            "density_peak",
            "formula_type",
            "F_mean",
            "F_tail",
            "F_cluster",
            "F_scale",
            "F_MAAT_multi",
        ],
        "raw_graph_plus_frustration": [
            "alpha",
            "density_peak",
            "formula_type",
            "degree_cv",
            "max_degree_norm",
            "literal_imbalance",
            "interaction_cycle_density",
            "backdoor_density",
            "propagation_depth_mean",
            "contradiction_rate_mean",
            "covariance_condition_log",
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
            "R_sat",
            "F_frust",
            "F_mean",
            "F_tail",
            "F_cluster",
            "F_scale",
            "F_MAAT_multi",
            "F_MAAT",
            "R_rob",
        ],
    }

    model_rows = []
    for estimator in ("ridge", "random_forest"):
        for name, features in models.items():
            model_rows.append(evaluate_model(rows, name, features, estimator))
    write_csv(OUTDIR / "paper50b_model_comparison.csv", model_rows)

    importance_features = models["raw_graph_plus_frustration"]
    importances = feature_importance(rows, importance_features)
    write_csv(OUTDIR / "paper50b_feature_importance.csv", importances)

    y = np.asarray([r["log_nodes"] for r in rows])
    scalar_scores = {
        "F_MAAT": np.asarray([r["F_MAAT"] for r in rows]),
        "R_rob_inverse": 1.0 - np.asarray([r["R_rob"] for r in rows]),
        "F_frust": np.asarray([r["F_frust"] for r in rows]),
        "F_MAAT_multi": np.asarray([r["F_MAAT_multi"] for r in rows]),
        "R_sat_inverse": 1.0 - np.asarray([r["R_sat"] for r in rows]),
        "frust_q95": np.asarray([r["frust_q95"] for r in rows]),
        "hotspot_cluster_frac": np.asarray([r["hotspot_cluster_frac"] for r in rows]),
        "covariance_condition_log": np.asarray([r["covariance_condition_log"] for r in rows]),
    }
    scalar_rows = [
        {"score": name, "pearson": corr(score, y), "spearman": spearman(score, y)}
        for name, score in scalar_scores.items()
    ]
    write_csv(OUTDIR / "paper50b_scalar_correlations.csv", scalar_rows)

    def row(model: str, estimator: str) -> dict:
        return next(r for r in model_rows if r["model"] == model and r["estimator"] == estimator)

    best_row = max(model_rows, key=lambda r: r["cv_r2"])
    summary = {
        "model": "Paper 50b SAT Frustration-Field Compression Test",
        "status": "Synthetic multi-family SAT frustration-field benchmark; not a proof of P vs NP and not a modern CDCL benchmark.",
        "random_seed": SEED,
        "n_instances": len(records),
        "formula_types": {t: sum(1 for r in records if r.formula_type == t) for t in sorted({r.formula_type for r in records})},
        "n_vars_range": [min(r.n_vars for r in records), max(r.n_vars for r in records)],
        "alpha_range": [float(min(r.alpha for r in records)), float(max(r.alpha for r in records))],
        "sat_fraction": float(np.mean([r.satisfiable for r in records])),
        "timeout_fraction": float(np.mean([r.dpll_timed_out for r in records])),
        "best_model": best_row,
        "ridge_density_r2": row("density_only", "ridge")["cv_r2"],
        "ridge_global_maat_r2": row("global_MAAT_scalar", "ridge")["cv_r2"],
        "ridge_local_frustration_density_r2": row("local_frustration_plus_density", "ridge")["cv_r2"],
        "rf_density_r2": row("density_only", "random_forest")["cv_r2"],
        "rf_global_maat_r2": row("global_MAAT_scalar", "random_forest")["cv_r2"],
        "rf_multi_scalar_50c_r2": row("multi_scalar_50c", "random_forest")["cv_r2"],
        "rf_multi_components_50c_r2": row("multi_components_50c", "random_forest")["cv_r2"],
        "rf_multi_plus_density_50c_r2": row("multi_plus_density_50c", "random_forest")["cv_r2"],
        "rf_local_frustration_density_r2": row("local_frustration_plus_density", "random_forest")["cv_r2"],
        "rf_raw_graph_frustration_r2": row("raw_graph_plus_frustration", "random_forest")["cv_r2"],
        "rf_multi_plus_density_minus_global_maat_r2": float(
            row("multi_plus_density_50c", "random_forest")["cv_r2"]
            - row("global_MAAT_scalar", "random_forest")["cv_r2"]
        ),
        "rf_local_minus_global_maat_r2": float(
            row("local_frustration_plus_density", "random_forest")["cv_r2"]
            - row("global_MAAT_scalar", "random_forest")["cv_r2"]
        ),
        "rf_local_minus_density_type_r2": float(
            row("local_frustration_plus_density", "random_forest")["cv_r2"]
            - row("density_plus_type", "random_forest")["cv_r2"]
        ),
        "key_interpretation": (
            "If local frustration fields outperform global MAAT scalar compression, "
            "SAT hardness is better interpreted as multi-scale frustration geometry "
            "than as macroscopic robustness."
        ),
        "outputs": [
            "paper50b_sat_instances.csv",
            "paper50b_model_comparison.csv",
            "paper50b_scalar_correlations.csv",
            "paper50b_feature_importance.csv",
            "paper50b_summary.json",
            "fig1_compression_ladder_rf.png",
            "fig2_tail_q95_vs_nodes.png",
            "fig3_hotspot_cluster_vs_nodes.png",
            "fig4_feature_importance.png",
        ],
    }
    with (OUTDIR / "paper50b_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    plot_results(rows, model_rows, importances)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
