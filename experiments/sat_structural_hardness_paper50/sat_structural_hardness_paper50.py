#!/usr/bin/env python3
"""
Paper 50: SAT Structural Hardness II.

This experiment directly addresses the SAT weakness identified in Paper 49.
It generates a deterministic random-3-SAT ensemble, computes graph-local and
propagation-local structural diagnostics, and compares MAAT-style structural
hardness against simple baselines.

The target is a bounded DPLL node-count proxy.  The MAAT features themselves
are graph/proxy diagnostics, not learned from the solver trace.

Run:
    python3 sat_structural_hardness_paper50.py
"""

from __future__ import annotations

import csv
import json
import math
import time
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np


SEED = 50050
OUTDIR = Path("outputs_paper50")
EPS = 1.0e-9


@dataclass
class InstanceRecord:
    instance_id: int
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
    expansion_defect: float
    interaction_cycle_density: float
    backdoor_density: float
    propagation_depth: float
    contradiction_rate: float
    covariance_condition_log: float
    H: float
    B: float
    S: float
    V: float
    R_resp: float
    R_rob: float
    F_MAAT: float
    F_MAAT_primary: float
    density_peak: float


def generate_random_3sat(rng: np.random.Generator, n_vars: int, n_clauses: int) -> list[tuple[int, int, int]]:
    clauses: set[tuple[int, int, int]] = set()
    attempts = 0
    while len(clauses) < n_clauses and attempts < 50 * n_clauses:
        attempts += 1
        vars_ = rng.choice(np.arange(1, n_vars + 1), size=3, replace=False)
        lits = []
        for v in vars_:
            sign = -1 if rng.random() < 0.5 else 1
            lits.append(int(sign * v))
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
    # Degree first, polarity imbalance as deterministic tie-breaker.
    return max(counts, key=lambda v: (counts[v], abs(polarity.get(v, 0)), -v))


def dpll_solve(clauses: list[tuple[int, ...]], n_vars: int, node_limit: int = 60000) -> tuple[bool, int, bool]:
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
        # Try dominant polarity first.
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
    timed_out = result is None
    return bool(result), nodes, timed_out


def graph_features(clauses: list[tuple[int, int, int]], n_vars: int, rng: np.random.Generator) -> dict[str, float]:
    m = len(clauses)
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

    degree_cv = float(np.std(deg) / (np.mean(deg) + EPS))
    max_degree_norm = float(np.max(deg) / (np.mean(deg) + EPS))
    literal_imbalance = float(np.mean(np.abs(pos - neg) / (deg + EPS)))

    # Expansion defect: sample variable subsets and measure clause-neighborhood expansion.
    expansion_vals = []
    sizes = [3, 5, 8]
    active_vars = np.where(deg > 0)[0]
    for size in sizes:
        if len(active_vars) < size:
            continue
        for _ in range(24):
            subset = rng.choice(active_vars, size=size, replace=False)
            touched = set()
            for v in subset:
                touched.update(var_to_clauses[int(v)])
            expansion_vals.append(len(touched) / (3.0 * size))
    mean_expansion = float(np.mean(expansion_vals)) if expansion_vals else 1.0
    expansion_defect = float(max(0.0, 1.0 - mean_expansion))

    # Interaction cycles proxy: triangle density in the variable co-occurrence graph.
    tri = 0
    wedges = 0
    for v in range(n_vars):
        nb = list(neighbors[v])
        d = len(nb)
        wedges += d * (d - 1) / 2
        if d >= 2:
            nb_set = neighbors[v]
            for a, b in combinations(nb, 2):
                if b in neighbors[a]:
                    tri += 1
    interaction_cycle_density = float((tri / 3.0) / (wedges + EPS))

    # Backdoor-density proxy: greedy variable set needed to touch 80% of clauses.
    uncovered = set(range(m))
    chosen = 0
    while uncovered and (m - len(uncovered)) / max(1, m) < 0.80:
        best_v = max(range(n_vars), key=lambda v: len(var_to_clauses[v] & uncovered))
        gain = len(var_to_clauses[best_v] & uncovered)
        if gain == 0:
            break
        uncovered -= var_to_clauses[best_v]
        chosen += 1
    backdoor_density = float(chosen / n_vars)

    # Local covariance conditioning over variable neighborhoods.
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
            ]
        )
    local_arr = np.asarray(local, dtype=float)
    cov = np.cov(local_arr.T) + 1.0e-6 * np.eye(local_arr.shape[1])
    eig = np.linalg.eigvalsh(cov)
    covariance_condition_log = float(np.log10((float(np.max(eig)) + EPS) / (float(np.min(eig)) + EPS)))

    return {
        "degree_cv": degree_cv,
        "max_degree_norm": max_degree_norm,
        "literal_imbalance": literal_imbalance,
        "expansion_defect": expansion_defect,
        "interaction_cycle_density": interaction_cycle_density,
        "backdoor_density": backdoor_density,
        "covariance_condition_log": covariance_condition_log,
        "deg": deg,
    }


def propagation_features(clauses: list[tuple[int, int, int]], n_vars: int) -> dict[str, float]:
    forced_counts = []
    contradictions = 0
    probes = 0
    for var in range(1, n_vars + 1):
        for val in (True, False):
            status, _, _, forced = unit_propagate([tuple(c) for c in clauses], {var: val})
            probes += 1
            forced_counts.append(forced / max(1, n_vars))
            if status is False:
                contradictions += 1
    return {
        "propagation_depth": float(np.mean(forced_counts)),
        "contradiction_rate": float(contradictions / max(1, probes)),
    }


def support_from_defects(feat: dict[str, float], alpha: float) -> dict[str, float]:
    # H: local consistency support from propagation contradictions and covariance blow-up.
    d_H = 2.0 * feat["contradiction_rate"] + 0.25 * feat["covariance_condition_log"]
    # B: literal and degree balance.
    d_B = 1.3 * feat["literal_imbalance"] + 0.75 * feat["degree_cv"]
    # S: activity should sit in a controlled propagation window, not zero or explosive.
    activity = feat["propagation_depth"]
    s_eff = math.exp(-0.5 * ((activity - 0.18) / 0.10) ** 2)
    d_S = (1.0 - s_eff) + 0.20 * feat["max_degree_norm"] / 6.0
    # V: connectivity/constraint geometry via expansion loss and short-cycle excess.
    d_V = 1.5 * feat["expansion_defect"] + 2.0 * feat["interaction_cycle_density"] + 0.8 * feat["backdoor_density"]

    H = 1.0 / (1.0 + d_H)
    B = 1.0 / (1.0 + d_B)
    S = 1.0 / (1.0 + d_S)
    V = 1.0 / (1.0 + d_V)
    R_resp = (H * B * V) ** (1.0 / 3.0)
    R_rob = min(R_resp, (H * B * S * V) ** 0.25)
    F_primary = -sum(math.log(EPS + x) for x in (H, B, S, V))
    F = F_primary - math.log(EPS + R_rob)
    density_peak = math.exp(-0.5 * ((alpha - 4.26) / 0.38) ** 2)
    return {
        "H": H,
        "B": B,
        "S": S,
        "V": V,
        "R_resp": R_resp,
        "R_rob": R_rob,
        "F_MAAT": F,
        "F_MAAT_primary": F_primary,
        "density_peak": density_peak,
    }


def generate_dataset(n_instances: int = 520) -> list[InstanceRecord]:
    rng = np.random.default_rng(SEED)
    records: list[InstanceRecord] = []
    alpha_grid = np.linspace(3.25, 5.15, 10)
    for idx in range(n_instances):
        n_vars = int(rng.integers(24, 37))
        alpha = float(alpha_grid[idx % len(alpha_grid)] + rng.normal(0.0, 0.035))
        n_clauses = int(round(alpha * n_vars))
        clauses = generate_random_3sat(rng, n_vars, n_clauses)
        t0 = time.perf_counter()
        sat, nodes, timed_out = dpll_solve(clauses, n_vars)
        dt = time.perf_counter() - t0
        gf = graph_features(clauses, n_vars, rng)
        pf = propagation_features(clauses, n_vars)
        feat = {**gf, **pf}
        supports = support_from_defects(feat, alpha)
        records.append(
            InstanceRecord(
                instance_id=idx,
                n_vars=n_vars,
                n_clauses=n_clauses,
                alpha=alpha,
                satisfiable=sat,
                dpll_nodes=nodes,
                dpll_timed_out=timed_out,
                solve_time_s=dt,
                log_nodes=math.log1p(nodes),
                degree_cv=feat["degree_cv"],
                max_degree_norm=feat["max_degree_norm"],
                literal_imbalance=feat["literal_imbalance"],
                expansion_defect=feat["expansion_defect"],
                interaction_cycle_density=feat["interaction_cycle_density"],
                backdoor_density=feat["backdoor_density"],
                propagation_depth=feat["propagation_depth"],
                contradiction_rate=feat["contradiction_rate"],
                covariance_condition_log=feat["covariance_condition_log"],
                H=supports["H"],
                B=supports["B"],
                S=supports["S"],
                V=supports["V"],
                R_resp=supports["R_resp"],
                R_rob=supports["R_rob"],
                F_MAAT=supports["F_MAAT"],
                F_MAAT_primary=supports["F_MAAT_primary"],
                density_peak=supports["density_peak"],
            )
        )
    return records


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


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if np.std(x) < EPS or np.std(y) < EPS:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return pearson(rankdata(np.asarray(x)), rankdata(np.asarray(y)))


def standardize_train_test(X_train: np.ndarray, X_test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = X_train.mean(axis=0)
    std = X_train.std(axis=0)
    std[std < 1.0e-9] = 1.0
    return (X_train - mean) / std, (X_test - mean) / std


def ridge_fit_predict(X_train: np.ndarray, y_train: np.ndarray, X_test: np.ndarray, alpha: float = 1.0) -> np.ndarray:
    X_train, X_test = standardize_train_test(X_train, X_test)
    Xb = np.column_stack([np.ones(len(X_train)), X_train])
    Xt = np.column_stack([np.ones(len(X_test)), X_test])
    reg = alpha * np.eye(Xb.shape[1])
    reg[0, 0] = 0.0
    beta = np.linalg.solve(Xb.T @ Xb + reg, Xb.T @ y_train)
    return Xt @ beta


def cv_regression(records: list[InstanceRecord], feature_names: list[str], folds: int = 5) -> dict[str, float]:
    rows = [asdict(r) for r in records]
    X = np.asarray([[row[f] for f in feature_names] for row in rows], dtype=float)
    y = np.asarray([row["log_nodes"] for row in rows], dtype=float)
    indices = np.arange(len(records))
    preds = np.zeros_like(y)
    for fold in range(folds):
        test = indices[indices % folds == fold]
        train = indices[indices % folds != fold]
        preds[test] = ridge_fit_predict(X[train], y[train], X[test], alpha=1.0)
    ss_res = float(np.sum((y - preds) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / (ss_tot + EPS)
    return {
        "cv_r2": r2,
        "pearson": pearson(preds, y),
        "spearman": spearman(preds, y),
        "rmse": float(np.sqrt(np.mean((y - preds) ** 2))),
    }


def permutation_null(records: list[InstanceRecord], feature_names: list[str], n_perm: int = 120) -> dict[str, float]:
    rng = np.random.default_rng(SEED + 1)
    rows = [asdict(r) for r in records]
    X = np.asarray([[row[f] for f in feature_names] for row in rows], dtype=float)
    y0 = np.asarray([row["log_nodes"] for row in rows], dtype=float)
    values = []
    indices = np.arange(len(records))
    folds = 5
    for _ in range(n_perm):
        y = rng.permutation(y0)
        preds = np.zeros_like(y)
        for fold in range(folds):
            test = indices[indices % folds == fold]
            train = indices[indices % folds != fold]
            preds[test] = ridge_fit_predict(X[train], y[train], X[test], alpha=1.0)
        values.append(spearman(preds, y))
    return {
        "null_spearman_mean": float(np.mean(values)),
        "null_spearman_std": float(np.std(values)),
        "null_spearman_p95": float(np.percentile(values, 95)),
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_results(records: list[InstanceRecord], model_rows: list[dict]) -> None:
    rows = [asdict(r) for r in records]
    log_nodes = np.array([r["log_nodes"] for r in rows])
    fmaat = np.array([r["F_MAAT"] for r in rows])
    alpha = np.array([r["alpha"] for r in rows])
    cond = np.array([r["covariance_condition_log"] for r in rows])
    backdoor = np.array([r["backdoor_density"] for r in rows])

    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    sc = ax.scatter(fmaat, log_nodes, c=alpha, s=24, cmap="viridis", alpha=0.8)
    ax.set_xlabel(r"$F_{\mathrm{MAAT}}$ structural hardness")
    ax.set_ylabel(r"$\log(1+\mathrm{DPLL\ nodes})$")
    ax.set_title("SAT structural hardness vs DPLL node count")
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(r"Clause density $\alpha=m/n$")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig1_fmaat_vs_dpll_nodes.png", dpi=260)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    sc = ax.scatter(cond, log_nodes, c=backdoor, s=24, cmap="magma", alpha=0.8)
    ax.set_xlabel("local covariance condition log")
    ax.set_ylabel(r"$\log(1+\mathrm{DPLL\ nodes})$")
    ax.set_title("Graph-local conditioning and structural hardness")
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("backdoor-density proxy")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig2_conditioning_vs_nodes.png", dpi=260)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    names = [r["model"] for r in model_rows]
    r2 = [r["cv_r2"] for r in model_rows]
    colors = ["#1a9850" if "MAAT" in n else "#4575b4" if "graph" in n else "#d73027" for n in names]
    ax.bar(names, r2, color=colors, edgecolor="#222222")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel(r"5-fold CV $R^2$ on $\log(1+\mathrm{nodes})$")
    ax.set_title("SAT hardness predictors")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig3_model_comparison.png", dpi=260)
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    records = generate_dataset()
    rows = [asdict(r) for r in records]
    write_csv(OUTDIR / "paper50_sat_instances.csv", rows)

    y = np.array([r.log_nodes for r in records])
    scalar_scores = {
        "F_MAAT": np.array([r.F_MAAT for r in records]),
        "density_peak": np.array([r.density_peak for r in records]),
        "R_rob_inverse": 1.0 - np.array([r.R_rob for r in records]),
        "covariance_condition_log": np.array([r.covariance_condition_log for r in records]),
        "backdoor_density": np.array([r.backdoor_density for r in records]),
    }
    scalar_rows = []
    for name, score in scalar_scores.items():
        scalar_rows.append(
            {
                "score": name,
                "pearson": pearson(score, y),
                "spearman": spearman(score, y),
            }
        )
    write_csv(OUTDIR / "paper50_scalar_correlations.csv", scalar_rows)

    models = {
        "density_only": ["alpha", "density_peak"],
        "standard_graph": ["alpha", "degree_cv", "max_degree_norm", "literal_imbalance", "interaction_cycle_density"],
        "stability_only": ["R_rob"],
        "MAAT_scalar": ["F_MAAT"],
        "MAAT_graph_local": [
            "F_MAAT",
            "H",
            "B",
            "S",
            "V",
            "R_rob",
            "expansion_defect",
            "degree_cv",
            "max_degree_norm",
            "literal_imbalance",
            "interaction_cycle_density",
            "backdoor_density",
            "contradiction_rate",
            "propagation_depth",
            "covariance_condition_log",
        ],
        "MAAT_graph_local_plus_density": [
            "alpha",
            "density_peak",
            "F_MAAT",
            "H",
            "B",
            "S",
            "V",
            "R_rob",
            "expansion_defect",
            "degree_cv",
            "max_degree_norm",
            "literal_imbalance",
            "interaction_cycle_density",
            "backdoor_density",
            "contradiction_rate",
            "propagation_depth",
            "covariance_condition_log",
        ],
    }
    model_rows = []
    for name, features in models.items():
        metrics = cv_regression(records, features)
        model_rows.append({"model": name, "features": "+".join(features), **metrics})
    write_csv(OUTDIR / "paper50_model_comparison.csv", model_rows)

    null = permutation_null(records, models["MAAT_graph_local"], n_perm=120)
    maat_row = next(r for r in model_rows if r["model"] == "MAAT_graph_local")
    maat_plus_row = next(r for r in model_rows if r["model"] == "MAAT_graph_local_plus_density")
    summary = {
        "model": "Paper 50 SAT Structural Hardness II",
        "status": "Synthetic random-3-SAT graph-local benchmark; not a proof of P vs NP and not a solver competition.",
        "random_seed": SEED,
        "n_instances": len(records),
        "n_vars_range": [min(r.n_vars for r in records), max(r.n_vars for r in records)],
        "alpha_range": [float(min(r.alpha for r in records)), float(max(r.alpha for r in records))],
        "timeout_fraction": float(np.mean([r.dpll_timed_out for r in records])),
        "sat_fraction": float(np.mean([r.satisfiable for r in records])),
        "scalar_correlations": scalar_rows,
        "model_comparison": model_rows,
        "permutation_null": null,
        "maat_graph_local_minus_density_r2": float(maat_row["cv_r2"] - next(r for r in model_rows if r["model"] == "density_only")["cv_r2"]),
        "maat_graph_local_minus_standard_graph_r2": float(maat_row["cv_r2"] - next(r for r in model_rows if r["model"] == "standard_graph")["cv_r2"]),
        "maat_graph_local_plus_density_minus_density_r2": float(maat_plus_row["cv_r2"] - next(r for r in model_rows if r["model"] == "density_only")["cv_r2"]),
        "maat_graph_local_plus_density_minus_standard_graph_r2": float(maat_plus_row["cv_r2"] - next(r for r in model_rows if r["model"] == "standard_graph")["cv_r2"]),
        "maat_graph_local_spearman_above_null_p95": bool(maat_row["spearman"] > null["null_spearman_p95"]),
        "key_interpretation": (
            "Graph-local MAAT defects are evaluated against density-only, standard graph, "
            "stability-only, scalar MAAT, and shuffled-null baselines."
        ),
        "outputs": [
            "paper50_sat_instances.csv",
            "paper50_scalar_correlations.csv",
            "paper50_model_comparison.csv",
            "paper50_summary.json",
            "fig1_fmaat_vs_dpll_nodes.png",
            "fig2_conditioning_vs_nodes.png",
            "fig3_model_comparison.png",
        ],
    }
    with (OUTDIR / "paper50_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    plot_results(records, model_rows)

    print("Paper 50 SAT structural hardness benchmark complete.")
    print(f"Instances: {len(records)}")
    print(f"Timeout fraction: {summary['timeout_fraction']:.4f}")
    print(f"SAT fraction: {summary['sat_fraction']:.4f}")
    for row in model_rows:
        print(f"{row['model']:18s} CV R2={row['cv_r2']:.4f} Spearman={row['spearman']:.4f}")
    print(f"Null p95 Spearman: {null['null_spearman_p95']:.4f}")


if __name__ == "__main__":
    main()
