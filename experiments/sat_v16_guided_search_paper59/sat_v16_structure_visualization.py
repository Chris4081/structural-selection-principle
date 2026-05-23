#!/usr/bin/env python3
"""
SAT structure visualisation for Paper 59.

The goal is to make a SAT instance visually inspectable in the same spirit as
a dense mathematical graph plot: variables, clauses, connectivity, and
MAAT-v1.6 branching pressure are shown as graph geometry.

This is a visual diagnostic, not a solver result.
"""

from __future__ import annotations

import os
from itertools import combinations
from pathlib import Path

_CACHE_ROOT = Path(os.environ.get("TMPDIR", "/tmp")) / "maat_paper59_visual_cache"
(_CACHE_ROOT / "matplotlib").mkdir(parents=True, exist_ok=True)
(_CACHE_ROOT / "xdg").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE_ROOT / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE_ROOT / "xdg"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from sat_v16_guided_search_paper59 import (
    SEED,
    active_clauses,
    build_instances,
    generate_modular_3sat,
    score_candidates,
)


OUTDIR = Path(__file__).resolve().parent / "outputs_paper59" / "sat_graph_visuals"
OUTDIR.mkdir(parents=True, exist_ok=True)


def choose_visual_instance():
    rng = np.random.default_rng(SEED)
    # Use a larger modular instance for visual density. This is not part of the
    # Paper-59 solver table; it is a display instance chosen to make the SAT
    # geometry visually comparable to dense mathematical graph drawings.
    n_vars = 72
    clauses = generate_modular_3sat(rng, n_vars=n_vars, n_clauses=int(round(4.25 * n_vars)))
    return "dense_modular_3sat_visual", n_vars, clauses


def variable_scores(n_vars: int, clauses: list[tuple[int, ...]]) -> dict[int, dict[str, float]]:
    assignment: dict[int, bool] = {}
    active = active_clauses(clauses, assignment)
    return score_candidates(active, list(range(1, n_vars + 1)))


def normalise(vals: list[float]) -> list[float]:
    arr = np.asarray(vals, dtype=float)
    lo, hi = float(arr.min()), float(arr.max())
    if hi - lo < 1e-12:
        return [0.5 for _ in vals]
    return ((arr - lo) / (hi - lo)).tolist()


def plot_variable_cooccurrence(family: str, n_vars: int, clauses: list[tuple[int, ...]]) -> None:
    scores = variable_scores(n_vars, clauses)
    graph = nx.Graph()
    for v in range(1, n_vars + 1):
        graph.add_node(v)
    weights: dict[tuple[int, int], int] = {}
    for clause in clauses:
        vars_ = sorted({abs(lit) for lit in clause})
        for u, v in combinations(vars_, 2):
            key = (u, v)
            weights[key] = weights.get(key, 0) + 1
    for (u, v), w in weights.items():
        graph.add_edge(u, v, weight=w)

    pos = nx.spring_layout(graph, seed=SEED, k=0.55, iterations=250, weight="weight")
    maat = [scores[v]["maat"] for v in graph.nodes]
    sizes = [90 + 430 * scores[v]["V"] for v in graph.nodes]
    edge_widths = [0.35 + 0.65 * graph.edges[e]["weight"] for e in graph.edges]

    fig, ax = plt.subplots(figsize=(10, 9))
    nx.draw_networkx_edges(graph, pos, ax=ax, width=edge_widths, alpha=0.23, edge_color="#1e40ff")
    nodes = nx.draw_networkx_nodes(
        graph,
        pos,
        ax=ax,
        node_size=sizes,
        node_color=maat,
        cmap="plasma",
        linewidths=0.8,
        edgecolors="#111111",
        alpha=0.95,
    )
    nx.draw_networkx_labels(graph, pos, ax=ax, font_size=7, font_color="#111111")
    cbar = fig.colorbar(nodes, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("MAAT v1.6 branching priority")
    ax.set_title(
        f"SAT variable co-occurrence graph ({family}, n={n_vars}, clauses={len(clauses)})\n"
        "node colour = MAAT branching priority, node size = connectivity support V"
    )
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(OUTDIR / "sat_variable_cooccurrence_maat_priority.png", dpi=220)
    plt.close(fig)


def plot_clause_variable_bipartite(family: str, n_vars: int, clauses: list[tuple[int, ...]]) -> None:
    scores = variable_scores(n_vars, clauses)
    graph = nx.Graph()
    var_nodes = [f"x{v}" for v in range(1, n_vars + 1)]
    clause_nodes = [f"C{i}" for i in range(len(clauses))]
    graph.add_nodes_from(var_nodes, kind="variable")
    graph.add_nodes_from(clause_nodes, kind="clause")
    edge_signs = []
    for i, clause in enumerate(clauses):
        cnode = f"C{i}"
        for lit in clause:
            vnode = f"x{abs(lit)}"
            graph.add_edge(vnode, cnode, sign=1 if lit > 0 else -1)
            edge_signs.append(1 if lit > 0 else -1)

    pos = nx.spring_layout(graph, seed=SEED + 1, k=0.22, iterations=350)
    var_scores = [scores[int(node[1:])]["maat"] for node in var_nodes]
    var_sizes = [35 + 230 * scores[int(node[1:])]["V"] for node in var_nodes]
    clause_lengths = [len(c) for c in clauses]
    clause_sizes = [18 + 14 * length for length in clause_lengths]

    pos_edges = [(u, v) for u, v, d in graph.edges(data=True) if d["sign"] > 0]
    neg_edges = [(u, v) for u, v, d in graph.edges(data=True) if d["sign"] < 0]

    fig, ax = plt.subplots(figsize=(12, 10))
    nx.draw_networkx_edges(graph, pos, edgelist=pos_edges, ax=ax, width=0.55, alpha=0.18, edge_color="#2563eb")
    nx.draw_networkx_edges(graph, pos, edgelist=neg_edges, ax=ax, width=0.55, alpha=0.18, edge_color="#ef4444")
    nx.draw_networkx_nodes(
        graph,
        pos,
        nodelist=clause_nodes,
        node_size=clause_sizes,
        node_color="#e5e7eb",
        edgecolors="#64748b",
        linewidths=0.4,
        alpha=0.60,
        ax=ax,
    )
    nodes = nx.draw_networkx_nodes(
        graph,
        pos,
        nodelist=var_nodes,
        node_size=var_sizes,
        node_color=var_scores,
        cmap="plasma",
        edgecolors="#111111",
        linewidths=0.5,
        alpha=0.96,
        ax=ax,
    )
    cbar = fig.colorbar(nodes, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("MAAT v1.6 branching priority")
    ax.set_title(
        f"SAT clause-variable graph ({family}, n={n_vars}, clauses={len(clauses)})\n"
        "blue edges = positive literal, red edges = negative literal"
    )
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(OUTDIR / "sat_clause_variable_bipartite_maat_priority.png", dpi=220)
    plt.close(fig)


def plot_priority_components(family: str, n_vars: int, clauses: list[tuple[int, ...]]) -> None:
    scores = variable_scores(n_vars, clauses)
    vars_ = list(range(1, n_vars + 1))
    components = {
        "degree": [scores[v]["degree"] for v in vars_],
        "JW": [scores[v]["jw"] for v in vars_],
        "MOMS": [scores[v]["moms"] for v in vars_],
        "MAAT": [scores[v]["maat"] for v in vars_],
        "V": [scores[v]["V"] for v in vars_],
    }

    fig, axes = plt.subplots(5, 1, figsize=(11, 8), sharex=True)
    for ax, (name, vals) in zip(axes, components.items()):
        vals_norm = normalise(vals)
        ax.bar(vars_, vals_norm, color="#0f766e" if name == "MAAT" else "#64748b")
        ax.set_ylabel(name)
        ax.set_ylim(0, 1.05)
    axes[-1].set_xlabel("variable")
    fig.suptitle(f"Root-level branching signals ({family})", y=0.995)
    fig.tight_layout()
    fig.savefig(OUTDIR / "sat_branching_signal_components.png", dpi=220)
    plt.close(fig)


def main() -> None:
    family, n_vars, clauses = choose_visual_instance()
    plot_variable_cooccurrence(family, n_vars, clauses)
    plot_clause_variable_bipartite(family, n_vars, clauses)
    plot_priority_components(family, n_vars, clauses)
    print(f"Wrote SAT visualisations to {OUTDIR}")


if __name__ == "__main__":
    main()
