#!/usr/bin/env python3
"""
MAAT v1.6 structural search geometry toy benchmark.

This script tests whether a structural search score can navigate simple
synthetic search spaces with traps and bridges differently from simple
distance/energy/random baselines. It is a diagnostic toy benchmark, not a
claim of algorithmic optimality.
"""

from __future__ import annotations

import heapq
import json
import math
import os
import random
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

_CACHE_ROOT = Path(os.environ.get("TMPDIR", "/tmp")) / "maat_v16_search_cache"
(_CACHE_ROOT / "matplotlib").mkdir(parents=True, exist_ok=True)
(_CACHE_ROOT / "xdg").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE_ROOT / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE_ROOT / "xdg"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


EPS = 1e-9
GRID_N = 34
N_INSTANCES = 14
SEARCH_BUDGET = 850
RANDOM_REPEATS = 8
OUTDIR = Path(__file__).resolve().parent / "outputs_v16"


@dataclass(frozen=True)
class SearchInstance:
    family: str
    seed: int
    graph: nx.Graph
    start: tuple[int, int]
    target: tuple[int, int]
    energy: dict[tuple[int, int], float]
    base_distance: dict[tuple[int, int], float]
    supports: dict[str, dict[tuple[int, int], float]]
    score_maat_anchored: dict[tuple[int, int], float]
    score_maat_heavy: dict[tuple[int, int], float]


def normalise(values: dict[tuple[int, int], float], invert: bool = False) -> dict[tuple[int, int], float]:
    arr = np.array(list(values.values()), dtype=float)
    lo, hi = float(np.min(arr)), float(np.max(arr))
    if hi - lo < 1e-12:
        return {k: 1.0 for k in values}
    out = {k: (v - lo) / (hi - lo) for k, v in values.items()}
    if invert:
        out = {k: 1.0 - v for k, v in out.items()}
    return out


def build_graph(family: str, seed: int) -> tuple[nx.Graph, tuple[int, int], tuple[int, int], tuple[int, int]]:
    rng = random.Random(seed)
    graph = nx.grid_2d_graph(GRID_N, GRID_N)
    start = (2, 2)
    target = (GRID_N - 3, GRID_N - 3)
    decoy = (rng.randint(7, 12), rng.randint(GRID_N - 13, GRID_N - 7))

    if family in {"bridge", "trap_bridge"}:
        wall_x = GRID_N // 2 + rng.choice([-2, -1, 0, 1, 2])
        # Put the gate away from the direct Euclidean route. This makes the
        # benchmark a search problem rather than a nearly straight-line task.
        gate_y = rng.choice([5, 7, GRID_N - 8, GRID_N - 6])
        gates = {gate_y, max(1, gate_y - 1), min(GRID_N - 2, gate_y + 1)}
        obstacles = [(wall_x, y) for y in range(1, GRID_N - 1) if y not in gates]
        graph.remove_nodes_from(obstacles)

        # Add a small cul-de-sac on the near side to make pure energy attractive
        # but structurally fragile.
        if family == "trap_bridge":
            block_y = min(GRID_N - 4, gate_y + 5)
            graph.remove_nodes_from([(wall_x - 3, y) for y in range(block_y, min(GRID_N - 1, block_y + 4))])

    return graph, start, target, decoy


def compute_energy(
    graph: nx.Graph,
    family: str,
    target: tuple[int, int],
    decoy: tuple[int, int],
    seed: int,
) -> tuple[dict[tuple[int, int], float], dict[tuple[int, int], float]]:
    rng = np.random.default_rng(seed)
    max_dist = math.sqrt(2) * (GRID_N - 1)
    base_distance = {}
    raw_energy = {}

    for x, y in graph.nodes:
        dist_target = math.hypot(target[0] - x, target[1] - y) / max_dist
        dist_decoy = math.hypot(decoy[0] - x, decoy[1] - y)
        trap = 0.0
        if family in {"trap", "trap_bridge"}:
            trap = -0.55 * math.exp(-(dist_decoy**2) / (2 * 4.5**2))
        ripple = 0.025 * math.sin(0.31 * x + 0.17 * y + 0.13 * seed)
        rough = 0.01 * float(rng.normal())
        base_distance[(x, y)] = dist_target
        raw_energy[(x, y)] = dist_target + trap + ripple + rough

    energy_norm = normalise(raw_energy)
    return energy_norm, base_distance


def compute_supports(
    graph: nx.Graph,
    energy: dict[tuple[int, int], float],
    base_distance: dict[tuple[int, int], float],
    target: tuple[int, int],
    seed: int,
) -> dict[str, dict[tuple[int, int], float]]:
    # Approximate bridge/connectivity centrality. This is intentionally a
    # graph-local structural cue, not a full shortest-path oracle.
    k = min(80, graph.number_of_nodes())
    centrality = nx.betweenness_centrality(graph, k=k, seed=seed, normalized=True)
    centrality = normalise(centrality)
    degree = {n: graph.degree[n] / 4.0 for n in graph.nodes}

    fragility = {}
    exit_balance = {}
    activity = {}
    reachability = {}

    for node in graph.nodes:
        neigh = list(graph.neighbors(node))
        if not neigh:
            fragility[node] = 1.0
            exit_balance[node] = 0.0
            activity[node] = 0.0
            reachability[node] = 0.0
            continue

        e0 = energy[node]
        neigh_e = np.array([energy[v] for v in neigh], dtype=float)
        improvements = np.maximum(0.0, e0 - neigh_e)
        worsening = np.maximum(0.0, neigh_e - e0)

        # Positive curvature catches attractive local wells: low local energy
        # surrounded by worse neighbours.
        curvature = max(0.0, float(np.mean(neigh_e) - e0))
        fragility[node] = curvature + 0.35 * float(np.std(neigh_e))

        improve_fraction = float(np.mean(neigh_e < e0))
        exit_balance[node] = 1.0 - abs(2.0 * improve_fraction - 1.0)

        # Activity should be useful novelty, not maximum jumpiness.
        improvement_scale = float(np.mean(improvements) - 0.4 * np.mean(worsening))
        activity[node] = math.exp(-0.5 * ((improvement_scale - 0.035) / 0.055) ** 2)

        # Local reachability plus approximate bridge centrality.
        reachability[node] = 0.45 * degree[node] + 0.55 * centrality[node]

    fragility_norm = normalise(fragility)
    reachability_norm = normalise(reachability)

    h_support = {n: 1.0 / (1.0 + 2.5 * fragility_norm[n]) for n in graph.nodes}
    b_support = {n: 0.15 + 0.85 * max(0.0, min(1.0, exit_balance[n])) for n in graph.nodes}
    s_support = {n: 0.10 + 0.90 * max(0.0, min(1.0, activity[n])) for n in graph.nodes}
    v_support = {n: 0.10 + 0.90 * max(0.0, min(1.0, reachability_norm[n])) for n in graph.nodes}

    r_resp = {n: (h_support[n] * b_support[n] * v_support[n]) ** (1.0 / 3.0) for n in graph.nodes}
    r_rob = {
        n: min(r_resp[n], (h_support[n] * b_support[n] * s_support[n] * v_support[n]) ** 0.25)
        for n in graph.nodes
    }

    return {
        "H": h_support,
        "B": b_support,
        "S": s_support,
        "V": v_support,
        "Rresp": r_resp,
        "Rrob": r_rob,
    }


def build_instance(family: str, seed: int) -> SearchInstance:
    graph, start, target, decoy = build_graph(family, seed)
    energy, base_distance = compute_energy(graph, family, target, decoy, seed)
    supports = compute_supports(graph, energy, base_distance, target, seed)

    score_maat_anchored = {}
    score_maat_heavy = {}
    for node in graph.nodes:
        h = supports["H"][node]
        b = supports["B"][node]
        s = supports["S"][node]
        v = supports["V"][node]
        r = supports["Rrob"][node]
        support_cost = -sum(math.log(EPS + q) for q in (h, b, s, v, r))
        verification_distance = base_distance[node]
        fragility_penalty = 1.0 - r
        # v1.6 search score: verification-dominant structural guidance.
        # The target anchor must remain primary; otherwise structural support can
        # trap the search in coherent but irrelevant regions.
        score_maat_heavy[node] = (
            0.70 * verification_distance
            + 0.18 * energy[node]
            + 0.06 * support_cost
            + 0.08 * fragility_penalty
            - 0.16 * v
        )
        score_maat_anchored[node] = (
            0.88 * verification_distance
            + 0.06 * energy[node]
            + 0.025 * support_cost
            + 0.030 * fragility_penalty
            - 0.040 * v
        )

    return SearchInstance(
        family,
        seed,
        graph,
        start,
        target,
        energy,
        base_distance,
        supports,
        score_maat_anchored,
        score_maat_heavy,
    )


def best_first_search(
    inst: SearchInstance,
    priority: dict[tuple[int, int], float] | None,
    algorithm: str,
    seed: int,
    budget: int = SEARCH_BUDGET,
) -> tuple[bool, int, list[tuple[int, int]]]:
    rng = random.Random(seed)
    visited = {inst.start}
    parent: dict[tuple[int, int], tuple[int, int] | None] = {inst.start: None}
    frontier: list[tuple[float, float, tuple[int, int]]] = []

    def push(node: tuple[int, int], jitter: float = 0.0) -> None:
        if algorithm == "random_frontier":
            key = rng.random()
        else:
            assert priority is not None
            key = priority[node] + jitter
        heapq.heappush(frontier, (key, rng.random(), node))

    push(inst.start)
    expansions = 0

    while frontier and expansions < budget:
        _, _, node = heapq.heappop(frontier)
        expansions += 1
        if node == inst.target:
            path = []
            cur: tuple[int, int] | None = node
            while cur is not None:
                path.append(cur)
                cur = parent[cur]
            return True, expansions, path[::-1]

        for nb in inst.graph.neighbors(node):
            if nb in visited:
                continue
            visited.add(nb)
            parent[nb] = node
            push(nb, jitter=1e-6 * rng.random())

    return False, expansions, []


def run_benchmark() -> tuple[pd.DataFrame, pd.DataFrame]:
    families = ["smooth", "trap", "bridge", "trap_bridge"]
    algorithms = {
        "random_frontier": None,
        "distance_only": "base_distance",
        "energy_only": "energy",
        "connectivity_only": "connectivity",
        "maat_v16_anchored": "maat_anchored",
        "maat_v16_heavy": "maat_heavy",
    }

    rows = []
    example_instances: dict[str, SearchInstance] = {}

    for family in families:
        for i in range(N_INSTANCES):
            seed = 1000 + 97 * i + 13 * families.index(family)
            inst = build_instance(family, seed)
            if i == 0:
                example_instances[family] = inst

            priority_maps = {
                "base_distance": inst.base_distance,
                "energy": inst.energy,
                "connectivity": {n: inst.base_distance[n] - 0.28 * inst.supports["V"][n] for n in inst.graph.nodes},
                "maat_anchored": inst.score_maat_anchored,
                "maat_heavy": inst.score_maat_heavy,
            }

            for algorithm, map_name in algorithms.items():
                repeats = RANDOM_REPEATS if algorithm == "random_frontier" else 1
                for rep in range(repeats):
                    priority = None if map_name is None else priority_maps[map_name]
                    found, expansions, path = best_first_search(inst, priority, algorithm, seed + 10000 * rep)
                    rows.append(
                        {
                            "family": family,
                            "instance": i,
                            "seed": seed,
                            "algorithm": algorithm,
                            "repeat": rep,
                            "found": int(found),
                            "expansions": expansions,
                            "path_length": len(path) if found else np.nan,
                        }
                    )

    results = pd.DataFrame(rows)
    summary = (
        results.groupby(["family", "algorithm"], as_index=False)
        .agg(
            success_rate=("found", "mean"),
            median_expansions=("expansions", "median"),
            mean_expansions=("expansions", "mean"),
            n_runs=("found", "size"),
        )
        .sort_values(["family", "median_expansions"])
    )

    plot_results(results, summary, example_instances)
    return results, summary


def plot_results(
    results: pd.DataFrame,
    summary: pd.DataFrame,
    examples: dict[str, SearchInstance],
) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    order = [
        "random_frontier",
        "distance_only",
        "energy_only",
        "connectivity_only",
        "maat_v16_anchored",
        "maat_v16_heavy",
    ]
    colors = {
        "random_frontier": "#8c8c8c",
        "distance_only": "#2b6cb0",
        "energy_only": "#dd6b20",
        "connectivity_only": "#805ad5",
        "maat_v16_anchored": "#00897b",
        "maat_v16_heavy": "#004d40",
    }

    fig, ax = plt.subplots(figsize=(10, 5.5))
    pivot = summary.pivot(index="family", columns="algorithm", values="success_rate").reindex(columns=order)
    pivot.plot(kind="bar", ax=ax, color=[colors[c] for c in pivot.columns], width=0.82)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("success rate under fixed budget")
    ax.set_xlabel("")
    ax.set_title("MAAT v1.6 structural search benchmark: success")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig1_success_rate.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5.5))
    pivot = summary.pivot(index="family", columns="algorithm", values="median_expansions").reindex(columns=order)
    pivot.plot(kind="bar", ax=ax, color=[colors[c] for c in pivot.columns], width=0.82)
    ax.set_ylabel("median node expansions")
    ax.set_xlabel("")
    ax.set_title("Search cost by environment family")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig2_median_expansions.png", dpi=180)
    plt.close(fig)

    inst = examples["trap_bridge"]
    priority_maps = {
        "distance_only": inst.base_distance,
        "energy_only": inst.energy,
        "maat_v16_anchored": inst.score_maat_anchored,
        "maat_v16_heavy": inst.score_maat_heavy,
    }
    paths = {}
    for alg, pmap in priority_maps.items():
        _, _, path = best_first_search(inst, pmap, alg, inst.seed)
        paths[alg] = path

    fig, ax = plt.subplots(figsize=(7, 7))
    grid = np.full((GRID_N, GRID_N), np.nan)
    for node, val in inst.energy.items():
        grid[node[1], node[0]] = val
    ax.imshow(grid, origin="lower", cmap="magma", alpha=0.9)
    for alg, path in paths.items():
        if not path:
            continue
        xs = [p[0] for p in path]
        ys = [p[1] for p in path]
        ax.plot(xs, ys, lw=2.1, label=alg, alpha=0.9)
    ax.scatter([inst.start[0]], [inst.start[1]], c="white", edgecolors="black", s=90, label="start", zorder=5)
    ax.scatter([inst.target[0]], [inst.target[1]], c="#00e676", edgecolors="black", s=90, label="target", zorder=5)
    ax.set_title("Example trap+bridge landscape and recovered paths")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig3_example_paths.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    for ax, key, title in zip(axes, ["H", "V", "Rrob"], ["H coherence", "V connectivity", "Rrob closure"]):
        grid = np.full((GRID_N, GRID_N), np.nan)
        for node, val in inst.supports[key].items():
            grid[node[1], node[0]] = val
        im = ax.imshow(grid, origin="lower", cmap="viridis", vmin=0, vmax=1)
        ax.scatter([inst.start[0]], [inst.start[1]], c="white", edgecolors="black", s=35)
        ax.scatter([inst.target[0]], [inst.target[1]], c="red", edgecolors="black", s=35)
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.suptitle("Structural supports on the same search space", y=1.03)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig4_support_maps.png", dpi=180)
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    results, summary = run_benchmark()

    results.to_csv(OUTDIR / "v16_search_results.csv", index=False)
    summary.to_csv(OUTDIR / "v16_search_summary.csv", index=False)

    aggregate = (
        results.groupby("algorithm", as_index=False)
        .agg(
            success_rate=("found", "mean"),
            median_expansions=("expansions", "median"),
            mean_expansions=("expansions", "mean"),
            n_runs=("found", "size"),
        )
        .sort_values("median_expansions")
    )

    best_by_family = {}
    for family, sub in summary.groupby("family"):
        best = sub.sort_values(["success_rate", "median_expansions"], ascending=[False, True]).iloc[0]
        best_by_family[family] = {
            "algorithm": str(best["algorithm"]),
            "success_rate": float(best["success_rate"]),
            "median_expansions": float(best["median_expansions"]),
        }

    payload = {
        "benchmark": "MAAT v1.6 structural search geometry toy benchmark",
        "grid_n": GRID_N,
        "instances_per_family": N_INSTANCES,
        "budget": SEARCH_BUDGET,
        "families": sorted(results["family"].unique().tolist()),
        "algorithms": sorted(results["algorithm"].unique().tolist()),
        "aggregate": aggregate.to_dict(orient="records"),
        "best_by_family": best_by_family,
        "scientific_status": (
            "Toy benchmark. It tests search-navigation behaviour on synthetic "
            "trap and bridge landscapes; it is not an optimality proof."
        ),
    }
    with open(OUTDIR / "v16_summary.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
