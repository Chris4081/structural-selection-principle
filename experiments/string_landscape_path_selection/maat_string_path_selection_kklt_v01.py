#!/usr/bin/env python3
"""
MAAT String Selection II — Structural Path Selection in a Reduced KKLT Landscape
================================================================================

This script turns the existing IIB/KKLT bridge vacuum-ranking output into a
directed transition graph.  Paper I ranked candidate vacua as endpoints.  This
benchmark asks a different question:

    which transitions between candidate vacua are structurally preferred?

The graph is deliberately phenomenological:

* nodes are KKLT bridge candidates classified as minima;
* edges connect near neighbours in reduced flux/moduli feature space;
* edge weights combine a tunnelling/barrier proxy with a MAAT structural
  uphill penalty, Tadpole/control penalties, and stability-loss penalty.

This is not a Coleman-De Luccia calculation, not a string decay-rate
derivation, and not a complete landscape dynamics model.  It is a reproducible
path-selection benchmark that upgrades static vacuum ranking into a structural
landscape-flow diagnostic.

Run:
    python3 maat_string_path_selection_kklt_v01.py
"""

from __future__ import annotations

import csv
import heapq
import json
import math
import os
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "maat_string_path_selection_mpl"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


BASE_DIR = Path(__file__).resolve().parent
OUT = BASE_DIR / "outputs_path_selection"
OUT.mkdir(exist_ok=True)

EPS = 1e-12


@dataclass(frozen=True)
class PathConfig:
    source_results: str = "../string_landscape_selection/structural_selection_iib_kklt_bridge_results.json"
    k_neighbours: int = 6
    top_path_sources: int = 8
    lambda_structural_uphill: float = 0.85
    lambda_control: float = 1.10
    lambda_tadpole: float = 0.65
    lambda_stability_loss: float = 0.75
    barrier_flux_weight: float = 0.45
    barrier_tau_weight: float = 0.20
    barrier_w0_weight: float = 0.15
    barrier_energy_weight: float = 0.10
    barrier_tadpole_step_weight: float = 0.10


CONFIG = PathConfig()


def resolve_source(path_text: str) -> Path:
    candidates = [
        BASE_DIR / path_text,
        Path.cwd() / path_text,
        Path("/Volumes/MAATSSD/structural-selection-principle/experiments/string_landscape_selection/structural_selection_iib_kklt_bridge_results.json"),
    ]
    for path in candidates:
        if path.exists():
            return path.resolve()
    raise FileNotFoundError(f"Could not find source results JSON. Tried: {candidates}")


def rankdata(values: list[float]) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    order = np.argsort(arr, kind="mergesort")
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(arr) + 1, dtype=float)
    return ranks


def spearman(x: list[float], y: list[float]) -> float:
    xr = rankdata(x)
    yr = rankdata(y)
    if np.std(xr) == 0.0 or np.std(yr) == 0.0:
        return float("nan")
    return float(np.corrcoef(xr, yr)[0, 1])


def support(defect: float) -> float:
    return 1.0 / (1.0 + max(0.0, float(defect)))


def load_minima(path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    with open(path, "r", encoding="utf-8") as f:
        payload = json.load(f)

    nodes: list[dict[str, Any]] = []
    for idx, candidate in enumerate(payload["candidates"]):
        if "minimum" not in str(candidate["vacuum_type"]):
            continue
        d = candidate["bridge_defects"]
        scores = candidate["kklt_scores"]
        flux_vec = np.array(candidate["F_flux"] + candidate["H_flux"], dtype=float)
        node = {
            "node_id": len(nodes),
            "candidate_id": idx,
            "background_id": int(candidate["background_id"]),
            "vacuum_type": candidate["vacuum_type"],
            "F_flux": candidate["F_flux"],
            "H_flux": candidate["H_flux"],
            "flux_vec": flux_vec,
            "N_flux": float(candidate["N_flux"]),
            "N_D3": float(candidate["N_D3"]),
            "Delta_Q_D3": float(candidate["Delta_Q_D3"]),
            "abs_Delta_Q_D3": float(candidate["abs_Delta_Q_D3"]),
            "W0": float(candidate["W0"]),
            "D": float(candidate["D"]),
            "tau": float(candidate["tau"]),
            "energy": float(candidate["energy"]),
            "abs_energy": float(candidate["abs_energy"]),
            "lambda_min": float(candidate["lambda_min"]),
            "F_bridge": float(candidate["F_bridge_10"]),
            "F_pheno": float(candidate["F_pheno_bridge"]),
            "kklt_f_master": float(candidate["kklt_f_master"]),
            "rank_energy_original": float(candidate["rank_energy"]),
            "rank_bridge_original": float(candidate["rank_bridge"]),
            "rank_pheno_original": float(candidate["rank_pheno"]),
            "Delta_Q_10": float(d["Delta_Q_10"]),
            "Delta_dim": float(d["Delta_dim_KKLT"]),
            "Delta_E": float(d["Delta_E_KKLT"]),
            "tachyon_defect": float(d["tachyon_defect"]),
            "R_support": float(scores["R"]),
            "B_support": float(scores["B"]),
            "Sigma_support": float(scores["Sigma"]),
        }
        nodes.append(node)

    bridge_ranks = rankdata([n["F_bridge"] for n in nodes])
    energy_ranks = rankdata([n["abs_energy"] for n in nodes])
    pheno_ranks = rankdata([n["F_pheno"] for n in nodes])
    for node, rb, re, rp in zip(nodes, bridge_ranks, energy_ranks, pheno_ranks):
        node["rank_bridge_minima"] = float(rb)
        node["rank_energy_minima"] = float(re)
        node["rank_pheno_minima"] = float(rp)
    return payload, nodes


def feature_matrix(nodes: list[dict[str, Any]]) -> np.ndarray:
    flux = np.vstack([n["flux_vec"] for n in nodes])
    extra = np.array(
        [
            [
                n["tau"],
                math.log10(abs(n["W0"]) + EPS),
                math.log10(abs(n["D"]) + EPS),
                math.log10(n["abs_energy"] + 1.0e-20),
                n["Delta_Q_10"],
            ]
            for n in nodes
        ],
        dtype=float,
    )
    x = np.hstack([flux / 4.0, extra])
    mu = np.mean(x, axis=0)
    sig = np.std(x, axis=0)
    sig[sig < EPS] = 1.0
    return (x - mu) / sig


def normalized_abs_delta(a: float, b: float, scale: float) -> float:
    return abs(float(a) - float(b)) / max(scale, EPS)


def make_edge(i: int, j: int, nodes: list[dict[str, Any]], feature_dist: float, config: PathConfig) -> dict[str, Any]:
    src = nodes[i]
    dst = nodes[j]

    tau_scale = max(n["tau"] for n in nodes) - min(n["tau"] for n in nodes)
    w0_scale = max(abs(n["W0"]) for n in nodes) - min(abs(n["W0"]) for n in nodes)
    energy_scale = max(math.log10(n["abs_energy"] + 1e-20) for n in nodes) - min(math.log10(n["abs_energy"] + 1e-20) for n in nodes)

    flux_step = min(1.0, feature_dist / 6.0)
    tau_step = normalized_abs_delta(src["tau"], dst["tau"], tau_scale)
    w0_step = normalized_abs_delta(abs(src["W0"]), abs(dst["W0"]), w0_scale)
    energy_step = normalized_abs_delta(
        math.log10(src["abs_energy"] + 1e-20),
        math.log10(dst["abs_energy"] + 1e-20),
        energy_scale,
    )
    tadpole_step = abs(dst["Delta_Q_10"] - src["Delta_Q_10"])

    tunnel_barrier = (
        config.barrier_flux_weight * flux_step
        + config.barrier_tau_weight * tau_step
        + config.barrier_w0_weight * w0_step
        + config.barrier_energy_weight * energy_step
        + config.barrier_tadpole_step_weight * tadpole_step
    )

    delta_F = dst["F_bridge"] - src["F_bridge"]
    delta_F_pheno = dst["F_pheno"] - src["F_pheno"]
    structural_uphill = max(0.0, delta_F)
    structural_downhill = max(0.0, -delta_F)
    avg_control = 0.5 * (src["Delta_dim"] + dst["Delta_dim"])
    avg_tadpole = 0.5 * (src["Delta_Q_10"] + dst["Delta_Q_10"])
    stability_loss = max(0.0, src["R_support"] - dst["R_support"])

    structural_action = (
        tunnel_barrier
        + config.lambda_structural_uphill * structural_uphill
        + config.lambda_control * avg_control
        + config.lambda_tadpole * avg_tadpole
        + config.lambda_stability_loss * stability_loss
    )
    energy_action = tunnel_barrier + 0.85 * max(0.0, dst["abs_energy"] - src["abs_energy"]) / max(src["abs_energy"], dst["abs_energy"], 1e-20)

    gamma_proxy = math.exp(-tunnel_barrier)
    transition_weight = gamma_proxy * math.exp(
        -(
            config.lambda_structural_uphill * structural_uphill
            + config.lambda_control * avg_control
            + config.lambda_tadpole * avg_tadpole
            + config.lambda_stability_loss * stability_loss
        )
    )

    return {
        "source": i,
        "target": j,
        "source_candidate_id": src["candidate_id"],
        "target_candidate_id": dst["candidate_id"],
        "source_background_id": src["background_id"],
        "target_background_id": dst["background_id"],
        "source_type": src["vacuum_type"],
        "target_type": dst["vacuum_type"],
        "feature_distance": feature_dist,
        "flux_step": flux_step,
        "tau_step": tau_step,
        "w0_step": w0_step,
        "energy_step": energy_step,
        "tadpole_step": tadpole_step,
        "tunnel_barrier_proxy": tunnel_barrier,
        "delta_F_bridge": delta_F,
        "delta_F_pheno": delta_F_pheno,
        "structural_uphill": structural_uphill,
        "structural_downhill": structural_downhill,
        "avg_control": avg_control,
        "avg_tadpole": avg_tadpole,
        "stability_loss": stability_loss,
        "gamma_proxy": gamma_proxy,
        "transition_weight": transition_weight,
        "structural_action": structural_action,
        "energy_action": energy_action,
        "source_F_bridge": src["F_bridge"],
        "target_F_bridge": dst["F_bridge"],
        "source_abs_energy": src["abs_energy"],
        "target_abs_energy": dst["abs_energy"],
    }


def build_edges(nodes: list[dict[str, Any]], config: PathConfig) -> list[dict[str, Any]]:
    x = feature_matrix(nodes)
    diff = x[:, None, :] - x[None, :, :]
    dist = np.sqrt(np.sum(diff * diff, axis=2))
    edges: list[dict[str, Any]] = []
    seen: set[tuple[int, int]] = set()
    for i in range(len(nodes)):
        nearest = np.argsort(dist[i])
        neighbours = [int(j) for j in nearest if int(j) != i][: config.k_neighbours]
        for j in neighbours:
            if (i, j) in seen:
                continue
            seen.add((i, j))
            edges.append(make_edge(i, j, nodes, float(dist[i, j]), config))
    return edges


def normalize_edge_probabilities(edges: list[dict[str, Any]], field: str, out_field: str) -> None:
    by_source: dict[int, list[dict[str, Any]]] = {}
    for edge in edges:
        by_source.setdefault(int(edge["source"]), []).append(edge)
    for source_edges in by_source.values():
        weights = np.array([math.exp(-float(e[field])) for e in source_edges], dtype=float)
        denom = float(np.sum(weights))
        if denom <= 0.0:
            for edge in source_edges:
                edge[out_field] = 0.0
        else:
            for edge, weight in zip(source_edges, weights):
                edge[out_field] = float(weight / denom)


def dijkstra(edges: list[dict[str, Any]], n_nodes: int, source: int, target: int, weight_field: str) -> tuple[float, list[int]]:
    adj: dict[int, list[tuple[int, float]]] = {}
    for edge in edges:
        adj.setdefault(int(edge["source"]), []).append((int(edge["target"]), float(edge[weight_field])))

    dist = [float("inf")] * n_nodes
    prev = [-1] * n_nodes
    dist[source] = 0.0
    heap = [(0.0, source)]
    while heap:
        cost, node = heapq.heappop(heap)
        if cost > dist[node]:
            continue
        if node == target:
            break
        for nxt, w in adj.get(node, []):
            new = cost + w
            if new < dist[nxt]:
                dist[nxt] = new
                prev[nxt] = node
                heapq.heappush(heap, (new, nxt))

    if not math.isfinite(dist[target]):
        return float("inf"), []
    path = [target]
    cur = target
    while cur != source:
        cur = prev[cur]
        if cur < 0:
            return float("inf"), []
        path.append(cur)
    path.reverse()
    return dist[target], path


def path_edges(path: list[int], edges: list[dict[str, Any]]) -> list[dict[str, Any]]:
    edge_map = {(int(e["source"]), int(e["target"])): e for e in edges}
    return [edge_map[(a, b)] for a, b in zip(path[:-1], path[1:])]


def follow_greedy(edges: list[dict[str, Any]], start: int, probability_field: str, max_steps: int = 30) -> list[int]:
    by_source: dict[int, list[dict[str, Any]]] = {}
    for edge in edges:
        by_source.setdefault(int(edge["source"]), []).append(edge)
    path = [start]
    seen = {start}
    current = start
    for _ in range(max_steps):
        candidates = by_source.get(current, [])
        if not candidates:
            break
        best = max(candidates, key=lambda e: float(e[probability_field]))
        nxt = int(best["target"])
        path.append(nxt)
        if nxt in seen:
            break
        seen.add(nxt)
        current = nxt
    return path


def save_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            clean = {}
            for field in fields:
                value = row.get(field, "")
                if isinstance(value, np.generic):
                    value = value.item()
                clean[field] = value
            writer.writerow(clean)


def main() -> None:
    source_path = resolve_source(CONFIG.source_results)
    payload, nodes = load_minima(source_path)
    edges = build_edges(nodes, CONFIG)
    normalize_edge_probabilities(edges, "structural_action", "P_structural")
    normalize_edge_probabilities(edges, "energy_action", "P_energy")

    best_energy_node = min(nodes, key=lambda n: n["abs_energy"])
    best_bridge_node = min(nodes, key=lambda n: n["F_bridge"])
    best_pheno_node = min(nodes, key=lambda n: n["F_pheno"])

    path_cost, path = dijkstra(edges, len(nodes), best_energy_node["node_id"], best_bridge_node["node_id"], "structural_action")
    e_path_cost, e_path = dijkstra(edges, len(nodes), best_energy_node["node_id"], best_bridge_node["node_id"], "energy_action")
    selected_path_edges = path_edges(path, edges) if path else []

    top_source_nodes = sorted(nodes, key=lambda n: n["abs_energy"])[: CONFIG.top_path_sources]
    greedy_rows = []
    for node in top_source_nodes:
        s_path = follow_greedy(edges, node["node_id"], "P_structural")
        e_path_greedy = follow_greedy(edges, node["node_id"], "P_energy")
        greedy_rows.append(
            {
                "source_node": node["node_id"],
                "source_candidate_id": node["candidate_id"],
                "source_rank_energy_minima": node["rank_energy_minima"],
                "structural_path": "->".join(map(str, s_path)),
                "structural_terminal": s_path[-1],
                "structural_terminal_F_bridge": nodes[s_path[-1]]["F_bridge"],
                "energy_path": "->".join(map(str, e_path_greedy)),
                "energy_terminal": e_path_greedy[-1],
                "energy_terminal_F_bridge": nodes[e_path_greedy[-1]]["F_bridge"],
            }
        )

    top_edges_struct = sorted(edges, key=lambda e: e["P_structural"], reverse=True)[:20]
    top_edges_energy = sorted(edges, key=lambda e: e["P_energy"], reverse=True)[:20]
    top_overlap = len({(e["source"], e["target"]) for e in top_edges_struct} & {(e["source"], e["target"]) for e in top_edges_energy})

    downhill_structural = sum(1 for e in edges if e["delta_F_bridge"] < 0)
    top_downhill_structural = sum(1 for e in top_edges_struct if e["delta_F_bridge"] < 0)

    node_fields = [
        "node_id",
        "candidate_id",
        "background_id",
        "vacuum_type",
        "tau",
        "W0",
        "D",
        "energy",
        "abs_energy",
        "lambda_min",
        "F_bridge",
        "F_pheno",
        "Delta_Q_10",
        "Delta_dim",
        "Delta_E",
        "R_support",
        "rank_energy_minima",
        "rank_bridge_minima",
        "rank_pheno_minima",
    ]
    save_csv(OUT / "path_selection_nodes.csv", nodes, node_fields)

    edge_fields = [
        "source",
        "target",
        "source_candidate_id",
        "target_candidate_id",
        "source_background_id",
        "target_background_id",
        "source_type",
        "target_type",
        "feature_distance",
        "tunnel_barrier_proxy",
        "delta_F_bridge",
        "structural_uphill",
        "structural_downhill",
        "avg_control",
        "avg_tadpole",
        "stability_loss",
        "gamma_proxy",
        "transition_weight",
        "structural_action",
        "energy_action",
        "P_structural",
        "P_energy",
        "source_F_bridge",
        "target_F_bridge",
        "source_abs_energy",
        "target_abs_energy",
    ]
    save_csv(OUT / "path_selection_edges.csv", edges, edge_fields)
    save_csv(
        OUT / "path_selection_greedy_paths.csv",
        greedy_rows,
        [
            "source_node",
            "source_candidate_id",
            "source_rank_energy_minima",
            "structural_path",
            "structural_terminal",
            "structural_terminal_F_bridge",
            "energy_path",
            "energy_terminal",
            "energy_terminal_F_bridge",
        ],
    )

    path_rows = []
    for step, edge in enumerate(selected_path_edges, start=1):
        path_rows.append(
            {
                "step": step,
                "source": edge["source"],
                "target": edge["target"],
                "source_F_bridge": edge["source_F_bridge"],
                "target_F_bridge": edge["target_F_bridge"],
                "delta_F_bridge": edge["delta_F_bridge"],
                "tunnel_barrier_proxy": edge["tunnel_barrier_proxy"],
                "structural_action": edge["structural_action"],
                "P_structural": edge["P_structural"],
            }
        )
    if path_rows:
        save_csv(
            OUT / "best_energy_to_best_bridge_path.csv",
            path_rows,
            [
                "step",
                "source",
                "target",
                "source_F_bridge",
                "target_F_bridge",
                "delta_F_bridge",
                "tunnel_barrier_proxy",
                "structural_action",
                "P_structural",
            ],
        )

    # Plots
    plt.rcParams.update({"font.size": 11, "figure.dpi": 150, "savefig.dpi": 170})

    fig, ax = plt.subplots(figsize=(9.2, 6.8))
    sc = ax.scatter(
        [n["rank_energy_minima"] for n in nodes],
        [n["rank_bridge_minima"] for n in nodes],
        c=[n["Delta_Q_10"] for n in nodes],
        s=55,
        cmap="viridis",
        alpha=0.82,
        edgecolors="white",
        linewidth=0.4,
    )
    ax.scatter(best_energy_node["rank_energy_minima"], best_energy_node["rank_bridge_minima"], marker="*", s=260, color="tab:orange", edgecolors="black", label="best energy")
    ax.scatter(best_bridge_node["rank_energy_minima"], best_bridge_node["rank_bridge_minima"], marker="*", s=260, color="tab:blue", edgecolors="black", label="best structural")
    ax.invert_yaxis()
    ax.set_xlabel("Energy rank among minima")
    ax.set_ylabel("Structural rank among minima")
    ax.set_title("Static ranking plane for reduced KKLT minima")
    ax.grid(alpha=0.25)
    ax.legend()
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label(r"$\Delta Q_{10}$")
    fig.tight_layout()
    fig.savefig(OUT / "fig1_static_ranking_plane.png", bbox_inches="tight")
    plt.close(fig)

    fig2, ax = plt.subplots(figsize=(10.5, 7.0))
    x_rank = np.array([n["rank_energy_minima"] for n in nodes])
    y_rank = np.array([n["rank_bridge_minima"] for n in nodes])
    ax.scatter(x_rank, y_rank, s=42, color="0.75", alpha=0.7)
    for edge in top_edges_struct[:35]:
        src = nodes[int(edge["source"])]
        dst = nodes[int(edge["target"])]
        ax.annotate(
            "",
            xy=(dst["rank_energy_minima"], dst["rank_bridge_minima"]),
            xytext=(src["rank_energy_minima"], src["rank_bridge_minima"]),
            arrowprops=dict(arrowstyle="->", lw=1.0, alpha=0.45, color="tab:blue"),
        )
    if path:
        for edge in selected_path_edges:
            src = nodes[int(edge["source"])]
            dst = nodes[int(edge["target"])]
            ax.annotate(
                "",
                xy=(dst["rank_energy_minima"], dst["rank_bridge_minima"]),
                xytext=(src["rank_energy_minima"], src["rank_bridge_minima"]),
                arrowprops=dict(arrowstyle="->", lw=2.8, alpha=0.85, color="tab:red"),
            )
    ax.scatter(best_energy_node["rank_energy_minima"], best_energy_node["rank_bridge_minima"], marker="*", s=260, color="tab:orange", edgecolors="black", label="source: best energy")
    ax.scatter(best_bridge_node["rank_energy_minima"], best_bridge_node["rank_bridge_minima"], marker="*", s=260, color="tab:blue", edgecolors="black", label="target: best structural")
    ax.invert_yaxis()
    ax.set_xlabel("Energy rank among minima")
    ax.set_ylabel("Structural rank among minima")
    ax.set_title("Structural transition flow on ranking plane")
    ax.grid(alpha=0.25)
    ax.legend()
    fig2.tight_layout()
    fig2.savefig(OUT / "fig2_structural_flow_graph.png", bbox_inches="tight")
    plt.close(fig2)

    fig3, ax = plt.subplots(figsize=(10.6, 6.2))
    if selected_path_edges:
        labels = [f"{e['source']}->{e['target']}" for e in selected_path_edges]
        x = np.arange(len(selected_path_edges))
        tunnel = np.array([e["tunnel_barrier_proxy"] for e in selected_path_edges])
        uphill = np.array([CONFIG.lambda_structural_uphill * e["structural_uphill"] for e in selected_path_edges])
        control = np.array([CONFIG.lambda_control * e["avg_control"] for e in selected_path_edges])
        tadpole = np.array([CONFIG.lambda_tadpole * e["avg_tadpole"] for e in selected_path_edges])
        stability = np.array([CONFIG.lambda_stability_loss * e["stability_loss"] for e in selected_path_edges])
        bottom = np.zeros_like(tunnel)
        for vals, name, color in [
            (tunnel, "barrier", "#455a64"),
            (uphill, "structural uphill", "#d81b60"),
            (control, "EFT control", "#43a047"),
            (tadpole, "tadpole", "#fb8c00"),
            (stability, "stability loss", "#5e35b1"),
        ]:
            ax.bar(x, vals, bottom=bottom, label=name, color=color)
            bottom = bottom + vals
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.set_ylabel("edge action contribution")
        ax.set_title("Selected structural path action decomposition")
        ax.grid(axis="y", alpha=0.22)
        ax.legend(ncol=3)
    else:
        ax.text(0.5, 0.5, "No path found", ha="center", va="center")
    fig3.tight_layout()
    fig3.savefig(OUT / "fig3_path_action_decomposition.png", bbox_inches="tight")
    plt.close(fig3)

    sorted_nodes = sorted(nodes, key=lambda n: n["F_bridge"])
    node_order = {n["node_id"]: i for i, n in enumerate(sorted_nodes)}
    mat = np.zeros((len(nodes), len(nodes)))
    for edge in edges:
        mat[node_order[int(edge["source"])], node_order[int(edge["target"])]] = max(mat[node_order[int(edge["source"])], node_order[int(edge["target"])]], edge["P_structural"])
    fig4, ax = plt.subplots(figsize=(8.5, 7.2))
    im = ax.imshow(mat, origin="lower", aspect="auto", cmap="magma")
    ax.set_xlabel("target node sorted by structural rank")
    ax.set_ylabel("source node sorted by structural rank")
    ax.set_title("Structural transition probability matrix")
    cbar = fig4.colorbar(im, ax=ax)
    cbar.set_label(r"$P_{ij}^{struct}$")
    fig4.tight_layout()
    fig4.savefig(OUT / "fig4_transition_probability_matrix.png", bbox_inches="tight")
    plt.close(fig4)

    summary = {
        "model": "MAAT String Selection II: Structural Path Selection in a Reduced KKLT Landscape",
        "status": (
            "Phenomenological transition-graph benchmark. Not a Coleman-De Luccia "
            "calculation, not a microscopic string decay-rate derivation, and not a "
            "complete landscape dynamics model."
        ),
        "source_results": str(source_path),
        "config": asdict(CONFIG),
        "source_summary": payload.get("summary", {}),
        "n_minima_nodes": len(nodes),
        "n_directed_edges": len(edges),
        "class_counts_minima": {
            name: sum(1 for n in nodes if n["vacuum_type"] == name)
            for name in sorted(set(n["vacuum_type"] for n in nodes))
        },
        "best_energy_node": {
            "node_id": best_energy_node["node_id"],
            "candidate_id": best_energy_node["candidate_id"],
            "background_id": best_energy_node["background_id"],
            "F_bridge": best_energy_node["F_bridge"],
            "abs_energy": best_energy_node["abs_energy"],
            "Delta_Q_10": best_energy_node["Delta_Q_10"],
            "rank_bridge_minima": best_energy_node["rank_bridge_minima"],
        },
        "best_bridge_node": {
            "node_id": best_bridge_node["node_id"],
            "candidate_id": best_bridge_node["candidate_id"],
            "background_id": best_bridge_node["background_id"],
            "F_bridge": best_bridge_node["F_bridge"],
            "abs_energy": best_bridge_node["abs_energy"],
            "Delta_Q_10": best_bridge_node["Delta_Q_10"],
            "rank_energy_minima": best_bridge_node["rank_energy_minima"],
        },
        "best_pheno_node": {
            "node_id": best_pheno_node["node_id"],
            "candidate_id": best_pheno_node["candidate_id"],
            "background_id": best_pheno_node["background_id"],
            "F_pheno": best_pheno_node["F_pheno"],
        },
        "spearman_abs_energy_vs_F_bridge_minima": spearman([n["abs_energy"] for n in nodes], [n["F_bridge"] for n in nodes]),
        "spearman_abs_delta_q_vs_F_bridge_minima": spearman([n["Delta_Q_10"] for n in nodes], [n["F_bridge"] for n in nodes]),
        "edge_statistics": {
            "mean_delta_F_bridge": float(np.mean([e["delta_F_bridge"] for e in edges])),
            "fraction_structural_downhill_edges": float(downhill_structural / len(edges)),
            "top20_structural_edges_downhill": int(top_downhill_structural),
            "top20_overlap_structural_energy": int(top_overlap),
            "mean_structural_action": float(np.mean([e["structural_action"] for e in edges])),
            "mean_energy_action": float(np.mean([e["energy_action"] for e in edges])),
        },
        "best_energy_to_best_bridge_structural_path": {
            "cost": path_cost,
            "nodes": path,
            "length_edges": max(0, len(path) - 1),
        },
        "best_energy_to_best_bridge_energy_action_path": {
            "cost": e_path_cost,
            "nodes": e_path,
            "length_edges": max(0, len(e_path) - 1),
        },
        "greedy_paths": greedy_rows,
        "outputs": {
            "nodes_csv": "outputs_path_selection/path_selection_nodes.csv",
            "edges_csv": "outputs_path_selection/path_selection_edges.csv",
            "greedy_paths_csv": "outputs_path_selection/path_selection_greedy_paths.csv",
            "selected_path_csv": "outputs_path_selection/best_energy_to_best_bridge_path.csv",
            "figures": [
                "outputs_path_selection/fig1_static_ranking_plane.png",
                "outputs_path_selection/fig2_structural_flow_graph.png",
                "outputs_path_selection/fig3_path_action_decomposition.png",
                "outputs_path_selection/fig4_transition_probability_matrix.png",
            ],
        },
    }
    with open(OUT / "path_selection_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("=" * 78)
    print("MAAT String Selection II — Structural Path Selection")
    print("=" * 78)
    print(f"source: {source_path}")
    print(f"nodes (minima): {len(nodes)}")
    print(f"directed edges: {len(edges)}")
    print(f"best energy node -> best structural node path: {path} cost={path_cost:.6g}")
    print(f"top20 overlap structural vs energy edges: {top_overlap}/20")
    print(f"fraction structural downhill edges: {downhill_structural / len(edges):.3f}")
    print(f"outputs: {OUT.resolve()}")


if __name__ == "__main__":
    main()
