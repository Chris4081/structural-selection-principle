#!/usr/bin/env python3
"""Paper 70 diagnostic analysis: local structural channels in SAT.

This script is intentionally diagnostic and post-hoc. It does not retune the
Paper 69 gate, does not alter CDCL policies, and does not re-run the solver.
It asks a narrower question:

    Which local SAT channels correlate with the Gate-vs-MOMS regret observed in
    the first external SATLIB smoke execution?

Raw CNF files are used only as local input. The emitted CSV/JSON/PNG outputs are
derived diagnostics and do not contain DIMACS clauses.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


EPS = 1.0e-12


Clause = Tuple[int, ...]


def read_csv_dict(path: Path) -> List[dict]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv_dict(path: Path, rows: Sequence[dict], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_dimacs(path: Path) -> Tuple[int, List[Clause]]:
    n_vars = 0
    clauses: List[Clause] = []
    current: List[int] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("c") or line.startswith("%"):
                continue
            if line.startswith("p"):
                parts = line.split()
                if len(parts) >= 4 and parts[1].lower() == "cnf":
                    n_vars = int(parts[2])
                continue
            for token in line.split():
                lit = int(token)
                if lit == 0:
                    if current:
                        clauses.append(tuple(current))
                    current = []
                else:
                    current.append(lit)
    if current:
        clauses.append(tuple(current))
    if n_vars <= 0:
        n_vars = max((abs(lit) for clause in clauses for lit in clause), default=0)
    return n_vars, clauses


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def entropy_from_counts(counts: Iterable[int]) -> float:
    vals = [c for c in counts if c > 0]
    total = sum(vals)
    if total <= 0:
        return 0.0
    probs = [c / total for c in vals]
    return -sum(p * math.log(p + EPS) for p in probs)


def normalized_entropy_from_counts(counts: Iterable[int]) -> float:
    vals = [c for c in counts if c > 0]
    if len(vals) <= 1:
        return 0.0
    return entropy_from_counts(vals) / math.log(len(vals))


def mean(vals: Sequence[float]) -> float:
    return sum(vals) / len(vals) if vals else 0.0


def std(vals: Sequence[float]) -> float:
    if len(vals) < 2:
        return 0.0
    return statistics.pstdev(vals)


def cv(vals: Sequence[float]) -> float:
    m = mean(vals)
    return std(vals) / (abs(m) + EPS)


def rankdata(values: Sequence[float]) -> List[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        avg_rank = (i + j - 1) / 2.0 + 1.0
        for k in range(i, j):
            ranks[order[k]] = avg_rank
        i = j
    return ranks


def pearson(x: Sequence[float], y: Sequence[float]) -> float:
    if len(x) != len(y) or len(x) < 2:
        return 0.0
    mx, my = mean(x), mean(y)
    dx = [v - mx for v in x]
    dy = [v - my for v in y]
    sx = math.sqrt(sum(v * v for v in dx))
    sy = math.sqrt(sum(v * v for v in dy))
    if sx <= EPS or sy <= EPS:
        return 0.0
    return sum(a * b for a, b in zip(dx, dy)) / (sx * sy)


def spearman(x: Sequence[float], y: Sequence[float]) -> float:
    return pearson(rankdata(x), rankdata(y))


def unit_propagate(clauses: Sequence[Clause], assumption: int) -> Tuple[int, bool]:
    """Return (number of assigned variables, conflict flag) after unit propagation."""
    assignment: Dict[int, bool] = {abs(assumption): assumption > 0}
    changed = True
    conflict = False
    while changed and not conflict:
        changed = False
        for clause in clauses:
            satisfied = False
            unassigned: List[int] = []
            for lit in clause:
                var = abs(lit)
                val = assignment.get(var)
                if val is None:
                    unassigned.append(lit)
                elif (lit > 0 and val) or (lit < 0 and not val):
                    satisfied = True
                    break
            if satisfied:
                continue
            if not unassigned:
                conflict = True
                break
            if len(unassigned) == 1:
                lit = unassigned[0]
                var = abs(lit)
                val = lit > 0
                old = assignment.get(var)
                if old is not None and old != val:
                    conflict = True
                    break
                if old is None:
                    assignment[var] = val
                    changed = True
    return len(assignment), conflict


def compute_local_channels(n_vars: int, clauses: Sequence[Clause]) -> dict:
    lengths = [len(c) for c in clauses]
    n_clauses = len(clauses)
    min_len = min(lengths) if lengths else 0
    short_clauses = [c for c in clauses if len(c) == min_len]

    var_occ = Counter(abs(lit) for clause in clauses for lit in clause)
    lit_occ = Counter(lit for clause in clauses for lit in clause)
    short_var_occ = Counter(abs(lit) for clause in short_clauses for lit in clause)
    short_total_lits = sum(short_var_occ.values())
    jw_var = Counter()
    for clause in clauses:
        w = 2.0 ** (-len(clause))
        for lit in clause:
            jw_var[abs(lit)] += w

    var_occ_values = [var_occ.get(i, 0) for i in range(1, n_vars + 1)]
    short_var_values = [short_var_occ.get(i, 0) for i in range(1, n_vars + 1)]
    jw_values = [jw_var.get(i, 0.0) for i in range(1, n_vars + 1)]

    pos_counts = [lit_occ.get(i, 0) for i in range(1, n_vars + 1)]
    neg_counts = [lit_occ.get(-i, 0) for i in range(1, n_vars + 1)]
    literal_imbalance_vals = [
        abs(p - q) / (p + q + EPS) for p, q in zip(pos_counts, neg_counts) if p + q > 0
    ]

    length_counts = Counter(lengths)
    length_entropy = normalized_entropy_from_counts(length_counts.values())
    weighted_clause_pressure = sum(2.0 ** (-l) for l in lengths) / (n_clauses + EPS)
    short_clause_fraction = len(short_clauses) / (n_clauses + EPS)
    short_clause_pressure = len(short_clauses) / (n_vars + EPS)
    short_variable_concentration = (
        max(short_var_values) / (short_total_lits + EPS) if short_total_lits else 0.0
    )
    short_variable_entropy = normalized_entropy_from_counts(short_var_values)
    moms_signal = short_clause_pressure * (1.0 + short_variable_concentration)
    jw_concentration = max(jw_values) / (sum(jw_values) + EPS) if jw_values else 0.0

    cascade_sizes: List[int] = []
    cascade_conflicts = 0
    # Keep this transparent rather than optimized; Paper 70 uses only 30 instances.
    for var in range(1, n_vars + 1):
        for lit in (var, -var):
            size, conflict = unit_propagate(clauses, lit)
            cascade_sizes.append(size)
            cascade_conflicts += int(conflict)

    return {
        "min_clause_len": min_len,
        "mean_clause_len": mean(lengths),
        "clause_length_cv": cv([float(v) for v in lengths]),
        "clause_length_entropy": length_entropy,
        "unit_binary_fraction": sum(1 for l in lengths if l <= 2) / (n_clauses + EPS),
        "short_clause_fraction": short_clause_fraction,
        "short_clause_pressure": short_clause_pressure,
        "short_variable_concentration": short_variable_concentration,
        "short_variable_entropy": short_variable_entropy,
        "moms_local_signal": moms_signal,
        "weighted_clause_pressure": weighted_clause_pressure,
        "jw_variable_concentration": jw_concentration,
        "variable_occurrence_cv": cv([float(v) for v in var_occ_values]),
        "literal_polarity_imbalance": mean(literal_imbalance_vals),
        "propagation_cascade_mean": mean([float(v) for v in cascade_sizes]),
        "propagation_cascade_max": max(cascade_sizes) if cascade_sizes else 0,
        "propagation_cascade_cv": cv([float(v) for v in cascade_sizes]),
        "propagation_conflict_fraction": cascade_conflicts / (len(cascade_sizes) + EPS),
    }


def load_manifest_map(paper69_dir: Path, outputs_dir: Path) -> Dict[str, dict]:
    candidates = [
        outputs_dir / "paper69_dataset_manifest_detected.csv",
        paper69_dir / "dataset_manifest.csv",
    ]
    for path in candidates:
        if path.exists():
            rows = read_csv_dict(path)
            return {r["dataset_id"]: r for r in rows}
    raise FileNotFoundError("No Paper-69 manifest CSV found")


def plot_bar_correlations(rows: Sequence[dict], output: Path) -> None:
    rows = sorted(rows, key=lambda r: abs(float(r["spearman_regret_gate_to_moms"])), reverse=True)
    names = [r["channel"] for r in rows]
    vals = [float(r["spearman_regret_gate_to_moms"]) for r in rows]
    plt.figure(figsize=(9, 6))
    colors = ["#b84a3a" if v > 0 else "#386fa4" for v in vals]
    plt.barh(range(len(names)), vals, color=colors)
    plt.axvline(0, color="black", linewidth=0.8)
    plt.yticks(range(len(names)), names, fontsize=8)
    plt.xlabel("Spearman correlation with gate-to-MOMS regret")
    plt.title("Paper 70: Local SAT channels vs missing MOMS advantage")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(output, dpi=180)
    plt.close()


def plot_scatter(rows: Sequence[dict], xkey: str, ykey: str, output: Path, title: str) -> None:
    families = sorted({r["family"] for r in rows})
    cmap = plt.get_cmap("tab10")
    color_map = {fam: cmap(i % 10) for i, fam in enumerate(families)}
    plt.figure(figsize=(7, 5))
    for fam in families:
        sub = [r for r in rows if r["family"] == fam]
        plt.scatter(
            [float(r[xkey]) for r in sub],
            [float(r[ykey]) for r in sub],
            label=fam.replace("satlib_", ""),
            s=42,
            alpha=0.85,
            color=color_map[fam],
            edgecolor="black",
            linewidth=0.4,
        )
    plt.xlabel(xkey.replace("_", " "))
    plt.ylabel(ykey.replace("_", " "))
    plt.title(title)
    plt.legend(fontsize=7, loc="best")
    plt.tight_layout()
    plt.savefig(output, dpi=180)
    plt.close()


def plot_family_summary(rows: Sequence[dict], output: Path) -> None:
    families = sorted({r["family"] for r in rows})
    channels = ["moms_local_signal", "short_clause_pressure", "propagation_cascade_mean", "variable_occurrence_cv"]
    family_means: Dict[str, Dict[str, float]] = {}
    for fam in families:
        sub = [r for r in rows if r["family"] == fam]
        family_means[fam] = {c: mean([float(r[c]) for r in sub]) for c in channels}

    # Min-max scale each channel for a readable diagnostic heatmap.
    scaled = []
    for fam in families:
        scaled_row = []
        for c in channels:
            vals = [family_means[f][c] for f in families]
            lo, hi = min(vals), max(vals)
            scaled_row.append((family_means[fam][c] - lo) / (hi - lo + EPS))
        scaled.append(scaled_row)

    plt.figure(figsize=(8, 4.5))
    plt.imshow(scaled, aspect="auto", cmap="viridis")
    plt.colorbar(label="family-normalized channel intensity")
    plt.xticks(range(len(channels)), [c.replace("_", " ") for c in channels], rotation=25, ha="right")
    plt.yticks(range(len(families)), [f.replace("satlib_", "") for f in families], fontsize=8)
    plt.title("Paper 70: Family-level local channel structure")
    plt.tight_layout()
    plt.savefig(output, dpi=180)
    plt.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--paper69-dir",
        type=Path,
        default=Path("../gate_challenge_sat_paper69"),
        help="Path to the Paper 69 experiment directory containing data_external and outputs.",
    )
    ap.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs_paper70_local_sat_channels"),
    )
    args = ap.parse_args()

    paper69_dir = args.paper69_dir.resolve()
    outputs69 = paper69_dir / "outputs_paper69_sat_gate_challenge"
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    instance_path = outputs69 / "paper69_moms_vs_gate_instance_diagnostics.csv"
    solve_path = outputs69 / "paper69_solve_records.csv"
    features_path = outputs69 / "paper69_gate_features.csv"
    if not instance_path.exists():
        raise FileNotFoundError(instance_path)

    instances = read_csv_dict(instance_path)
    gate_features = {r["dataset_id"]: r for r in read_csv_dict(features_path)}
    manifest = load_manifest_map(paper69_dir, outputs69)

    rows: List[dict] = []
    missing_cnf: List[str] = []
    for inst in instances:
        dataset_id = inst["dataset_id"]
        mrow = manifest.get(dataset_id)
        if not mrow:
            missing_cnf.append(dataset_id)
            continue
        cnf_path = paper69_dir / mrow["local_path"]
        if not cnf_path.exists():
            missing_cnf.append(dataset_id)
            continue
        expected_hash = mrow.get("sha256", "")
        actual_hash = sha256_file(cnf_path)
        if expected_hash and actual_hash != expected_hash:
            raise RuntimeError(f"SHA256 mismatch for {dataset_id}: {cnf_path}")
        n_vars, clauses = parse_dimacs(cnf_path)
        local = compute_local_channels(n_vars, clauses)
        gf = gate_features.get(dataset_id, {})
        row = {
            "dataset_id": dataset_id,
            "family": inst["family"],
            "split": inst["split"],
            "n_vars": n_vars,
            "n_clauses": len(clauses),
            "regret_gate_to_moms": float(inst["regret_gate_to_moms"]),
            "delta_gate_vs_moms": float(inst["delta_gate_vs_moms"]),
            "gate_v17_cost": float(inst["gate_v17"]),
            "moms_cost": float(inst["moms"]),
            "score_with_R_cost": float(inst["score_with_R"]),
            "G_gate": float(inst["G_gate"]),
            "R_rob": float(inst["R_rob"]),
            "H": float(inst["H"]),
            "B": float(inst["B"]),
            "S": float(inst["S"]),
            "V": float(inst["V"]),
            "paper69_degree_cv": float(gf.get("degree_cv", 0.0)),
            "paper69_literal_imbalance": float(gf.get("literal_imbalance", 0.0)),
            "paper69_clustering_mean": float(gf.get("clustering_mean", 0.0)),
            "paper69_largest_component_frac": float(gf.get("largest_component_frac", 0.0)),
        }
        row.update(local)
        rows.append(row)

    if missing_cnf:
        raise FileNotFoundError(
            "Missing local CNFs for Paper 70 analysis: " + ", ".join(missing_cnf[:10])
        )

    feature_fields = [
        "dataset_id",
        "family",
        "split",
        "n_vars",
        "n_clauses",
        "regret_gate_to_moms",
        "delta_gate_vs_moms",
        "gate_v17_cost",
        "moms_cost",
        "score_with_R_cost",
        "G_gate",
        "R_rob",
        "H",
        "B",
        "S",
        "V",
        "paper69_degree_cv",
        "paper69_literal_imbalance",
        "paper69_clustering_mean",
        "paper69_largest_component_frac",
        "min_clause_len",
        "mean_clause_len",
        "clause_length_cv",
        "clause_length_entropy",
        "unit_binary_fraction",
        "short_clause_fraction",
        "short_clause_pressure",
        "short_variable_concentration",
        "short_variable_entropy",
        "moms_local_signal",
        "weighted_clause_pressure",
        "jw_variable_concentration",
        "variable_occurrence_cv",
        "literal_polarity_imbalance",
        "propagation_cascade_mean",
        "propagation_cascade_max",
        "propagation_cascade_cv",
        "propagation_conflict_fraction",
    ]
    write_csv_dict(output_dir / "paper70_local_channel_features.csv", rows, feature_fields)

    target = [float(r["regret_gate_to_moms"]) for r in rows]
    channel_names = [
        "G_gate",
        "R_rob",
        "paper69_degree_cv",
        "paper69_literal_imbalance",
        "paper69_clustering_mean",
        "paper69_largest_component_frac",
        "min_clause_len",
        "mean_clause_len",
        "clause_length_cv",
        "clause_length_entropy",
        "unit_binary_fraction",
        "short_clause_fraction",
        "short_clause_pressure",
        "short_variable_concentration",
        "short_variable_entropy",
        "moms_local_signal",
        "weighted_clause_pressure",
        "jw_variable_concentration",
        "variable_occurrence_cv",
        "literal_polarity_imbalance",
        "propagation_cascade_mean",
        "propagation_cascade_max",
        "propagation_cascade_cv",
        "propagation_conflict_fraction",
    ]
    corr_rows = []
    for c in channel_names:
        vals = [float(r[c]) for r in rows]
        corr_rows.append(
            {
                "channel": c,
                "spearman_regret_gate_to_moms": spearman(vals, target),
                "pearson_regret_gate_to_moms": pearson(vals, target),
                "mean": mean(vals),
                "std": std(vals),
            }
        )
    write_csv_dict(
        output_dir / "paper70_channel_correlations.csv",
        corr_rows,
        ["channel", "spearman_regret_gate_to_moms", "pearson_regret_gate_to_moms", "mean", "std"],
    )

    family_rows = []
    families = sorted({r["family"] for r in rows})
    family_fields = [
        "family",
        "n_instances",
        "mean_regret_gate_to_moms",
        "median_regret_gate_to_moms",
        "mean_moms_local_signal",
        "mean_short_clause_pressure",
        "mean_short_variable_concentration",
        "mean_weighted_clause_pressure",
        "mean_propagation_cascade_mean",
        "mean_variable_occurrence_cv",
        "mean_G_gate",
    ]
    for fam in families:
        sub = [r for r in rows if r["family"] == fam]
        regrets = [float(r["regret_gate_to_moms"]) for r in sub]
        family_rows.append(
            {
                "family": fam,
                "n_instances": len(sub),
                "mean_regret_gate_to_moms": mean(regrets),
                "median_regret_gate_to_moms": statistics.median(regrets),
                "mean_moms_local_signal": mean([float(r["moms_local_signal"]) for r in sub]),
                "mean_short_clause_pressure": mean([float(r["short_clause_pressure"]) for r in sub]),
                "mean_short_variable_concentration": mean(
                    [float(r["short_variable_concentration"]) for r in sub]
                ),
                "mean_weighted_clause_pressure": mean(
                    [float(r["weighted_clause_pressure"]) for r in sub]
                ),
                "mean_propagation_cascade_mean": mean(
                    [float(r["propagation_cascade_mean"]) for r in sub]
                ),
                "mean_variable_occurrence_cv": mean([float(r["variable_occurrence_cv"]) for r in sub]),
                "mean_G_gate": mean([float(r["G_gate"]) for r in sub]),
            }
        )
    write_csv_dict(output_dir / "paper70_family_channel_summary.csv", family_rows, family_fields)

    plot_bar_correlations(corr_rows, output_dir / "fig1_channel_correlations.png")
    plot_scatter(
        rows,
        "moms_local_signal",
        "regret_gate_to_moms",
        output_dir / "fig2_moms_signal_vs_regret.png",
        "MOMS-local signal vs gate-to-MOMS regret",
    )
    plot_family_summary(rows, output_dir / "fig3_family_channel_heatmap.png")
    plot_scatter(
        rows,
        "propagation_cascade_mean",
        "regret_gate_to_moms",
        output_dir / "fig4_propagation_vs_regret.png",
        "Propagation cascade proxy vs gate-to-MOMS regret",
    )

    strongest = sorted(corr_rows, key=lambda r: abs(float(r["spearman_regret_gate_to_moms"])), reverse=True)[:8]
    summary = {
        "paper": 70,
        "title": "Local Structural Channels in SAT",
        "analysis_type": "post_hoc_diagnostic_not_retuning",
        "source": "Paper 69 SATLIB family-balanced smoke execution",
        "paper68_preregistration_doi": "10.5281/zenodo.20882852",
        "n_instances": len(rows),
        "families": families,
        "target": "regret_gate_to_moms (positive means gate costs more than MOMS)",
        "raw_cnf_policy": "Raw CNFs are local inputs only and are not emitted in outputs.",
        "top_abs_spearman_channels": strongest,
        "outputs": [
            "paper70_local_channel_features.csv",
            "paper70_channel_correlations.csv",
            "paper70_family_channel_summary.csv",
            "fig1_channel_correlations.png",
            "fig2_moms_signal_vs_regret.png",
            "fig3_family_channel_heatmap.png",
            "fig4_propagation_vs_regret.png",
        ],
        "interpretation_rule": (
            "This diagnostic identifies missing local channels. It must not be used "
            "to alter Paper 69 gate parameters."
        ),
    }
    (output_dir / "paper70_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
