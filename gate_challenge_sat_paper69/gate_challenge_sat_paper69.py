#!/usr/bin/env python3
"""Paper 69: Gate Challenge I -- SAT/CDCL.

First execution scaffold for the preregistered MAAT v1.7 Gate Challenge.

This script is intentionally conservative:

- it requires external DIMACS CNF files under ``data_external/``;
- it does not generate local CNFs for primary evidence;
- it freezes the SAT polarity p_D=+1 from Paper 68;
- it compares gate-aware MAAT against score-with-R and classical heuristics;
- if no external CNFs are present, it writes a "not executed" run log.

The transparent CDCL backend is imported from Paper 63 to keep the solving loop
auditable and comparable with prior MAAT SAT work.  Paper 69 adds external
DIMACS parsing, protocol-gate calibration, classical MOMS/Jeroslow-Wang
baselines, shuffled-R null control, paired bootstrap confidence intervals, and
a strict dataset manifest.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import lzma
import math
import os
import random
import sys
import time
from bz2 import open as bz2_open
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
PAPER63_DIR = REPO / "experiments" / "sat_cdcl_structure_gated_mode_a_paper63"
sys.path.insert(0, str(PAPER63_DIR))

from sat_cdcl_structure_gated_mode_a_paper63 import (  # type: ignore
    CDCLSolver,
    EPS,
    canonical_clause,
)

_CACHE_ROOT = Path(os.environ.get("TMPDIR", "/tmp")) / "maat_paper69_cache"
(_CACHE_ROOT / "matplotlib").mkdir(parents=True, exist_ok=True)
(_CACHE_ROOT / "xdg").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE_ROOT / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE_ROOT / "xdg"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SEED = 69069
OUTDIR = ROOT / "outputs_paper69_sat_gate_challenge"
DATA_DIR = ROOT / "data_external"

POLICIES = [
    "vsids",
    "moms",
    "jeroslow_wang",
    "progress_only",
    "score_with_R",
    "gate_v17",
    "gate_shuffled_R",
]

PAPER68_DOI = "10.5281/zenodo.20882852"
SAT_POLARITY = +1
TAU_STRUCTURAL = 0.38
SOFTMAX_BETA = 8.0


@dataclass
class Instance:
    dataset_id: str
    family: str
    source_name: str
    source_url: str
    license_or_terms: str
    download_tls_mode: str
    archive_sha256: str
    local_path: str
    sha256: str
    n_vars: int
    n_clauses: int
    clauses: list[tuple[int, ...]]


@dataclass
class InstanceFeatures:
    dataset_id: str
    family: str
    split: str
    n_vars: int
    n_clauses: int
    alpha: float
    H: float
    B: float
    S: float
    V: float
    R_resp: float
    R_rob: float
    cal_center_R: float
    cal_mad_R: float
    z_R: float
    z_t: float
    z_x: float
    gate_signal: float
    G_gate: float
    G_hard: int
    G_shuffled: float
    degree_cv: float
    literal_imbalance: float
    clustering_mean: float
    largest_component_frac: float
    component_frac: float


@dataclass
class SolveRecord:
    dataset_id: str
    family: str
    split: str
    policy: str
    n_vars: int
    n_clauses: int
    alpha: float
    solved: int
    sat: int
    timeout: int
    decisions: int
    propagations: int
    conflicts: int
    learned_clauses: int
    max_level: int
    gate_value: float
    runtime_s: float
    compute_cost: float


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def stable_fold(key: str, n_folds: int = 5) -> int:
    digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
    return int(digest[:8], 16) % n_folds


def open_text(path: Path):
    suffixes = "".join(path.suffixes)
    if suffixes.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    if suffixes.endswith(".bz2"):
        return bz2_open(path, "rt", encoding="utf-8", errors="ignore")
    if suffixes.endswith(".xz") or suffixes.endswith(".lzma"):
        return lzma.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("rt", encoding="utf-8", errors="ignore")


def parse_dimacs(path: Path) -> tuple[int, list[tuple[int, ...]]]:
    n_vars = 0
    clauses: list[tuple[int, ...]] = []
    current: list[int] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("c"):
                continue
            if line.startswith("%"):
                break
            if line.startswith("p"):
                parts = line.split()
                if len(parts) >= 4 and parts[1].lower() == "cnf":
                    n_vars = int(parts[2])
                continue
            for tok in line.split():
                if tok.startswith("%"):
                    break
                lit = int(tok)
                if lit == 0:
                    clause = canonical_clause(current)
                    if clause is not None:
                        clauses.append(tuple(clause))
                    current = []
                else:
                    current.append(lit)
    if current:
        clause = canonical_clause(current)
        if clause is not None:
            clauses.append(tuple(clause))
    if n_vars <= 0:
        n_vars = max((abs(lit) for clause in clauses for lit in clause), default=0)
    return n_vars, clauses


def infer_family(path: Path) -> str:
    parts = [p.lower() for p in path.parts]
    name = path.name.lower()
    joined = "/".join(parts)
    if "uf" in name or "uuf" in name or "random" in joined:
        return "random_3sat"
    if "aim" in joined:
        return "aim"
    if "dubois" in joined:
        return "dubois"
    if "hole" in joined or "pigeon" in joined or "php" in joined:
        return "pigeonhole"
    if "color" in joined or "flat" in joined:
        return "graph_coloring"
    if "satcomp" in joined or "competition" in joined:
        return "competition"
    return path.parent.name or "external"


def load_manifest(manifest_path: Path) -> dict[str, dict[str, str]]:
    if not manifest_path.exists():
        raise FileNotFoundError(f"dataset manifest is required for external CNF execution: {manifest_path}")
    with manifest_path.open("r", newline="", encoding="utf-8") as handle:
        return {row["local_path"]: row for row in csv.DictReader(handle)}


def cnf_paths(data_dir: Path) -> list[Path]:
    patterns = ["*.cnf", "*.dimacs", "*.cnf.gz", "*.cnf.bz2", "*.cnf.xz", "*.cnf.lzma"]
    paths: list[Path] = []
    for pattern in patterns:
        paths.extend(data_dir.rglob(pattern))
    return sorted(set(paths))


def validate_dataset(data_dir: Path, manifest_path: Path, output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    required_columns = {
        "dataset_id",
        "family",
        "source_name",
        "source_url",
        "license_or_terms",
        "local_path",
        "sha256",
    }
    validation: dict[str, object] = {
        "paper": 69,
        "title": "Gate Challenge I: SAT/CDCL",
        "validation_type": "dataset_manifest_sha256_gatekeeper",
        "paper68_preregistration_doi": PAPER68_DOI,
        "manifest_path": str(manifest_path),
        "data_dir": str(data_dir),
        "checks": {},
        "counts": {},
        "failures": {},
        "dataset_ready": False,
    }

    if not manifest_path.exists():
        validation["status"] = "VALIDATION FAIL"
        validation["checks"] = {
            "Manifest": "FAIL",
            "CNFs": "NOT_RUN",
            "SHA256": "NOT_RUN",
            "Unexpected Files": "NOT_RUN",
            "Dataset Ready": "NO",
        }
        validation["failures"] = {"manifest": [f"missing manifest: {manifest_path}"]}
        (output_dir / "paper69_dataset_validation.json").write_text(json.dumps(validation, indent=2), encoding="utf-8")
        return validation

    with manifest_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = set(reader.fieldnames or [])
        rows = list(reader)

    missing_columns = sorted(required_columns - fieldnames)
    duplicate_paths = sorted(
        path
        for path in {row.get("local_path", "") for row in rows}
        if path and sum(1 for row in rows if row.get("local_path", "") == path) > 1
    )
    manifest_paths = {row.get("local_path", "") for row in rows if row.get("local_path", "")}
    actual_paths = {str(path.relative_to(ROOT)) for path in cnf_paths(data_dir)}

    missing_files: list[str] = []
    missing_hashes: list[str] = []
    sha_mismatches: list[dict[str, str]] = []
    parse_failures: list[dict[str, str]] = []
    verified = 0
    parseable = 0
    for row in rows:
        rel = row.get("local_path", "").strip()
        if not rel:
            missing_files.append("<empty local_path>")
            continue
        path = ROOT / rel
        if not path.exists():
            missing_files.append(rel)
            continue
        expected = row.get("sha256", "").strip().lower()
        if not expected:
            missing_hashes.append(rel)
            continue
        observed = sha256_file(path).lower()
        if observed != expected:
            sha_mismatches.append({"local_path": rel, "manifest_sha256": expected, "computed_sha256": observed})
        else:
            verified += 1
        try:
            n_vars, clauses = parse_dimacs(path)
            if n_vars <= 0 or not clauses:
                parse_failures.append({"local_path": rel, "error": "no_variables_or_clauses_detected"})
            else:
                parseable += 1
        except Exception as exc:  # noqa: BLE001 - validator should report all parse failures.
            parse_failures.append({"local_path": rel, "error": str(exc)})

    unexpected_files = sorted(actual_paths - manifest_paths)
    untracked_manifest_rows = sorted(manifest_paths - actual_paths)
    manifest_ok = bool(rows) and not missing_columns and not duplicate_paths
    cnfs_ok = bool(actual_paths) and not missing_files and not untracked_manifest_rows and not parse_failures
    sha_ok = not missing_hashes and not sha_mismatches and verified == len(rows)
    unexpected_ok = not unexpected_files
    ready = manifest_ok and cnfs_ok and sha_ok and unexpected_ok

    validation["status"] = "VALIDATION PASS" if ready else "VALIDATION FAIL"
    validation["dataset_ready"] = ready
    validation["checks"] = {
        "Manifest": "PASS" if manifest_ok else "FAIL",
        "CNFs": "PASS" if cnfs_ok else "FAIL",
        "SHA256": "PASS" if sha_ok else "FAIL",
        "Unexpected Files": "PASS" if unexpected_ok else "FAIL",
        "Dataset Ready": "YES" if ready else "NO",
    }
    validation["counts"] = {
        "manifest_rows": len(rows),
        "actual_cnf_files": len(actual_paths),
        "verified_sha256_files": verified,
        "parseable_dimacs_files": parseable,
        "missing_manifest_columns": len(missing_columns),
        "duplicate_local_paths": len(duplicate_paths),
        "missing_files": len(missing_files),
        "manifest_rows_without_hash": len(missing_hashes),
        "sha256_mismatches": len(sha_mismatches),
        "parse_failures": len(parse_failures),
        "unexpected_cnf_files": len(unexpected_files),
    }
    validation["failures"] = {
        "missing_manifest_columns": missing_columns,
        "duplicate_local_paths": duplicate_paths,
        "missing_files": missing_files,
        "manifest_rows_without_hash": missing_hashes,
        "sha256_mismatches": sha_mismatches,
        "parse_failures": parse_failures,
        "unexpected_cnf_files": unexpected_files,
    }
    validation["protocol_note"] = (
        "Dataset validation checks only data provenance and integrity. "
        "It does not alter gate parameters, calibration rules, or solver policies."
    )
    (output_dir / "paper69_dataset_validation.json").write_text(json.dumps(validation, indent=2), encoding="utf-8")
    return validation


def print_validation(validation: dict[str, object]) -> None:
    checks = validation.get("checks", {})
    print("Dataset Validation")
    for label in ["Manifest", "CNFs", "SHA256", "Unexpected Files", "Dataset Ready"]:
        value = checks.get(label, "UNKNOWN") if isinstance(checks, dict) else "UNKNOWN"
        print(f"{label + ' ':<18} {value}")
    print(validation.get("status", "VALIDATION UNKNOWN"))


def select_family_balanced_paths(
    paths: list[Path],
    manifest: dict[str, dict[str, str]],
    per_family: int,
) -> list[Path]:
    grouped: dict[str, list[Path]] = {}
    for path in paths:
        rel = str(path.relative_to(ROOT))
        meta = manifest.get(rel, {})
        family = meta.get("family") or infer_family(path)
        grouped.setdefault(family, []).append(path)

    selected: list[Path] = []
    for family in sorted(grouped):
        selected.extend(sorted(grouped[family])[:per_family])
    return sorted(selected)


def discover_instances(
    data_dir: Path,
    manifest_path: Path,
    max_instances: int | None,
    family_balanced: int | None,
) -> list[Instance]:
    paths = cnf_paths(data_dir)
    if not paths:
        return []
    manifest = load_manifest(manifest_path)
    if family_balanced is not None:
        if family_balanced <= 0:
            raise RuntimeError("--family-balanced must be a positive integer.")
        paths = select_family_balanced_paths(paths, manifest, family_balanced)
    if max_instances is not None:
        paths = paths[:max_instances]

    instances: list[Instance] = []
    for idx, path in enumerate(paths):
        rel = str(path.relative_to(ROOT))
        if rel not in manifest:
            raise RuntimeError(
                f"external CNF is not documented in dataset_manifest.csv: {rel}. "
                "Paper 69 requires source, terms, and SHA256 metadata before execution."
            )
        meta = manifest[rel]
        n_vars, clauses = parse_dimacs(path)
        if n_vars <= 0 or not clauses:
            continue
        digest = sha256_file(path)
        manifest_sha = meta.get("sha256", "").strip().lower()
        if not manifest_sha:
            raise RuntimeError(
                f"manifest row for {rel} has no SHA256. "
                "Run the SATLIB stager or fill the hash before executing the benchmark."
            )
        if manifest_sha != digest.lower():
            raise RuntimeError(
                f"SHA256 mismatch for {rel}: manifest={manifest_sha}, computed={digest}. "
                "The benchmark is refused until the manifest and raw file agree."
            )
        dataset_id = meta.get("dataset_id") or f"external_{idx:04d}_{path.stem}"
        instances.append(
            Instance(
                dataset_id=dataset_id,
                family=meta.get("family") or infer_family(path),
                source_name=meta.get("source_name") or "external_user_supplied",
                source_url=meta.get("source_url") or "",
                license_or_terms=meta.get("license_or_terms") or "user_must_verify_original_terms",
                download_tls_mode=meta.get("download_tls_mode") or "not_recorded",
                archive_sha256=meta.get("archive_sha256") or "",
                local_path=rel,
                sha256=digest,
                n_vars=n_vars,
                n_clauses=len(clauses),
                clauses=clauses,
            )
        )
    return instances


def family_counts(instances: list[Instance]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for inst in instances:
        counts[inst.family] = counts.get(inst.family, 0) + 1
    return dict(sorted(counts.items()))


def graph_components(neighbors: list[set[int]]) -> list[int]:
    seen = set()
    sizes = []
    for start in range(len(neighbors)):
        if start in seen:
            continue
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
        sizes.append(size)
    return sizes


def root_supports(instance: Instance) -> dict[str, float]:
    n_vars = instance.n_vars
    clauses = instance.clauses
    deg = np.zeros(n_vars, dtype=float)
    pos = np.zeros(n_vars, dtype=float)
    neg = np.zeros(n_vars, dtype=float)
    neighbors: list[set[int]] = [set() for _ in range(n_vars)]
    for clause in clauses:
        vars_ = sorted({abs(lit) - 1 for lit in clause if 1 <= abs(lit) <= n_vars})
        for lit in clause:
            v = abs(lit) - 1
            if 0 <= v < n_vars:
                deg[v] += 1.0
                if lit > 0:
                    pos[v] += 1.0
                else:
                    neg[v] += 1.0
        for a, b in combinations(vars_, 2):
            neighbors[a].add(b)
            neighbors[b].add(a)

    graph_deg = np.array([len(x) for x in neighbors], dtype=float)
    degree_cv = float(np.std(deg) / (np.mean(deg) + EPS))
    literal_imbalance = float(np.mean(np.abs(pos - neg) / (pos + neg + EPS)))

    clustering = []
    for v, nb_set in enumerate(neighbors):
        nb = list(nb_set)
        if len(nb) < 2:
            clustering.append(0.0)
            continue
        possible = len(nb) * (len(nb) - 1) / 2.0
        actual = sum(1 for a, b in combinations(nb, 2) if b in neighbors[a])
        clustering.append(actual / (possible + EPS))
    clustering_mean = float(np.mean(clustering)) if clustering else 0.0
    components = graph_components(neighbors)
    largest_component_frac = max(components, default=0) / max(1, n_vars)
    component_frac = len(components) / max(1, n_vars)

    alpha = instance.n_clauses / max(1, n_vars)
    # SAT-specific operational supports, frozen for Paper 69.
    H = 1.0 / (1.0 + literal_imbalance)
    B = 1.0 / (1.0 + degree_cv)
    S = math.exp(-0.5 * ((alpha - 4.26) / 1.25) ** 2)
    V = max(0.0, min(1.0, 0.65 * largest_component_frac + 0.35 * clustering_mean))
    R_resp = (max(EPS, H) * max(EPS, B) * max(EPS, V)) ** (1.0 / 3.0)
    R_rob = min(R_resp, (max(EPS, H) * max(EPS, B) * max(EPS, S) * max(EPS, V)) ** 0.25)
    return {
        "alpha": float(alpha),
        "H": float(H),
        "B": float(B),
        "S": float(S),
        "V": float(V),
        "R_resp": float(R_resp),
        "R_rob": float(R_rob),
        "degree_cv": degree_cv,
        "literal_imbalance": literal_imbalance,
        "clustering_mean": clustering_mean,
        "largest_component_frac": float(largest_component_frac),
        "component_frac": float(component_frac),
    }


def median_abs_dev(values: np.ndarray) -> float:
    med = float(np.median(values))
    return float(np.median(np.abs(values - med)))


def build_features(instances: list[Instance]) -> list[InstanceFeatures]:
    base = []
    for inst in instances:
        supports = root_supports(inst)
        split = "calibration" if stable_fold(inst.dataset_id) == 0 else "test"
        base.append((inst, split, supports))

    cal_values = np.array([x[2]["R_rob"] for x in base if x[1] == "calibration"], dtype=float)
    if len(cal_values) < 3:
        cal_values = np.array([x[2]["R_rob"] for x in base], dtype=float)
    center = float(np.median(cal_values)) if len(cal_values) else 0.0
    mad = median_abs_dev(cal_values) if len(cal_values) else 1.0
    scale = mad + EPS

    rng = random.Random(SEED + 991)
    raw_gates = []
    for inst, split, supports in base:
        z_R = SAT_POLARITY * (supports["R_rob"] - center) / scale
        gate_signal = z_R  # z_t and z_x unavailable for static CNF; frozen to zero.
        G_gate = 1.0 / (1.0 + math.exp(-gate_signal))
        raw_gates.append(G_gate)
    shuffled = raw_gates[:]
    rng.shuffle(shuffled)

    features: list[InstanceFeatures] = []
    for (inst, split, supports), G_shuf in zip(base, shuffled):
        z_R = SAT_POLARITY * (supports["R_rob"] - center) / scale
        gate_signal = z_R
        G_gate = 1.0 / (1.0 + math.exp(-gate_signal))
        features.append(
            InstanceFeatures(
                dataset_id=inst.dataset_id,
                family=inst.family,
                split=split,
                n_vars=inst.n_vars,
                n_clauses=inst.n_clauses,
                alpha=supports["alpha"],
                H=supports["H"],
                B=supports["B"],
                S=supports["S"],
                V=supports["V"],
                R_resp=supports["R_resp"],
                R_rob=supports["R_rob"],
                cal_center_R=center,
                cal_mad_R=mad,
                z_R=float(z_R),
                z_t=0.0,
                z_x=0.0,
                gate_signal=float(gate_signal),
                G_gate=float(G_gate),
                G_hard=int(G_gate >= 0.5),
                G_shuffled=float(G_shuf),
                degree_cv=supports["degree_cv"],
                literal_imbalance=supports["literal_imbalance"],
                clustering_mean=supports["clustering_mean"],
                largest_component_frac=supports["largest_component_frac"],
                component_frac=supports["component_frac"],
            )
        )
    return features


class Paper69Solver(CDCLSolver):
    def literal_scores_jw_moms(self) -> dict[int, dict[str, float]]:
        unassigned = [v for v in range(1, self.n_vars + 1) if v not in self.assignment]
        active = [c for c in self.active_clauses() if c]
        out = {
            v: {
                "jw_pos": 0.0,
                "jw_neg": 0.0,
                "moms_pos": 0.0,
                "moms_neg": 0.0,
            }
            for v in unassigned
        }
        if not active:
            return out
        min_len = min(len(c) for c in active)
        for clause in active:
            jw_weight = 2.0 ** (-len(clause))
            moms_weight = 1.0 if len(clause) == min_len else 0.0
            for lit in clause:
                v = abs(lit)
                if v not in out:
                    continue
                if lit > 0:
                    out[v]["jw_pos"] += jw_weight
                    out[v]["moms_pos"] += moms_weight
                else:
                    out[v]["jw_neg"] += jw_weight
                    out[v]["moms_neg"] += moms_weight
        return out

    def choose_decision(self) -> int | None:  # noqa: C901 - policy dispatch is explicit.
        unassigned = [v for v in range(1, self.n_vars + 1) if v not in self.assignment]
        if not unassigned:
            return None
        scores = self.score_candidates()
        lit_scores = self.literal_scores_jw_moms()

        if self.policy == "vsids":
            var = max(unassigned, key=lambda v: (self.vsids[v], scores.get(v, {}).get("degree", 0.0), -v))
            polarity = self.phase.get(var, True)
        elif self.policy == "progress_only":
            var = max(unassigned, key=lambda v: (scores[v]["progress"], self.vsids[v], -v))
            polarity = scores[var]["polarity_true_score"] >= scores[var]["polarity_false_score"]
        elif self.policy == "score_with_R":
            var = self.softmax_choice({v: scores[v]["mode_a"] for v in unassigned})
            polarity = scores[var]["polarity_true_score"] >= scores[var]["polarity_false_score"]
        elif self.policy in {"gate_v17", "gate_shuffled_R"}:
            var = self.softmax_choice({v: scores[v]["progress"] - self.structure_gate * self.tau * scores[v]["h_maat"] for v in unassigned})
            polarity = scores[var]["polarity_true_score"] >= scores[var]["polarity_false_score"]
        elif self.policy == "jeroslow_wang":
            var = max(
                unassigned,
                key=lambda v: (
                    lit_scores[v]["jw_pos"] + lit_scores[v]["jw_neg"],
                    self.vsids[v],
                    -v,
                ),
            )
            polarity = lit_scores[var]["jw_pos"] >= lit_scores[var]["jw_neg"]
        elif self.policy == "moms":
            def moms_score(v: int) -> tuple[float, float, float, int]:
                p = lit_scores[v]["moms_pos"]
                n = lit_scores[v]["moms_neg"]
                return (p + n + 2.0 * p * n, p + n, self.vsids[v], -v)

            var = max(unassigned, key=moms_score)
            polarity = lit_scores[var]["moms_pos"] >= lit_scores[var]["moms_neg"]
        else:
            raise ValueError(f"unknown Paper 69 policy: {self.policy}")

        if var in self.phase and self.rng.random() < 0.25 and self.policy in {"vsids", "progress_only", "score_with_R", "gate_v17", "gate_shuffled_R"}:
            polarity = self.phase[var]
        return var if polarity else -var


def solve_instance(inst: Instance, feat: InstanceFeatures, policy: str, conflict_budget: int) -> SolveRecord:
    gate_value = feat.G_gate if policy != "gate_shuffled_R" else feat.G_shuffled
    start = time.perf_counter()
    solver = Paper69Solver(
        inst.clauses,
        inst.n_vars,
        policy=policy,
        seed=SEED + 1009 * stable_fold(inst.dataset_id, 1_000_003) + 37 * POLICIES.index(policy),
        structure_gate=gate_value,
        conflict_budget=conflict_budget,
        beta=SOFTMAX_BETA,
        tau=TAU_STRUCTURAL,
    )
    sat, solved = solver.solve()
    runtime_s = time.perf_counter() - start
    timeout = int(not solved)
    compute_cost = (
        solver.conflicts
        + 0.30 * solver.decisions
        + 0.004 * solver.propagations
        + timeout * conflict_budget
    )
    return SolveRecord(
        dataset_id=inst.dataset_id,
        family=inst.family,
        split=feat.split,
        policy=policy,
        n_vars=inst.n_vars,
        n_clauses=inst.n_clauses,
        alpha=feat.alpha,
        solved=int(solved),
        sat=int(bool(sat)),
        timeout=timeout,
        decisions=int(solver.decisions),
        propagations=int(solver.propagations),
        conflicts=int(solver.conflicts),
        learned_clauses=int(solver.learned_clauses),
        max_level=int(solver.max_level),
        gate_value=float(gate_value),
        runtime_s=float(runtime_s),
        compute_cost=float(compute_cost),
    )


def paired_bootstrap(values: np.ndarray, n_boot: int = 1000) -> dict[str, float]:
    if len(values) == 0:
        return {"mean": math.nan, "ci95_low": math.nan, "ci95_high": math.nan}
    rng = np.random.default_rng(SEED + 222)
    means = []
    for _ in range(n_boot):
        sample = rng.choice(values, size=len(values), replace=True)
        means.append(float(np.mean(sample)))
    lo, hi = np.quantile(means, [0.025, 0.975])
    return {"mean": float(np.mean(values)), "ci95_low": float(lo), "ci95_high": float(hi)}


def summarize_results(df: pd.DataFrame, features: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    test = df[df["split"] == "test"].copy()
    if test.empty:
        return pd.DataFrame(), pd.DataFrame()
    pivot = test.pivot_table(index=["dataset_id", "family"], columns="policy", values="compute_cost")
    rows = []
    for policy in POLICIES:
        if policy not in pivot.columns:
            continue
        vals = pivot[policy].dropna()
        rows.append(
            {
                "policy": policy,
                "n_test": int(vals.shape[0]),
                "mean_compute_cost": float(vals.mean()),
                "median_compute_cost": float(vals.median()),
                "solve_rate": float(test[test["policy"] == policy]["solved"].mean()),
                "timeout_rate": float(test[test["policy"] == policy]["timeout"].mean()),
            }
        )

    comparisons = []
    baselines = ["score_with_R", "vsids", "moms", "jeroslow_wang", "progress_only", "gate_shuffled_R"]
    for baseline in baselines:
        if "gate_v17" not in pivot.columns or baseline not in pivot.columns:
            continue
        # Positive delta means gate improves by reducing compute cost.
        delta = (pivot[baseline] - pivot["gate_v17"]).dropna().to_numpy(dtype=float)
        boot = paired_bootstrap(delta)
        comparisons.append(
            {
                "comparison": f"gate_v17_vs_{baseline}",
                "n_pairs": int(len(delta)),
                "delta_utility_mean": boot["mean"],
                "ci95_low": boot["ci95_low"],
                "ci95_high": boot["ci95_high"],
                "wins": int((delta > 0).sum()),
                "losses": int((delta < 0).sum()),
                "ties": int((delta == 0).sum()),
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(comparisons)


def write_empty_outputs(reason: str, outdir: Path, data_dir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    empty_manifest = pd.DataFrame(
        columns=[
            "dataset_id",
            "family",
            "source_name",
            "source_url",
            "license_or_terms",
            "download_tls_mode",
            "archive_sha256",
            "local_path",
            "sha256",
            "n_vars",
            "n_clauses",
        ]
    )
    empty_manifest.to_csv(outdir / "paper69_dataset_manifest_detected.csv", index=False)
    pd.DataFrame(
        columns=[
            "dataset_id",
            "family",
            "split",
            "policy",
            "n_vars",
            "n_clauses",
            "alpha",
            "solved",
            "sat",
            "timeout",
            "decisions",
            "propagations",
            "conflicts",
            "learned_clauses",
            "max_level",
            "gate_value",
            "runtime_s",
            "compute_cost",
        ]
    ).to_csv(outdir / "paper69_solve_records.csv", index=False)
    pd.DataFrame(
        columns=[
            "policy",
            "n_test",
            "mean_compute_cost",
            "median_compute_cost",
            "solve_rate",
            "timeout_rate",
        ]
    ).to_csv(outdir / "paper69_policy_summary.csv", index=False)
    pd.DataFrame(
        columns=[
            "comparison",
            "n_pairs",
            "delta_utility_mean",
            "ci95_low",
            "ci95_high",
            "wins",
            "losses",
            "ties",
        ]
    ).to_csv(outdir / "paper69_gate_comparisons.csv", index=False)
    summary = {
        "paper": 69,
        "title": "Gate Challenge I: SAT/CDCL",
        "status": "not_executed",
        "reason": reason,
        "data_dir": str(data_dir),
        "paper68_preregistration_doi": PAPER68_DOI,
        "required_action": "Place public external DIMACS CNF files under data_external/ and document sources in dataset_manifest.csv.",
    }
    (outdir / "paper69_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (outdir / "paper69_run_log.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


def plot_policy_summary(summary_df: pd.DataFrame, outdir: Path) -> None:
    if summary_df.empty:
        return
    order = summary_df.sort_values("median_compute_cost")["policy"].tolist()
    data = summary_df.set_index("policy").loc[order]
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    ax.bar(data.index, data["median_compute_cost"], color="#2a9d8f", edgecolor="#1f2933")
    ax.set_ylabel("median compute cost (test split)")
    ax.set_title("Paper 69 SAT Gate Challenge: policy comparison")
    ax.tick_params(axis="x", rotation=25)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir / "fig1_policy_compute_cost.png", dpi=180)
    plt.close(fig)


def plot_gate_comparisons(comp_df: pd.DataFrame, outdir: Path) -> None:
    if comp_df.empty:
        return
    fig, ax = plt.subplots(figsize=(8.6, 4.7))
    x = np.arange(len(comp_df))
    means = comp_df["delta_utility_mean"].to_numpy(dtype=float)
    lows = comp_df["ci95_low"].to_numpy(dtype=float)
    highs = comp_df["ci95_high"].to_numpy(dtype=float)
    yerr = np.vstack([means - lows, highs - means])
    ax.axhline(0, color="#555", lw=1)
    ax.bar(x, means, color="#264653", edgecolor="#1f2933")
    ax.errorbar(x, means, yerr=yerr, fmt="none", ecolor="#e76f51", capsize=4, lw=1.4)
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("gate_v17_vs_", "vs ") for c in comp_df["comparison"]], rotation=25, ha="right")
    ax.set_ylabel("paired utility delta (baseline cost - gate cost)")
    ax.set_title("Gate-vs-score and gate-vs-baseline comparisons")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir / "fig2_gate_comparisons.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--manifest", type=Path, default=ROOT / "dataset_manifest.csv")
    parser.add_argument("--max-instances", type=int, default=None)
    parser.add_argument(
        "--family-balanced",
        type=int,
        default=None,
        metavar="N",
        help="Select up to N deterministic CNFs per benchmark family before solving.",
    )
    parser.add_argument("--conflict-budget", type=int, default=5000)
    parser.add_argument("--output-dir", type=Path, default=OUTDIR)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate dataset_manifest.csv, CNF existence, SHA256 hashes, and unexpected CNFs, then exit without solving.",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / ".mplconfig").mkdir(parents=True, exist_ok=True)

    if args.validate_only:
        validation = validate_dataset(args.data_dir, args.manifest, args.output_dir)
        print_validation(validation)
        return

    try:
        instances = discover_instances(args.data_dir, args.manifest, args.max_instances, args.family_balanced)
    except (FileNotFoundError, RuntimeError) as exc:
        write_empty_outputs(f"refused_manifest_verification: {exc}", args.output_dir, args.data_dir)
        print(json.dumps(json.loads((args.output_dir / "paper69_summary.json").read_text()), indent=2))
        return
    if not instances:
        write_empty_outputs("no_external_dimacs_cnf_files_found", args.output_dir, args.data_dir)
        print(json.dumps(json.loads((args.output_dir / "paper69_summary.json").read_text()), indent=2))
        return

    features = build_features(instances)
    by_id = {f.dataset_id: f for f in features}

    manifest_rows = []
    for inst in instances:
        manifest_rows.append(
            {
                "dataset_id": inst.dataset_id,
                "family": inst.family,
                "source_name": inst.source_name,
                "source_url": inst.source_url,
                "license_or_terms": inst.license_or_terms,
                "download_tls_mode": inst.download_tls_mode,
                "archive_sha256": inst.archive_sha256,
                "local_path": inst.local_path,
                "sha256": inst.sha256,
                "n_vars": inst.n_vars,
                "n_clauses": inst.n_clauses,
            }
        )
    pd.DataFrame(manifest_rows).to_csv(args.output_dir / "paper69_dataset_manifest_detected.csv", index=False)
    pd.DataFrame([asdict(f) for f in features]).to_csv(args.output_dir / "paper69_gate_features.csv", index=False)

    records: list[SolveRecord] = []
    for idx, inst in enumerate(instances):
        feat = by_id[inst.dataset_id]
        for policy in POLICIES:
            records.append(solve_instance(inst, feat, policy, args.conflict_budget))
        print(f"[{idx+1}/{len(instances)}] solved {inst.dataset_id}")

    solve_df = pd.DataFrame([asdict(r) for r in records])
    solve_df.to_csv(args.output_dir / "paper69_solve_records.csv", index=False)

    summary_df, comp_df = summarize_results(solve_df, pd.DataFrame([asdict(f) for f in features]))
    summary_df.to_csv(args.output_dir / "paper69_policy_summary.csv", index=False)
    comp_df.to_csv(args.output_dir / "paper69_gate_comparisons.csv", index=False)

    plot_policy_summary(summary_df, args.output_dir)
    plot_gate_comparisons(comp_df, args.output_dir)

    fair = comp_df[comp_df["comparison"] == "gate_v17_vs_score_with_R"]
    fair_result = fair.iloc[0].to_dict() if not fair.empty else {}
    best_classical = {}
    classical_rows = comp_df[comp_df["comparison"].isin([f"gate_v17_vs_{b}" for b in ["vsids", "moms", "jeroslow_wang"]])]
    if not classical_rows.empty:
        best_classical = classical_rows.sort_values("delta_utility_mean").iloc[0].to_dict()
    execution_scope = (
        "family_balanced_smoke_test"
        if args.family_balanced is not None
        else "smoke_test"
        if args.max_instances is not None or args.conflict_budget != 5000
        else "primary_candidate"
    )

    summary = {
        "paper": 69,
        "title": "Gate Challenge I: SAT/CDCL",
        "status": "smoke_executed" if "smoke_test" in execution_scope else "executed",
        "execution_scope": execution_scope,
        "paper68_preregistration_doi": PAPER68_DOI,
        "n_instances": len(instances),
        "n_test_instances": int(sum(1 for f in features if f.split == "test")),
        "n_calibration_instances": int(sum(1 for f in features if f.split == "calibration")),
        "max_instances": args.max_instances,
        "family_balanced_per_family": args.family_balanced,
        "family_counts": family_counts(instances),
        "conflict_budget": args.conflict_budget,
        "policies": POLICIES,
        "manifest_verification": "required_sha256_match_before_solver_execution",
        "download_tls_mode_usage": "audit_metadata_only_not_used_by_gate_or_solver_policy",
        "primary_fair_comparison": fair_result,
        "best_classical_comparison_proxy": best_classical,
        "interpretation_rule": "Positive support requires gate_v17_vs_score_with_R CI lower bound > 0; classical defeat requires positive CI against strongest classical baseline.",
    }
    (args.output_dir / "paper69_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (args.output_dir / "paper69_run_log.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
