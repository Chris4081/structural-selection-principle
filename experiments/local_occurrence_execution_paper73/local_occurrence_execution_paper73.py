#!/usr/bin/env python3
"""Paper 73: Local Occurrence Hypothesis execution.

This script executes the Paper-72 preregistered SAT-specific protocol. It does
not tune parameters, alter thresholds, redefine L_star, or modify the frozen
Paper-69 gate policy. The only policy extension is the Paper-72 Gate+L wrapper:

    G_72 = sigmoid(logit(G_gate) + L_star)

where L_star is computed by calibration-only residualization of occurrence and
degree channels against H,B,S,V,R_rob,G_gate.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import signal
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
_repo_from_layout = ROOT.parents[1]
if (_repo_from_layout / "experiments" / "gate_challenge_sat_paper69").exists():
    REPO = _repo_from_layout
elif os.environ.get("MAAT_SSL_REPO"):
    REPO = Path(os.environ["MAAT_SSL_REPO"])
else:
    raise RuntimeError(
        "Cannot locate the structural-selection-principle repo. Run this script "
        "from the in-repo experiment folder or set MAAT_SSL_REPO."
    )
P69_DIR = REPO / "experiments" / "gate_challenge_sat_paper69"
P72_PROTOCOL_PORTABLE = (
    "../local_occurrence_hypothesis_paper72/"
    "outputs_paper72_local_occurrence_protocol/paper72_protocol.json"
)
P72_PROTOCOL = (
    REPO
    / "experiments"
    / "local_occurrence_hypothesis_paper72"
    / "outputs_paper72_local_occurrence_protocol"
    / "paper72_protocol.json"
)

sys.path.insert(0, str(P69_DIR))
import gate_challenge_sat_paper69 as p69  # type: ignore  # noqa: E402


_CACHE_ROOT = Path(os.environ.get("TMPDIR", "/tmp")) / "maat_paper73_cache"
(_CACHE_ROOT / "matplotlib").mkdir(parents=True, exist_ok=True)
(_CACHE_ROOT / "xdg").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE_ROOT / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE_ROOT / "xdg"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


EPS = 1.0e-12
SEED = 72072
OUTDIR = ROOT / "outputs_paper73_local_occurrence_execution"
DATA_DIR = ROOT / "data_external"
MANIFEST = ROOT / "dataset_manifest.csv"

POLICIES = [
    "vsids",
    "moms",
    "jeroslow_wang",
    "score_with_R",
    "frozen_gate_v17",
    "gate_plus_L",
    "gate_shuffled_L",
]

COMPARATORS = ["frozen_gate_v17", "score_with_R", "moms", "jeroslow_wang", "vsids", "gate_shuffled_L"]


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
class Feature73:
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
    G_gate: float
    G_72: float
    G_shuffled_L: float
    L_star: float
    L_star_shuffled: float
    O_raw: float
    D_raw: float
    O_hat: float
    D_hat: float
    O_res: float
    D_res: float
    z_O_res: float
    z_D_res: float
    cal_center_R: float
    cal_mad_R: float
    z_R: float
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
    runtime_budget_exceeded: int
    decisions: int
    propagations: int
    conflicts: int
    learned_clauses: int
    max_level: int
    gate_value: float
    runtime_s: float
    compute_cost: float


class RuntimeBudgetExceeded(Exception):
    pass


def load_protocol() -> dict:
    return json.loads(P72_PROTOCOL.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def local_rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def cnf_paths(data_dir: Path) -> list[Path]:
    patterns = ["*.cnf", "*.dimacs", "*.cnf.gz", "*.cnf.bz2", "*.cnf.xz", "*.cnf.lzma"]
    paths: list[Path] = []
    for pattern in patterns:
        paths.extend(data_dir.rglob(pattern))
    return sorted(set(paths))


def load_manifest(manifest_path: Path) -> dict[str, dict[str, str]]:
    if not manifest_path.exists():
        raise FileNotFoundError(f"dataset manifest is required: {manifest_path}")
    with manifest_path.open("r", newline="", encoding="utf-8") as handle:
        return {row["local_path"]: row for row in csv.DictReader(handle)}


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
        "paper": 73,
        "title": "Local Occurrence Hypothesis Execution",
        "validation_type": "dataset_manifest_sha256_gatekeeper",
        "paper72_protocol": P72_PROTOCOL_PORTABLE,
        "series_ii_doi": "10.5281/zenodo.21062386",
        "manifest_path": "./dataset_manifest.csv",
        "data_dir": "./data_external",
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
            "Parseability": "NOT_RUN",
            "Unexpected Files": "NOT_RUN",
            "Dataset Ready": "NO",
        }
        validation["failures"] = {"manifest": [f"missing manifest: {manifest_path.name}"]}
        (output_dir / "paper73_dataset_validation.json").write_text(json.dumps(validation, indent=2), encoding="utf-8")
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
    actual_paths = {local_rel(path) for path in cnf_paths(data_dir)}

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
            n_vars, clauses = p69.parse_dimacs(path)
            if n_vars <= 0 or not clauses:
                parse_failures.append({"local_path": rel, "error": "no_variables_or_clauses_detected"})
            else:
                parseable += 1
        except Exception as exc:  # noqa: BLE001
            parse_failures.append({"local_path": rel, "error": str(exc)})

    unexpected_files = sorted(actual_paths - manifest_paths)
    untracked_manifest_rows = sorted(manifest_paths - actual_paths)
    manifest_ok = bool(rows) and not missing_columns and not duplicate_paths
    cnfs_ok = bool(actual_paths) and not missing_files and not untracked_manifest_rows
    sha_ok = not missing_hashes and not sha_mismatches and verified == len(rows)
    parse_ok = not parse_failures and parseable == len(rows)
    unexpected_ok = not unexpected_files
    ready = manifest_ok and cnfs_ok and sha_ok and parse_ok and unexpected_ok

    validation["status"] = "VALIDATION PASS" if ready else "VALIDATION FAIL"
    validation["dataset_ready"] = ready
    validation["checks"] = {
        "Manifest": "PASS" if manifest_ok else "FAIL",
        "CNFs": "PASS" if cnfs_ok else "FAIL",
        "SHA256": "PASS" if sha_ok else "FAIL",
        "Parseability": "PASS" if parse_ok else "FAIL",
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
        "Dataset validation checks provenance and integrity only. It does not "
        "alter Gate+L, residualization, calibration, budgets, baselines, or interpretation."
    )
    (output_dir / "paper73_dataset_validation.json").write_text(json.dumps(validation, indent=2), encoding="utf-8")
    return validation


def print_validation(validation: dict[str, object]) -> None:
    checks = validation.get("checks", {})
    print("Dataset Validation")
    for label in ["Manifest", "CNFs", "SHA256", "Parseability", "Unexpected Files", "Dataset Ready"]:
        value = checks.get(label, "UNKNOWN") if isinstance(checks, dict) else "UNKNOWN"
        print(f"{label + ' ':<18} {value}")
    print(validation.get("status", "VALIDATION UNKNOWN"))


def select_family_balanced_paths(paths: list[Path], manifest: dict[str, dict[str, str]], per_family: int) -> list[Path]:
    grouped: dict[str, list[Path]] = {}
    for path in paths:
        rel = local_rel(path)
        family = manifest.get(rel, {}).get("family", path.parent.name)
        grouped.setdefault(family, []).append(path)
    selected: list[Path] = []
    for family in sorted(grouped):
        selected.extend(sorted(grouped[family])[:per_family])
    return sorted(selected)


def discover_instances(data_dir: Path, manifest_path: Path, per_family: int) -> list[Instance]:
    manifest = load_manifest(manifest_path)
    paths = select_family_balanced_paths(cnf_paths(data_dir), manifest, per_family)
    instances: list[Instance] = []
    for path in paths:
        rel = local_rel(path)
        if rel not in manifest:
            raise RuntimeError(f"CNF not documented in manifest: {rel}")
        meta = manifest[rel]
        digest = sha256_file(path)
        expected = meta.get("sha256", "").strip().lower()
        if not expected or expected != digest.lower():
            raise RuntimeError(f"SHA256 mismatch or missing hash for {rel}")
        n_vars, clauses = p69.parse_dimacs(path)
        if n_vars <= 0 or not clauses:
            continue
        instances.append(
            Instance(
                dataset_id=meta.get("dataset_id") or path.stem,
                family=meta.get("family") or path.parent.name,
                source_name=meta.get("source_name") or "external",
                source_url=meta.get("source_url") or "",
                license_or_terms=meta.get("license_or_terms") or "check_original_terms",
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
    out: dict[str, int] = {}
    for inst in instances:
        out[inst.family] = out.get(inst.family, 0) + 1
    return dict(sorted(out.items()))


def stable_int(text: str) -> int:
    return int(hashlib.sha256(text.encode("utf-8")).hexdigest()[:12], 16)


def split_by_family(instances: list[Instance], calibration_fraction: float) -> dict[str, str]:
    split: dict[str, str] = {}
    grouped: dict[str, list[Instance]] = {}
    for inst in instances:
        grouped.setdefault(inst.family, []).append(inst)
    for family, rows in grouped.items():
        rows = sorted(rows, key=lambda x: x.dataset_id)
        rng = random.Random(SEED + stable_int(family))
        rng.shuffle(rows)
        n_cal = int(round(len(rows) * calibration_fraction))
        if len(rows) >= 2:
            n_cal = max(1, min(len(rows) - 1, n_cal))
        cal_ids = {inst.dataset_id for inst in rows[:n_cal]}
        for inst in rows:
            split[inst.dataset_id] = "calibration" if inst.dataset_id in cal_ids else "test"
    return split


def median_abs_dev(values: np.ndarray) -> float:
    med = float(np.median(values))
    return float(np.median(np.abs(values - med)))


def z_stats(values: np.ndarray) -> tuple[float, float]:
    mean = float(np.mean(values)) if len(values) else 0.0
    std = float(np.std(values)) if len(values) else 1.0
    return mean, std if std > EPS else 1.0


def raw_occurrence_channel(inst: Instance) -> float:
    occ = np.zeros(inst.n_vars, dtype=float)
    for clause in inst.clauses:
        for lit in clause:
            v = abs(lit) - 1
            if 0 <= v < inst.n_vars:
                occ[v] += 1.0
    return float(np.std(occ) / (np.mean(occ) + EPS))


def sigmoid(x: float) -> float:
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def logit(p: float) -> float:
    q = float(np.clip(p, 1.0e-6, 1.0 - 1.0e-6))
    return math.log(q / (1.0 - q))


def build_features(instances: list[Instance], protocol: dict) -> list[Feature73]:
    split_map = split_by_family(instances, protocol["dataset_rules"]["calibration_fraction_per_family"])
    support_rows = []
    for inst in instances:
        supports = p69.root_supports(
            p69.Instance(
                dataset_id=inst.dataset_id,
                family=inst.family,
                source_name=inst.source_name,
                source_url=inst.source_url,
                license_or_terms=inst.license_or_terms,
                download_tls_mode=inst.download_tls_mode,
                archive_sha256=inst.archive_sha256,
                local_path=inst.local_path,
                sha256=inst.sha256,
                n_vars=inst.n_vars,
                n_clauses=inst.n_clauses,
                clauses=inst.clauses,
            )
        )
        support_rows.append((inst, split_map[inst.dataset_id], supports))

    cal_support_rows = [row for row in support_rows if row[1] == "calibration"]
    cal_R = np.array([row[2]["R_rob"] for row in cal_support_rows], dtype=float)
    center_R = float(np.median(cal_R)) if len(cal_R) else 0.0
    mad_R = median_abs_dev(cal_R) if len(cal_R) else 1.0
    scale_R = mad_R + EPS

    support_names = protocol["frozen_support_basis"]
    raw_matrix = []
    O_raw_all = []
    D_raw_all = []
    for inst, split, supports in support_rows:
        z_R = (supports["R_rob"] - center_R) / scale_R
        G_gate = sigmoid(z_R)
        row_values = {
            "H": supports["H"],
            "B": supports["B"],
            "S": supports["S"],
            "V": supports["V"],
            "R_rob": supports["R_rob"],
            "G_gate": G_gate,
        }
        raw_matrix.append([row_values[name] for name in support_names])
        O_raw_all.append(raw_occurrence_channel(inst))
        D_raw_all.append(float(supports["degree_cv"]))

    raw_matrix_np = np.array(raw_matrix, dtype=float)
    O_raw_np = np.array(O_raw_all, dtype=float)
    D_raw_np = np.array(D_raw_all, dtype=float)
    cal_mask = np.array([split == "calibration" for _, split, _ in support_rows], dtype=bool)
    cal_X_raw = raw_matrix_np[cal_mask]
    means = np.mean(cal_X_raw, axis=0)
    stds = np.std(cal_X_raw, axis=0)
    stds = np.where(stds <= EPS, 1.0, stds)
    X_z = (raw_matrix_np - means) / stds
    X_cal = X_z[cal_mask]
    design_cal = np.column_stack([np.ones(X_cal.shape[0]), X_cal])
    design_all = np.column_stack([np.ones(X_z.shape[0]), X_z])

    coef_O, *_ = np.linalg.lstsq(design_cal, O_raw_np[cal_mask], rcond=None)
    coef_D, *_ = np.linalg.lstsq(design_cal, D_raw_np[cal_mask], rcond=None)
    O_hat = design_all @ coef_O
    D_hat = design_all @ coef_D
    O_res = O_raw_np - O_hat
    D_res = D_raw_np - D_hat

    O_res_mean, O_res_std = z_stats(O_res[cal_mask])
    D_res_mean, D_res_std = z_stats(D_res[cal_mask])
    z_O_res = (O_res - O_res_mean) / O_res_std
    z_D_res = (D_res - D_res_mean) / D_res_std
    L_star = 0.5 * (z_D_res - z_O_res)

    shuffled_L = L_star.copy()
    grouped_idx: dict[str, list[int]] = {}
    for idx, (inst, _, _) in enumerate(support_rows):
        grouped_idx.setdefault(inst.family, []).append(idx)
    rng = random.Random(SEED + 773)
    for family, idxs in grouped_idx.items():
        vals = [float(shuffled_L[i]) for i in idxs]
        rng.shuffle(vals)
        for i, val in zip(idxs, vals):
            shuffled_L[i] = val

    features: list[Feature73] = []
    for idx, (inst, split, supports) in enumerate(support_rows):
        z_R = (supports["R_rob"] - center_R) / scale_R
        G_gate = sigmoid(z_R)
        G_72 = sigmoid(logit(G_gate) + float(L_star[idx]))
        G_shuf = sigmoid(logit(G_gate) + float(shuffled_L[idx]))
        features.append(
            Feature73(
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
                G_gate=float(G_gate),
                G_72=float(G_72),
                G_shuffled_L=float(G_shuf),
                L_star=float(L_star[idx]),
                L_star_shuffled=float(shuffled_L[idx]),
                O_raw=float(O_raw_np[idx]),
                D_raw=float(D_raw_np[idx]),
                O_hat=float(O_hat[idx]),
                D_hat=float(D_hat[idx]),
                O_res=float(O_res[idx]),
                D_res=float(D_res[idx]),
                z_O_res=float(z_O_res[idx]),
                z_D_res=float(z_D_res[idx]),
                cal_center_R=center_R,
                cal_mad_R=mad_R,
                z_R=float(z_R),
                degree_cv=supports["degree_cv"],
                literal_imbalance=supports["literal_imbalance"],
                clustering_mean=supports["clustering_mean"],
                largest_component_frac=supports["largest_component_frac"],
                component_frac=supports["component_frac"],
            )
        )
    return features


class Paper73Solver(p69.Paper69Solver):
    def choose_decision(self) -> int | None:  # noqa: C901
        unassigned = [v for v in range(1, self.n_vars + 1) if v not in self.assignment]
        if not unassigned:
            return None
        scores = self.score_candidates()
        lit_scores = self.literal_scores_jw_moms()

        if self.policy == "vsids":
            var = max(unassigned, key=lambda v: (self.vsids[v], scores.get(v, {}).get("degree", 0.0), -v))
            polarity = self.phase.get(var, True)
        elif self.policy == "score_with_R":
            var = self.softmax_choice({v: scores[v]["mode_a"] for v in unassigned})
            polarity = scores[var]["polarity_true_score"] >= scores[var]["polarity_false_score"]
        elif self.policy in {"frozen_gate_v17", "gate_plus_L", "gate_shuffled_L"}:
            var = self.softmax_choice({v: scores[v]["progress"] - self.structure_gate * self.tau * scores[v]["h_maat"] for v in unassigned})
            polarity = scores[var]["polarity_true_score"] >= scores[var]["polarity_false_score"]
        elif self.policy == "jeroslow_wang":
            var = max(
                unassigned,
                key=lambda v: (lit_scores[v]["jw_pos"] + lit_scores[v]["jw_neg"], self.vsids[v], -v),
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
            raise ValueError(f"unknown Paper 73 policy: {self.policy}")

        if var in self.phase and self.rng.random() < 0.25 and self.policy in {"vsids", "score_with_R", "frozen_gate_v17", "gate_plus_L", "gate_shuffled_L"}:
            polarity = self.phase[var]
        return var if polarity else -var


def gate_for_policy(feat: Feature73, policy: str) -> float:
    if policy == "frozen_gate_v17":
        return feat.G_gate
    if policy == "gate_plus_L":
        return feat.G_72
    if policy == "gate_shuffled_L":
        return feat.G_shuffled_L
    return feat.G_gate


def solve_instance(inst: Instance, feat: Feature73, policy: str, conflict_budget: int, runtime_budget: int) -> SolveRecord:
    gate_value = gate_for_policy(feat, policy)
    start = time.perf_counter()
    solver = Paper73Solver(
        inst.clauses,
        inst.n_vars,
        policy=policy,
        seed=SEED + 1009 * p69.stable_fold(inst.dataset_id, 1_000_003) + 37 * POLICIES.index(policy),
        structure_gate=gate_value,
        conflict_budget=conflict_budget,
        beta=p69.SOFTMAX_BETA,
        tau=p69.TAU_STRUCTURAL,
    )
    runtime_budget_exceeded = 0
    old_handler = signal.getsignal(signal.SIGALRM)

    def _timeout_handler(signum, frame):  # noqa: ANN001
        raise RuntimeBudgetExceeded()

    signal.signal(signal.SIGALRM, _timeout_handler)
    signal.setitimer(signal.ITIMER_REAL, runtime_budget)
    try:
        sat, solved = solver.solve()
    except RuntimeBudgetExceeded:
        sat, solved = False, False
        runtime_budget_exceeded = 1
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, old_handler)
    runtime_s = time.perf_counter() - start
    timeout = int(not solved or runtime_budget_exceeded)
    compute_cost = solver.conflicts + 0.30 * solver.decisions + 0.004 * solver.propagations + timeout * conflict_budget
    return SolveRecord(
        dataset_id=inst.dataset_id,
        family=inst.family,
        split=feat.split,
        policy=policy,
        n_vars=inst.n_vars,
        n_clauses=inst.n_clauses,
        alpha=feat.alpha,
        solved=int(bool(solved)),
        sat=int(bool(sat)),
        timeout=timeout,
        runtime_budget_exceeded=runtime_budget_exceeded,
        decisions=int(solver.decisions),
        propagations=int(solver.propagations),
        conflicts=int(solver.conflicts),
        learned_clauses=int(solver.learned_clauses),
        max_level=int(solver.max_level),
        gate_value=float(gate_value),
        runtime_s=float(runtime_s),
        compute_cost=float(compute_cost),
    )


def paired_bootstrap(values: np.ndarray, n_boot: int, seed: int) -> dict[str, float]:
    if len(values) == 0:
        return {"mean": math.nan, "ci95_low": math.nan, "ci95_high": math.nan}
    rng = np.random.default_rng(seed)
    means = np.empty(n_boot, dtype=float)
    for idx in range(n_boot):
        sample = rng.choice(values, size=len(values), replace=True)
        means[idx] = float(np.mean(sample))
    lo, hi = np.quantile(means, [0.025, 0.975])
    return {"mean": float(np.mean(values)), "ci95_low": float(lo), "ci95_high": float(hi)}


def summarize_results(solve_df: pd.DataFrame, protocol: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    test = solve_df[solve_df["split"] == "test"].copy()
    pivot = test.pivot_table(index=["dataset_id", "family"], columns="policy", values="compute_cost")

    policy_rows = []
    for policy in POLICIES:
        vals = pivot[policy].dropna() if policy in pivot.columns else pd.Series(dtype=float)
        sub = test[test["policy"] == policy]
        policy_rows.append(
            {
                "policy": policy,
                "n_test": int(vals.shape[0]),
                "mean_compute_cost": float(vals.mean()) if len(vals) else math.nan,
                "median_compute_cost": float(vals.median()) if len(vals) else math.nan,
                "solve_rate": float(sub["solved"].mean()) if not sub.empty else math.nan,
                "timeout_rate": float(sub["timeout"].mean()) if not sub.empty else math.nan,
            }
        )

    n_boot = int(protocol["metrics"]["confidence_interval"]["resamples"])
    seed = int(protocol["metrics"]["confidence_interval"]["seed"])
    comparison_rows = []
    family_rows = []
    for comp in COMPARATORS:
        if "gate_plus_L" not in pivot.columns or comp not in pivot.columns:
            continue
        delta_series = (pivot[comp] - pivot["gate_plus_L"]).dropna()
        delta = delta_series.to_numpy(dtype=float)
        boot = paired_bootstrap(delta, n_boot, seed + stable_int(comp) % 100000)
        comparison_rows.append(
            {
                "comparison": f"gate_plus_L_vs_{comp}",
                "n_pairs": int(len(delta)),
                "delta_utility_mean": boot["mean"],
                "ci95_low": boot["ci95_low"],
                "ci95_high": boot["ci95_high"],
                "wins": int((delta > 0).sum()),
                "losses": int((delta < 0).sum()),
                "ties": int((delta == 0).sum()),
            }
        )
        family_delta = delta_series.reset_index(name="delta")
        for family, sub in family_delta.groupby("family"):
            vals = sub["delta"].to_numpy(dtype=float)
            family_rows.append(
                {
                    "comparison": f"gate_plus_L_vs_{comp}",
                    "family": family,
                    "n_pairs": int(len(vals)),
                    "delta_utility_mean": float(np.mean(vals)),
                    "wins": int((vals > 0).sum()),
                    "losses": int((vals < 0).sum()),
                    "ties": int((vals == 0).sum()),
                }
            )

    null_row = [r for r in comparison_rows if r["comparison"] == "gate_plus_L_vs_gate_shuffled_L"]
    shuffled_rows = null_row if null_row else []
    return (
        pd.DataFrame(policy_rows),
        pd.DataFrame(comparison_rows),
        pd.DataFrame(family_rows),
        pd.DataFrame(shuffled_rows),
    )


def classify_result(comparisons: pd.DataFrame) -> dict[str, object]:
    rows = {row["comparison"]: row for row in comparisons.to_dict(orient="records")}

    def ci_low(name: str) -> float:
        return float(rows.get(name, {}).get("ci95_low", math.nan))

    min_gate = ci_low("gate_plus_L_vs_frozen_gate_v17") > 0.0
    min_score = ci_low("gate_plus_L_vs_score_with_R") > 0.0
    null_survives = ci_low("gate_plus_L_vs_gate_shuffled_L") > 0.0
    strong_moms = ci_low("gate_plus_L_vs_moms") > 0.0
    strong_jw = ci_low("gate_plus_L_vs_jeroslow_wang") > 0.0

    minimum_support = bool(min_gate and min_score and null_survives)
    strong_sat_support = bool(minimum_support and strong_moms and strong_jw)
    if strong_sat_support:
        interpretation = "strong_sat_support"
        conclusion = "The preregistered Local Occurrence Hypothesis was supported at the strong SAT level under the predefined protocol."
    elif minimum_support:
        interpretation = "minimum_support_only"
        conclusion = "The preregistered Local Occurrence Hypothesis received minimum support only under the predefined protocol."
    else:
        interpretation = "no_minimum_support"
        conclusion = "The preregistered Local Occurrence Hypothesis was not supported under the predefined protocol."

    return {
        "interpretation_class": interpretation,
        "minimum_support_conditions": {
            "beats_frozen_gate_v17": min_gate,
            "beats_score_with_R": min_score,
            "survives_shuffled_L_null": null_survives,
        },
        "strong_sat_support_conditions": {
            "beats_MOMS": strong_moms,
            "beats_Jeroslow_Wang": strong_jw,
        },
        "conclusion": conclusion,
    }


def plot_outputs(policy_df: pd.DataFrame, comp_df: pd.DataFrame, family_df: pd.DataFrame, features_df: pd.DataFrame, solve_df: pd.DataFrame, outdir: Path) -> None:
    if not policy_df.empty:
        order = policy_df.sort_values("median_compute_cost")["policy"].tolist()
        data = policy_df.set_index("policy").loc[order]
        fig, ax = plt.subplots(figsize=(8.8, 4.8))
        ax.bar(data.index, data["median_compute_cost"], color="#2a9d8f", edgecolor="#1f2933")
        ax.set_ylabel("median compute cost (test split)")
        ax.set_title("Paper 73 frozen-policy comparison")
        ax.tick_params(axis="x", rotation=25)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(outdir / "fig1_policy_compute_cost.png", dpi=180)
        plt.close(fig)

    if not comp_df.empty:
        fig, ax = plt.subplots(figsize=(9.2, 4.9))
        x = np.arange(len(comp_df))
        means = comp_df["delta_utility_mean"].to_numpy(dtype=float)
        lows = comp_df["ci95_low"].to_numpy(dtype=float)
        highs = comp_df["ci95_high"].to_numpy(dtype=float)
        yerr = np.vstack([means - lows, highs - means])
        ax.axhline(0, color="#555", lw=1)
        ax.bar(x, means, color="#264653", edgecolor="#1f2933")
        ax.errorbar(x, means, yerr=yerr, fmt="none", ecolor="#e76f51", capsize=4, lw=1.4)
        labels = [c.replace("gate_plus_L_vs_", "vs ") for c in comp_df["comparison"]]
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=25, ha="right")
        ax.set_ylabel("paired utility delta (comparator cost - Gate+L cost)")
        ax.set_title("Gate+L preregistered comparisons")
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(outdir / "fig2_gate_l_comparisons.png", dpi=180)
        plt.close(fig)

    if not family_df.empty:
        pivot = family_df.pivot_table(index="family", columns="comparison", values="delta_utility_mean")
        fig, ax = plt.subplots(figsize=(10.5, 5.4))
        image = ax.imshow(pivot.to_numpy(dtype=float), aspect="auto", cmap="coolwarm")
        ax.set_yticks(np.arange(len(pivot.index)))
        ax.set_yticklabels(pivot.index)
        ax.set_xticks(np.arange(len(pivot.columns)))
        ax.set_xticklabels([c.replace("gate_plus_L_vs_", "vs ") for c in pivot.columns], rotation=30, ha="right")
        ax.set_title("Family-wise mean utility deltas")
        fig.colorbar(image, ax=ax, label="delta")
        fig.tight_layout()
        fig.savefig(outdir / "fig3_family_delta_heatmap.png", dpi=180)
        plt.close(fig)

    if not features_df.empty and not solve_df.empty:
        test_features = features_df[features_df["split"] == "test"][["dataset_id", "family", "L_star"]]
        pivot = solve_df[solve_df["split"] == "test"].pivot_table(index=["dataset_id", "family"], columns="policy", values="compute_cost").reset_index()
        if {"gate_plus_L", "frozen_gate_v17"}.issubset(pivot.columns):
            pivot["delta_vs_frozen_gate"] = pivot["frozen_gate_v17"] - pivot["gate_plus_L"]
            merged = test_features.merge(pivot[["dataset_id", "family", "delta_vs_frozen_gate"]], on=["dataset_id", "family"], how="inner")
            fig, ax = plt.subplots(figsize=(6.6, 4.8))
            for family, sub in merged.groupby("family"):
                ax.scatter(sub["L_star"], sub["delta_vs_frozen_gate"], label=family, s=28)
            ax.axhline(0, color="#555", lw=1)
            ax.set_xlabel("L_star")
            ax.set_ylabel("frozen gate cost - Gate+L cost")
            ax.set_title("Local occurrence channel vs Gate+L gain")
            ax.legend(fontsize=7, loc="best")
            ax.spines[["top", "right"]].set_visible(False)
            fig.tight_layout()
            fig.savefig(outdir / "fig4_lstar_vs_gate_gain.png", dpi=180)
            plt.close(fig)


def write_not_executed(reason: str, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    summary = {
        "paper": 73,
        "title": "Local Occurrence Hypothesis Execution",
        "status": "not_executed",
        "reason": reason,
        "interpretation_class": "not_applicable",
    }
    (outdir / "paper73_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (outdir / "paper73_run_log.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--output-dir", type=Path, default=OUTDIR)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    protocol = load_protocol()
    per_family = int(protocol["dataset_rules"]["primary_family_balanced_instances_per_family"])
    conflict_budget = int(protocol["budgets"]["primary_conflict_budget"])
    runtime_budget = int(protocol["budgets"]["runtime_budget_seconds_per_policy_instance"])

    validation = validate_dataset(args.data_dir, args.manifest, args.output_dir)
    if args.validate_only:
        print_validation(validation)
        return
    if not validation.get("dataset_ready"):
        write_not_executed("dataset_gatekeeper_failed", args.output_dir)
        print_validation(validation)
        return

    instances = discover_instances(args.data_dir, args.manifest, per_family)
    if not instances:
        write_not_executed("no_validated_instances_discovered", args.output_dir)
        print(json.dumps(json.loads((args.output_dir / "paper73_summary.json").read_text()), indent=2))
        return

    features = build_features(instances, protocol)
    by_id = {f.dataset_id: f for f in features}

    manifest_rows = [
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
        for inst in instances
    ]
    pd.DataFrame(manifest_rows).to_csv(args.output_dir / "paper73_dataset_manifest_detected.csv", index=False)
    features_df = pd.DataFrame([asdict(f) for f in features])
    features_df.to_csv(args.output_dir / "paper73_gate_l_features.csv", index=False)

    records: list[SolveRecord] = []
    for idx, inst in enumerate(instances):
        feat = by_id[inst.dataset_id]
        for policy in POLICIES:
            records.append(solve_instance(inst, feat, policy, conflict_budget, runtime_budget))
        print(f"[{idx + 1}/{len(instances)}] solved {inst.dataset_id}")

    solve_df = pd.DataFrame([asdict(r) for r in records])
    solve_df.to_csv(args.output_dir / "paper73_solve_records.csv", index=False)
    policy_df, comp_df, family_df, shuffled_df = summarize_results(solve_df, protocol)
    policy_df.to_csv(args.output_dir / "paper73_policy_summary.csv", index=False)
    comp_df.to_csv(args.output_dir / "paper73_gate_l_comparisons.csv", index=False)
    family_df.to_csv(args.output_dir / "paper73_family_results.csv", index=False)
    shuffled_df.to_csv(args.output_dir / "paper73_shuffled_l_null.csv", index=False)
    plot_outputs(policy_df, comp_df, family_df, features_df, solve_df, args.output_dir)

    classification = classify_result(comp_df)
    summary = {
        "paper": 73,
        "title": "Local Occurrence Hypothesis Execution",
        "status": classification["interpretation_class"],
        "execution_status": "executed",
        "paper72_protocol": P72_PROTOCOL_PORTABLE,
        "series_ii_doi": "10.5281/zenodo.21062386",
        "n_instances": len(instances),
        "n_test_instances": int(sum(1 for f in features if f.split == "test")),
        "n_calibration_instances": int(sum(1 for f in features if f.split == "calibration")),
        "family_counts": family_counts(instances),
        "policies": POLICIES,
        "conflict_budget": conflict_budget,
        "runtime_budget_seconds": runtime_budget,
        "bootstrap_resamples": int(protocol["metrics"]["confidence_interval"]["resamples"]),
        "bootstrap_seed": int(protocol["metrics"]["confidence_interval"]["seed"]),
        "dataset_gatekeeper": validation.get("status"),
        "interpretation_rule": "Only no_minimum_support, minimum_support_only, or strong_sat_support are allowed.",
        **classification,
    }
    (args.output_dir / "paper73_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (args.output_dir / "paper73_run_log.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
