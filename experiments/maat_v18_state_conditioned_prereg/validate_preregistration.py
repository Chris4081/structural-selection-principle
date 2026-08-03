#!/usr/bin/env python3
"""Dependency-free gatekeeper for the Paper-75 release artifact."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PROTOCOL = ROOT / "outputs_paper75_preregistration" / "paper75_preregistration.json"
SCHEMA = ROOT / "paper75_preregistration.schema.json"
MANIFEST = ROOT / "dataset_manifest_paper75.csv"
FORBIDDEN_SUFFIXES = {".cnf", ".gz", ".bz2", ".xz", ".zip", ".tar", ".part"}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def validate() -> list[str]:
    errors: list[str] = []
    for path in (PROTOCOL, SCHEMA, MANIFEST):
        if not path.exists():
            errors.append(f"missing required artifact: {path.name}")
    if errors:
        return errors
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    for key in schema["required"]:
        if key not in protocol:
            errors.append(f"protocol missing required key: {key}")
    expected = {
        "paper": 75, "year": 2026,
        "status": "preregistration_only_no_external_execution",
    }
    for key, value in expected.items():
        if protocol.get(key) != value:
            errors.append(f"frozen value mismatch: {key}")
    constants = protocol.get("frozen_constants", {})
    frozen = {
        "activation_budget_q": 0.25, "no_harm_delta_relative": 0.05,
        "distinctness_epsilon": 0.01, "reconstruction_rho_star": 0.90,
        "conflict_budget": 5000, "runtime_budget_seconds_per_policy_instance": 30,
        "bootstrap_resamples": 10000, "bootstrap_seed": 75075,
        "minimum_test_activation": 0.125,
    }
    for key, value in frozen.items():
        if constants.get(key) != value:
            errors.append(f"frozen constant mismatch: {key}")
    if protocol.get("dataset", {}).get("sha256") != sha256_file(MANIFEST):
        errors.append("frozen manifest SHA256 mismatch")
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if protocol.get("dataset", {}).get("rows") != len(rows):
        errors.append("manifest row count mismatch")
    if any(r.get("paper75_split") not in {"calibration", "test"} for r in rows):
        errors.append("manifest contains invalid split")
    for path in ROOT.rglob("*"):
        if path.is_file() and (path.suffix.lower() in FORBIDDEN_SUFFIXES or path.name == ".DS_Store"):
            errors.append(f"forbidden release artifact: {path.relative_to(ROOT)}")
    if any("result" in p.name.lower() or "solve_record" in p.name.lower() for p in ROOT.rglob("*")):
        errors.append("possible execution result artifact present")
    return errors


def main() -> None:
    errors = validate()
    if errors:
        print("VALIDATION FAIL")
        for error in errors:
            print(f"- {error}")
        raise SystemExit(1)
    print("Dataset metadata .... PASS")
    print("Manifest SHA256 ..... PASS")
    print("Frozen constants .... PASS")
    print("No raw CNFs ......... PASS")
    print("No execution outputs  PASS")
    print("VALIDATION PASS")


if __name__ == "__main__":
    main()
