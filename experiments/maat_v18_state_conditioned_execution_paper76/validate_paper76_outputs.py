#!/usr/bin/env python3
"""Validate completeness and frozen interpretation of Paper-76 outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs_paper76_state_conditioned_execution"
EXPECTED_STATUS = {
    "execution_status": "executed",
    "construct_status": "distinct",
    "activation_status": "adequate_activation",
    "safety_status": "harmful",
    "utility_status": "negative_result",
}


def validate() -> list[str]:
    failures: list[str] = []
    required = [
        "paper76_dataset_validation.json", "paper76_calibration.json",
        "paper76_solve_records.csv", "paper76_policy_summary.csv",
        "paper76_primary_comparisons.csv", "paper76_family_results.csv",
        "paper76_summary.json", "paper76_run_log.json",
    ]
    for name in required:
        if not (OUT / name).exists():
            failures.append(f"missing output: {name}")
    if failures:
        return failures

    validation = json.loads((OUT / "paper76_dataset_validation.json").read_text(encoding="utf-8"))
    summary = json.loads((OUT / "paper76_summary.json").read_text(encoding="utf-8"))
    records = pd.read_csv(OUT / "paper76_solve_records.csv")
    if validation.get("status") != "VALIDATION PASS":
        failures.append("dataset gatekeeper did not pass")
    if len(records) != 1283:
        failures.append(f"expected 1283 solver records, found {len(records)}")
    calibration = records[records.split == "calibration"]
    test = records[records.split == "test"]
    if len(calibration) != 19 or set(calibration.policy) != {"moms"}:
        failures.append("calibration records are not exactly 19 MOMS-only runs")
    if test.dataset_id.nunique() != 79 or test.policy.nunique() != 16:
        failures.append("test matrix is not 79 instances by 16 policies")
    counts = test.groupby("dataset_id").policy.nunique()
    if not counts.eq(16).all():
        failures.append("at least one test instance lacks a policy result")
    if records.audit_max_abs_error.max() > 1.0e-12:
        failures.append("incremental support audit exceeded 1e-12")
    statuses = summary.get("status_axes", {})
    for key, expected in EXPECTED_STATUS.items():
        if statuses.get(key) != expected:
            failures.append(f"status mismatch for {key}: {statuses.get(key)!r}")
    text_outputs = "\n".join(
        path.read_text(encoding="utf-8", errors="ignore")
        for path in OUT.iterdir() if path.suffix in {".json", ".csv"}
    )
    if "/Users/" in text_outputs or "/Volumes/" in text_outputs:
        failures.append("personal absolute path found in derived outputs")
    if '"clauses"' in text_outputs or "original_clauses" in text_outputs:
        failures.append("raw clause content field found in derived outputs")
    return failures


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-only", action="store_true", help="validate without rewriting the PASS artifact")
    args = parser.parse_args()
    failures = validate()
    checks = {
        "dataset_gatekeeper": "PASS" if not any("gatekeeper" in f for f in failures) else "FAIL",
        "record_matrix": "PASS" if not any("record" in f or "matrix" in f or "policy result" in f for f in failures) else "FAIL",
        "incremental_audit": "PASS" if not any("audit" in f for f in failures) else "FAIL",
        "status_axes": "PASS" if not any("status mismatch" in f for f in failures) else "FAIL",
        "portable_derived_outputs": "PASS" if not any("path" in f or "clause" in f for f in failures) else "FAIL",
    }
    result = {
        "paper": 76, "validation": checks, "failures": failures,
        "status": "VALIDATION PASS" if not failures else "VALIDATION FAIL",
    }
    if not args.check_only:
        (OUT / "paper76_output_validation.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if not failures else 1)


if __name__ == "__main__":
    main()
