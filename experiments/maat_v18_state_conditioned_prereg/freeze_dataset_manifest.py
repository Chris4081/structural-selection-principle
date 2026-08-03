#!/usr/bin/env python3
"""Freeze Paper-75 metadata, hashes and splits without loading raw CNFs."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

from state_conditioned_spec import deterministic_split


FIELDS = [
    "dataset_id", "family", "source_name", "source_url", "license_or_terms",
    "local_path", "sha256", "download_tls_mode", "archive_sha256", "notes",
    "paper75_split", "paper75_split_key", "paper75_instance_weight",
]


def freeze(source: Path, destination: Path) -> None:
    with source.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("source manifest is empty")
    required = {"dataset_id", "family", "sha256", "local_path"}
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"missing source columns: {sorted(missing)}")
    families: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        digest = row["sha256"].strip().lower()
        if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError(f"invalid SHA256 for {row['dataset_id']}")
        row["paper75_split_key"] = deterministic_split(row["dataset_id"], digest)
        families[row["family"]].append(row)
    n_families = len(families)
    for family_rows in families.values():
        family_rows.sort(key=lambda r: (r["paper75_split_key"], r["dataset_id"]))
        n_cal = max(1, math.floor(0.20 * len(family_rows)))
        for index, row in enumerate(family_rows):
            row["paper75_split"] = "calibration" if index < n_cal else "test"
            row["paper75_instance_weight"] = f"{1.0 / (n_families * len(family_rows)):.17g}"
    frozen = sorted((row for group in families.values() for row in group), key=lambda r: (r["family"], r["dataset_id"]))
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(frozen)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_manifest", type=Path)
    parser.add_argument("--output", type=Path, default=Path("dataset_manifest_paper75.csv"))
    args = parser.parse_args()
    freeze(args.source_manifest, args.output)
    print(f"FROZEN MANIFEST WRITTEN: {args.output}")


if __name__ == "__main__":
    main()
