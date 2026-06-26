#!/usr/bin/env python3
"""Download and stage a documented SATLIB subset for Paper 69.

The script downloads public SATLIB tarballs, extracts DIMACS CNF files under
``data_external/satlib/``, and writes ``dataset_manifest.csv`` for the solver
script. Raw SATLIB files remain ignored by git; only the manifest and derived
outputs should be committed after checking the original benchmark terms.

Default mode is deliberately small and deterministic: it extracts at most
``--limit-per-archive`` CNF files from each archive after lexicographic sorting.
Use ``--limit-per-archive 0`` to extract all CNFs from the selected archives.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import ssl
import tarfile
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data_external" / "satlib"
ARCHIVE_DIR = DATA_DIR / "_archives"
MANIFEST = ROOT / "dataset_manifest.csv"

SATLIB_SOURCE = "SATLIB"
SATLIB_BENCHMARK_PAGE = "https://www.cs.ubc.ca/~hoos/SATLIB/benchm.html"
SATLIB_TERMS_NOTE = (
    "SATLIB public benchmark archive; terms are not stated as an explicit "
    "software/data license on the benchmark page. Use for reproducible research "
    "with attribution; check original terms before redistributing raw CNFs."
)


@dataclass(frozen=True)
class SatlibArchive:
    family: str
    label: str
    url: str
    expected_status: str
    notes: str


ARCHIVES: dict[str, SatlibArchive] = {
    "uf20": SatlibArchive(
        family="satlib_random_3sat_uf20",
        label="uf20-91",
        url="https://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/RND3SAT/uf20-91.tar.gz",
        expected_status="sat",
        notes="Uniform random 3-SAT, 20 variables, 91 clauses; SATLIB lists 1000 satisfiable instances.",
    ),
    "uf50": SatlibArchive(
        family="satlib_random_3sat_uf50",
        label="uf50-218",
        url="https://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/RND3SAT/uf50-218.tar.gz",
        expected_status="sat",
        notes="Uniform random 3-SAT, 50 variables, 218 clauses; SATLIB lists 1000 satisfiable instances.",
    ),
    "uuf50": SatlibArchive(
        family="satlib_random_3sat_uuf50",
        label="uuf50-218",
        url="https://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/RND3SAT/uuf50-218.tar.gz",
        expected_status="unsat",
        notes="Uniform random 3-SAT, 50 variables, 218 clauses; SATLIB lists 1000 unsatisfiable instances.",
    ),
    "flat30": SatlibArchive(
        family="satlib_graph_coloring_flat30",
        label="flat30-60",
        url="https://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/GCP/flat30-60.tar.gz",
        expected_status="sat",
        notes="SAT-encoded flat graph-colouring instances; SATLIB lists 100 satisfiable instances.",
    ),
    "aim": SatlibArchive(
        family="satlib_dimacs_aim",
        label="aim",
        url="https://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/DIMACS/AIM/aim.tar.gz",
        expected_status="mixed",
        notes="DIMACS AIM artificial random 3-SAT family mirrored by SATLIB.",
    ),
    "dubois": SatlibArchive(
        family="satlib_dimacs_dubois",
        label="dubois",
        url="https://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/DIMACS/DUBOIS/dubois.tar.gz",
        expected_status="unsat",
        notes="DIMACS Dubois randomly generated unsatisfiable instances mirrored by SATLIB.",
    ),
    "phole": SatlibArchive(
        family="satlib_dimacs_pigeonhole",
        label="pigeon-hole",
        url="https://www.cs.ubc.ca/~hoos/SATLIB/Benchmarks/SAT/DIMACS/PHOLE/pigeon-hole.tar.gz",
        expected_status="unsat",
        notes="DIMACS pigeon-hole instances mirrored by SATLIB.",
    ),
}

DEFAULT_ARCHIVES = ["uf20", "uuf50", "flat30", "aim", "dubois", "phole"]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def safe_member_name(name: str) -> Path:
    path = Path(name)
    clean_parts = [p for p in path.parts if p not in {"", ".", ".."}]
    return Path(*clean_parts)


def download_archive(archive: SatlibArchive, archive_dir: Path, allow_insecure_ssl: bool) -> Path:
    archive_dir.mkdir(parents=True, exist_ok=True)
    out = archive_dir / f"{archive.label}.tar.gz"
    if out.exists() and out.stat().st_size > 0:
        return out
    tmp = out.with_name(out.name + ".part")
    request = urllib.request.Request(archive.url, headers={"User-Agent": "MAAT-Paper69-SATLIB-Stager/1.0"})
    context = ssl._create_unverified_context() if allow_insecure_ssl else None
    try:
        with urllib.request.urlopen(request, timeout=60, context=context) as response:
            data = response.read()
        tmp.write_bytes(data)
        tmp.replace(out)
    except urllib.error.URLError as exc:
        tmp.unlink(missing_ok=True)
        if "CERTIFICATE_VERIFY_FAILED" in str(exc) and not allow_insecure_ssl:
            raise RuntimeError(
                "TLS certificate verification failed for a SATLIB download. "
                "If the source URL has been verified separately, rerun with "
                "--allow-insecure-ssl. In that mode, the manifest will mark the "
                "TLS exception and every extracted CNF will still be fixed by SHA256."
            ) from exc
        raise
    return out


def extract_cnfs(
    archive: SatlibArchive,
    tar_path: Path,
    limit_per_archive: int,
    allow_insecure_ssl: bool,
) -> list[dict[str, str]]:
    target_dir = DATA_DIR / archive.label
    target_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    archive_sha256 = sha256_file(tar_path)
    tls_mode = "tls_verification_disabled_opt_in" if allow_insecure_ssl else "tls_verified"

    with tarfile.open(tar_path, "r:gz") as tf:
        members = [m for m in tf.getmembers() if m.isfile() and m.name.lower().endswith((".cnf", ".dimacs"))]
        members = sorted(members, key=lambda m: m.name)
        if limit_per_archive > 0:
            members = members[:limit_per_archive]

        for member in members:
            file_obj = tf.extractfile(member)
            if file_obj is None:
                continue
            data = file_obj.read()
            rel_member = safe_member_name(member.name)
            if rel_member.name == "":
                continue
            out_path = target_dir / rel_member.name
            # Avoid accidental path traversal and keep archive-internal
            # directory noise out of the run folder.
            out_path.write_bytes(data)
            rows.append(
                {
                    "dataset_id": f"{archive.label}_{rel_member.stem}",
                    "family": archive.family,
                    "source_name": SATLIB_SOURCE,
                    "source_url": archive.url,
                    "license_or_terms": SATLIB_TERMS_NOTE,
                    "local_path": str(out_path.relative_to(ROOT)),
                    "sha256": sha256_file(out_path),
                    "download_tls_mode": tls_mode,
                    "archive_sha256": archive_sha256,
                    "notes": f"{archive.notes} Expected status: {archive.expected_status}.",
                }
            )
    return rows


def write_manifest(rows: list[dict[str, str]], manifest_path: Path) -> None:
    fieldnames = [
        "dataset_id",
        "family",
        "source_name",
        "source_url",
        "license_or_terms",
        "local_path",
        "sha256",
        "download_tls_mode",
        "archive_sha256",
        "notes",
    ]
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--archives",
        nargs="+",
        default=DEFAULT_ARCHIVES,
        choices=sorted(ARCHIVES),
        help="SATLIB archive keys to stage.",
    )
    parser.add_argument(
        "--limit-per-archive",
        type=int,
        default=20,
        help="Deterministic CNF extraction limit per archive. Use 0 for all files.",
    )
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument(
        "--allow-insecure-ssl",
        action="store_true",
        help=(
            "Opt-in only: disable TLS certificate verification for legacy SATLIB "
            "downloads if the local Python certificate store fails. The manifest "
            "will record this mode and all extracted CNFs are SHA256-fixed."
        ),
    )
    args = parser.parse_args()

    rows: list[dict[str, str]] = []
    for key in args.archives:
        archive = ARCHIVES[key]
        tar_path = download_archive(archive, ARCHIVE_DIR, args.allow_insecure_ssl)
        extracted = extract_cnfs(archive, tar_path, args.limit_per_archive, args.allow_insecure_ssl)
        rows.extend(extracted)
        # Write after every successful archive so staged raw data never remain
        # undocumented if a later archive fails.
        write_manifest(rows, args.manifest)
        print(f"{archive.label}: extracted {len(extracted)} CNF files")

    write_manifest(rows, args.manifest)
    print(f"Wrote {len(rows)} manifest rows to {args.manifest}")
    print("SATLIB benchmark page:", SATLIB_BENCHMARK_PAGE)


if __name__ == "__main__":
    main()
