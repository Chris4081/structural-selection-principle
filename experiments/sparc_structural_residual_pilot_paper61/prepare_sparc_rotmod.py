#!/usr/bin/env python3
"""
Prepare SPARC Rotmod_LTG.zip for the Paper 61 pipeline.

Input:
    data/raw/Rotmod_LTG.zip

Output:
    data/sparc_rotation_curves.csv

The raw SPARC Rotmod files use the columns:

    Rad Vobs errV Vgas Vdisk Vbul SBdisk SBbul

This converter keeps the velocity columns needed by the MAAT/NFW/RAR pilot.
It does not modify the original measurements except for assigning canonical
column names and adding fixed default mass-to-light factors used by the pilot.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from zipfile import ZipFile

import pandas as pd


def parse_rotmod_text(name: str, text: str, upsilon_disk: float, upsilon_bulge: float) -> pd.DataFrame:
    galaxy = re.sub(r"_rotmod\.dat$", "", Path(name).name)
    distance_mpc = None
    rows = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            if "Distance" in stripped:
                match = re.search(r"Distance\s*=\s*([0-9.+-Ee]+)", stripped)
                if match:
                    distance_mpc = float(match.group(1))
            continue
        parts = stripped.split()
        if len(parts) < 6:
            continue
        rad, vobs, errv, vgas, vdisk, vbul = map(float, parts[:6])
        rows.append(
            {
                "galaxy": galaxy,
                "distance_mpc": distance_mpc,
                "r_kpc": rad,
                "v_obs_km_s": vobs,
                "e_v_obs_km_s": errv,
                "v_gas_km_s": vgas,
                "v_disk_km_s": vdisk,
                "v_bulge_km_s": vbul,
                "upsilon_disk": upsilon_disk,
                "upsilon_bulge": upsilon_bulge,
            }
        )
    return pd.DataFrame(rows)


def run(args: argparse.Namespace) -> dict:
    in_zip = Path(args.input)
    out_csv = Path(args.output)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    if not in_zip.exists():
        raise FileNotFoundError(f"Missing SPARC Rotmod zip: {in_zip}")

    frames = []
    with ZipFile(in_zip) as zf:
        names = sorted(name for name in zf.namelist() if name.endswith("_rotmod.dat"))
        for name in names:
            text = zf.read(name).decode("utf-8", errors="replace")
            frame = parse_rotmod_text(name, text, args.upsilon_disk, args.upsilon_bulge)
            if not frame.empty:
                frames.append(frame)

    df = pd.concat(frames, ignore_index=True)
    df.to_csv(out_csv, index=False)
    summary = {
        "input": str(in_zip),
        "output": str(out_csv),
        "n_galaxies": int(df["galaxy"].nunique()),
        "n_rows": int(len(df)),
        "upsilon_disk": float(args.upsilon_disk),
        "upsilon_bulge": float(args.upsilon_bulge),
        "license_note": (
            "Prepared from SPARC Rotmod_LTG.zip. Cite Lelli, McGaugh & Schombert "
            "(2016). If using the Zenodo record, cite DOI 10.5281/zenodo.16284118 "
            "and respect CC-BY-4.0 attribution."
        ),
    }
    summary_path = out_csv.with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare SPARC Rotmod_LTG.zip for Paper 61")
    parser.add_argument("--input", default="data/raw/Rotmod_LTG.zip")
    parser.add_argument("--output", default="data/sparc_rotation_curves.csv")
    parser.add_argument("--upsilon-disk", type=float, default=0.5)
    parser.add_argument("--upsilon-bulge", type=float, default=0.7)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())

