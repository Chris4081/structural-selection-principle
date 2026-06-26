# Paper 69 -- Gate Challenge I: SAT/CDCL

This folder contains the execution scaffold for the first external SAT/CDCL
run of the preregistered MAAT v1.7 Gate Challenge protocol.

Scientific status: **execution scaffold and reproducible SATLIB smoke execution**.
Raw external CNF files are not committed; published outputs are derived from
the validated SATLIB subset documented in the SHA256 manifest. The script still
writes a `not_executed` summary when no external DIMACS instances are present.

Preregistration archive:

```text
Krieg, C. (2026). MAAT Research Series I -- Papers 1--68.
Zenodo. DOI: 10.5281/zenodo.20882852
https://doi.org/10.5281/zenodo.20882852
```

## What Is Frozen

The run follows the Paper 68 gate protocol:

- SAT polarity: `p_D = +1`
- gate-vs-score comparison: `gate_v17` versus `score_with_R`
- gate-null comparison: deterministic shuffled-`R_rob`
- primary metrics: conflicts, decisions, runtime, compute cost, and paired
  utility deltas
- confidence intervals: paired bootstrap on instance-level utility deltas
- failure rule: a negative or non-positive lower confidence bound against
  score-with-`R` is not support for the gate hypothesis

## Data Policy

Raw benchmark files are intentionally not committed here.

Raw SATLIB benchmark instances are intentionally excluded from this repository.
They must be obtained from the original SATLIB source and verified against the
provided SHA256 manifest before running the benchmarks.

Place public external DIMACS CNF files under:

```text
data_external/
```

Then copy the manifest template:

```bash
cp dataset_manifest_template.csv dataset_manifest.csv
```

and fill in source, license/terms, local path, and optional SHA256 information
for each benchmark file. The script also computes SHA256 hashes automatically
for detected files and writes a detected manifest into the output folder.

Examples of suitable external sources include SATLIB, SAT Competition benchmark
families, and DIMACS-style public CNF archives. Before redistribution or upload,
the original benchmark terms must be checked. No endorsement by any dataset
provider is implied.

For the recommended SATLIB smoke execution, use:

```bash
python3 download_satlib_paper69.py --limit-per-archive 20
```

This stages a deterministic SATLIB subset, writes `dataset_manifest.csv`, and
keeps raw CNFs under the git-ignored `data_external/` folder. Use
`--limit-per-archive 0` only when you intentionally want to stage full selected
archives.

If the legacy SATLIB HTTPS endpoint fails because the local Python certificate
store cannot verify the chain, the stager may be rerun with the explicit opt-in
flag:

```bash
python3 download_satlib_paper69.py --limit-per-archive 20 --allow-insecure-ssl
```

This does not change any gate parameter. The manifest marks
`download_tls_mode=tls_verification_disabled_opt_in`, records each archive
SHA256, and records every extracted CNF SHA256. The stager downloads archives
via temporary `.part` files and updates `dataset_manifest.csv` after every
successfully staged archive, so completed raw files are not left undocumented.
Do not commit raw CNF files.

## Run

Validate the staged external dataset before any solver execution:

```bash
python3 gate_challenge_sat_paper69.py --validate-only
```

This writes:

```text
outputs_paper69_sat_gate_challenge/paper69_dataset_validation.json
```

Only after `VALIDATION PASS`, run:

```bash
python3 gate_challenge_sat_paper69.py
```

Optional:

```bash
python3 gate_challenge_sat_paper69.py --max-instances 50 --conflict-budget 5000
```

Family-balanced smoke test:

```bash
python3 gate_challenge_sat_paper69.py --family-balanced 5 --conflict-budget 1000
```

This selects up to five CNFs per benchmark family after manifest/SHA256
validation. It is preferred over `--max-instances` for infrastructure tests
because it avoids lexicographic family bias.

Outputs are written to:

```text
outputs_paper69_sat_gate_challenge/
```

Main outputs:

- `paper69_dataset_manifest_detected.csv`
- `paper69_dataset_validation.json`
- `paper69_gate_features.csv` (only after execution with CNFs)
- `paper69_solve_records.csv`
- `paper69_policy_summary.csv`
- `paper69_gate_comparisons.csv`
- `paper69_summary.json`
- `paper69_run_log.json`
- `fig1_policy_compute_cost.png` (only after execution with CNFs)
- `fig2_gate_comparisons.png` (only after execution with CNFs)

Post-hoc diagnostic, not retuning:

```bash
python3 analyze_moms_vs_gate_paper69.py
```

This writes family-level and instance-level diagnostics explaining where MOMS
beats the frozen gate. It must not be used to change gate parameters for Paper
69.

## Baseline Scope

The transparent CDCL backend currently supports VSIDS, MOMS, Jeroslow--Wang,
progress-only, score-with-`R`, gate-aware MAAT, and shuffled-`R` gate nulls.
LRB/CHB are listed in the Paper 68 protocol as stronger future solver baselines;
they require either extending this auditable backend or wrapping a modern
external CDCL solver.
