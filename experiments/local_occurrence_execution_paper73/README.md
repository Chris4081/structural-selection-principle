# Paper 73 -- Local Occurrence Hypothesis Execution

This folder executes the Paper-72 preregistered SAT-specific protocol.

Scientific status: execution paper only. It does not change the Paper-72
Gate+L definition, residualization, calibration rule, budgets, bootstrap rule,
success criteria, failure criteria, baselines, or family-balanced split.

Execution result: `no_minimum_support`. The frozen Local Occurrence Hypothesis
was tested under the predefined protocol and was not supported. Paper 72 remains
the preregistered protocol document; Paper 73 reports only that the frozen
hypothesis did not pass its first execution test.

Paper 72 protocol:

```text
experiments/local_occurrence_hypothesis_paper72/outputs_paper72_local_occurrence_protocol/paper72_protocol.json
```

Research Series II / Gate Challenge Cycle I archive:

Krieg, C. (2026). MAAT Research Series II -- Gate Challenge Cycle I.
Zenodo. DOI: 10.5281/zenodo.21062386
https://doi.org/10.5281/zenodo.21062386

## Data Policy

Raw CNFs are not committed here. Stage public SATLIB/DIMACS instances locally:

```bash
python3 download_satlib_paper73.py
```

If the legacy SATLIB HTTPS endpoint cannot be verified by the local Python
certificate store, use the explicit opt-in mode only after checking the source:

```bash
python3 download_satlib_paper73.py --allow-insecure-ssl
```

The manifest records download TLS mode, archive SHA256, and CNF SHA256. Raw CNFs
and SATLIB archives stay under `data_external/`, which is git-ignored and must
not be committed.

SATLIB benchmark instances are not redistributed in this repository or release.
The included manifests contain only source URLs, hashes, and metadata required
for reproducibility. Users must obtain the original benchmark instances from
SATLIB and cite the SATLIB source.

## Run

Validate first:

```bash
python3 local_occurrence_execution_paper73.py --validate-only
```

Execute the frozen Paper-72 protocol:

```bash
python3 local_occurrence_execution_paper73.py
```

No parameter optimization, threshold change, weight adjustment, or retuning is
allowed.

## Outputs

Outputs are written to:

```text
outputs_paper73_local_occurrence_execution/
```

Main files:

- `paper73_local_occurrence_hypothesis_execution.tex`
- `paper73_dataset_validation.json`
- `paper73_dataset_manifest_detected.csv`
- `paper73_gate_l_features.csv`
- `paper73_solve_records.csv`
- `paper73_policy_summary.csv`
- `paper73_gate_l_comparisons.csv`
- `paper73_family_results.csv`
- `paper73_shuffled_l_null.csv`
- `paper73_summary.json`
- `paper73_run_log.json`
- `fig1_policy_compute_cost.png`
- `fig2_gate_l_comparisons.png`
- `fig3_family_delta_heatmap.png`
- `fig4_lstar_vs_gate_gain.png`

The committed CSV/JSON files are derived metadata and solver outputs. They do
not contain raw DIMACS clauses or redistributed SATLIB benchmark files.

## Interpretation Classes

Paper 73 may report only one of the Paper-72 classes:

- `no_minimum_support`
- `minimum_support_only`
- `strong_sat_support`

If Gate+L loses, the conclusion is:

```text
The preregistered Local Occurrence Hypothesis was not supported under the predefined protocol.
```
