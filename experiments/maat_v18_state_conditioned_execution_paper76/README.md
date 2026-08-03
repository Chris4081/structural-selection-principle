# Paper 76 - State-Conditioned Gate Arbitration Execution

This folder executes, without retuning, the protocol frozen by Paper 75 and
publicly archived before execution at:

`https://doi.org/10.5281/zenodo.21767033`

## Frozen execution result

The complete 98-instance run has been executed under the released protocol.
Its independent status axes are:

```text
execution_status  = executed
construct_status  = distinct
activation_status = adequate_activation
safety_status     = harmful
utility_status    = negative_result
```

The family-balanced primary comparison against MOMS is `delta_U=-642.96`
with 95% bootstrap interval `[-841.63,-295.82]`, where positive values would
favor state v1.8. Relative regret against MOMS is `0.9135` with interval
`[0.6956,1.1756]`. The preregistered state-conditioned gate is therefore not
supported as a practical SAT/CDCL improvement under this protocol.

No parameter or interpretation rule was changed after observing this result.

## Non-negotiable boundary

Paper 76 does not modify the Paper-75 state definition, supports, closure,
normalization, calibration policy, activation budget, tie handling, fallback,
Mode-A policy, splits, budgets, bootstrap, baselines, nulls, ablations, or
status axes. The runner refuses execution when the release, protocol,
manifest, raw-file hashes, DIMACS parser, or unexpected-file checks fail.

## Data policy

Raw SATLIB CNFs and archives are excluded from this folder and from repository
releases. Obtain them from SATLIB using the source URLs in
`dataset_manifest_paper75.csv`. The manifest records the exact SHA256 required
for every instance. SATLIB authors and maintainers do not endorse MAAT.

## Licensing

Paper 76, its figures, and the derived result tables are released under the
[Creative Commons Attribution 4.0 International license](https://creativecommons.org/licenses/by/4.0/)
(CC BY 4.0). Source code remains under the repository MIT License. External
SATLIB benchmark materials are not included and retain their original terms.
The manifest provides attribution, source URLs, and integrity metadata; it does
not grant additional rights in the external benchmark instances.

For an existing Paper-69 staging directory containing the 98 frozen files:

```bash
export MAAT_SSL_REPO="/path/to/structural-selection-principle"
python3 paper76_execution.py \
  --data-dir "../gate_challenge_sat_paper69/data_external" \
  --validate-only
```

Only after `VALIDATION PASS`:

```bash
python3 paper76_execution.py \
  --data-dir "../gate_challenge_sat_paper69/data_external"
python3 paper76_report.py
python3 validate_paper76_outputs.py
```

For a read-only verification of an archived result folder:

```bash
python3 validate_paper76_outputs.py --check-only
```

## Reproduced artifacts

- `paper76_dataset_validation.json`
- `paper76_calibration.json`
- `paper76_calibration_states.csv`
- `paper76_solve_records.csv`
- `paper76_policy_summary.csv`
- `paper76_primary_comparisons.csv`
- `paper76_family_results.csv`
- `paper76_cost_decomposition_by_family.csv`
- `paper76_cost_decomposition.json`
- `paper76_summary.json`
- `paper76_run_log.json`
- `paper76_output_validation.json`
- publication figures `fig1`--`fig4`

The CSV and JSON outputs are derived solver measurements and contain no clause
lists or redistributed CNF contents.

## Development tests

The tests use only explicitly excluded synthetic fixtures:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
```
