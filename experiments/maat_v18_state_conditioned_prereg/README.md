# Paper 75 - State-Conditioned Gate Arbitration Preregistration

This folder freezes the SAT/CDCL protocol for MAAT v1.8 Gate Challenge Cycle II. It is a preregistration only. It contains no external benchmark execution and no Paper-76 result.

Paper 76 must not run until the exact Paper-75 release has a public Zenodo DOI.

## Frozen scientific scope

The protocol fixes the decision state `X_t`, code-level `H/B/S/V` supports, normalization, temporal response, causal update and rollback semantics, MOMS-only calibration traces, MOMS fallback, state-bound common random numbers, activation budget `q=0.25`, minimum test activation `q_test_min=0.125`, no-harm bound `delta=0.05`, distinctness `epsilon=0.01`, reconstruction threshold `rho_star=0.90`, quantile-tie handling, data splits, budgets, overhead accounting, bootstrap, baselines, ablations, and independent status axes.

Uncertainty is a secondary report only. It cannot alter the primary gate, threshold, tie fraction, activation budget, or status axes. `no_harm_not_established` is explicitly separated from `harmful`.

## Preregistration artifacts

- `outputs_paper75_preregistration/paper75_preregistration.json`: machine-readable frozen protocol.
- `paper75_preregistration.schema.json`: JSON Schema for the protocol.
- `dataset_manifest_paper75.csv`: frozen SATLIB metadata, source URLs, SHA256 hashes, family-balanced weights, and calibration/test assignments.
- `state_conditioned_spec.py`: pure state semantics used by synthetic unit tests; it is not a solver runner.
- `tests/test_paper75.py`: synthetic/excluded-fixture tests only.

Raw SATLIB CNFs and archives are intentionally excluded. Users must obtain benchmark files from the original SATLIB source and verify them against the frozen hashes. No endorsement by SATLIB or its authors is implied.

## Rebuild and validate

The manifest was derived from the already documented Paper-69 metadata without loading or solving any CNF. To write the protocol and validate the release:

```bash
python3 paper75_preregistration.py
python3 validate_preregistration.py
python3 -m unittest discover -s tests -v
```

Expected final validator line:

```text
VALIDATION PASS
```

## Development rule

Before the Paper-75 Zenodo release, implementation checks may use only small synthetic fixtures in `tests/` or fixtures explicitly excluded from primary evidence. External state-conditioned calibration, threshold estimation, solver comparison, or utility measurement is forbidden before the DOI exists.

## Paper 76 boundary

Paper 76 may execute only the released protocol. It may not change support equations, update semantics, fallback, common-random-number rule, `q`, `q_test_min`, `delta`, `epsilon`, `rho_star`, tie handling, splits, budgets, baselines, bootstrap, or status rules. Any change requires a new preregistration and a new public timestamp.
