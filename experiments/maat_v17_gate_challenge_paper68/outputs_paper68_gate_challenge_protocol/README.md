# Paper 68 Gate Challenge Protocol

This folder contains the preregistered protocol artifacts for Paper 68.
It does not run the external benchmarks.  It freezes the gate equation,
calibration rule, domains, baselines, metrics, null controls, and failure
criteria before future external validation.

Run:

```bash
python3 gate_challenge_protocol.py
```

Outputs are written to:

```text
outputs_paper68_gate_challenge_protocol/
```

Key files:

- `gate_challenge_preregistration.json`
- `gate_challenge_domain_matrix.csv`
- `gate_challenge_baseline_matrix.csv`
- `gate_challenge_metric_registry.csv`
- `fig1_gate_challenge_protocol.png`

Scientific status: protocol only, no external validation result.
