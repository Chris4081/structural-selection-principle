# Paper 68 -- The Gate Challenge

This folder contains the preregistration artifacts for:

**The Gate Challenge: A Predeclared Falsification Protocol for MAAT v1.7 Gate-Aware Structural Selection**

Scientific status: protocol only. This folder does not contain external
validation results and should not be cited as evidence that MAAT v1.7 is
supported.

The protocol freezes:

- gate equation,
- calibration rule,
- domain polarity,
- metrics,
- baselines,
- null controls,
- bootstrap requirements,
- success and failure criteria.

Run:

```bash
python3 gate_challenge_protocol.py
```

Generated files are written to:

```text
outputs_paper68_gate_challenge_protocol/
```

Main artifacts:

- `gate_challenge_preregistration.json`
- `gate_challenge_domain_matrix.csv`
- `gate_challenge_baseline_matrix.csv`
- `gate_challenge_metric_registry.csv`
- `fig1_gate_challenge_protocol.png`

Future external validation runs must not alter the gate equation or threshold
calibration rule after seeing test-set outcomes.
