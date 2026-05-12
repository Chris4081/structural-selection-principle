# Paper 52 — Cross-Domain Transfer in Structural Selection

This experiment tests whether a frozen MAAT defect architecture transfers
between two previously generated domains without target-domain retuning.

## Domains

| Domain | Source artifact | Target convention |
|---|---|---|
| SAT | `../sat_frustration_fields_paper50b/outputs_paper50b/paper50b_sat_instances.csv` | higher quality = lower normalised DPLL effort |
| Quantum pointer states | `../quantum_pointer_state_selection_paper51/outputs_paper51/paper51_pointer_instances.csv` | higher quality = higher pointer robustness / fidelity target |

The common transfer architecture uses only the shared supports:

```text
H, B, S, V, R_rob
```

## Frozen Transfer Protocol

For each source domain, support defects are computed as:

```text
D_a = -log(epsilon + Gamma_a)
```

Non-negative weights are fitted on the source domain only:

```text
lambda_source = argmin ||D_source lambda - (1 - quality_source)||^2
```

The weights are normalised and frozen.  The target-domain score is then:

```text
S_frozen(X) = - sum_a lambda_a D_a(X)
```

Higher score is predeclared to mean higher structural quality.  No target-domain
sign flips, retuning, or metric switching are allowed.

## Null Models

Two null families are generated:

- `shuffled defects`: each target defect column is independently shuffled,
  preserving marginal distributions but destroying sample-level sector
  coherence.
- `lambda permutation`: the frozen source weights are permuted across sectors.

## Headline Results

| Transfer | Frozen Spearman | Equal/Scalar Spearman | Interpretation |
|---|---:|---:|---|
| SAT -> Quantum | 0.6034 | 0.3771 | strong asymmetric transfer, above shuffled-defect null |
| Quantum -> SAT | -0.0408 | -0.2683 | weak/failed positive transfer, but less bad than global scalar baselines |

The result is therefore not a universal transfer proof.  It is a partial,
asymmetric transfer result: the SAT-trained architecture collapses onto a
connectivity-dominated mode that transfers well to the quantum pointer-state
benchmark, while the quantum-trained mixture does not solve SAT hardness.

## Reproduction

Run from the repository root:

```bash
python3 experiments/cross_domain_transfer_paper52/cross_domain_transfer_paper52.py
```

The script writes CSV, JSON, and PNG outputs to:

```text
experiments/cross_domain_transfer_paper52/outputs_paper52/
```

## Data and License Note

This experiment reuses repository-internal synthetic/derived artifacts from
Paper 50b and Paper 51.  No new external dataset is introduced here.  The
outputs are reproducibility artifacts for the structural-selection benchmark.

