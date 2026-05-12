# Paper 53 — Structural Mode Transfer Robustness

**Robustness of Structural Mode Transfer in MAAT v1.4**  
*Alternative Frozen Fits, Regularisation Tests, and Multi-Domain Null Validation*

This experiment extends Paper 52 by testing whether the observed
`V`/connectivity transfer is stable under alternative source-only fitting
rules and additional internal benchmark domains.

## Domains

| Domain | Samples | Quality convention |
|---|---:|---|
| SAT-Frustration | 550 | `1 - minmax(log DPLL nodes)` |
| Quantum-Pointer | 3500 | `minmax(pointer robustness target)` |
| Active-Significance | 6000 | `minmax(R_sig)` |
| Societal-CCI | 5000 | `minmax(ASI_soc)` |

All domains are internal synthetic or derived artifacts from the repository.
This is not an independent external validation test.

## Fit Rules

The shared support architecture uses:

```text
H, B, S, V, R_rob
```

The following source-only fitting rules are compared:

- `equal`
- `nnls`
- `ridge_pos`
- `lasso_pos`
- `response_top20`
- `v_only`

No target-domain retuning or sign flips are allowed.

## Headline Results

Paper 53 confirms that the Paper 52 SAT -> Quantum `V` result is not merely an
NNLS artifact:

| Rule | SAT -> Quantum Spearman | `lambda_V` |
|---|---:|---:|
| `nnls` | 0.6034 | 1.0000 |
| `ridge_pos` | 0.5981 | 0.8394 |
| `lasso_pos` | 0.6034 | 1.0000 |
| `response_top20` | 0.5975 | 0.8264 |
| `v_only` | 0.6034 | 1.0000 |
| `equal` | 0.3771 | 0.2000 |

However, `V` is not a universal transfer law.  Across all sources and fitting
rules, the dominant sector is `V` in only about 58% of frozen architectures.
Equal weights achieve the largest mean cross-domain Spearman because two
additional internal domains use support-composite quality targets.

## Reproduction

Run from the repository root:

```bash
python3 experiments/structural_mode_transfer_paper53/structural_mode_transfer_robustness_paper53.py
```

Outputs are written to:

```text
experiments/structural_mode_transfer_paper53/outputs_paper53/
```

## Data and License Note

This experiment reuses repository-internal synthetic/derived artifacts from
earlier papers.  No new external dataset is redistributed.

