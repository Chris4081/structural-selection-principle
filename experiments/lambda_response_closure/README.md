# Response-Based Derivation of MAAT Structural Weights

This directory contains the reproducibility bundle for a response-theoretic
closure of MAAT structural weights, used in

**Response-Based Derivation of MAAT Structural Weights: From Covariance
Geometry to Selection Pressure**.

The purpose is to move one step beyond a purely fitted `lambda_a` hierarchy.
Here the weights are derived as effective response coefficients of a concrete
defect ensemble.

## Scientific Status

This is an effective physical/statistical derivation, not a proof that
`lambda_a` are universal constants of nature.

The result should be read as:

```text
lambda_a = sector response pressure induced by defect fluctuations under mu_0
```

not as:

```text
lambda_a = unique microscopic fundamental constant
```

## Method

Let

```text
d(X) = (d_H, d_B, d_S, d_V, d_R)
```

be the primitive MAAT defect vector of an ensemble with empirical reference
measure `mu_0`.  Around `mu_0`, exponential tilting gives

```text
P_lambda(X) ∝ mu_0(X) exp[-lambda · d(X)].
```

Linear response implies

```text
delta <d_a> = - sum_b Cov_mu0(d_a,d_b) lambda_b.
```

Therefore, for a target selected sub-ensemble,

```text
lambda = (C + eta tr(C)/A I)^(-1) (<d>_mu0 - <d>_target).
```

The ridge term stabilises inversion of the covariance matrix.  Negative
components are projected to zero because `lambda_a` is interpreted as a
non-negative penalty weight; a negative raw value means that the corresponding
sector is not acting as a penalty for the chosen target.

## Input Data

The script uses the fused boundary-aware dataset from:

```text
../boundary_aware_lambda_calibration/maat_defects_fused.csv
```

The dataset contains 3400 samples:

| Source | Samples |
| --- | ---: |
| SAT hardness atlas | 2000 |
| MAAT-Core | 500 |
| MAAT-Core boundary | 900 |

## Run

From this directory:

```bash
python3 lambda_response_closure.py
```

This generates:

- `lambda_response_closure_results.json`
- `plots/lambda_response_shares.png`
- `plots/lambda_response_correlation.png`
- `plots/lambda_response_target_match.png`
- `plots/lambda_response_vs_closed_fit.png`

## Main Interpretation

Different target definitions correspond to different physical selection
questions:

| Target | Meaning |
| --- | --- |
| `low_defect_20pct` | generic low-defect structural selection |
| `safe_and_core_safe` | safety / robustness dominated selected ensemble |
| `safe_boundary_only` | explicit safe boundary states |
| `not_violated` | exclusion of violated configurations |

The important point is not that all targets produce the same weights.  The
important point is that the weights are no longer arbitrary: they are computed
from the covariance geometry of the defect ensemble and a specified target
selection condition.

## Results

Using ridge `eta = 1e-3`, the response-derived shares are:

| Target | H | B | S | V | R | Interpretation |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `low_defect_20pct` | 0.1587 | 0.2099 | 0.1808 | 0.1534 | 0.2972 | Generic low-defect selection; robustness is the largest sector but not exclusive. |
| `safe_and_core_safe` | 0.0945 | 0.0701 | 0.0000 | 0.0023 | 0.8331 | Explicit safe targets are overwhelmingly robustness-driven. |
| `safe_boundary_only` | 0.1093 | 0.1830 | 0.1911 | 0.1597 | 0.3569 | Boundary-safe states reproduce a balanced hierarchy with dominant robustness. |
| `not_violated` | 0.0000 | 0.0000 | 0.2602 | 0.0000 | 0.7398 | Excluding violations activates mainly robustness and activity. |

The strongest physical conclusion is:

```text
R dominance emerges as a covariance-response effect of safety/boundary target
selection, not as a manually imposed assumption.
```

## Relation to Previous Lambda Fits

The previous boundary-aware MaxEnt fit solved for a closed lambda hierarchy by
optimisation.  This experiment supplies the missing response-theoretic
interpretation:

```text
lambda_a is the response pressure needed to move the reference ensemble
toward a selected lower-defect ensemble.
```

In statistical-mechanics language, this is analogous to deriving inverse
temperature as a Lagrange multiplier conjugate to an energy constraint, but
with multiple structural defect constraints.
