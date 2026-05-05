# Paper 42 — Blind Projection Test

This folder reproduces **Paper 42: Blind Projection Test**.

The benchmark removes the projection-shape tuning used in earlier residual
diagnostics. The projection parameters are derived from response-closed MAAT
v1.2.1 sector weights:

```text
defects -> covariance C -> lambda
        -> response shares pi_a
        -> gamma_lambda, Bstar_lambda, alpha_lambda
        -> C_proj_lambda(z)
        -> residual / CCI / permutation tests
```

The test is blind with respect to the projection shape: there is no epsilon
fit, no gamma tuning, no Bstar scan, and no alpha scan. The regulator
`EPS = 1e-12` is fixed only to avoid numerical division by zero.

## Run

```bash
cd experiments/maat_paper42_blind_projection_test
python3 maat_paper42_blind_projection_test.py
```

The script uses a fixed random seed (`42`) and `20,000` permutations for the
main null tests.

## Outputs

The script writes into `outputs_paper42/`:

- `paper42_summary.json`
- `paper42_correlation_table.csv`
- `paper42_blind_signature_table.csv`
- `paper42_lambda_projection_curve.csv`
- `fig1_blind_projection_curve.png`
- `fig2_lambda_response_shares.png`
- `fig3_blind_cci_vs_abs_residual.png`
- `fig4_blind_cproj_vs_abs_residual.png`
- `fig5_null_test_blind_cci.png`

## Main Results

| Quantity | Result |
| --- | --- |
| Growth comparison points | `13` |
| Response shares | `pi_H=0.0960`, `pi_B=0.1257`, `pi_S=0.3879`, `pi_V=0.3903` |
| Derived projection parameters | `gamma_lambda=1.3305`, `Bstar_lambda=4.5091`, `alpha_lambda=1.9937` |
| Transition marker | `z_tr = 0.6986` |
| Spearman `CCI_diag` vs `abs(residual_sigma)` | `0.4286`, permutation `p = 0.1456` |
| Spearman `C_proj_lambda` vs `abs(residual_sigma)` | `0.4286`, permutation `p = 0.1432` |

The result is positive but non-significant. It is a stress test of intrinsic
projection robustness, not a detection claim.

## Data Attribution and License Notes

The Planck-normalised reference parameters and compact `f sigma_8` comparison
points are external scientific data and should be cited to the original
publications/collaborations when reused or discussed. The CSV tables and
figures in this folder are derived analysis artifacts generated for
reproducibility of the MAAT diagnostic calculation.

No endorsement by the Planck Collaboration, survey collaborations, or original
data authors is implied. Repository code is released under the repository
license.
