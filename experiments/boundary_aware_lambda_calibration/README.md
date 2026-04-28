# Boundary-Aware MAAT Lambda Calibration

This directory contains the reproducibility bundle for

**Boundary-Aware Calibration of MAAT Structural Weights: A Reproducible
Benchmark for Constraint-Dominated Structural Selection**.

The benchmark calibrates structural selection weights over a fused defect
dataset containing SAT hardness instances, unconstrained MAAT-Core states, and
explicit MAAT-Core boundary regimes.

## Scientific Status

This is a boundary-aware calibration benchmark, not a proof that the fitted
weights are universal constants. The goal is narrower: to test whether explicit
constraint-boundary information changes the inferred hierarchy of structural
selection weights.

## Dataset

The fused dataset contains 3400 samples:

| Source | Samples | File |
| --- | ---: | --- |
| SAT hardness atlas | 2000 | `maat_defects_sat.csv` |
| MAAT-Core | 500 | `maat_defects_core.csv` |
| MAAT-Core boundary | 900 | `maat_defects_core_boundary.csv` |
| Fused dataset | 3400 | `maat_defects_fused.csv` |

The common defect columns are:

```text
d_H, d_B, d_S, d_V, d_R, label, source
```

## Run

From this directory:

```bash
python3 merge_defects.py
python3 fit_closed_maat_lambda_v1.py
python3 plot_closed_maat_lambda_v2.py
```

This generates:

- `maat_defects_fused.csv`
- `closed_maat_lambda_fit_results.json`
- `plots/fig1_lambda_distribution.png`
- `plots/fig2_lambda_shares.png`
- `plots/fig3_sat_vs_core_defects.png`
- `plots/fig4_structural_energy_distribution.png`
- `plots/fig5_c_hat_vs_flambda_sat.png`

## Result

The fitted weights are:

| Sector | lambda_a | Share |
| --- | ---: | ---: |
| H | 3.218288 | 0.1561 |
| B | 2.469024 | 0.1197 |
| S | 2.526500 | 0.1225 |
| V | 4.328860 | 0.2099 |
| R | 8.078183 | 0.3917 |

The main qualitative hierarchy is:

```text
R > V > H > S ~= B
```

This supports the interpretation that robustness / respect becomes the
dominant selector once explicit boundary regimes are included.

The weak prior used in the regularised fit is the previous MaxEnt calibration
from the constants benchmark:

```text
lambda^(0) = (1.610873, 2.428297, 4.623655, 4.644912, 8.832344)
```

for `(H, B, S, V, R)`.

## Optional Provenance

The file `export_core_boundary_defects.py` documents how the boundary dataset
was generated from MAAT-Core. It requires the MAAT-Core package and is included
for provenance. The reproducible analysis in this folder starts from the
exported CSV files above.

## Paper

The compiled PDF is stored in:

```text
documentation/27_Boundary_Aware_Calibration_of_MAAT_Structural_Weights.pdf
```
