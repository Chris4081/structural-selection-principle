# Paper 40 — MAAT v1.2.1 Structural Signature Test in Growth Data

This folder contains the reproducibility package for:

**MAAT v1.2.1 Structural Signature Test in Growth Data**  
Projection CCI, Residual Magnitudes, and Exploratory Null Tests

## Purpose

The benchmark tests whether MAAT structural diagnostics correlate with the
residual structure of a compact `f sigma_8(z)` growth comparison set relative
to a Planck-normalised LCDM baseline.

This is a diagnostic signature test only. It is not a full cosmological
likelihood, not a Boltzmann-code result, and not evidence for modified growth.
Because the balance support `B` is constructed from residual information, the
CCI-residual correlation is semi-supervised rather than a fully independent
blind prediction.

## Run

```bash
python3 maat_paper40_structural_signature_test.py
```

## Main Results

| Quantity | Result |
| --- | --- |
| Growth comparison points | `13` |
| Projection transition estimate | `z_tr = 0.6508` |
| Spearman `R_proj` vs signed residual | `0.6319`, `p = 0.0228` |
| Spearman `V` vs signed residual | `-0.6319`, `p = 0.0228` |
| Spearman `CCI_diag` vs `|residual_sigma|` | `0.5934`, `p = 0.0338` |
| Random-field null for `CCI_diag` | `p = 0.0363` |
| Redshift-shuffle null for `CCI_diag` | `p = 0.0359` |
| Spearman `B` vs `|residual_sigma|` | `-0.7912`, `p = 0.0017` |

The `B` result is expected because `B` is residual-sensitive in this benchmark.
The scientific reading is therefore cautious: Paper 40 is a structural
diagnostic consistency test, not a detection claim.

## Outputs

The script writes:

- `outputs_paper40/paper40_summary.json`
- `outputs_paper40/paper40_signature_table.csv`
- `outputs_paper40/fig1_cci_diag_vs_signed_residual.png`
- `outputs_paper40/fig2_scatter_cci_diag_signed_residual.png`
- `outputs_paper40/fig3_scatter_cci_diag_abs_residual.png`
- `outputs_paper40/fig4_spearman_signed_residuals.png`
- `outputs_paper40/fig5_spearman_abs_residuals.png`
- `outputs_paper40/fig6_breakthrough_null_test_cci_diag.png`

## Data Attribution and License Notes

Planck-normalised reference parameters and the compact `f sigma_8` comparison
points are external scientific data and should be cited to the original
publications/collaborations when reused or discussed. The CSV tables and
figures in this folder are derived analysis artifacts generated for
reproducibility of the MAAT diagnostic calculation. No endorsement by the
Planck Collaboration, survey collaborations, or original data authors is
implied.

