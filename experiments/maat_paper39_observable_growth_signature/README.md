# Paper 39 — MAAT v1.2.1 Observable Growth Signature Proxy

This folder contains the reproducibility package for:

**MAAT v1.2.1 Observable Growth Signature Proxy**  
Projection-Modulated `f sigma_8`, Emergent Robustness, and a Boundary-Limited Diagnostic Scan

## Purpose

The benchmark tests whether a bounded MAAT projection template can be inserted
as a small observable signature in the growth observable `f sigma_8(z)` while
preserving the MAAT v1.2.1 robustness closure:

```text
R_resp = (H B V)^(1/3)
R_rob  = min(R_resp, (H B S V)^(1/4))
Stability = R_rob
```

The proxy is diagnostic only. It is not a full perturbation calculation, not a
Boltzmann-code result, and not a precision cosmological likelihood.

## Run

```bash
python3 maat_paper39_observable_signature_v121.py
```

## Main Results

| Quantity | Result |
| --- | --- |
| Growth comparison points | `13` |
| Projection transition estimate | `z_tr = 0.6508` |
| Scan range | `epsilon in [-0.01, 0.01]` |
| Best epsilon | `-0.0100` |
| LCDM chi2 | `12.4373` |
| MAAT proxy chi2 | `12.3772` |
| Delta chi2 | `-0.0601` |
| Max `|Delta f sigma_8 / f sigma_8|` | `0.9891%` |
| Mean `|Delta f sigma_8 / f sigma_8|` | `0.5725%` |
| Mean `R_resp` | `0.7247` |
| Mean `R_rob` | `0.6673` |
| Mean `CCI_min` | `0.2662` |
| Mean `CCI_diag` | `0.2036` |

The best `epsilon` value lies at the lower boundary of the scan interval, so it
is interpreted as a weak directional diagnostic rather than a measured
parameter.

## Outputs

The script writes:

- `outputs_paper39/paper39_summary.json`
- `outputs_paper39/paper39_growth_signature_results.csv`
- `outputs_paper39/paper39_epsilon_chi2_scan.csv`
- `outputs_paper39/paper39_projection_template_v121.csv`
- `outputs_paper39/fig1_projection_template.png`
- `outputs_paper39/fig2_growth_comparison.png`
- `outputs_paper39/fig3_residuals.png`
- `outputs_paper39/fig4_chi2_scan.png`
- `outputs_paper39/fig5_v121_robustness_closure.png`
- `outputs_paper39/fig6_v121_cci.png`
- `outputs_paper39/fig7_v121_support_fields.png`

## Data Attribution and License Notes

Planck-normalised reference parameters and the compact `f sigma_8` comparison
points are external scientific data and should be cited to the original
publications/collaborations when reused or discussed. The CSV tables and
figures in this folder are derived analysis artifacts generated for
reproducibility of the MAAT diagnostic proxy. No endorsement by the Planck
Collaboration, survey collaborations, or original data authors is implied.

