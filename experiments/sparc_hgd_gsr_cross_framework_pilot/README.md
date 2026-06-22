# SPARC MAAT x HGD-GSR Cross-Framework Pilot

This experiment is a collaboration-ready pilot for comparing MAAT Paper-65
galaxy-level structural residuals with externally supplied HGD-GSR `Q_bar`
values on the same SPARC galaxy basis.

It is not a claim that either framework replaces the other. The test asks
whether the MAAT structural residual signal and the HGD-GSR structural signal
are related, independent, complementary, or redundant.

## Required External Input

Place a collaborator-provided file here:

```text
data/hgd_gsr_qbar.csv
```

Required columns:

```text
galaxy,Q_bar
NGC2403,0.123
...
```

Accepted aliases:

- Galaxy column: `galaxy`, `Galaxy`, `name`, `Name`, `SPARC`, `galaxy_name`
- Q-bar column: `Q_bar`, `q_bar`, `Qbar`, `qbar`, `hgd_qbar`, `HGD_Q_bar`

If the file is absent, the script creates:

```text
data/hgd_gsr_qbar_template.csv
outputs_cross_framework/DATA_CONTRACT.md
```

## Run

```bash
python3 sparc_hgd_gsr_cross_framework_pilot.py
```

## Tests

The pilot computes:

- galaxy-level correlation between `mean_D_struct`, `peak_D_struct`, `mean_R_rob`,
  and `Q_bar`;
- predictive comparison against Paper-65 mismatch targets:
  `nfw_like_chi2`, `nfw_like_rmse_v2`, `rar_mean_sigma_residual`,
  and `mean_residual_fraction`;
- bootstrap confidence intervals for correlation signals;
- cross-validated OLS comparisons:
  baseline baryonic proxies, baseline plus `D_struct`, baseline plus `Q_bar`,
  and baseline plus both.

## Outputs

```text
outputs_cross_framework/cross_framework_joined.csv
outputs_cross_framework/cross_framework_correlations.csv
outputs_cross_framework/cross_framework_bootstrap.csv
outputs_cross_framework/cross_framework_cv_results.csv
outputs_cross_framework/cross_framework_summary.json
outputs_cross_framework/fig1_dstruct_vs_qbar.png
outputs_cross_framework/fig2_cv_model_comparison.png
outputs_cross_framework/fig3_bootstrap_correlations.png
```

## Licence and Attribution

MAAT/SPARC-derived inputs follow the Paper-65 SPARC attribution and license
notes. SPARC rotation-curve data are external scientific data listed under
CC-BY-4.0 at Zenodo DOI `10.5281/zenodo.16284118`.

HGD-GSR `Q_bar` values are not redistributed by this repository unless
explicitly supplied with permission and clear citation/licence terms by the
HGD-GSR author.

No endorsement by SPARC, HGD-GSR, Ali Alhawarat, Zenodo, or any original data
provider is implied.
