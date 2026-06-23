# External Q_bar comparison

This folder contains a neutral comparison between the externally supplied HGD-GSR Q_bar table and the existing Paper-65 SPARC galaxy-level diagnostics.

## Coverage

- Paper-65 summary rows: `175`
- Q_bar rows: `175`
- Matched galaxies: `175`

## Strongest rank associations

| target | Spearman rho | permutation p |
|---|---:|---:|
| `vobs_max` | 0.6175 | 0.0002 |
| `median_bulge_fraction_proxy` | 0.5584 | 0.0002 |
| `mean_V` | -0.5469 | 0.0002 |
| `nfw_like_rmse_v2` | 0.5459 | 0.0002 |
| `peak_D_struct` | -0.4783 | 0.0002 |
| `mean_D_struct` | -0.4345 | 0.0002 |

## Independence-test highlights

- Direct `Q_bar`--`mean_D_struct` Spearman correlation is reported in `qbar_correlations.csv`.
- Partial rank checks include controls for `log_mbar_proxy_msun` and `log_r_max_kpc` in `qbar_partial_rank_checks.csv`.
- Mismatch prediction comparisons for `D_struct only`, `Q_bar only`, and `D_struct + Q_bar` are reported in `qbar_independence_model_comparison.csv`.

The baryonic mass control is a simple enclosed baryonic mass-scale proxy derived from `V_bar(R_max)^2 R_max / G`; it is not a full photometric stellar-mass estimate.

## Interpretation note

These results are diagnostic correlations only. They do not establish a physical model or causal relation. Q_bar is treated as an external input and compared against existing MAAT/SPARC structural summaries without changing the Paper-61 or Paper-65 pipelines.

## Data attribution and license note

SPARC-derived MAAT inputs follow the Paper-65 SPARC attribution and CC-BY-4.0 notes. SPARC rotation-curve data should be cited to Lelli, McGaugh, and Schombert and to the Zenodo-hosted SPARC record when reused.

The `Q_bar` table is a collaborator-supplied derived HGD-GSR descriptor from Ali Alhawarat. It is not an original SPARC measurement and not a MAAT-derived quantity. Unless broader redistribution terms are separately supplied by the HGD-GSR author, treat the table as a neutral comparison input for this pilot and cite/acknowledge Ali Alhawarat and HGD-GSR when discussing or reusing the comparison.

No endorsement by SPARC, Zenodo, VizieR, CDS, HGD-GSR, Ali Alhawarat, the original SPARC data providers, or any catalogue authors is implied.
