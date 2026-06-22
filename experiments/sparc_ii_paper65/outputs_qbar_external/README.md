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

## Interpretation note

These results are diagnostic correlations only. They do not establish a physical model or causal relation. Q_bar is treated as an external input and compared against existing MAAT/SPARC structural summaries without changing the Paper-61 or Paper-65 pipelines.

No endorsement by SPARC, Zenodo, VizieR, CDS, the original data providers, or the external Q_bar author is implied.
