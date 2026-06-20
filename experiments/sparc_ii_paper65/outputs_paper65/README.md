# Paper 65: SPARC II

SPARC II extends Paper 61 with shuffled-radius nulls, morphology/proxy splits,
and an optional published halo-catalogue cross-check.

Headline radius-null result:

- D_struct vs NFW-like pooled radial Spearman: `0.3930`
- D_struct vs RAR pooled radial Spearman: `0.0123`
- NFW-like exceeds radius-shuffle |rho|95: `True`
- RAR exceeds radius-shuffle |rho|95: `False`

Halo-catalogue status: `not_run_no_local_catalogue`.

Outputs:

- `paper65_radius_shuffle_nulls.csv`
- `paper65_radius_null_summary.json`
- `paper65_galaxywise_radial_correlations.csv`
- `paper65_morphology_proxy_splits.csv`
- `paper65_halo_catalogue_join.csv` if a local halo catalogue is supplied
- `paper65_summary.json`
- four PNG figures

No endorsement by SPARC, Li et al., Zenodo, VizieR, CDS, or any original data
provider is implied.
