# Paper 61 - SPARC Structural Residual Pilot

This folder contains the reproducibility pipeline for:

**Paper 61 - SPARC Pilot Test of Structural Dark-Matter Residuals**  
*A MAAT v1.6 Real-Data Bridge from Rotation-Curve Residuals to NFW/MOND Baselines*

## Status

This folder now contains a SPARC data preparation path and a real-data pilot
run based on the Zenodo-hosted `Rotmod_LTG.zip` file.  The result is a pilot
diagnostic benchmark, not a precision dark-matter halo analysis.

## Expected local data format

The prepared CSV is:

```text
data/sparc_rotation_curves.csv
```

Required columns, with common aliases accepted:

| Canonical column | Meaning |
|------------------|---------|
| `galaxy` | galaxy identifier |
| `r_kpc` | radius in kpc |
| `v_obs_km_s` | observed circular velocity |
| `e_v_obs_km_s` | velocity uncertainty; optional |
| `v_gas_km_s` | gas velocity contribution |
| `v_disk_km_s` | disk velocity contribution |
| `v_bulge_km_s` | bulge velocity contribution |
| `upsilon_disk` | optional disk mass-to-light factor |
| `upsilon_bulge` | optional bulge mass-to-light factor |

## Run

Prepare the SPARC Rotmod files:

```bash
python3 prepare_sparc_rotmod.py
```

Run the SPARC pilot:

```bash
python3 sparc_structural_residual_pilot.py --input data/sparc_rotation_curves.csv --output outputs_sparc_real
```

## Outputs

The SPARC output folder is:

```text
outputs_sparc_real/
```

| File | Role |
|------|------|
| `paper61_radius_diagnostics.csv` | Radius-level baryonic, residual, NFW-like, RAR, and MAAT diagnostics |
| `paper61_galaxy_summary.csv` | Galaxy-level summaries and baseline mismatch values |
| `paper61_summary.json` | Run metadata, mode, summary statistics, and license note |
| `fig1_rotation_decomposition_examples.png` | Example rotation decompositions |
| `fig2_population_structural_vs_baselines.png` | Population-level structural support versus NFW/RAR mismatch |
| `fig3_shuffled_galaxy_nulls.png` | Predeclared shuffled-galaxy null comparison |

## Data attribution and license note

SPARC data were obtained from the Zenodo record `10.5281/zenodo.16284118`,
which lists CC-BY-4.0.  When using or redistributing these data:

- Cite Lelli, McGaugh, and Schombert (2016), *The Astronomical Journal* 152, 157.
- If using the Zenodo-hosted SPARC record, cite DOI `10.5281/zenodo.16284118`.
- Respect the Zenodo record's CC-BY-4.0 attribution requirement.
- If using VizieR, include the standard VizieR/CDS acknowledgment.
- Do not imply endorsement by the SPARC authors, Zenodo, VizieR, CDS, or the original data providers.
- Clearly distinguish original SPARC measurements from derived MAAT artifacts.

## Scope checklist

- The reported Paper-61 numbers are pilot correlations, not a dark-matter model fit.
- The NFW-like and MOND/RAR-like channels are comparison diagnostics only.
- They are not full professional halo fits or precision modified-dynamics analyses.
- No endorsement by SPARC authors, Zenodo, VizieR, CDS, or original data providers is claimed or implied.
