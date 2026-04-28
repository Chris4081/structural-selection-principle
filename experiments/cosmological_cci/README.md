# Cosmological Critical Coherence Index

This directory contains the reproducibility bundle for

**Cosmological Critical Coherence Index: A Structural-Stress Observable for
Cosmic Evolution**.

This is Paper 28 in the structural-selection series. It defines a compact
cosmological CCI diagnostic that combines expansion stress, redshift activity,
linear growth coherence, and a structural imbalance term. The benchmark is an
observable proposal, not a new cosmological model or a replacement for
precision parameter inference.

## Contents

```text
cosmological_cci/
├── README.md
├── maat_cci_cosmology_v02.py
├── maat_cci_cosmology_v02_chronometers.csv
├── maat_cci_cosmology_v02.csv
├── maat_cci_cosmology_v02_data_comparison.csv
└── plots/
    ├── maat_cci_cosmology_v02_plot.png
    ├── maat_cci_cosmology_v02_Hz_comparison.png
    ├── maat_cci_cosmology_v02_data_comparison.png
    └── maat_cci_cosmology_v02_residuals.png
```

## Reproduce

Run from this directory:

```bash
python3 maat_cci_cosmology_v02.py
```

This regenerates:

- `maat_cci_cosmology_v02.csv`
- `maat_cci_cosmology_v02_data_comparison.csv`
- `plots/maat_cci_cosmology_v02_plot.png`
- `plots/maat_cci_cosmology_v02_Hz_comparison.png`
- `plots/maat_cci_cosmology_v02_data_comparison.png`
- `plots/maat_cci_cosmology_v02_residuals.png`

## Inputs

The script uses Planck-2018 flat LCDM reference parameters:

- `H0 = 67.4 km s^-1 Mpc^-1`
- `Omega_m = 0.315`
- `Omega_Lambda = 0.685`
- `sigma8_0 = 0.811`

The chronometer table `maat_cci_cosmology_v02_chronometers.csv` is a compact
literature-derived H(z) table used for the observational projection.

## Data Provenance

The original Planck parameter values and Cosmic Chronometer measurements belong
to the cited publications and collaborations. The CSV files and plots in this
folder are derived analysis artifacts generated for reproducibility of the
Cosmological CCI diagnostic.

No endorsement by the Planck Collaboration or the chronometer-data authors is
implied.

## License and Reuse

The code in this directory follows the repository license. The derived CSV
tables and PNG figures are provided as reproducibility artifacts for the
Cosmological CCI diagnostic.

The underlying Planck-2018 cosmological parameters and Cosmic Chronometer H(z)
measurements are external scientific data. Reuse of those original data remains
subject to the terms of the original publications, journals, and
collaborations. When reusing or discussing the numerical inputs, cite the
Planck-2018 cosmological-parameter paper and the Cosmic Chronometer sources
listed in the manuscript.

## Scientific Status

This benchmark:

- is a diagnostic projection over Planck LCDM plus chronometer H(z) data;
- uses an approximate growth factor rather than a Boltzmann-code pipeline;
- does not include BAO, DESI, supernovae, weak lensing, or RSD data;
- does not perform parameter inference or model comparison;
- should be read as a structural-stress observable proposal.

## Paper

The compiled paper is stored at:

```text
documentation/28_Cosmological_Critical_Coherence_Index.pdf
```
