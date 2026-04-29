# Cosmological CCI v0.3: Growth Connectivity and Robustness

This folder contains the reproducibility bundle for Paper 29:

**Cosmological Critical Coherence Index with Growth Connectivity and Robustness Margins**  
*A Multi-Sector Structural-Stress Benchmark Using H(z) and f sigma_8(z)*

The benchmark extends the chronometer-only cosmological CCI projection from
Paper 28 by making the previously open connectivity and robustness sectors
operational:

- `H`: expansion consistency from Cosmic Chronometer H(z) residuals
- `S`: structural activity / imbalance from the E(z)-D(z) mismatch
- `V`: growth connectivity from BOSS DR12 f sigma_8 residuals
- `R`: robustness from joint expansion-growth consistency

It also derives:

- a companion maximum-entropy sector-weight fit `lambda_a`
- a data-driven transition proxy `z_c` from the curvature of
  `log CCI_v03`

## Scientific Status

This is a diagnostic benchmark, not a precision cosmology likelihood. It does
not fit cosmological parameters, replace LambdaCDM, or claim evidence for
modified gravity. Its purpose is to close the operational definitions of the
multi-sector cosmological CCI and provide a small reproducible test pipeline.

## Contents

```text
cosmological_cci_v03/
├── README.md
├── maat_cci_cosmology_v03_growth.py
├── maat_cci_cosmology_v03_chronometers.csv
├── maat_cci_cosmology_v03_fsigma8.csv
├── maat_cci_cosmology_v03_grid.csv
├── maat_cci_cosmology_v03_results.json
└── plots/
    ├── cci_v02_v03_comparison.png
    ├── fsigma8_growth_connectivity.png
    ├── lambda_calibration.png
    ├── transition_curvature.png
    └── v_r_supports.png
```

## Reproducing the Results

Run from this directory:

```bash
python3 maat_cci_cosmology_v03_growth.py
```

This regenerates:

- `maat_cci_cosmology_v03_grid.csv`
- `maat_cci_cosmology_v03_results.json`
- all PNG files in `plots/`

The script uses fixed model settings and deterministic calculations, so the
outputs should be reproducible for the same Python/package versions.

## Main Numerical Outputs

The current run gives:

| Quantity | Value |
| --- | ---: |
| Transition proxy `z_c` | `1.1140311804` |
| H(z) pull RMS | `0.6927` |
| f sigma_8 pull RMS | `0.6975` |
| `CCI_v03_norm(z=1)` | `6.6006` |
| `CCI_v03_norm(z=2)` | `15.4152` |
| `min V_corr` | `0.5096` |
| `min R_robust` | `0.4599` |

The fitted companion weights are:

| Sector | lambda | Share |
| --- | ---: | ---: |
| H | `1.7739` | `0.2930` |
| S | `1.9607` | `0.3239` |
| V | `1.0509` | `0.1736` |
| R | `1.2677` | `0.2094` |

Hierarchy:

```text
lambda_S > lambda_H > lambda_R > lambda_V
```

## Data Provenance and License Note

The Planck-2018 parameter values are quoted from the Planck cosmological
parameter paper. The Cosmic Chronometer H(z) values are a compact
literature-derived table from the chronometer sources cited in Paper 29. The
f sigma_8 values are BOSS DR12 consensus measurements reported in Alam et al.
(2017).

The original observational measurements remain subject to the terms of their
original publications, journals, and collaborations. The CSV and PNG files in
this folder are derived analysis artifacts provided only for reproducibility of
this diagnostic benchmark. No endorsement by the Planck Collaboration,
SDSS/BOSS Collaboration, or the chronometer-data authors is implied.

## Repository Link

Public folder:

<https://github.com/Chris4081/structural-selection-principle/tree/main/experiments/cosmological_cci_v03>

