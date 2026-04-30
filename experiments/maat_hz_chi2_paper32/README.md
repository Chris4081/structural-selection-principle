# Paper 32: MAAT H(z) Chronometer Comparison

This folder contains the reproducibility package for Paper 32:

**First H(z) Data Comparison of MAAT Structural-Selection Cosmology**  
*Cosmic Chronometer Diagnostic and Chi-Square Scan*

The benchmark compares MAAT expansion histories with a compact Cosmic
Chronometer `H(z)` table and a Planck-like LCDM reference. It is a diagnostic
comparison, not a full cosmological inference.

## Scripts

```bash
python3 maat_hz_data_comparison_v01.py
python3 maat_hz_chi2_fit_v01.py
```

## Inputs

`maat_hz_data_comparison_v01.py` reads the v0.10 observable projection from:

```text
../maat_observable_predictions_v10/outputs/maat_v10_observables.csv
```

Both scripts include the compact Cosmic Chronometer table internally for
reproducibility.

## Outputs

`maat_hz_data_comparison_v01.py` writes to:

```text
maat_hz_data_comparison_v01/
```

Main outputs:

- `maat_vs_chronometer_Hz_table.csv`
- `maat_vs_chronometer_summary.json`
- `fig1_MAAT_vs_Hz_data.png`
- `fig2_Hz_pulls.png`
- `fig3_MAAT_relative_H_deviation.png`

`maat_hz_chi2_fit_v01.py` writes to:

```text
maat_hz_chi2_fit_v01/
```

Main outputs:

- `maat_hz_chi2_scan_results.csv`
- `maat_hz_chi2_fit_summary.json`
- `best_maat_trajectory.csv`
- `fig1_best_MAAT_vs_Hz_data.png`
- `fig2_chi2_heatmap.png`
- `fig3_best_relative_H_deviation.png`

## Key Numbers

- Chronometer points: `31`
- LCDM chi-square: `14.8759`
- Fixed v0.10 MAAT chi-square: `15.3661`
- Best scan MAAT chi-square: `15.3674`
- Best scan reduced chi-square with `k=2`: `0.5299`
- Stable scan models: `581 / 616`
- Best scan parameters: `mu = 15.8489`, `phidot0 = 2.7889`
- Best scan maximum `Omega_MAAT`: `0.02668`
- Best scan maximum relative Hubble deviation: `0.02347`

## Interpretation

The MAAT branch is not favoured over LCDM in this first diagnostic comparison.
However, it remains dynamically stable, subdominant in the energy budget, and
close to the LCDM chronometer fit. This makes the result useful as a
falsifiability bridge rather than a detection claim.

## Data and License Note

The scripts embed a compact Cosmic Chronometer `H(z)` table compiled from the
published literature. No Hubble Space Telescope image/data products, DESI/BAO
catalogues, CMB likelihoods, or weak-lensing datasets are bundled here.

The repository license applies only to the code, generated tables, and
generated figures produced by the scripts. It does not relicense third-party
observational measurements. The observational `H(z)` values are treated as
published scientific data and should be cited to the original Cosmic
Chronometer publications and compilations when reused or discussed.
