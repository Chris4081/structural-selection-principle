# MAAT v0.10 Observable Predictions Layer

This folder contains the reproducibility package for the MAAT v0.10 observable
projection benchmark.

The script converts the representative stable v0.9 FLRW trajectory into
observable cosmology proxies:

- expansion history `E(z)=H(z)/H0`,
- relative Hubble deviation from an internally normalised LCDM reference,
- MAAT equation of state `w_MAAT(z)`,
- MAAT density fraction `Omega_MAAT(z)`,
- approximate growth proxy `f sigma_8(z)`.

This is a diagnostic projection, not a fit to cosmological data and not a full
linear-perturbation calculation.

## Reproduce

From this folder:

```bash
python3 maat_observable_predictions_v10.py
```

The script reads:

```text
../maat_dynamic_fields_v05_v09/v09_flrw_stability_scan/outputs/representative_flrw_trajectory.csv
```

and writes all outputs to:

```text
outputs/
```

## Outputs

| File | Description |
|------|-------------|
| `maat_v10_observables.csv` | redshift grid and observable columns |
| `maat_v10_summary.json` | numerical summary of the v0.10 projection |
| `observable_predictions_summary.png` | four-panel summary figure |
| `hubble_history_vs_lcdm.png` | expansion history comparison |
| `relative_hubble_deviation.png` | relative Hubble deviation |
| `maat_equation_of_state.png` | MAAT equation-of-state history |
| `maat_density_fraction.png` | MAAT density fraction |
| `growth_proxy_fsigma8.png` | approximate growth proxy |

## Key Numbers

- Redshift range: `0 <= z <= 2.329`
- Reference LCDM fractions: `Omega_m0 = 0.3007`, `Omega_Lambda0 = 0.6992`
- Present MAAT density fraction: `Omega_MAAT0 = 0.00331`
- Present MAAT equation of state: `w_MAAT0 = -0.801`
- Maximum `|Delta H/H|`: `0.0399`
- Mean `|Delta H/H|`: `0.00481`
- Maximum `Omega_MAAT`: `0.0291`
- Maximum `|Delta f sigma_8 / f sigma_8|`: `0.0336`

## Status

The v0.10 layer is the first observable-signature test of the v0.9 stable
branch. It should be read as a falsifiability scaffold: it defines measurable
quantities that later versions can confront with real expansion and growth
data.

## Data and License Note

This benchmark does not bundle or use external Hubble Space Telescope data,
observational `H(z)` catalogues, DESI/BAO measurements, CMB data, or
weak-lensing data. The quantity `H(z)` is the Hubble parameter reconstructed
from the synthetic v0.9 toy FLRW trajectory, and the LCDM curve is an
internally normalised reference model.

The repository license applies to the code, generated toy-model outputs, and
figures in this folder. If future versions compare against real observational
datasets, the corresponding survey/data-release citations and license terms
must be added next to the imported data files.
