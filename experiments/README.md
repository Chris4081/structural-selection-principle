# Phenomenological String-Landscape Selection Experiments

This directory contains the reproducibility bundle for the paper

**A Phenomenological Structural Selection Measure for String Backgrounds**.

The experiments are exploratory benchmark tests for extending the structural
selection framework to string-background ranking. They are not claimed to be a
complete string-theoretic derivation, a unique vacuum measure, a Standard Model
compactification, or a final Theory of Everything. Their purpose is narrower:
to test whether layered structural ranking carries information beyond
energy-only ordering in progressively more structured toy and bridge ensembles.

## Layout

```text
experiments/
├── phenomenological_structural_selection_string_measure.pdf
├── cosmology_structural_selection/
│   ├── README.md
│   ├── maat_cosmology_toy_v2.py
│   ├── maat_cosmology_toy_v2_results.json
│   └── maat_cosmology_toy_v2_plots/
├── fixed_energy_structural_selection/
│   ├── README.md
│   ├── structural_selection_fixed_energy_benchmarks.py
│   ├── fixed_energy_structural_selection_results.json
│   └── fixed_energy_structural_selection_plots/
├── boundary_aware_lambda_calibration/
│   ├── README.md
│   ├── merge_defects.py
│   ├── fit_closed_maat_lambda_v1.py
│   ├── plot_closed_maat_lambda_v2.py
│   ├── closed_maat_lambda_fit_results.json
│   ├── maat_defects_*.csv
│   └── plots/
├── cosmological_cci/
│   ├── README.md
│   ├── maat_cci_cosmology_v02.py
│   ├── maat_cci_cosmology_v02_chronometers.csv
│   ├── maat_cci_cosmology_v02.csv
│   ├── maat_cci_cosmology_v02_data_comparison.csv
│   └── plots/
├── natural_constants_selection/
│   ├── README.md
│   ├── naturkonstante.py
│   ├── naturkonstante_v2_structure_scan.py
│   ├── ...
│   ├── naturkonstante_v10_multiseed.py
│   ├── naturkonstante_v12_maxent_lambda.py
│   ├── naturkonstante_v13_maxent_sm_bridge.py
│   ├── naturkonstante_v12_maxent_lambda_results.json
│   ├── naturkonstante_v13_maxent_sm_bridge_results.json
│   └── plots/
├── standard_model_bridge/
│   ├── README.md
│   ├── standard_model_rg_maat_bridge.py
│   ├── standard_model_rg_maat_summary_figure.py
│   ├── standard_model_rg_maat_v11_holdout.py
│   ├── sm_bridge_result_section.tex
│   ├── sm_bridge_v11_holdout_section.tex
│   ├── standard_model_rg_maat_results.json
│   ├── standard_model_rg_maat_v11_holdout_results.json
│   └── standard_model_rg_maat_plots/
└── string_landscape_selection/
    ├── structural_selection_10d_tadpole_toy.py
    ├── structural_selection_iib_kklt_scan.py
    ├── structural_selection_iib_kklt_bridge.py
    ├── structural_selection_iib_period_kklt_bridge.py
    ├── structural_selection_iib_exact_period_kklt_bridge.py
    ├── structural_selection_iib_backreaction_sm_bridge.py
    ├── *_results.json
    └── *_plots/
```

## Requirements

The scripts require Python 3 and the following packages:

```bash
python3 -m pip install numpy scipy matplotlib mpmath
```

No network access is required once the Python packages are installed.

## Reproducing the String-Landscape Results

Run the scripts from the `string_landscape_selection/` directory:

```bash
cd experiments/string_landscape_selection

python3 structural_selection_10d_tadpole_toy.py
python3 structural_selection_iib_kklt_bridge.py
python3 structural_selection_iib_period_kklt_bridge.py
python3 structural_selection_iib_exact_period_kklt_bridge.py
python3 structural_selection_iib_backreaction_sm_bridge.py
```

Each script writes a JSON result file and a plot directory in the same folder.
The default seeds are fixed inside the scripts, so the benchmark outputs are
deterministic for the same Python/package versions.

## Reproducing the Fixed-Energy Field Tests

The `fixed_energy_structural_selection/` directory contains a separate
field-theory benchmark bundle testing whether structural ranking differs from
energy-only ranking in `phi^4` and Sine-Gordon systems.

Run:

```bash
cd experiments/fixed_energy_structural_selection
python3 structural_selection_fixed_energy_benchmarks.py
```

This generates:

- `fixed_energy_structural_selection_results.json`
- `fixed_energy_structural_selection_plots/`

The 2D plots use a 2x2 constrained layout with separated labels to avoid
overlap in the domain-wall panels and energy-vs-structure comparison.

## Reproducing the Cosmology Toy Benchmark

The `cosmology_structural_selection/` directory contains the FLRW scalar-field
toy benchmark extending structural selection from static field configurations
to simple cosmological histories.

Run:

```bash
cd experiments/cosmology_structural_selection
python3 maat_cosmology_toy_v2.py
```

This generates:

- `maat_cosmology_toy_v2_results.json`
- `maat_cosmology_toy_v2_plots/`

## Reproducing the Boundary-Aware Lambda Calibration

The `boundary_aware_lambda_calibration/` directory contains the fused
SAT/MAAT-Core/boundary defect benchmark used to calibrate closed MAAT
structural weights.

Run:

```bash
cd experiments/boundary_aware_lambda_calibration
python3 merge_defects.py
python3 fit_closed_maat_lambda_v1.py
python3 plot_closed_maat_lambda_v2.py
```

This generates:

- `maat_defects_fused.csv`
- `closed_maat_lambda_fit_results.json`
- `plots/fig1_lambda_distribution.png`
- `plots/fig2_lambda_shares.png`
- `plots/fig3_sat_vs_core_defects.png`
- `plots/fig4_structural_energy_distribution.png`
- `plots/fig5_c_hat_vs_flambda_sat.png`

## Reproducing the Cosmological CCI Observable

The `cosmological_cci/` directory contains the Paper 28 diagnostic pipeline
for the Cosmological Critical Coherence Index.

Run:

```bash
cd experiments/cosmological_cci
python3 maat_cci_cosmology_v02.py
```

This generates:

- `maat_cci_cosmology_v02.csv`
- `maat_cci_cosmology_v02_data_comparison.csv`
- `plots/maat_cci_cosmology_v02_plot.png`
- `plots/maat_cci_cosmology_v02_Hz_comparison.png`
- `plots/maat_cci_cosmology_v02_data_comparison.png`
- `plots/maat_cci_cosmology_v02_residuals.png`

**Data attribution and license note:** The Planck-2018 parameter values and
Cosmic Chronometer H(z) measurements are external scientific data from the
cited literature. The repository contains a compact chronometer table and
derived CSV/PNG artifacts only for reproducibility of this diagnostic. Reuse of
the original measurements remains subject to the terms of the original
publications, journals, and collaborations. No endorsement by the Planck
Collaboration or the chronometer-data authors is implied.

## Reproducing the Natural-Constants Benchmark

The `natural_constants_selection/` directory contains the v1--v13
natural-constants benchmark sequence leading into the Standard-Model bridge.
The scripts test whether broad MAAT-type structural diagnostics define stable
low-defect basins in effective-constant space. Versions v12 and v13 add a
maximum-entropy calibration of the sector weights `lambda_a`.

Run:

```bash
cd experiments/natural_constants_selection
python3 naturkonstante.py
python3 naturkonstante_v2_structure_scan.py
python3 naturkonstante_v3_physics_constraints.py
python3 naturkonstante_v4_stellar_chemistry.py
python3 naturkonstante_v5_robustness_scan.py
python3 naturkonstante_v6_ablation_scan.py
python3 naturkonstante_v7_landscape_heatmap.py
python3 naturkonstante_v8_gradient_flow.py
python3 naturkonstante_v9_rg_maat.py
python3 naturkonstante_v10_multiseed.py
python3 naturkonstante_v12_maxent_lambda.py
python3 naturkonstante_v13_maxent_sm_bridge.py
```

Publication figures copied from the local workspace are stored in:

- `natural_constants_selection/plots/maat_constants_v7_heatmap.png`
- `natural_constants_selection/plots/maat_constants_v8_gradient_flow.png`

## Reproducing the Standard-Model Bridge Benchmark

The `standard_model_bridge/` directory contains a phenomenological bridge test
between one-loop Standard-Model-like RG flow and MAAT structural selection.
Observed Standard Model values are used only as comparison markers, not as
optimisation targets in the score.

Run:

```bash
cd experiments/standard_model_bridge
python3 standard_model_rg_maat_bridge.py
python3 standard_model_rg_maat_summary_figure.py
python3 standard_model_rg_maat_v11_holdout.py
```

This generates:

- `standard_model_rg_maat_results.json`
- `standard_model_rg_maat_v11_holdout_results.json`
- `standard_model_rg_maat_plots/`
including the publication-style summary figure
`standard_model_rg_maat_plots/sm_bridge_nature_summary.png` and the v11
holdout plots.

## String-Landscape Script Overview

| Script | Purpose | Main outputs |
| --- | --- | --- |
| `structural_selection_10d_tadpole_toy.py` | Type-IIB-inspired flux toy ensemble with an explicit D3 tadpole proxy. | `structural_selection_10d_tadpole_toy_results.json`, `structural_selection_10d_tadpole_toy_plots/` |
| `structural_selection_iib_kklt_scan.py` | Shared reduced KKLT potential and structural diagnostic utilities. | Helper module used by later scripts. |
| `structural_selection_iib_kklt_bridge.py` | Flux-to-KKLT bridge ensemble with stationarity and stability tests. | `structural_selection_iib_kklt_bridge_results.json`, `structural_selection_iib_kklt_bridge_plots/` |
| `structural_selection_iib_period_kklt_bridge.py` | Large-complex-structure period proxy connected to KKLT ranking. | `structural_selection_iib_period_kklt_bridge_results.json`, `structural_selection_iib_period_kklt_bridge_plots/` |
| `structural_selection_iib_exact_period_kklt_bridge.py` | Mirror-quintic Picard-Fuchs period benchmark connected to KKLT ranking. | `structural_selection_iib_exact_period_kklt_bridge_results.json`, `structural_selection_iib_exact_period_kklt_bridge_plots/` |
| `structural_selection_iib_backreaction_sm_bridge.py` | Adds phenomenological backreaction and Standard-Model-sector proxy layers. | `structural_selection_iib_backreaction_sm_bridge_results.json`, `structural_selection_iib_backreaction_sm_bridge_plots/` |

## Fixed-Energy Benchmark Overview

| Script | Purpose | Main outputs |
| --- | --- | --- |
| `fixed_energy_structural_selection/structural_selection_fixed_energy_benchmarks.py` | Consolidated 1D `phi^4`, 2D `phi^4`, equal-energy `phi^4`, and Sine-Gordon structural-selection tests. | `fixed_energy_structural_selection_results.json`, `fixed_energy_structural_selection_plots/` |

## Cosmology Benchmark Overview

| Script | Purpose | Main outputs |
| --- | --- | --- |
| `cosmology_structural_selection/maat_cosmology_toy_v2.py` | Flat-FLRW scalar-field toy histories ranked by MAAT structural diagnostics. | `maat_cosmology_toy_v2_results.json`, `maat_cosmology_toy_v2_plots/` |

## Natural-Constants Benchmark Overview

| Script | Purpose | Main outputs |
| --- | --- | --- |
| `natural_constants_selection/naturkonstante.py` | Initial direct structural optimisation over dimensionless effective constants. | Console output / script-defined outputs |
| `natural_constants_selection/naturkonstante_v2_structure_scan.py` | Log-space structure scan. | Script-defined outputs |
| `natural_constants_selection/naturkonstante_v3_physics_constraints.py` | Broad physics viability bands. | Script-defined outputs |
| `natural_constants_selection/naturkonstante_v4_stellar_chemistry.py` | Chemistry and fusion proxies. | Script-defined outputs |
| `natural_constants_selection/naturkonstante_v5_robustness_scan.py` | Multi-seed robustness. | Script-defined outputs |
| `natural_constants_selection/naturkonstante_v6_ablation_scan.py` | Sector ablation tests. | Script-defined outputs |
| `natural_constants_selection/naturkonstante_v7_landscape_heatmap.py` | Basin visualisation. | `plots/maat_constants_v7_heatmap.png` |
| `natural_constants_selection/naturkonstante_v8_gradient_flow.py` | Gradient-flow attractor test. | `plots/maat_constants_v8_gradient_flow.png` |
| `natural_constants_selection/naturkonstante_v9_rg_maat.py` | Toy RG plus MAAT score. | Script-defined outputs |
| `natural_constants_selection/naturkonstante_v10_multiseed.py` | RG multi-seed robustness. | Script-defined outputs |
| `natural_constants_selection/naturkonstante_v12_maxent_lambda.py` | Maximum-entropy calibration of sector weights. | `naturkonstante_v12_maxent_lambda_results.json` |
| `natural_constants_selection/naturkonstante_v13_maxent_sm_bridge.py` | MaxEnt-weighted constants bridge. | `naturkonstante_v13_maxent_sm_bridge_results.json` |

## Boundary-Aware Lambda Calibration Overview

| Script | Purpose | Main outputs |
| --- | --- | --- |
| `boundary_aware_lambda_calibration/merge_defects.py` | Merges SAT, MAAT-Core, and boundary defect CSVs into a fused calibration ensemble. | `maat_defects_fused.csv` |
| `boundary_aware_lambda_calibration/fit_closed_maat_lambda_v1.py` | Fits closed MAAT sector weights over the fused defect ensemble. | `closed_maat_lambda_fit_results.json` |
| `boundary_aware_lambda_calibration/plot_closed_maat_lambda_v2.py` | Generates lambda, share, defect-comparison, and structural-energy plots from the fitted result. | `plots/fig*.png` |

## Cosmological CCI Overview

| Script | Purpose | Main outputs |
| --- | --- | --- |
| `cosmological_cci/maat_cci_cosmology_v02.py` | Generates the cosmological CCI model grid, chronometer H(z) projection, residual table, and plots. | `maat_cci_cosmology_v02.csv`, `maat_cci_cosmology_v02_data_comparison.csv`, `plots/*.png` |

## Standard-Model Bridge Overview

| Script | Purpose | Main outputs |
| --- | --- | --- |
| `standard_model_bridge/standard_model_rg_maat_bridge.py` | One-loop Standard-Model-like RG bridge from UV parameters to IR effective observables ranked by MAAT structural diagnostics. | `standard_model_rg_maat_results.json`, `standard_model_rg_maat_plots/` |
| `standard_model_bridge/standard_model_rg_maat_summary_figure.py` | Builds the four-panel publication-style summary figure from the benchmark outputs. | `standard_model_rg_maat_plots/sm_bridge_nature_summary.png`, `standard_model_rg_maat_plots/sm_bridge_nature_summary.pdf` |
| `standard_model_bridge/standard_model_rg_maat_v11_holdout.py` | Direct-term holdout benchmark testing cross-sector predictivity for selected SM-like observables. | `standard_model_rg_maat_v11_holdout_results.json`, `standard_model_rg_maat_plots/sm_bridge_v11_holdout_*.png` |

## Scientific Status

These scripts implement toy, bridge, and proxy models. In particular:

- The tadpole and backreaction terms are operational proxies, not full 10D
  supergravity solutions.
- The Standard-Model-sector layer is a phenomenological diagnostic proxy, not a
  constructed chiral compactification.
- The Standard-Model bridge uses broad viability bands and one-loop toy
  matching; it is not a precision electroweak fit or a derivation of the
  observed constants.
- The natural-constants benchmark uses broad phenomenological viability bands;
  it identifies structural basins rather than deriving constants from first
  principles.
- The v12/v13 MaxEnt weights are effective benchmark weights calibrated on a
  synthetic defect ensemble; they are not unique microscopic constants.
- The boundary-aware lambda calibration is a fused benchmark over SAT and
  MAAT-Core-derived data. It shows constraint dominance in that closed system,
  but does not establish universal MAAT constants.
- The cosmological CCI benchmark is a diagnostic projection over Planck LCDM
  and Cosmic Chronometer H(z) data. It is not a new cosmological model and does
  not perform parameter inference.
- The v11 holdout benchmark removes direct score terms only; indirect
  cross-sector appearances of held-out observables remain active by design.
- The Picard-Fuchs benchmark uses a controlled mirror-quintic period model, but
  does not constitute a full global compactification scan.
- The fixed-energy field tests use operational structural diagnostics, not a
  first-principles microscopic derivation.
- The cosmology benchmark is a flat-FLRW scalar-field toy model, not a complete
  cosmological theory.
- The structural measure is intended as a testable ranking architecture, not as
  a final microscopic string measure.

This status is intentional. The experiments provide a compact reproducibility
loop:

```text
model assumptions -> generated ensemble -> structural ranking -> JSON results -> plots
```

## Citation / Paper Reference

If using or discussing these experiments, cite the accompanying paper:

Christof Krieg,
*A Phenomenological Structural Selection Measure for String Backgrounds*,
2026.

For the natural-constants and Standard-Model bridge benchmarks, cite:

Christof Krieg,
*Structural Selection of Effective Constants: From MAAT Basins to
MaxEnt-Weighted RG Bridge Tests*,
2026.

For the boundary-aware lambda calibration benchmark, cite:

Christof Krieg,
*Boundary-Aware Calibration of MAAT Structural Weights: A Reproducible
Benchmark for Constraint-Dominated Structural Selection*,
2026.

For the cosmological CCI benchmark, cite:

Christof Krieg,
*Cosmological Critical Coherence Index: A Structural-Stress Observable for
Cosmic Evolution*,
2026.

For the external cosmology inputs used in Paper 28, cite the Planck-2018
cosmological-parameter paper and the Cosmic Chronometer literature referenced
in the manuscript.

Reference comparison values for fundamental constants and Standard-Model
inputs should cite CODATA/NIST and PDG 2024:

- NIST/CODATA values of the fundamental physical constants:
  <https://www.nist.gov/programs-projects/codata-values-fundamental-physical-constants>
- Particle Data Group, Review of Particle Physics 2024:
  <https://pdg.lbl.gov/index-2024.html>
