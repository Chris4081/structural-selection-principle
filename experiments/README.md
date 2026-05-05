# Phenomenological String-Landscape Selection Experiments

This directory contains the reproducibility bundle for the phenomenological
string-landscape papers

**A Phenomenological Structural Selection Measure for String Backgrounds**.

and

**Structural Selection in the String Landscape: A MAAT-Based Phenomenological
Framework for Vacuum Ranking**.

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
├── ../documentation/Structural_Selection_in_the_String_Landscape_MAAT_Framework.pdf
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
├── lambda_response_closure/
│   ├── README.md
│   ├── lambda_response_closure.py
│   ├── lambda_response_closure_results.json
│   └── plots/
├── cosmological_cci/
│   ├── README.md
│   ├── maat_cci_cosmology_v02.py
│   ├── maat_cci_cosmology_v02_chronometers.csv
│   ├── maat_cci_cosmology_v02.csv
│   ├── maat_cci_cosmology_v02_data_comparison.csv
│   └── plots/
├── cosmological_cci_v03/
│   ├── README.md
│   ├── maat_cci_cosmology_v03_growth.py
│   ├── maat_cci_cosmology_v03_chronometers.csv
│   ├── maat_cci_cosmology_v03_fsigma8.csv
│   ├── maat_cci_cosmology_v03_grid.csv
│   ├── maat_cci_cosmology_v03_results.json
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
├── maat_v121_observables_stability_paper37/
│   ├── README.md
│   ├── paper37_observables_emergent_robustness.py
│   ├── paper37_stability_landscape.py
│   ├── observable_outputs/
│   ├── stability_landscape_outputs/
│   └── sat_validation/
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
python3 -m pip install numpy scipy pandas matplotlib mpmath
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

## Reproducing the Lambda Response Closure

The `lambda_response_closure/` directory contains the response-theoretic
closure of MAAT sector weights.  It derives effective `lambda_a` values from
the covariance geometry of the fused defect ensemble:

```text
lambda = (Cov_mu0[d] + eta tr(C)/A I)^(-1)(<d>_mu0 - <d>_target)
```

Run:

```bash
cd experiments/lambda_response_closure
python3 lambda_response_closure.py
```

This generates:

- `lambda_response_closure_results.json`
- `plots/lambda_response_shares.png`
- `plots/lambda_response_correlation.png`
- `plots/lambda_response_target_match.png`
- `plots/lambda_response_vs_closed_fit.png`

Main result: `R` dominance emerges as a covariance-response effect when the
target ensemble encodes safety or boundary stability, rather than being imposed
as a fixed assumption.

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

## Reproducing the Cosmological CCI v0.3 Multi-Sector Benchmark

The `cosmological_cci_v03/` directory contains the Paper 29 benchmark that
extends the cosmological CCI with measured growth connectivity, robustness
margins, maximum-entropy companion weights, and a data-driven transition proxy.

Run:

```bash
cd experiments/cosmological_cci_v03
python3 maat_cci_cosmology_v03_growth.py
```

This generates:

- `maat_cci_cosmology_v03_grid.csv`
- `maat_cci_cosmology_v03_results.json`
- `plots/fsigma8_growth_connectivity.png`
- `plots/v_r_supports.png`
- `plots/cci_v02_v03_comparison.png`
- `plots/transition_curvature.png`
- `plots/lambda_calibration.png`

**Data attribution and license note:** The benchmark uses Planck-2018
reference parameters, a compact Cosmic Chronometer H(z) table, and BOSS DR12
consensus f sigma_8 values from the cited literature. The repository contains
derived CSV/PNG artifacts only for reproducibility. Reuse of the original
measurements remains subject to the terms of the original publications,
journals, and collaborations. No endorsement by the Planck Collaboration,
SDSS/BOSS Collaboration, or the chronometer-data authors is implied.

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

## Reproducing Paper 37 v1.2.1 Observable and Stability Benchmarks

The `maat_v121_observables_stability_paper37/` directory contains the Paper 37
v1.2.1 companion benchmark. It reproduces the baseline observable proxy run,
the two-parameter stability landscape, and the SAT correlation table used in
the empirical discussion of emergent robustness.

Run:

```bash
cd experiments/maat_v121_observables_stability_paper37
python3 paper37_observables_emergent_robustness.py
python3 paper37_stability_landscape.py
cd sat_validation
python3 maat_v121_sat_validation.py maat_cosmos_full_results.csv
```

This generates:

- `observable_outputs/maat_v121_observable_predictions.csv`
- `observable_outputs/maat_v121_observable_summary.json`
- `stability_landscape_outputs/paper37_maat_v121_landscape_scan.csv`
- `stability_landscape_outputs/paper37_landscape_summary.json`
- `sat_validation/maat_v121_sat_validation_results.csv`
- `sat_validation/maat_v121_sat_correlations.csv`
- all figures used in the Paper 37 proxy and scan sections

## Reproducing Paper 38 v1.2.1 Linear-Growth Closure Benchmark

The `maat_paper38_v121_robustness_closure/` directory contains the Paper 38
linear-growth closure benchmark. It updates the Paper 35 growth pipeline with
the v1.2.1 closure convention, tests selection-field perturbations, checks
positivity of the lambda-relaxation equation, and exports the closure
diagnostics.

Run:

```bash
cd experiments/maat_paper38_v121_robustness_closure
python3 maat_paper38_v121_robustness_closure.py
```

This generates:

- `paper38_v121_outputs/paper38_summary.json`
- `paper38_v121_outputs/paper38_growth_curves.csv`
- `paper38_v121_outputs/paper38_v121_structural_closure.csv`
- `paper38_v121_outputs/paper38_dlambda_evolution.csv`
- `paper38_v121_outputs/fig_paper38_growth.png`
- `paper38_v121_outputs/fig_paper38_v121_closure.png`
- `paper38_v121_outputs/fig_paper38_perturbations.png`

## Reproducing Paper 39 v1.2.1 Observable Growth Signature Proxy

The `maat_paper39_observable_growth_signature/` directory contains the Paper 39
observable-signature proxy. It modulates a Planck-normalised `f sigma_8`
baseline by a bounded MAAT projection template, scans a small epsilon interval,
and exports the v1.2.1 robustness closure diagnostics.

Run:

```bash
cd experiments/maat_paper39_observable_growth_signature
python3 maat_paper39_observable_signature_v121.py
```

This generates:

- `outputs_paper39/paper39_summary.json`
- `outputs_paper39/paper39_growth_signature_results.csv`
- `outputs_paper39/paper39_epsilon_chi2_scan.csv`
- `outputs_paper39/paper39_projection_template_v121.csv`
- `outputs_paper39/fig1_projection_template.png`
- `outputs_paper39/fig2_growth_comparison.png`
- `outputs_paper39/fig3_residuals.png`
- `outputs_paper39/fig4_chi2_scan.png`
- `outputs_paper39/fig5_v121_robustness_closure.png`
- `outputs_paper39/fig6_v121_cci.png`
- `outputs_paper39/fig7_v121_support_fields.png`

Main diagnostic results:

- Growth comparison points: `13`
- Best epsilon: `-0.0100`
- LCDM chi2: `12.4373`
- MAAT proxy chi2: `12.3772`
- Delta chi2: `-0.0601`
- Max `|Delta f sigma_8 / f sigma_8|`: `0.9891%`
- Mean `R_rob`: `0.6673`

The best epsilon value lies at the scan boundary. The benchmark should
therefore be read as a compatibility/signature test, not as a parameter
measurement.

## Reproducing Paper 40 v1.2.1 Structural Signature Test

The `maat_paper40_structural_signature_test/` directory contains the Paper 40
residual-structure diagnostic. It compares MAAT v1.2.1 projection and CCI
diagnostics against signed and absolute `f sigma_8` residuals relative to a
Planck-normalised LCDM baseline, then runs permutation null tests.

Run:

```bash
cd experiments/maat_paper40_structural_signature_test
python3 maat_paper40_structural_signature_test.py
```

This generates:

- `outputs_paper40/paper40_summary.json`
- `outputs_paper40/paper40_signature_table.csv`
- `outputs_paper40/fig1_cci_diag_vs_signed_residual.png`
- `outputs_paper40/fig2_scatter_cci_diag_signed_residual.png`
- `outputs_paper40/fig3_scatter_cci_diag_abs_residual.png`
- `outputs_paper40/fig4_spearman_signed_residuals.png`
- `outputs_paper40/fig5_spearman_abs_residuals.png`
- `outputs_paper40/fig6_breakthrough_null_test_cci_diag.png`

Main diagnostic results:

- Growth comparison points: `13`
- Spearman `CCI_diag` vs `|residual_sigma|`: `0.5934`, `p = 0.0338`
- Random-field null for `CCI_diag`: `p = 0.0363`
- Redshift-shuffle null for `CCI_diag`: `p = 0.0359`
- Spearman `R_proj` vs signed residual: `0.6319`, `p = 0.0228`
- Spearman `V` vs signed residual: `-0.6319`, `p = 0.0228`

Because the balance support `B` is residual-sensitive, this benchmark is a
semi-supervised structural consistency test, not a blind prediction or
detection claim.

## Reproducing the SO(10)-Motivated Structural Selection Benchmark

The `maat_so10_structural_selection/` directory contains the extra
phenomenological SO(10)-motivated bridge benchmark. It tests MAAT v1.2.1 as a
structural-selection layer over one-loop gauge running and a third-generation
Yukawa benchmark with an imposed high-scale `b-tau` boundary condition.

Run:

```bash
cd experiments/maat_so10_structural_selection
python3 maat_so10_gauge_2d_fit.py
python3 maat_so10_yukawa_coupled_rg_v04.py
```

This generates:

- `outputs/gauge_2d_summary.json`
- `outputs/yukawa_v04_summary.json`
- `outputs/maat_selection_landscape.png`

Main diagnostic results:

- Gauge benchmark: `alpha_GUT = 0.022676468338`,
  `M_GUT = 6.639198e15 GeV`, `SM chi2 = 100.789861`,
  `R_rob = 0.60301490`.
- Yukawa benchmark: `M_GUT = 1.858006e16 GeV`,
  `Delta_b = 0.0506451697`, `chi2_yukawa = 0.00024341`,
  `R_rob = 0.99916912`.

The gauge benchmark is a one-loop diagnostic and does not achieve precision
unification. The Yukawa benchmark is an SO(10)-motivated compatibility test
with a fitted bottom-threshold correction, not a first-principles derivation
of fermion masses.

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

## Lambda Response Closure Overview

| Script | Purpose | Main outputs |
| --- | --- | --- |
| `lambda_response_closure/lambda_response_closure.py` | Derives effective MAAT sector weights from defect covariance and target-response geometry. | `lambda_response_closure_results.json`, `plots/lambda_response_*.png` |

## Cosmological CCI Overview

| Script | Purpose | Main outputs |
| --- | --- | --- |
| `cosmological_cci/maat_cci_cosmology_v02.py` | Generates the cosmological CCI model grid, chronometer H(z) projection, residual table, and plots. | `maat_cci_cosmology_v02.csv`, `maat_cci_cosmology_v02_data_comparison.csv`, `plots/*.png` |
| `cosmological_cci_v03/maat_cci_cosmology_v03_growth.py` | Extends the cosmological CCI with f sigma_8 growth connectivity, robustness margins, MaxEnt companion weights, and a curvature-derived transition proxy. | `maat_cci_cosmology_v03_grid.csv`, `maat_cci_cosmology_v03_results.json`, `plots/*.png` |

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
- The lambda response closure gives an effective physical/statistical
  derivation of `lambda_a` as response coefficients of a specified defect
  ensemble. It does not prove that the resulting weights are universal
  microscopic constants.
- The cosmological CCI benchmark is a diagnostic projection over Planck LCDM
  and Cosmic Chronometer H(z) data. It is not a new cosmological model and does
  not perform parameter inference.
- The cosmological CCI v0.3 benchmark adds BOSS DR12 f sigma_8 growth
  connectivity and robustness margins, but remains a diagnostic benchmark, not
  a full cosmological likelihood or model comparison.
- The v11 holdout benchmark removes direct score terms only; indirect
  cross-sector appearances of held-out observables remain active by design.
- The Picard-Fuchs benchmark uses a controlled mirror-quintic period model, but
  does not constitute a full global compactification scan.
- The fixed-energy field tests use operational structural diagnostics, not a
  first-principles microscopic derivation.
- The cosmology benchmark is a flat-FLRW scalar-field toy model, not a complete
  cosmological theory.
- The Paper 37 v1.2.1 benchmark uses proxy observable and SAT-correlation
  diagnostics. It validates internal consistency of the closure convention; it
  is not a precision cosmological likelihood or a proof of universal SAT
  hardness laws.
- The Paper 38 v1.2.1 benchmark is a linear-growth and selection-field
  perturbation consistency check. It is not a Boltzmann-code calculation and
  not an observational detection.
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

For the updated MAAT string-landscape ranking framework, cite:

Christof Krieg,
*Structural Selection in the String Landscape: A MAAT-Based Phenomenological
Framework for Vacuum Ranking*,
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

For the response-theoretic lambda closure, cite:

Christof Krieg,
*Response-Based Derivation of MAAT Structural Weights: From Covariance Geometry
to Selection Pressure*,
2026.

For the cosmological CCI benchmark, cite:

Christof Krieg,
*Cosmological Critical Coherence Index: A Structural-Stress Observable for
Cosmic Evolution*,
2026.

For the cosmological CCI v0.3 multi-sector benchmark, cite:

Christof Krieg,
*Cosmological Critical Coherence Index with Growth Connectivity and Robustness
Margins: A Multi-Sector Structural-Stress Benchmark Using H(z) and f sigma_8(z)*,
2026.

For the external cosmology inputs used in Papers 28 and 29, cite the
Planck-2018 cosmological-parameter paper, the Cosmic Chronometer literature,
and the BOSS DR12 consensus analysis referenced in the manuscripts.

For the SO(10)-motivated structural selection benchmark, cite:

Christof Krieg,
*Structural Selection in SO(10)-Motivated Unified Field Theories:
A Phenomenological MAAT Layer for Gauge and Yukawa Benchmarks*,
2026.

Reference comparison values for fundamental constants and Standard-Model
inputs should cite CODATA/NIST and PDG 2024:

- NIST/CODATA values of the fundamental physical constants:
  <https://www.nist.gov/programs-projects/codata-values-fundamental-physical-constants>
- Particle Data Group, Review of Particle Physics 2024:
  <https://pdg.lbl.gov/index-2024.html>
