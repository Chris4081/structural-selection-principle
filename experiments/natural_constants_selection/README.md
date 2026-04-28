# Natural-Constants Structural Selection Benchmark

This directory contains the natural-constants benchmark sequence used in

**Structural Selection of Effective Constants: From MAAT Basins to
MaxEnt-Weighted RG Bridge Tests**.

The scripts are phenomenological toy and bridge tests. They are not a
first-principles derivation of the constants, not a precision Standard Model
fit, and not a solution of the cosmological-constant problem. Their purpose is
narrower: to test whether MAAT-type structural scores define stable low-defect
basins in effective-constant space.

## Pipeline

```text
dimensionless constants
-> structural viability sectors
-> MaxEnt sector-weight calibration
-> MAAT score
-> basin / robustness / gradient-flow diagnostics
-> generated plots
```

Later versions connect this natural-constants scan to the
`standard_model_bridge/` benchmark, where one-loop Standard-Model-like RG flow
maps UV parameters to IR effective observables.

## Run

Run the scripts from this directory:

```bash
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

The publication figures copied from the local workspace are stored in:

- `plots/maat_constants_v7_heatmap.png`
- `plots/maat_constants_v8_gradient_flow.png`

The MaxEnt scripts generate:

- `naturkonstante_v12_maxent_lambda_results.json`
- `naturkonstante_v13_maxent_sm_bridge_results.json`

## Version Overview

| Version | Script | Purpose |
| --- | --- | --- |
| v1 | `naturkonstante.py` | Initial direct structural optimisation. |
| v2 | `naturkonstante_v2_structure_scan.py` | Log-space structure scan. |
| v3 | `naturkonstante_v3_physics_constraints.py` | Broad physics viability bands. |
| v4 | `naturkonstante_v4_stellar_chemistry.py` | Chemistry and fusion proxies. |
| v5 | `naturkonstante_v5_robustness_scan.py` | Multi-seed robustness. |
| v6 | `naturkonstante_v6_ablation_scan.py` | Sector ablation tests. |
| v7 | `naturkonstante_v7_landscape_heatmap.py` | Basin visualisation. |
| v8 | `naturkonstante_v8_gradient_flow.py` | Gradient-flow attractor test. |
| v9 | `naturkonstante_v9_rg_maat.py` | Toy RG plus MAAT score. |
| v10 | `naturkonstante_v10_multiseed.py` | RG multi-seed robustness. |
| v12 | `naturkonstante_v12_maxent_lambda.py` | Maximum-entropy calibration of sector weights `lambda_a`. |
| v13 | `naturkonstante_v13_maxent_sm_bridge.py` | MaxEnt-weighted constants bridge using the calibrated sector weights. |

## MaxEnt Sector Weights

The v12 script calibrates effective sector weights as positive MaxEnt
Lagrange multipliers over a synthetic defect ensemble. The resulting weights
are:

| Sector | lambda_a | Normalized share |
| --- | ---: | ---: |
| H | 1.610873 | 0.0728 |
| B | 2.428297 | 0.1097 |
| S | 4.623655 | 0.2088 |
| V | 4.644912 | 0.2098 |
| R | 8.832344 | 0.3989 |

The v13 script uses these weights in a normalized weighted MAAT score. The
MaxEnt-calibrated weights are effective benchmark weights, not fundamental
microscopic constants.

## Data Sources

Reference comparison values for fundamental constants should be cited from
official sources. In the accompanying paper, low-energy constants follow
CODATA/NIST recommended values, and particle-physics comparison values follow
PDG-style 2024 inputs:

- NIST/CODATA values of the fundamental physical constants:
  <https://www.nist.gov/programs-projects/codata-values-fundamental-physical-constants>
- CODATA recommended values of the fundamental physical constants 2022:
  <https://www.nist.gov/publications/codata-recommended-values-fundamental-physical-constants-2022>
- Particle Data Group, Review of Particle Physics 2024:
  <https://pdg.lbl.gov/index-2024.html>
- PDG 2024 journal citation:
  S. Navas et al. (Particle Data Group), Phys. Rev. D 110, 030001 (2024),
  <https://doi.org/10.1103/PhysRevD.110.030001>

## Data Attribution and License Note

PDG 2024 content is licensed under Creative Commons Attribution 4.0
International (CC BY 4.0), except where otherwise noted. NIST web information
is generally public information unless marked as copyrighted, while NIST
Standard Reference Data may carry separate copyright and licensing conditions.

This benchmark uses only published reference comparison values and does not
redistribute NIST or PDG databases. When reusing or discussing the numerical
inputs, cite CODATA/NIST and PDG as the source of the reference values.

## Scientific Status

The benchmark supports only a cautious claim: structural scoring can define
low-defect basins near observed effective constants, and MaxEnt calibration can
replace the earlier equal-weight assumption by effective sector weights. It
does not show that the constants have been uniquely derived. The strongest next
step is a precision Standard-Model bridge with full RG running, threshold
corrections, uncertainty propagation, real defect ensembles, and stricter
held-out prediction tests.
