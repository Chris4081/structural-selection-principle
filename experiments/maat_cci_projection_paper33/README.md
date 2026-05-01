# MAAT CCI Projection Benchmark (Paper 33)

This folder contains the reproducibility scripts and generated artifacts for:

**Paper 33 -- The Critical Coherence Index as a Projection Observable**  
*Breadth--Depth Compression and Transition Robustness in a Minimal Cosmological Benchmark*

## Purpose

The scripts test a projection interpretation of the Critical Coherence Index (CCI).
Expansion defines a growing configuration breadth, linear growth and latent depth
define accessible coherent structure, and the CCI is interpreted as a stabilized
projection stress between the two.

This is a diagnostic benchmark, not a cosmological parameter fit.

## Main Scripts

Run from this folder:

```bash
python3 maat_cci_projection_test_v03.py
python3 maat_cci_projection_sensitivity_v04.py
```

The v0.2 script is included only as a stabilization predecessor:

```bash
python3 maat_cci_projection_test_v02.py
```

## Outputs

| Folder | Role |
|---|---|
| `maat_cci_projection_test_v02/` | Stabilized predecessor model with bounded compression and residual curvature transition. |
| `maat_cci_projection_test_v03/` | LCDM-growth proxy projection benchmark used as the main baseline. |
| `maat_cci_projection_sensitivity_v04/` | Parameter-sensitivity scan over projection strength, latent depth scale, sharpness, and depth floor. |

## Key Results

| Quantity | Result |
|---|---:|
| v0.3 transition estimate | `z_tr = 0.84985` |
| v0.3 final normalized CCI at z=3 | `474.184` |
| v0.4 parameter points | `6930` |
| v0.4 mean transition redshift | `0.84491` |
| v0.4 median transition redshift | `0.88889` |
| v0.4 fraction in 0.5 <= z_tr <= 1.1 | `0.55657` |

## Data and Attribution Note

The benchmark uses Planck-like LCDM parameter values and an analytic
Carroll--Press--Turner growth proxy. No observational likelihood or external
measurement table is fitted in this experiment.

The CSV/JSON/PNG files in this folder are derived reproducibility artifacts.
When discussing the background cosmology or growth approximation, cite the
original Planck and Carroll--Press--Turner sources. No endorsement by the
Planck Collaboration or by the original authors is implied.

## Scientific Status

This experiment supports an operational interpretation of the CCI as a
projection observable. It does not claim evidence for modified gravity,
dark-energy dynamics, or a measured cosmological transition.
