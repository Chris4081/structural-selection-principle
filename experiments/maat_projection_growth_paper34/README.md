# MAAT Projection Observable vs Growth and Expansion Data (Paper 34)

This folder contains the reproduction script and generated artifacts for:

**Paper 34 -- MAAT Projection Observable vs Growth and Expansion Data**  
*Transition Marker and Parameter Sensitivity in v0.11*

## Purpose

Paper 34 tests whether the MAAT projection observable `C_proj` collapses onto
standard cosmological growth data. It does not. This is the central diagnostic
result: `C_proj` is not a matter-growth observable but a projection-level
diagnostic measuring structural projection stress.

## Main Command

Run from this folder:

```bash
python3 maat_paper34_reproduction_v2.py
```

The script creates the `paper34_outputs/` folder and writes all CSV, JSON, and
PNG artifacts used in the paper.

## Key Results

| Quantity | Result |
|---|---:|
| Baseline transition marker | `z_tr = 1.04356` |
| `C_norm(z_tr)` | `22.64463` |
| `R_proj(z_tr)` | `1.11151e-03` |
| `C_norm(z=3)` | `647.49295` |
| LCDM `f sigma_8` chi-square/dof | `1.04813` |
| Rescaled `C_proj` vs `f sigma_8` chi-square/dof | `42.76835` |
| LCDM `H(z)` chi-square/dof | `0.47944` |
| Gamma-alpha scan points | `625` |
| Scan median transition | `1.18715` |
| Scan transition range | `[0.28073, 1.58353]` |

## Outputs

| File | Role |
|---|---|
| `paper34_summary.json` | Machine-readable summary of all key numerical results. |
| `paper34_core_curves.csv` | Redshift curves for expansion, growth, projection, residual, and curvature. |
| `paper34_fsigma8_comparison.csv` | Growth-data comparison table. |
| `paper34_scan_gamma_alpha.csv` | 25x25 scan over gamma and alpha. |
| `fig_paper34_main.png` | Main comparison and transition-marker figure. |
| `fig_paper34_components.png` | Projection component curves. |
| `fig_paper34_scan.png` | Transition-redshift heatmap over gamma and alpha. |

## Data and Attribution Note

The script includes small fixed tables of Cosmic Chronometer and `f sigma_8`
values used only for reproducible diagnostic comparisons. These measurements
come from external observational literature and should be cited to the
original surveys and compilations when reused or discussed. The repository
contains derived analysis artifacts only; no endorsement by the original
collaborations is implied.

## Scientific Status

This is a diagnostic preprint-level benchmark. It does not fit cosmological
parameters, replace LCDM, or claim evidence for modified gravity. Its purpose
is to show that the projection observable probes a different informational
layer than standard growth observables.
