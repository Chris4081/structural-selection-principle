# Dark Matter as Structural Coherence Residual in MAAT v1.6

This folder contains the reproducibility script for the extra phenomenological
paper:

**Dark Matter as Structural Coherence Residual in MAAT v1.6**  
*A Phenomenological Interpretation of Galactic Rotation Residuals*

## Purpose

The benchmark is a toy diagnostic, not a galaxy fit and not a replacement for
particle dark matter.  It asks how MAAT v1.6 reads the standard rotation-curve
residual:

```text
v_dark^2(r) = max(v_obs^2(r) - v_baryon^2(r), 0)
```

The interpretation is deliberately conservative:

> Dark matter is treated as the gravitationally inferred structural residual
> required to keep the observed dynamical trajectory coherent after visible
> contributions are accounted for.

## Run

```bash
python3 dark_matter_structural_coherence.py
```

The script writes all outputs to:

```text
outputs_phenomenological/
```

## Outputs

| File | Role |
|------|------|
| `dark_matter_structural_coherence_profile.csv` | Radial toy profile, inferred dark residual, and MAAT supports |
| `dark_matter_structural_coherence_summary.json` | Summary statistics used in the paper |
| `fig1_rotation_curve_residual.png` | Toy rotation curve, baryonic contribution, and dark residual |
| `fig2_structural_supports.png` | MAAT v1.6 supports and structural dark support |
| `fig3_residual_robustness_phase.png` | Residual fraction versus robustness phase portrait |

## Main summary values

| Quantity | Value |
|----------|------:|
| max observed velocity | `204.9972 km/s` |
| outer mean residual fraction | `0.8802` |
| outer mean robustness | `0.7806` |
| transition peak structural dark support | `0.3486` |
| radius of peak structural dark support | `20.3782 kpc` |
| mean H | `1.0000` |
| mean B | `0.9060` |
| mean S_eff | `0.7506` |
| mean V | `0.5919` |
| mean R_rob | `0.7611` |

## Data attribution and license note

No external astronomical dataset is redistributed.  The rotation curve,
baryonic profile, residual profile, CSV tables, JSON summary, and figures are
synthetic reproducibility artifacts generated locally by the script.

