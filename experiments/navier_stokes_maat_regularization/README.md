# Extra Phenomenological Paper - BKM-Aware MAAT Diagnostics for Navier-Stokes

This experiment operationalises a symbolic MAAT-inspired functional for 3D
incompressible Navier-Stokes diagnostics:

```text
ToE_MAAT = integral [ (H+B+S+V+R) * Z * cascade_resistance ]
                    / [ DeltaE + DeltaQ + D0 ] dt
```

It is deliberately not a proof of Navier-Stokes regularity. It is a numerical
diagnostic programme that asks whether MAAT v1.6 structural coordinates can
track risk channels in a 3D Taylor-Green vortex benchmark.

## Run

```bash
python3 navier_stokes_maat_regularization.py
```

## Operational Interpretation

- `H`: equation residual plus incompressibility consistency
- `B`: energy balance consistency
- `S_eff`: bounded activity pressure from palinstrophy
- `V`: low-mode spectral coherence
- `R_rob`: v1.2.1 robustness closure
- `Z_partition`: structural coherence partition factor
- `cascade_resistance`: spectral-tail and vortex-stretching resistance
- `DeltaE`: energy-balance defect
- `DeltaQ`: spectral-tail/numerical quality defect
- `D0`: fixed dimensionless baseline penalty set to `11`

## Core Results

| Scenario | min R_rob | max warning | final ToE_MAAT action | warning vs stretching |
|---|---:|---:|---:|---:|
| moderate viscosity | `0.9029` | `0.1677` | `0.9567` | `0.9912` |
| low viscosity / high stress | `0.7578` | `0.5692` | `0.6786` | `1.0000` |

## Main Finding

The MAAT warning does not merely track the vorticity infinity norm. In this toy
3D setting it aligns much more strongly with vortex-stretching pressure, while
the integrated ToE_MAAT action is larger for the more coherent moderate-viscosity
trajectory.
This is the structural value added by the diagnostic: it points to a genuinely
three-dimensional stress channel rather than merely renaming a standard
vorticity monitor.
The observed correlations should not be interpreted as evidence for singularity
prediction. They indicate only that the structural diagnostic appears sensitive
to the dynamically relevant stretching channel in these controlled toy
trajectories.

The benchmark is still narrow: only two Taylor-Green trajectories are tested,
so the correlations should be read as a proof-of-concept signal rather than
broad evidence.

This supports a narrow interpretation:

> MAAT v1.6 can be used as a phenomenological structural-risk coordinate for
> Navier-Stokes trajectories. It does not solve the regularity problem.

## Outputs

- `outputs_phenomenological/navier_stokes_maat_diagnostics.csv`
- `outputs_phenomenological/summary_by_scenario.csv`
- `outputs_phenomenological/navier_stokes_maat_summary.json`
- `outputs_phenomenological/fig1_bkm_warning_timeseries.png`
- `outputs_phenomenological/fig2_toe_maat_action.png`
- `outputs_phenomenological/fig3_warning_vs_vorticity.png`

## Data and License Note

No external fluid dataset is redistributed. All velocity fields are generated
locally from Taylor-Green initial conditions. CSV/JSON/PNG files are derived
reproducibility artifacts.
