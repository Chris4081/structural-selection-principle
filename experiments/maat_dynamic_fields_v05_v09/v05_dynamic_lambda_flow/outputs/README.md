# MAAT v0.5 Dynamic Lambda Flow

This folder contains a toy simulation of dynamic structural selection
pressures.

The simulated equation is:

```text
tau d lambda / dt = -lambda + lambda_star(t)
lambda_star(t) = (C(t) + eta tr(C(t))/A I)^(-1) (mean(d)(t) - d_target(t))
```

The interval `45 <= t < 95` introduces a boundary-pressure event by tightening
the target robustness defect `d_R*`.

## Main Outputs

| File | Meaning |
|------|---------|
| `lambda_dynamic_flow_timeseries.csv` | Full time series |
| `lambda_dynamic_flow_summary.json` | Parameters and final diagnostics |
| `lambda_trajectory.png` | Evolution of all lambda sectors |
| `loss_and_tracking.png` | Loss and tracking error |
| `mean_vs_target_defects.png` | Selected mean defects vs target defects |
| `lambda_phase_portrait.png` | Phase portrait in lambda space |
| `stability_and_conditioning.png` | Stability and covariance conditioning |

## Final Diagnostics

| Quantity | Value |
|----------|------:|
| final loss | 0.00332624 |
| min loss | 0.00250204 |
| final lambda norm | 1.98115 |
| max lambda_R | 2.87546 |
| mean tracking error last 100 steps | 2.75463e-05 |

Interpretation: v0.4 response closure appears as the moving fixed point of a
v0.5 dynamical relaxation system. This is an effective toy model, not yet a
microscopic physical derivation.
