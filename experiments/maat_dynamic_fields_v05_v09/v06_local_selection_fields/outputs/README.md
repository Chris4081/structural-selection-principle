# MAAT v0.6 Local Selection Fields

This folder contains a minimal 1D phi^4 benchmark for local MAAT
selection-pressure fields.

The experiment promotes global weights to local fields:

```text
lambda_a -> lambda_a(x,t)
```

and evolves them with:

```text
tau_lambda partial_t lambda_a =
    D_lambda partial_xx lambda_a - lambda_a + lambda_a_star(x,t)
```

The local fixed point `lambda_star(x,t)` is computed from coarse-grained local
defect covariance.

## Outputs

| File | Meaning |
|------|---------|
| `local_selection_phi4_timeseries.csv` | Time-series metrics |
| `local_selection_phi4_summary.json` | Parameters and final diagnostics |
| `field_profiles.png` | Initial, clean, baseline final, and MAAT final fields |
| `metric_comparison.png` | Baseline vs MAAT metrics |
| `lambda_field_snapshots.png` | Local lambda field snapshots |
| `lambda_field_heatmap.png` | Compact lambda heatmap |
| `final_defect_profiles.png` | Final defects vs clean-kink target |

## Final Ratios

| Metric | MAAT / baseline |
|--------|----------------:|
| distance to kink | 1.0060 |
| residual RMS | 0.1648 |
| roughness | 1.0000 |
| energy | 1.0000 |

## Post-Perturbation Mean Ratios

| Metric | MAAT / baseline |
|--------|----------------:|
| distance to kink | 0.9682 |
| residual RMS | 0.7576 |
| roughness | 0.9995 |
| energy | 0.9975 |

This is an effective local-field benchmark, not a first-principles microscopic
derivation.
