# MAAT v0.9 FLRW Stability Scan

This folder contains a toy FLRW stability scan for the stable sigma=+1 branch
of the v0.8 scalar worked example.

The scanned model is:

```text
P(X, lambda) = X + mu lambda U(X)
d_X = X / Lambda_X^4
U(X) = -log((epsilon + Gamma)/(epsilon + 1))
Gamma = 1/(1+d_X)
```

The scan checks:

- no ghost: `P_X > 0`
- gradient stability: `c_s^2 > 0`
- positive MAAT energy density
- background safety: `max Omega_MAAT < 0.05`

## Outputs

| File | Meaning |
|------|---------|
| `flrw_stability_scan_results.csv` | Full parameter scan |
| `flrw_stability_scan_summary.json` | Aggregate diagnostics |
| `representative_flrw_trajectory.csv` | Representative stable trajectory |
| `flrw_stability_scan_heatmaps.png` | Scan heatmaps |
| `representative_flrw_trajectory.png` | Stable trajectory observables |
| `lambda_and_kinetic_evolution.png` | Lambda and kinetic invariant |
| `kinetic_branch_phase_space.png` | Analytic stability map in `(X, mu lambda)` |

## Summary

| Quantity | Value |
|----------|------:|
| total cases | 784 |
| stable cases | 666 |
| background-safe cases | 618 |
| stable fraction | 0.8495 |
| background-safe fraction | 0.7883 |

This is not a cosmological data fit. It is a stability and consistency scan
for the kinetic structural-selection branch.
