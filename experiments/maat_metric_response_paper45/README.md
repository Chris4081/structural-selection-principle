# Paper 45 — Effective Metric Response in MAAT Structural Cosmology

This experiment extends the Paper-43 linear-growth benchmark from a
growth-only effective coupling to a minimal metric-response parameterization.

## Status

This is an effective benchmark, not a full relativistic perturbation solver.
It does not implement CLASS/CAMB, weak-lensing likelihoods, CMB anisotropies,
or nonlinear structure formation. The goal is to test whether the existing
MAAT projection kernel can be inserted into the standard metric-response
language in a bounded and reproducible way.

## Definitions

Growth coupling:

```text
mu(z) = G_eff/G = 1 + eta_g * C_hat_proj(z)
```

Gravitational slip:

```text
eta_slip(z) = Phi/Psi = 1 + beta_s * C_hat_proj(z)
```

Weyl/lensing response:

```text
Sigma(z) = mu(z) * [1 + eta_slip(z)] / 2
```

Two diagnostic observables are reported:

```text
Weyl_proxy(z) = Sigma(z) * D_MAAT(z) / D_LCDM(z)
EG_proxy(z)   = Sigma(z) * f_LCDM(z) / f_MAAT(z)
```

## Inputs

The experiment reuses the Paper-43 setup:

- Planck-like flat LCDM reference parameters
- response-derived bounded projection kernel `C_hat_proj(z)`
- compact `f sigma_8` comparison set used in Papers 40, 42, and 43

No lensing data are fitted in this benchmark. Consequently, the compact
growth-only comparison constrains `eta_g` but does not constrain `beta_s`.

## Reproduce

Run from this folder:

```bash
python3 maat_metric_response_solver_v01.py
```

## Key Results

| Quantity | Result |
| --- | --- |
| `eta_g` scan | `[0.00, 0.08]`, 41 points |
| `beta_s` scan | `[-0.06, 0.06]`, 61 points |
| Total branches | `2501` |
| Stable / positive branches | `2501 / 2501` |
| Growth-only best `eta_g` | `0.0000` |
| Growth-only `beta_s` | unconstrained without lensing data |
| Representative branch | `eta_g = 0.02`, `beta_s = 0.03` |
| Max `|mu - 1|`, representative | `1.9966%` |
| Max `|eta_slip - 1|`, representative | `2.9949%` |
| Max `|Sigma - 1|`, representative | `3.5240%` |
| Max Weyl-proxy deviation | `2.7314%` |
| Max `E_G`-proxy deviation | `2.3019%` |

## Outputs

| File | Role |
| --- | --- |
| `outputs_paper45/paper45_summary.json` | Summary of definitions, scan, and representative branch |
| `outputs_paper45/paper45_metric_curves.csv` | Redshift curves for growth and metric-response quantities |
| `outputs_paper45/paper45_metric_scan.csv` | Full `eta_g` / `beta_s` scan |
| `outputs_paper45/paper45_growth_metric_comparison.csv` | Compact growth-data comparison table with metric response |
| `outputs_paper45/fig1_metric_response_summary.png` | Main four-panel summary |
| `outputs_paper45/fig2_metric_response_scan.png` | Metric-response scan heatmaps |
| `outputs_paper45/fig3_growth_vs_metric_response.png` | Growth response vs metric response scatter |

## Data Attribution and License Note

The Planck-normalised reference parameters and compact `f sigma_8` comparison
points are external scientific data and should be cited to the original
publications/collaborations. Repository CSV/PNG files are derived
reproducibility artifacts only. No endorsement by the Planck Collaboration,
survey collaborations, or original data authors is implied.
