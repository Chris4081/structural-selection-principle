# Paper 46 — Environmental Screening in MAAT Structural Cosmology

This experiment extends the Paper-45 metric-response benchmark by making the
bounded MAAT projection kernel environment-dependent.

## Status

This is an effective screening benchmark. It is not a microscopic derivation
of chameleon, Vainshtein, symmetron, or other screening mechanisms, and it is
not a local-gravity test. The goal is to show how the MAAT projection response
can be suppressed in dense or structurally stabilized environments.

## Core Definitions

Environmental projection kernel:

```text
C_env(z, Delta, Sigma_env) = C_hat_proj(z) * S_env(Delta, Sigma_env)
```

Screening factor:

```text
S_env = [1 + alpha_rho Delta_+^n + alpha_sigma Sigma_env^m]^(-1)
Delta_+ = max(Delta - 1, 0)
```

Metric-response channels:

```text
mu_env(z)       = 1 + eta_g  C_env(z)
eta_slip_env(z)= 1 + beta_s C_env(z)
Sigma_lens_env = mu_env(z) * [1 + eta_slip_env(z)] / 2
```

## Parameters

| Parameter | Value |
| --- | --- |
| `eta_g` | `0.02` |
| `beta_s` | `0.03` |
| `alpha_rho` | `0.15` |
| `alpha_sigma` | `1.0` |
| `n_rho` | `0.75` |
| `m_sigma` | `2.0` |

## Environmental Archetypes

| Environment | `Delta` | `Sigma_env` | Interpretation |
| --- | ---: | ---: | --- |
| `void` | `0.20` | `0.05` | low-density, weakly screened |
| `sheet` | `0.80` | `0.15` | mildly underdense sheet |
| `field` | `1.00` | `0.20` | mean-density field |
| `filament` | `5.00` | `0.45` | moderately overdense structure |
| `cluster` | `100.00` | `0.85` | dense cluster-like environment |
| `local_dense` | `1e6` | `0.95` | local high-density proxy |

## Key Results

| Quantity | Result |
| --- | --- |
| Stable environments | `6 / 6` |
| `S_env(void)` | `0.9975` |
| `S_env(cluster)` | `0.1555` |
| `S_env(local_dense)` | `0.0002107` |
| Void max `|Sigma_lens - 1|` | `3.5151%` |
| Cluster max `|Sigma_lens - 1|` | `0.5441%` |
| Local-dense max `|Sigma_lens - 1|` | `0.000736%` |
| Screening transition, `Sigma_env=0.2` | `S_env ~= 0.5` at `Delta ~= 12.35` |

## Interpretation

The same MAAT projection response that can be percent-level in void-like
regions becomes strongly suppressed in cluster-like regions and practically
invisible in a high-density local-test proxy. This provides an effective
structural answer to the standard modified-gravity question: why a
cosmological metric response might coexist with near-GR local behaviour.

## Reproduce

Run from this folder:

```bash
python3 maat_environmental_screening_solver_v01.py
```

## Outputs

| File | Role |
| --- | --- |
| `outputs_paper46/paper46_summary.json` | Main summary of definitions, parameters, and key results |
| `outputs_paper46/paper46_environment_archetypes.csv` | Environment table and response metrics |
| `outputs_paper46/paper46_environment_curves.csv` | Redshift curves for all environments |
| `outputs_paper46/paper46_screening_grid.csv` | Environmental phase-space grid |
| `outputs_paper46/fig1_environmental_screening_summary.png` | Main four-panel summary |
| `outputs_paper46/fig2_screening_phase_space.png` | Screening phase-space heatmaps |
| `outputs_paper46/fig3_environment_response_bars.png` | Environment response hierarchy |
| `outputs_paper46/fig4_screening_transition.png` | Transition with overdensity |

## Data Attribution and License Note

The Planck-normalised reference parameters and compact `f sigma_8` comparison
points are external scientific data and should be cited to the original
publications/collaborations. Repository CSV/PNG files are derived
reproducibility artifacts only. No endorsement by the Planck Collaboration,
survey collaborations, or original data authors is implied.
