# Paper 58 Experiment — Stationarity Balance Beats Population Balance

**Stationarity Balance Beats Population Balance**  
*A Defect-Field Benchmark on Open Quantum Pointer States*

This experiment is fully classical. It does **not** require a quantum computer,
Qiskit runtime, IBM Quantum access, or paid hardware time.

The benchmark simulates one-qubit Lindblad dynamics under six open-system
channel families:

- `z_dephasing`
- `x_dephasing`
- `amplitude_damping`
- `thermal_relaxation`
- `mixed_zx`
- `depolarizing`

For each initial pure state and generator, the script computes a pointer
robustness target from the evolved trajectory:

```text
pointer_robustness =
0.50 fidelity_retention
+ 0.35 probability_stability
+ 0.15 purity_final
```

All predictive features are computed from the initial state and the Lindblad
generator, not from the final target.

## Main Question

Does a stationarity-sensitive balance support

```text
B_stat = 1 / (1 + 2 |dp/dt|_0 + 0.35 |d Tr(rho^2)/dt|_0)
```

predict pointer robustness better than raw population balance

```text
B_pop = 1 - |p0 - p1|?
```

## Run

From this directory:

```bash
python3 open_quantum_pointer_paper58.py
```

Outputs are written to:

```text
outputs_paper58/
```

## Main Outputs

- `paper58_pointer_instances.csv`
- `paper58_model_comparison.csv`
- `paper58_leave_channel_out.csv`
- `paper58_scalar_correlations.csv`
- `paper58_feature_importance.csv`
- `paper58_predictions.csv`
- `paper58_summary.json`
- `fig1_model_comparison.png`
- `fig2_leave_family_out.png`
- `fig3_feature_importance.png`
- `fig4_scalar_correlations.png`
- `fig5_prediction_scatter.png`
- `fig6_bpop_vs_bstat.png`

## Core Results

| Feature set | 5-fold RF R2 | Spearman |
|-------------|-------------:|---------:|
| Standard quantum baseline | `0.9011` | `0.9241` |
| Defect-field stationary | `0.8527` | `0.8856` |
| Scalar stationarity | `0.6457` | `0.7719` |
| Scalar population balance | `0.4561` | `0.6221` |
| Rates only | `0.1641` | `0.4646` |
| State geometry only | `0.0008` | `0.2945` |
| Shuffled defect null | `-0.0261` | `-0.0003` |

Leave-channel-family-out:

| Feature set | LFO R2 | LFO Spearman |
|-------------|------:|-------------:|
| Defect-field stationary | `0.7075` | `0.8286` |
| Standard quantum baseline | `0.6637` | `0.8105` |
| Scalar stationarity | `0.4928` | `0.7284` |
| Scalar population balance | `0.3108` | `0.5604` |

Key result:

```text
Delta R2(scalar stationarity - scalar population) = +0.1896
```

The stationarity-sensitive balance definition is therefore much better aligned
with pointer robustness than raw population balance in this benchmark.

## Scientific Status

This is a reproducible classical open-system simulation. It is not a proof of
quantum measurement, not a collapse theory, and not a hardware experiment. The
result supports a narrower claim: pointer robustness is better captured by
stationarity-sensitive defect fields than by raw population balance alone.

## Data and License Note

No external quantum dataset is redistributed. All CSV/JSON/PNG files are
derived reproducibility artifacts generated locally by the script.
