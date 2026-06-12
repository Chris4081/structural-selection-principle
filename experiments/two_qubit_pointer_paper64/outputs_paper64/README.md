# Paper 64: Two-Qubit Pointer Selection

This experiment extends the one-qubit pointer benchmarks of Papers 51 and 58
to two qubits. It asks which entangled states survive open-system Lindblad
evolution and whether correlated stationarity is a better balance coordinate
than raw local population balance.

Headline result:

- best feature set: `defect_field_stationary` with RF R2 `0.9160`
- scalar correlated stationarity RF R2: `0.7831`
- scalar population RF R2: `0.7232`
- stationarity gain over population: `+0.0599`

Outputs:

- `paper64_two_qubit_instances.csv`
- `paper64_model_comparison.csv`
- `paper64_scalar_correlations.csv`
- `paper64_feature_importance.csv`
- `paper64_summary.json`
- `fig1_model_comparison.png`
- `fig2_balance_scatter.png`
- `fig3_channel_family_summary.png`
- `fig4_feature_importance.png`
- `fig5_scalar_correlations.png`

No external quantum dataset is redistributed. All CSV/JSON/PNG files are
derived reproducibility artifacts generated locally by the script.
