# Paper 71 -- Decision-Local Structural Selection in SAT

This folder contains the post-hoc hypothesis-generation analysis following
Papers 69--70.

Paper 71 does **not** modify the preregistered Paper-69 gate, does **not**
introduce a new gate, and does **not** define MAAT v1.8. It asks whether the
local SAT channels identified in Paper 70 are independent structural
information or mostly reparameterizations of existing supports.

## Status

Post-hoc hypothesis generation only.

Input:

- `../sat_local_structural_channels_paper70/outputs_paper70_local_sat_channels/paper70_local_channel_features.csv`

No raw CNF files are required for Paper 71.

## Run

```bash
python3 analyze_local_channel_independence_paper71.py
```

## Outputs

Outputs are written to:

```text
outputs_paper71_local_channel_independence/
```

Main outputs:

- `paper71_local_support_correlations.csv`
- `paper71_local_independence_summary.csv`
- `paper71_regret_model_comparison.csv`
- `paper71_residual_shuffle_nulls.csv`
- `paper71_summary.json`
- `fig1_residual_local_channels.png`
- `fig2_model_comparison.png`
- `fig3_short_entropy_residual_vs_regret.png`

## Interpretation Rule

Candidate local channels may motivate a future preregistered test, but Paper 71
does not change Paper 69 and does not validate a new gate.

The useful output is a hypothesis boundary:

```text
local occurrence/degree residuals are weak but testable candidates;
short-clause entropy is not independent after existing supports are removed.
```
