# Paper 70 -- Local Structural Channels in SAT

This folder contains the diagnostic follow-up to Paper 69.

Paper 69 executed the preregistered MAAT v1.7 Gate Challenge protocol on a
family-balanced SATLIB smoke subset. The frozen gate did not receive positive
support and lost clearly against the classical MOMS heuristic.

Paper 70 does **not** retune the gate. It asks a diagnostic question:

```text
Which local SAT information channels does MOMS exploit that the frozen
root-level gate does not see?
```

## Status

Post-hoc diagnostic analysis only. The outputs are derived from the Paper 69
SATLIB smoke execution and local CNF parsing. Raw CNF files are used only as
local inputs and are not emitted or redistributed.

## Run

From this folder:

```bash
python3 analyze_sat_local_channels_paper70.py \
  --paper69-dir ../gate_challenge_sat_paper69
```

The script expects the Paper 69 directory to contain:

- `outputs_paper69_sat_gate_challenge/paper69_moms_vs_gate_instance_diagnostics.csv`
- `outputs_paper69_sat_gate_challenge/paper69_gate_features.csv`
- `outputs_paper69_sat_gate_challenge/paper69_dataset_manifest_detected.csv`
- local raw CNFs under `data_external/`

The raw CNFs must be staged and verified according to the Paper 69 SHA256
manifest. They are not committed to the repository.

## Outputs

Outputs are written to:

```text
outputs_paper70_local_sat_channels/
```

Main outputs:

- `paper70_local_channel_features.csv`
- `paper70_channel_correlations.csv`
- `paper70_family_channel_summary.csv`
- `paper70_summary.json`
- `fig1_channel_correlations.png`
- `fig2_moms_signal_vs_regret.png`
- `fig3_family_channel_heatmap.png`
- `fig4_propagation_vs_regret.png`

## Interpretation Rule

This analysis may identify missing local channels. It must not be used to
change Paper 69 gate parameters retroactively.

Any future SAT-specific gate or local-channel extension must be preregistered
as a new test.
