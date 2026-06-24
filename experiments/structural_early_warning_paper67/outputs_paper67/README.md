# Paper 67 -- Structural Early Warning as a Refinement Trigger

Utility benchmark for adaptive-refinement triggers in forced 2D turbulence.

## Headline

- Best lead-coverage utility monitor: `W_MAAT`.
- Utility: `1.2750` at false-alarm target `0.050`.
- Event count: `4`.
- Observed utility margin over high-k: `0.0500`.
- Event-bootstrap 95% CI for that margin: `[-0.7000, 0.8000]`.

## Important status

This first implementation uses a local forced-2D turbulence simulation.
JHTDB is specified as an external replication route but no JHTDB data are downloaded or redistributed.

## Outputs

- `paper67_forced2d_timeseries.csv`
- `paper67_trigger_results.csv`
- `paper67_ablation_results.csv`
- `paper67_bootstrap_ci.csv`
- `paper67_event_leads.csv`
- `paper67_summary.json`
- `paper67_jhtdb_protocol.json`
- `fig1_forced2d_warning_timeseries.png`
- `fig2_trigger_leadtime_comparison.png`
- `fig3_false_alarm_calibration.png`
- `fig4_ablation_utility.png`
- `fig5_bootstrap_delta_ci.png`
