# Paper 67 -- Structural Early Warning as a Refinement Trigger

This experiment tests whether a MAAT structural warning functional can act as
an adaptive-refinement trigger in forced 2D turbulence.

## Status

This is a local forced-2D utility pilot. It is not a Navier--Stokes regularity
claim and it does not use JHTDB data. JHTDB is included only as a predeclared
external replication route.

## Run

```bash
python3 structural_early_warning_paper67.py
```

The script writes all generated artifacts to `outputs_paper67/`.

## Primary Comparison

All monitors are thresholded at the same false-alarm target. The primary
utility score is:

```text
lead_coverage_utility = detection_rate * median_lead_time
```

The current local benchmark result is:

```text
W_MAAT utility          = 1.275
high-k monitor utility  = 1.225
max-vorticity utility   = 0.825
RMS-vorticity utility   = 0.5125
```

The observed `W_MAAT - high-k` utility margin is `0.050`, but the event-bootstrap
95% interval is `[-0.700, 0.800]`. This is therefore not yet a robust advantage
over direct spectral-tail monitoring.

The component ablation is also informative: `W_MAAT` without the robustness-loss
factor reaches utility `1.400`, while robustness-loss alone fails. In this pilot,
the useful trigger signal is tail/activity dominated.

## Data Attribution

All numerical results in this folder are generated locally. No external raw
data are redistributed. Future JHTDB validation should cite the Johns Hopkins
Turbulence Databases and the relevant JHTDB publications, including Li et al.
(2008), arXiv:0804.1703. No endorsement by JHTDB, Johns Hopkins University,
IDIES, SciServer, NSF, or the JHTDB dataset authors is implied.
