# Active Respect / Structural Significance Toy Simulation

This experiment tests a possible extension of the MAAT v1.2.1 interpretation
of `R`.

The working hypothesis is:

```text
H, B, S, V = primitive support sectors
R          = emergent meta-field / order field
```

The experiment does **not** return to a fifth primitive `R` sector. Instead it
separates three derived levels:

```text
R_resp = (H B V)^(1/3)
```

Passive structural adherence.

```text
R_rob = min(R_resp, (H B S_eff V)^(1/4))
```

Activity-sensitive robustness.

```text
R_sig = R_resp^(1-alpha) S_eff^alpha
```

Active structural significance.

Here `S_eff` is not raw activity. It is controlled activity in an optimal
window:

```text
S_eff(A) = exp[-0.5 ((A - A_star) / sigma_A)^2]
```

The motivating example is star-like: a system is not structurally significant
because it is merely stable, but because coherent order, balance, controlled
activity, and environmental coupling jointly produce persistent influence.

## Run

```bash
cd experiments/active_respect_significance
python3 active_respect_significance.py
```

## Outputs

The script writes into `outputs/`:

- `active_respect_significance_results.json`
- `active_respect_curves.csv`
- `active_respect_ensemble.csv`
- `fig1_respect_hierarchy_vs_activity.png`
- `fig2_rsig_phase_map.png`
- `fig3_passive_vs_active_significance.png`
- `fig4_archetype_supports.png`

## Scientific Status

This is a toy simulation and conceptual sanity check. It is not a physical
stellar model, not a cosmological fit, and not a new primitive-field proposal.

Its purpose is narrower: to test whether the v1.2.1 idea can be extended from
passive robustness to active significance while preserving the four-primary
sector architecture.
