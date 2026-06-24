# Paper 66 -- Frozen Transfer Triangle

This experiment extends the Paper 52/53 frozen-transfer protocol from two
domains to three:

- SAT frustration fields from Paper 50b,
- quantum pointer-state robustness from Paper 51,
- fluid warning/coherence diagnostics from Paper 60.

The central question is whether the connectivity sector `V` behaves like a
shared low-frequency target-alignment channel across all three internal
domains.

## Run

```bash
python3 frozen_transfer_triangle_paper66.py
```

## Protocol

For each source domain, the script learns source-only non-negative weights over
`H,B,S,V,R_rob`, freezes them, and evaluates them on the other two domains
without retuning or sign flips. It also evaluates equal weights, `V`-only
target alignment, and shuffled-defect nulls.

`R_rob` is used here as a derived diagnostic coordinate, not as a fifth
primitive support sector.

All targets are oriented as:

```text
higher score = higher structural quality
```

## Outputs

```text
outputs_paper66/paper66_transfer_results.csv
outputs_paper66/paper66_frozen_weights.csv
outputs_paper66/paper66_shuffle_nulls.csv
outputs_paper66/paper66_summary.json
outputs_paper66/fig1_transfer_triangle_nnls.png
outputs_paper66/fig2_vonly_target_alignment.png
outputs_paper66/fig3_rule_comparison.png
outputs_paper66/fig4_v_weights.png
outputs_paper66/fig5_null_margins.png
```

## Status

Internal benchmark only. This is not external validation and not a universality
proof. A failure is informative: it would mean the previously observed
connectivity transfer is not a stable three-domain channel.

Important caveats:

- `V`-only rows are source-independent and should be read as target alignment,
  not fair architecture-level transfer.
- The Fluid target is support-derived from Paper 60, so transfers into Fluid
  are partly circular. Transfer into SAT is the cleaner stress direction.
