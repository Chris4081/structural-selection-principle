# Paper 66 -- Frozen Transfer Triangle

Internal transfer benchmark across SAT, Quantum, and Fluid domains.

The protocol freezes source-domain weights and evaluates them on the other two domains without target-domain retuning.

## Headline

- Best cross-domain direction: `Quantum->Fluid` using `equal` with Spearman `0.9909`.
- Row-duplicated `v_only` target-alignment mean: `0.4735`.
- Mean cross-domain Spearman by `nnls`: `0.3814`.
- `v_only` is a target-mode diagnostic, not source-trained transfer.

## Outputs

- `paper66_transfer_results.csv`
- `paper66_frozen_weights.csv`
- `paper66_shuffle_nulls.csv`
- `paper66_summary.json`
- `fig1_transfer_triangle_nnls.png`
- `fig2_vonly_target_alignment.png`
- `fig3_rule_comparison.png`
- `fig4_v_weights.png`
- `fig5_null_margins.png`

## Status

This is an internal benchmark only. It tests whether the connectivity mode `V` behaves like a shared low-frequency target-alignment channel across existing MAAT domains. The Fluid target is support-derived, so transfers into Fluid are partly circular. It is not external validation and not a universality proof.
