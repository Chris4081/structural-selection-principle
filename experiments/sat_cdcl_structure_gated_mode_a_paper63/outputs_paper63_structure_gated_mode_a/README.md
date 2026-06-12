# Paper 63: Structure-Gated Mode-A CDCL Branching

This experiment extends Paper 62 by gating the MAAT v1.6-T Mode-A structural
penalty with a root-level instance diagnostic.

Global Mode A is:

```text
A(c) = prog(c) - tau * h_MAAT(c)
```

Structure-gated Mode A is:

```text
A_g(c) = prog(c) - g(F) * tau * h_MAAT(c)
```

Predeclared parameters:

- seed: `63063`
- conflict budget: `3500`
- beta: `8.0`
- tau: `0.38`
- rollout horizon h: `0`
- policies: `vsids, progress_only, maat_only, mode_a, mode_a_gated, random`

Important scope note: this is not an industrial SAT solver and not a comparison
against CaDiCaL, Kissat, Glucose, or MiniSat. It is a controlled internal policy
test in which VSIDS, global Mode A, and gated Mode A share the same propagation
and conflict-learning machinery.

Headline gated Mode-A result vs VSIDS: median regret 0, mean cost ratio 1.163, wins/losses/ties 47/48/4.

Gated Mode A vs global Mode A: median regret -0.008, mean regret -0.5703.

Median gate by family: `{'graph_coloring': 0.5219, 'modular_3sat': 0.5903, 'pigeonhole': 0.0, 'planted_3sat': 0.5629, 'random_3sat': 0.581, 'xor3_cnf': 0.7488}`.

Outputs:

- `paper63_cdcl_runs.csv`: all per-instance runs
- `paper63_structure_gates.csv`: root-level structure-gate diagnostics
- `paper63_summary_by_policy.csv`: aggregate policy statistics
- `paper63_regret_vs_vsids.csv`: paired regret against VSIDS
- `paper63_regret_vs_mode_a.csv`: paired regret against global Mode A
- `paper63_family_regret_vs_vsids.csv`: family-level paired regret
- `paper63_summary.json`: metadata and compact result table
- `fig1_policy_compute_cost.png`
- `fig2_regret_vs_vsids.png`
- `fig3_family_regret.png`
- `fig4_structure_gate_by_family.png`
- `fig5_gated_vs_global_mode_a.png`
- `fig6_gate_vs_regret.png`

Interpretation rule:

- Negative regret means the policy used less compute than VSIDS on the paired
  instance.
- Positive regret means VSIDS was better.
- A useful gate should reduce global Mode-A tail regret without merely
  collapsing to `progress_only` everywhere.
