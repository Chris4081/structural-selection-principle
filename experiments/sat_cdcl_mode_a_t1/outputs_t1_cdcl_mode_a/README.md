# MAAT v1.6-T T1: Mode-A CDCL Branching Pilot

This experiment inserts the MAAT v1.6-T Mode-A acquisition rule into a small
transparent CDCL solver and measures compute-to-solution regret against VSIDS.

Mode A is:

```text
A(c) = prog(c) - tau * h_MAAT(c)
```

Predeclared parameters:

- seed: `62062`
- conflict budget: `3500`
- beta: `8.0`
- tau: `0.38`
- rollout horizon h: `0`
- policies: `vsids, progress_only, maat_only, mode_a, random`

Important scope note: this is not an industrial SAT solver and not a comparison
against CaDiCaL, Kissat, Glucose, or MiniSat. It is a controlled internal policy
test in which VSIDS and Mode A share the same propagation and conflict-learning
machinery.

Headline Mode-A result vs VSIDS: median regret 0, mean cost ratio 1.295, wins/losses/ties 44/47/8.

Outputs:

- `t1_cdcl_runs.csv`: all per-instance runs
- `t1_cdcl_summary_by_policy.csv`: aggregate policy statistics
- `t1_cdcl_regret_vs_vsids.csv`: paired regret against VSIDS
- `t1_cdcl_family_regret_vs_vsids.csv`: family-level paired regret
- `t1_cdcl_summary.json`: metadata and compact result table
- `fig1_policy_compute_cost.png`
- `fig2_regret_vs_vsids.png`
- `fig3_family_regret.png`

Interpretation rule:

- Negative regret means the policy used less compute than VSIDS on the paired
  instance.
- Positive regret means VSIDS was better.
- A useful Mode-A signal must beat `progress_only`; otherwise the structural
  penalty did not earn its compute.
