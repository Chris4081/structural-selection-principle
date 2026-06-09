# MAAT v1.6-T T1: Mode-A CDCL Branching Pilot

This experiment is the first policy-level benchmark for the MAAT v1.6 search
geometry layer. It inserts the Mode-A acquisition rule into a small transparent
CDCL solver and measures compute-to-solution regret against VSIDS.

The tested acquisition rule is:

```text
A(c) = prog(c) - tau * h_MAAT(c)
```

where `prog(c)` is a normalized CDCL-local progress score and `h_MAAT(c)` is a
bounded structural-risk score computed from local SAT support coordinates.

## Scope

This is not an industrial SAT solver and not a benchmark against CaDiCaL,
Kissat, Glucose, or MiniSat. It is a controlled internal test in which VSIDS,
progress-only, MAAT-only, random, and Mode-A policies share the same propagation
and conflict-learning machinery.

The benchmark is deliberately small so that policy differences are auditable.

## Predeclared Parameters

- seed: `62062`
- conflict budget: `3500`
- beta: `8.0`
- tau: `0.38`
- rollout horizon h: `0`
- policies: `vsids`, `progress_only`, `maat_only`, `mode_a`, `random`

## Run

```bash
python3 sat_cdcl_mode_a_t1.py
```

## Outputs

The script writes all artifacts to:

```text
outputs_t1_cdcl_mode_a/
```

Main files:

- `t1_cdcl_runs.csv`
- `t1_cdcl_summary_by_policy.csv`
- `t1_cdcl_regret_vs_vsids.csv`
- `t1_cdcl_family_regret_vs_vsids.csv`
- `t1_cdcl_summary.json`
- `fig1_policy_compute_cost.png`
- `fig2_regret_vs_vsids.png`
- `fig3_family_regret.png`

## Initial Result

In the first fixed-seed run, all 99 instances are solved. Mode A ties VSIDS in
median cost but has positive mean regret and does not beat `progress_only`
overall. The family-level result is mixed: Mode A is helpful on modular SAT and
some structured families, but loses clearly on pigeonhole instances.

Interpretation:

- Mode A is a usable active branching signal in this transparent solver.
- The structural penalty does not yet earn compute globally.
- The result supports a family-/mode-sensitive view of v1.6 search geometry.
- The next test should add stronger CDCL baselines or make `tau` family-adaptive
  under a predeclared tuning protocol.
