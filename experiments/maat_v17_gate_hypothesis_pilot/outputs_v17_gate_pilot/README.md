# MAAT v1.7 Gate-Hypothesis Pilot

This experiment re-analyses existing repository outputs to test a narrow
hypothesis:

> R often works better as a conditional gate than as another global scoring
> factor.

Domains:

- SAT: Paper 63 global Mode A vs structure-gated Mode A.
- Quantum: Paper 64 scalar R score vs hard R gate for entanglement robustness.
- Fluid: Paper 67 W_MAAT score vs hard R gate for refinement triggering.

This is a pilot, not a new final MAAT version.  It uses no external dataset
beyond artifacts already present in this repository.

## Headline

The evidence supports the gate direction, but not yet a completed v2.0-style
framework.  SAT shows a small tail-damage reduction from gating.  The quantum
scalar test shows a large in-distribution ranking improvement.  The fluid test
shows that the best hard R gate improves lead-coverage utility, but does so by
sacrificing event coverage and is sensitive to the gate threshold.

Conclusion: v1.7 is justified as a gate-hypothesis / design-direction paper.
v2.0 would require a predeclared, calibrated gate that generalizes under harder
external tests.

## Summary Table

```text
 domain                                metric               score_variant           gate_variant  score_value  gate_value  delta_gate_minus_score  relative_improvement median_delta gate_win_rate    n                                                                                 interpretation
    SAT          compute_cost_lower_is_better               global Mode A structure-gated Mode A    20.556808   19.986465               -0.570343              0.027745       -0.008      0.505051   99                   gate slightly reduces mean/tail damage relative to ungated structural policy
Quantum Spearman with entanglement robustness           base * R_rob_stat          gate_R_ge_q75     0.251318    0.514968                0.263650              1.049072           --            -- 1800 hard R gate improves scalar ranking in-distribution, but this is not a field-model replacement
  Fluid                 lead_coverage_utility W_MAAT = base * (1 - R_rob)          gate_R_le_q75     1.275000    1.500000                0.225000              0.176471           --            --    4     best hard R gate improves utility but sacrifices event coverage and is threshold-sensitive
```
