# Societal Critical Coherence Index

Synthetic toy benchmark for structural stress, controlled activity, and active
significance in social-system archetypes.

## Status

This is an ethically bounded toy model. It does **not** rank real people,
parties, countries, organisations, or communities. All inputs are synthetic
archetypes chosen to test whether the MAAT/CCI logic can separate raw
mobilisation from constructive transformation.

## Core Idea

Social health is not defined as the absence of conflict. In this toy
framework, a social system is structurally healthier when it can transform
conflict into coherent, balanced, and connected activity.

The primitive support sectors are:

| Symbol | Meaning in the toy model |
| --- | --- |
| `H` | Coherent orientation, shared sense-making, institutional consistency |
| `B` | Balance between tensions and conflict integration |
| `A_raw` | Raw mobilisation, reform pressure, social activity |
| `V` | Connectedness, trust, communication, social cohesion |

Raw activity is converted into controlled activity:

```text
S_eff(A) = exp[-0.5 ((A_raw - A_star) / sigma_A)^2]
```

This encodes the central assumption: too little activity produces stagnation,
while too much activity can become overheated or chaotic. Constructive
transformation occurs in an activity window around `A_star`.

## Derived Quantities

Passive structural adherence:

```text
R_resp,soc = (H B V)^(1/3)
```

Activity-sensitive robustness:

```text
R_rob,soc = min(R_resp,soc, (H B S_eff V)^(1/4))
```

Active societal significance:

```text
R_sig,soc = R_resp,soc^(1-alpha) S_eff^alpha
```

Societal structural-stress index:

```text
CCI_soc = A_raw / (H + B + V + R_rob,soc + epsilon)
```

Constructive-significance index:

```text
ASI_soc = R_sig,soc / (1 + CCI_soc)
```

## Synthetic Archetypes

The benchmark evaluates six synthetic archetypes:

- `stagnant_order`
- `constructive_reform`
- `polarized_mobilization`
- `fragmented_activism`
- `authoritarian_stability`
- `creative_democratic_renewal`

It also generates a fixed-seed synthetic ensemble of 5000 random social-system
states for phase-space visualisation.

## Reproduce

Run from this folder:

```bash
python3 societal_cci_toy.py
```

The script writes JSON, CSV, and PNG outputs to `outputs/`.

## Outputs

| File | Role |
| --- | --- |
| `outputs/societal_cci_results.json` | Summary metrics, definitions, key results |
| `outputs/societal_cci_archetypes.csv` | Synthetic archetype table |
| `outputs/societal_cci_ensemble.csv` | Fixed-seed synthetic ensemble |
| `outputs/fig1_societal_indices_by_archetype.png` | Index comparison across archetypes |
| `outputs/fig2_stress_vs_significance_phase_space.png` | Stress vs constructive-significance phase space |
| `outputs/fig3_controlled_activity_window.png` | Controlled activity window |
| `outputs/fig4_synthetic_regime_counts.png` | Synthetic regime counts |

## Key Result

The toy benchmark separates three concepts that are often conflated:

- Low stress can still be stagnant.
- High activity can still be polarising or fragmenting.
- Constructive transformation requires controlled activity together with
  coherence, balance, and connectedness.

## Ethical Boundary

This model is a reflection and diagnostics toy framework. It is not an
instrument for political scoring, population ranking, surveillance, or
normative judgement of real communities.
