# Paper 59 - MAAT v1.6 Guided SAT Search

This experiment moves from passive SAT-hardness prediction to active SAT
search steering. It tests whether MAAT v1.6 structural-search coordinates can
serve as a branching heuristic in a transparent DPLL solver with unit
propagation.

It is deliberately modest:

- not a proof of `P != NP`
- not a CDCL solver competition
- not a replacement for Glucose, MiniSat, CaDiCaL, or Kissat
- not a claim that MAAT beats classical SAT heuristics

## Tested Heuristics

- `first`
- `random`
- `degree`
- `jeroslow_wang`
- `moms`
- `maat_v16_branch`
- `maat_v16_heavy`

## Benchmark Families

The script deterministically generates small standard SAT-family instances:

- random 3-SAT
- planted 3-SAT
- modular 3-SAT
- 3-XOR CNF
- graph coloring
- pigeonhole

## Run

```bash
python3 sat_v16_guided_search_paper59.py
```

## Core Aggregate Results

| Heuristic | Solved rate | Median decisions | Median conflicts | Median log cost |
|---|---:|---:|---:|---:|
| MOMS | 1.000 | 9.0 | 2.0 | 2.6283 |
| Jeroslow-Wang | 1.000 | 11.0 | 2.0 | 2.7788 |
| MAAT v1.6 branch | 1.000 | 12.0 | 2.0 | 2.9178 |
| MAAT v1.6 heavy | 1.000 | 12.0 | 2.0 | 2.9232 |
| Degree | 1.000 | 14.0 | 3.0 | 3.0587 |
| First | 0.993 | 17.0 | 6.0 | 3.2696 |
| Random | 1.000 | 18.0 | 5.0 | 3.3087 |

## Main Finding

The active MAAT v1.6 brancher is a usable structural search heuristic, but not
the strongest one in this benchmark. It improves over degree, first, and random
baselines, but it does not beat MOMS or Jeroslow-Wang.

This is the expected honest first result:

> MAAT v1.6 provides a usable structural branching signal, but in its present
> form it does not yet compete with the strongest SAT-specific heuristics.

The benchmark is deliberately light-to-moderate rather than adversarially hard.
Most instances are solved by almost every reasonable heuristic, so the
experiment mainly measures branching-cost ordering on solvable generated
instances, not hard timeout separation.
The improvement over degree, first, and random is therefore a moderate but real
signal, not a strong scalability claim.

The heavy MAAT variant is a useful negative control: increasing structural
weight without enough direct verification pressure does not improve the result
and is slightly worse than the anchored MAAT brancher.

## Outputs

The script writes:

- `outputs_paper59/paper59_solver_runs.csv`
- `outputs_paper59/paper59_summary_by_family.csv`
- `outputs_paper59/paper59_summary_aggregate.csv`
- `outputs_paper59/paper59_summary.json`
- `outputs_paper59/fig1_median_search_cost.png`
- `outputs_paper59/fig2_timeout_rate.png`
- `outputs_paper59/fig3_aggregate_cost.png`
- `outputs_paper59/fig4_pairwise_difference.png`

The optional visualisation script:

```bash
python3 sat_v16_structure_visualization.py
```

writes SAT graph drawings to:

- `outputs_paper59/sat_graph_visuals/sat_variable_cooccurrence_maat_priority.png`
- `outputs_paper59/sat_graph_visuals/sat_clause_variable_bipartite_maat_priority.png`
- `outputs_paper59/sat_graph_visuals/sat_branching_signal_components.png`

These figures show SAT formulas as graph objects:

- variable co-occurrence graph: variables connected when they share clauses
- clause-variable bipartite graph: variables and clauses as separate node types
- node colour: MAAT v1.6 branching priority
- node size: connectivity support `V`

The graph figures are root-state feature maps. They show which variables the
structural score treats as pressure points, and the component plot connects
those visual pressure points to verification gain, short-clause pressure,
connectivity, and robustness closure. The performance tables then show the
limit: this visible centrality signal is useful but still weaker than
short-clause-specialised MOMS and Jeroslow-Wang.

## Data and License Note

No external CNF dataset is redistributed. Instances are generated locally from
standard synthetic SAT-family generators. CSV/JSON/PNG outputs are derived
reproducibility artifacts.
