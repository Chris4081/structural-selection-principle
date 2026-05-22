# MAAT v1.6 Structural Search Geometry Toy Benchmark

This experiment tests a first MAAT v1.6 idea:

> Structural selection can be interpreted as priority-modulating search
> geometry over complex candidate spaces, but it must remain anchored to
> verification.

The benchmark uses synthetic two-dimensional search spaces with four families:

- `smooth`
- `trap`
- `bridge`
- `trap_bridge`

For each family, the script compares:

- `random_frontier`
- `distance_only`
- `energy_only`
- `connectivity_only`
- `maat_v16_anchored`
- `maat_v16_heavy`

The key lesson is deliberately modest. The anchored v1.6 score is robust and
beats energy-only search in trap-like cases, but simple distance/connectivity
baselines dominate these clean grid tasks. That dominance is expected because
the target is visible, Euclidean distance is cheap, and the bridges are simple
graph bottlenecks. Heavy structural weighting slows search, showing that
structural support must modulate the search rather than replace the validation
target.

## Run

```bash
python3 structural_search_geometry_v16.py
```

## Outputs

The script writes reproducibility artifacts to `outputs_v16/`:

- `v16_search_results.csv`
- `v16_search_summary.csv`
- `v16_summary.json`
- `fig1_success_rate.png`
- `fig2_median_expansions.png`
- `fig3_example_paths.png`
- `fig4_support_maps.png`

## Current Core Results

| Algorithm | Success rate | Median expansions |
|---|---:|---:|
| distance_only | 1.000 | 59.0 |
| connectivity_only | 1.000 | 59.0 |
| maat_v16_anchored | 1.000 | 92.0 |
| energy_only | 1.000 | 106.0 |
| maat_v16_heavy | 1.000 | 234.0 |
| random_frontier | 0.741 | 702.5 |

In trap-like families, `maat_v16_anchored` improves over `energy_only`:

| Family | MAAT anchored | Energy-only |
|---|---:|---:|
| trap | 84.5 | 150.5 |
| trap_bridge | 92.0 | 183.5 |

This corresponds to roughly 44% fewer median expansions in `trap` and 50%
fewer in `trap_bridge` relative to `energy_only`.

## Scientific Status

This is a toy search-geometry benchmark, not an algorithmic optimality proof.
It is intended to test the v1.6 design principle:

> Search geometry should modulate verification-guided search, not replace the
> target, proof obligation, or validation criterion.

The next meaningful tests should use richer spaces where local evaluation is
costly or noisy and where bridges are not visible as simple geometric
bottlenecks: theorem-search graphs, CDCL/SAT state spaces, graph-construction
problems, model-configuration spaces, or optimisation tasks with expensive
validation.

## License and Data

No external dataset is redistributed. All CSV/JSON/PNG outputs are generated
locally from synthetic search-space instances.
