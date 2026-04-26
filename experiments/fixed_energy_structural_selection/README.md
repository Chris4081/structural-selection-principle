# Fixed-Energy Structural-Selection Benchmarks

This directory contains a compact reproducibility bundle for testing whether a
structural selection functional distinguishes coherent field configurations
from rough or instability-dominated perturbations beyond ordinary energy
ranking.

The benchmarks are exploratory and operational. They are not a proof of a
complete theory of nature. Their purpose is to provide a clean numerical stress
test for the claim:

```text
structural quality is not reducible to total energy alone
```

## Contents

```text
fixed_energy_structural_selection/
├── README.md
├── structural_selection_fixed_energy_benchmarks.py
├── fixed_energy_structural_selection_results.json
└── fixed_energy_structural_selection_plots/
```

## Requirements

```bash
python3 -m pip install numpy matplotlib
```

## Reproduce

Run from this directory:

```bash
python3 structural_selection_fixed_energy_benchmarks.py
```

The script generates:

- `fixed_energy_structural_selection_results.json`
- `fixed_energy_structural_selection_plots/phi4_1d_fields.png`
- `fixed_energy_structural_selection_plots/phi4_2d_fields.png`
- `fixed_energy_structural_selection_plots/energy_vs_structural.png`
- `fixed_energy_structural_selection_plots/equal_energy_fields.png`
- `fixed_energy_structural_selection_plots/equal_energy_structural_score.png`
- `fixed_energy_structural_selection_plots/sine_gordon_fields.png`
- `fixed_energy_structural_selection_plots/sine_gordon_energy_vs_structural.png`

## Benchmarks

| Benchmark | Question tested |
| --- | --- |
| 1D `phi^4` kink ranking | Does the structural score prefer coherent kink-like configurations over rough chaos? |
| 2D `phi^4` domain-wall ranking | Does the same ordering persist in a two-dimensional domain-wall setting? |
| Equal-energy `phi^4` test | Can the structural score separate smooth and rough perturbations at approximately fixed total energy? |
| Sine-Gordon soliton ranking | Does the structural ordering generalise beyond the `phi^4` model? |

## Notes on the Plots

The original 2D plot layout used four panels in a single row, which could cause
titles, labels, and the colorbar to overlap. The current version uses a 2x2
`constrained_layout` grid with a shared colorbar, so the 2D figure should remain
readable in both the repository and paper exports.

## Scientific Status

The functional used here is an operational proxy:

```text
F_MAAT = - sum_a log(epsilon + Gamma_a),   Gamma_a = 1 / (1 + d_a)
```

with sectors for equation consistency, conservation consistency, activity,
connectivity, and robustness. The diagnostics are intentionally simple so the
tests remain transparent and reproducible.

The strongest result is the equal-energy test: configurations can have nearly
identical total energy while receiving very different structural scores.

