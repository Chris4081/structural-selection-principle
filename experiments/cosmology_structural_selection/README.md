# Cosmology Structural-Selection Toy Benchmark

This directory contains the reproducibility bundle for the FLRW scalar-field
toy benchmark used in the paper

**Structural Selection Across Field Theory and Cosmology: A MAAT Functional
Benchmark in Nonlinear and FLRW Toy Models**.

The benchmark extends the structural-selection diagnostics from static field
configurations to simple cosmological histories. It is an operational toy
model, not a complete cosmological theory.

## Requirements

```bash
python3 -m pip install numpy matplotlib
```

## Reproduce

Run from this directory:

```bash
python3 maat_cosmology_toy_v2.py
```

The script generates:

- `maat_cosmology_toy_v2_results.json`
- `maat_cosmology_toy_v2_plots/fig_cosmo_v2_expansion.png`
- `maat_cosmology_toy_v2_plots/fig_cosmo_v2_w.png`
- `maat_cosmology_toy_v2_plots/fig_cosmo_v2_maat_ranking.png`

## Tested Histories

| History | Interpretation |
| --- | --- |
| `slow_roll_coherent` | Coherent slow-roll-like expansion on a plateau potential. |
| `kinetic_relaxing` | Initially kinetic-dominated field relaxing under Hubble friction. |
| `static_dead` | Weakly generative expansion with low dynamical activity. |
| `ghost_like` | Rapid expansion with ghost-like kinetic sign reversal. |
| `collapse_negative_potential` | Negative-potential toy case with invalid-density/collapse behavior. |

## Scientific Status

The diagnostics are proxies:

- `H`: Friedmann consistency,
- `B`: positive density and finite equation-of-state behavior,
- `S`: controlled expansion activity,
- `V`: temporal smoothness/coherence,
- `R`: robustness against ghost, collapse, kinetic domination, or invalid density.

The core point is not that this toy model solves cosmology. The point is that
the same structural-selection architecture used in field configurations can be
applied to dynamical histories and can penalize highly productive but unstable
histories such as the ghost-like case.

