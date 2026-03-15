# Structural Selection Principle — FLRW Benchmark

Companion code for:

> Krieg, C. (2025). *A Structural Selection Principle for Dynamically Consistent Field Configurations*. Preprint, available from the author.

## Overview

This repository provides a reproducible numerical benchmark for the **Structural Selection Principle** proposed in the paper above.

The core idea: among all solutions of the fundamental field equations `δS[Φ] = 0`, physically preferred configurations are those that minimise a covariant structural energy functional `E[Φ]`, derived from a Maximum-Entropy principle on the solution manifold.

## The Three FLRW Benchmarks

The benchmark evaluates `E[Φ]` for three cosmological backgrounds:

| Case | Description | Expected |
|------|-------------|----------|
| A | Static-field de Sitter (KG-compatible) | Reference |
| B | Slow-roll scalar-tensor FLRW | **Preferred** (lowest E) |
| C | Ghost-near stress test | Rejected (highest E) |

**Result:** `E_B < E_A < E_C` — the functional correctly identifies the dynamically active slow-roll background as structurally preferred.

> **Disclaimer:** This benchmark demonstrates the *discriminative capacity* of the functional. It does not validate the selection principle as a law of nature, nor does it derive cosmological predictions from first principles.

## The Five Structural Diagnostics

| Symbol | Name | Physical content |
|--------|------|-----------------|
| `H` | Residual Stability | Distance from solution manifold |
| `B` | Conservation | Energy-momentum conservation |
| `S_C` | Dynamical Activity | Regulated non-equilibrium activity |
| `N` | Connectivity | Inter-field coupling (= 0 single field) |
| `Σ` | Consistency | Ghost-freeness / EFT integrity |

## Usage

```bash
pip install numpy
python structural_selection_flrw_v4.py
```

## Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `XI` | 1e-3 | Non-minimal coupling ξ |
| `M_C` | 1.0 | Coherence field scale (Planck units) |
| `H0` | 1e-3 | Hubble constant (normalised) |
| `W_M` | -1.0 | Matter equation of state |
| `LAM_*` | 1.0 | Lagrange multipliers (equal weights) |

## Output

```
Case                             F(C)        H        B       S      Sigma          E
A: de Sitter                   1.0000  -0.0000  -0.0000  0.0000    -1.0000     1.0000
B: slow-roll FLRW              1.0001  -0.0000  -0.0000  0.1111    -0.9998     0.8887
C: ghost-near (stress test)    0.0500  -0.2036  -0.0000  0.1111  -399.8401   399.9326

Ranking E_B < E_A < E_C: YES
```

## Related Work

- MAAT-Core framework: [Zenodo DOI 10.5281/zenodo.18489336](https://doi.org/10.5281/zenodo.18489336)
- Companion paper (CEST model): available from the author

## Author

Christof Krieg  
Independent Researcher  
[maat-research.com](https://maat-research.com) · [GitHub: Chris4081](https://github.com/Chris4081)

## License

MIT License — free to use, modify, and share with attribution.
