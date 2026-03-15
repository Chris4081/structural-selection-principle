# Structural Selection Principle

Companion code and papers for an ongoing research series on
structural selection in nonlinear field theories.

> *This work is part of an ongoing research programme on structural
> selection principles for dynamically consistent field configurations.*

---

## Papers in this Series

### Paper 1 — General Framework
**A Structural Selection Principle for Dynamically Consistent Field Configurations**
*Maximum-Entropy Weighting on Physical Solution Spaces*

- [PDF on GitHub](https://github.com/Chris4081/structural-selection-principle/blob/main/documentation/A_Structural_Selection_Principle_for_Dynamically_Consistent_Field_Configurations_large_Maximum_Entropy_Weighting_on_Physical_Solution_Spaces.pdf)
- [Academia.edu](https://independent.academia.edu/KriegChristof)

**Core idea:** Among all solutions of `δS[Φ] = 0`, physically preferred
configurations minimise a structural energy functional `E[Φ]`, derived
from a Maximum-Entropy principle on the solution manifold.

---

### Paper 2 — 1D Topological Selection
**Topological Selection in Nonlinear Field Theories: Robust Dominance of Kink Configurations under Structural Energy Weighting**

- [Academia.edu](https://independent.academia.edu/KriegChristof)

**Core result:** In a 1D φ⁴ field theory, kink configurations dominate
the structural Boltzmann weight across **all 81 tested parameter
combinations** (λ, β sweep).

| Class | Ē | Weight |
|-------|---|--------|
| Vacuum | 4.479 | ≈ 0 |
| **Kink** | **0.015** | **99.96%** ← preferred |
| Localised | 1.941 | 0.04% |
| Chaotic | 9.210 | ≈ 0 |

---

### Paper 3 — 2D Selection Transition
**Activity-Target Dependence in the Structural Selection of Two-Dimensional φ⁴ Field Configurations**
*Dimensional Dependence and Selection Transition in the Activity Sector*

- [Academia.edu](https://independent.academia.edu/KriegChristof)

**Core result:** In 2D, a clear selection transition occurs as a function
of the activity target s★:

| s★ | Dominant class | Domain-wall weight |
|----|---------------|-------------------|
| 0.01 | **Domain wall** | 87% |
| 0.05 | **Domain wall** | 87% |
| 0.10 | **Domain wall** | 86% |
| 0.20 | **Domain wall** | 72% |
| 0.42 | Chaotic | 0% |

The activity target s★ is a **dimension-dependent calibration parameter**,
not a universal constant of the framework.

---

## Repository Structure

```
structural-selection-principle/
├── documentation/
│   ├── A_Structural_Selection_Principle_for_Dynamically_Consistent_Field_Configurations_large_Maximum_Entropy_Weighting_on_Physical_Solution_Spaces.pdf   ← Paper 1: General framework
│   └── Topological_Selection_in_Nonlinear_Field_Theories__Robust_Dominance_of_Kink_Configurations_under_Structural_Energy_Weighting.pdf   ← Paper 2: Topological selection
│   └──  Activity_Target_Dependence_in_the_Structural_Selection_of_Two_Dimensional _Field_Configurations.pdf  ← 2D selection transition
├── structural_selection_flrw_v4.py      ← FLRW benchmark (Paper 1)
├── maat_structural_selection_study_v2.py ← 1D lattice study (Paper 2)
├── maat_structural_selection_2d_sstar_sweep.py ← 2D sweep (Paper 3)
├── requirements.txt
└── README.md
```

---

## Five Structural Diagnostics

| Symbol | Name | Physical content |
|--------|------|-----------------|
| `H` | Residual Stability | Distance from solution manifold |
| `B` | Conservation | Energy-momentum integrity |
| `S_C` | Dynamical Activity | Regulated non-equilibrium activity |
| `N` | Connectivity | Inter-field coupling |
| `Σ` | Consistency | Ghost-freeness / EFT integrity |

---

## Run the Benchmarks

```bash
pip install numpy pandas matplotlib
```

**FLRW benchmark (Paper 1):**
```bash
python structural_selection_flrw_v4.py
```

**1D lattice study (Paper 2):**
```bash
python maat_structural_selection_study_v2.py
```

**2D activity-target sweep (Paper 3):**
```bash
python maat_structural_selection_2d_sstar_sweep.py
```

---

## Key Results at a Glance

```
1D (N=128):  kink       wins 81/81 parameter settings
2D (48×48):  domainwall wins majority for s★ ≤ 0.20
             chaotic    wins for s★ = 0.42
             → selection transition at s★_c ∈ (0.20, 0.42)
```

---

## Author

Christof Krieg — Independent Researcher
[Academia.edu](https://independent.academia.edu/KriegChristof) ·
[GitHub: Chris4081](https://github.com/Chris4081)

## License

MIT License — free to use, modify, and share with attribution.
