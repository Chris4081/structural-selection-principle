# Structural Selection Principle

Companion code and papers for the research programme:

> Krieg, C. (2026). *A Structural Selection Principle for Dynamically Consistent Field Configurations — Maximum-Entropy Weighting on Physical Solution Spaces.*

> Krieg, C. (2026). *Topological Selection in Nonlinear Field Theories: Robust Dominance of Kink Configurations under Structural Energy Weighting.*

---

## Repository Structure

```
structural-selection-principle/
├── documentation/
│   ├── A_Structural_Selection_Principle_for_Dynamically_Consistent_Field_Configurations_large_Maximum_Entropy_Weighting_on_Physical_Solution_Spaces.pdf   ← Paper 1: General framework
│   └── Topological_Selection_in_Nonlinear_Field_Theories__Robust_Dominance_of_Kink_Configurations_under_Structural_Energy_Weighting.pdf   ← Paper 2: Topological selection
├── structural_selection_flrw_v4.py   ← FLRW benchmark (Paper 1)
├── maat_structural_selection_study_v2.py  ← Lattice study (Paper 2)
├── requirements.txt
└── README.md
```

---

## Paper 1 — Structural Selection Principle (General Framework)

**Full paper:** [PDF on GitHub](https://github.com/Chris4081/structural-selection-principle/blob/main/documentation/A_Structural_Selection_Principle_for_Dynamically_Consistent_Field_Configurations_large_Maximum_Entropy_Weighting_on_Physical_Solution_Spaces.pdf) · [Academia.edu](https://independent.academia.edu/KriegChristof)

### Core Idea

Among all solutions of `δS[Φ] = 0`, physically preferred configurations
are those that minimise a structural energy functional `E[Φ]`, derived
from a Maximum-Entropy principle on the solution manifold:

```
Φ_phys = arg min_{δS=0} E[Φ]
```

### Five Structural Diagnostics

| Symbol | Name | Physical content |
|--------|------|-----------------|
| `H` | Residual Stability | Distance from solution manifold |
| `B` | Conservation | Energy-momentum integrity |
| `S_C` | Dynamical Activity | Regulated non-equilibrium activity |
| `N` | Connectivity | Inter-field coupling |
| `Σ` | Consistency | Ghost-freeness / EFT integrity |

### FLRW Benchmark

```bash
python structural_selection_flrw_v4.py
```

| Case | Description | E |
|------|-------------|---|
| A: de Sitter | Static field | 1.000 |
| B: slow-roll | Dynamically active | **0.889** ← preferred |
| C: ghost-near | Stress test | 399.933 |

---

## Paper 2 — Topological Selection (Lattice Study)

### Core Result

Four classes of `φ⁴` field configurations compete under structural
energy weighting across **81 parameter combinations**:

| Class | E_mean | Weight | Result |
|-------|--------|--------|--------|
| Vacuum | 4.479 | ≈ 0 | Rejected |
| **Kink** | **0.015** | **99.96%** | **← Preferred** |
| Localised | 1.941 | 0.04% | Rejected |
| Chaotic | 9.210 | ≈ 0 | Rejected |

**Kink configurations win in all 81 parameter settings tested.**

### Run the Study

```bash
python maat_structural_selection_study_v2.py
```

Outputs saved to `maat_results_v2/`:
- `robustness.csv` — wins per class across all parameter settings
- `class_sweep.csv` — full sweep results
- `class_dominance.png` — bar chart
- `heatmap_kink_weight_beta4.png` — kink weight heatmap

---

## Installation

```bash
pip install numpy pandas matplotlib
```

---

## Parameters (Paper 2)

| Parameter | Value | Description |
|-----------|-------|-------------|
| `N` | 128 | Lattice sites |
| `RUNS_PER_CLASS` | 12 | Runs per class |
| `ACTIVITY_TARGET` | 0.42 | Preferred activity level |
| `BETA_LIST` | [2,4,8] | Inverse temperatures |
| `SEED` | 42 | Random seed |

---

## Author

Christof Krieg — Independent Researcher  
[Academia.edu](https://independent.academia.edu/KriegChristof) · [GitHub: Chris4081](https://github.com/Chris4081)

## License

MIT License — free to use, modify, and share with attribution.
