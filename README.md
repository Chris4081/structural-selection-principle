# Structural Selection Principle

Companion code and papers for an ongoing research series on
structural selection, coherence, and chaos in nonlinear field theories.

> *This work is part of an ongoing research programme on structural
> selection principles for dynamically consistent field configurations.*

---

## Papers in this Series

### Paper 1 — General Framework
**A Structural Selection Principle for Dynamically Consistent Field Configurations**

- [PDF on GitHub](https://github.com/Chris4081/structural-selection-principle/blob/main/documentation/)
- [Academia.edu](https://independent.academia.edu/KriegChristof)

**Core idea:** δS[Φ]=0 defines possible universes; structural energy E[Φ] selects the physical one.

---

### Paper 2 — 1D Topological Selection
**Topological Selection in Nonlinear Field Theories: Robust Dominance of Kink Configurations**

**Core result:** Kinks dominate across **all 81** tested parameter combinations (λ, β sweep).

| Class | Ē | Weight |
|-------|---|--------|
| Vacuum | 4.479 | ≈ 0 |
| **Kink** | **0.015** | **99.96%** |
| Localised | 1.941 | 0.04% |
| Chaotic | 9.210 | ≈ 0 |

---

### Paper 3 — 2D Selection Transition
**Activity-Target Dependence in the Structural Selection of 2D φ⁴ Field Configurations**

**Core result:** A clear selection transition at s★_c ∈ (0.20, 0.42):

| s★ | Dominant class | Domain-wall weight |
|----|---------------|--------------------|
| ≤ 0.20 | **Domain wall** | 72–87% |
| 0.42 | Chaotic | 100% |

---

### Paper 4 — CCI Framework
**A Structural Free Energy on the Solution Manifold**
*Coherence Diagnostics, Critical Coherence Index, and Regime Classification*

**Core idea:** CCI = destabilising / stabilising — a Reynolds-type number for field systems.

```
CCI < τ₁  →  Ordered  (solitons, kinks)
τ₁ ≤ CCI < τ₂  →  Critical  (intermittency)
CCI ≥ τ₂  →  Chaotic  (instability-dominated)
```

---

### Paper 5 — Physical Grounding of the CCI
**Physical Grounding of the Critical Coherence Index: Entropy Production and Structural Information**

**Core result:** CCI ≈ Ṡ_prod / (I_nn + ε) — numerically confirmed.

| Quantity | Spearman r |
|----------|-----------|
| **Ṡ/I ratio** | **0.760** ← strong |
| Ṡ alone | 0.530 |
| I alone | 0.246 |

---

## Repository Structure

```
structural-selection-principle/
├── documentation/              ← Papers (PDF)
├── structural_selection_flrw_v4.py      ← FLRW benchmark (Paper 1)
├── maat_structural_selection_study_v2.py ← 1D lattice study (Paper 2)
├── maat_structural_selection_2d_sstar_sweep.py ← 2D sweep (Paper 3)
├── cci_entropy_test.py                  ← CCI entropy test (Paper 5)
├── requirements.txt
└── README.md
```

---

## Run the Code

```bash
pip install numpy pandas matplotlib scipy
```

| Script | Paper | What it does |
|--------|-------|-------------|
| `structural_selection_flrw_v4.py` | Paper 1 | FLRW 3-case benchmark |
| `maat_structural_selection_study_v2.py` | Paper 2 | 1D kink dominance sweep |
| `maat_structural_selection_2d_sstar_sweep.py` | Paper 3 | 2D s★ phase diagram |
| `cci_entropy_test.py` | Paper 5 | CCI vs Ṡ/I correlation |

---

## Key Results at a Glance

```
FLRW (Paper 1):   slow-roll FLRW preferred  E=0.889 < E_deSitter=1.00
1D   (Paper 2):   kink wins 81/81 settings  weight 99.96%
2D   (Paper 3):   selection transition       s★_c ∈ (0.20, 0.42)
CCI  (Paper 5):   CCI ≈ Ṡ/I                Spearman r = 0.76
```

---

## Author

Christof Krieg — Independent Researcher
[Academia.edu](https://independent.academia.edu/KriegChristof) ·
[GitHub: Chris4081](https://github.com/Chris4081)

## License

MIT License — free to use, modify, and share with attribution.
