# Structural Selection Principle

Companion code and papers for an ongoing research series on
structural selection, coherence, and chaos in nonlinear field theories.

> *This work is part of an ongoing research programme on structural
> selection principles for dynamically consistent field configurations.*

**Author:** Christof Krieg — Independent Researcher
[Academia.edu](https://independent.academia.edu/KriegChristof) ·
[GitHub: Chris4081](https://github.com/Chris4081)

---

## Papers in this Series

### Paper 1 — General Framework
**A Structural Selection Principle for Dynamically Consistent Field Configurations**
*Maximum-Entropy Weighting on Physical Solution Spaces*

**Core idea:** δS[Φ]=0 defines possible configurations; a structural energy functional E[Φ] selects the physically preferred ones via a MaxEnt Boltzmann weight.

---

### Paper 2 — 1D Topological Selection
**Topological Selection in Nonlinear Field Theories: Robust Dominance of Kink Configurations**

**Core result:** Kinks dominate across **all 81** tested (λ, β) parameter combinations.

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

| s★ | Dominant class | Weight |
|----|---------------|--------|
| ≤ 0.20 | **Domain wall** | 72–87% |
| 0.42 | Chaotic | 100% |

---

### Paper 4 — CCI Framework
**A Structural Free Energy on the Solution Manifold**
*Coherence Diagnostics, Critical Coherence Index, and Regime Classification*

**Core idea:** CCI = destabilising / stabilising — a Reynolds-type number for field systems.

```
CCI < τ₁       →  Ordered   (solitons, kinks)
τ₁ ≤ CCI < τ₂  →  Critical  (intermittency)
CCI ≥ τ₂       →  Chaotic   (instability-dominated)
```

---

### Paper 5 — Physical Grounding of the CCI
**Physical Grounding of the Critical Coherence Index: Entropy Production and Structural Information**

**Core result:** CCI correlates with Ṡ/I ratio (discrete ensemble):

| Quantity | Spearman r |
|----------|-----------|
| **Ṡ/I ratio** | **0.760** |
| Ṡ alone | 0.530 |
| I alone | 0.246 |

---

### Paper 6 — Continuous Ensemble Test
**Structural Free Energy and the CCI in a Continuous Ensemble of φ⁴ Field Configurations**

**Core result:** CCI is primarily a proxy for structural free energy (160-trajectory ensemble):

| Quantity | Spearman r |
|----------|-----------|
| **F_struct** | **0.878** ← dominant |
| Ṡ/I ratio | 0.288 |
| I_nn alone | −0.341 |

---

### Paper 7 — CCI Entropy Test (Discrete)
**Physical Grounding of the CCI: Entropy Production and Structural Information in Nonlinear Field Dynamics**

Discrete 4-class ensemble (vacuum, kink, localised, chaotic).
Establishes the Ṡ/I hypothesis; limitations acknowledged.

---

### Paper 8 — Entropy–Information Scaling
**Entropy–Information Scaling of the Critical Coherence Index in Nonlinear Field Dynamics**

**Core result:** CCI is consistent with a power-law scaling relation:

```
CCI ∝ Ṡ_cg⁺ / I_nn^α,   α ≈ 2.5
Spearman r_s = 0.734,  p ≈ 1.6×10⁻²¹
```

Phase diagrams over (σ_π, σ_noise) reveal a clear ordered → chaotic transition. The exponent α is empirical; its theoretical origin is an open problem.

---

## Repository Structure

```
structural-selection-principle/
│
├── structural_selection_flrw_v4.py        ← Paper 1: FLRW benchmark
├── maat_structural_selection_study_v2.py  ← Paper 2: 1D kink sweep
├── maat_structural_selection_2d_sstar.py  ← Paper 3: 2D s★ sweep
├── cci_entropy_test.py                    ← Paper 5: CCI vs Ṡ/I
├── cci_continuous_ensemble.py             ← Paper 6: 160-run ensemble
├── cci_alpha_scaling.py                   ← Paper 8: α-scan + log-fit
├── cci_phase_diagram.py                   ← Paper 8: phase diagrams
│
├── documentation/                         ← PDFs of all papers
└── README.md
```

---

## Run the Code

```bash
pip install numpy pandas matplotlib scipy scikit-learn
```

| Script | Paper | What it does |
|--------|-------|-------------|
| `structural_selection_flrw_v4.py` | 1 | FLRW 3-case MaxEnt benchmark |
| `maat_structural_selection_study_v2.py` | 2 | 1D kink dominance, 81 settings |
| `maat_structural_selection_2d_sstar.py` | 3 | 2D s★ phase diagram |
| `cci_entropy_test.py` | 5 | CCI vs Ṡ/I, discrete classes |
| `cci_continuous_ensemble.py` | 6 | CCI vs F_struct, 160 runs |
| `cci_alpha_scaling.py` | 8 | α-scan, log-fit, scatter plots |
| `cci_phase_diagram.py` | 8 | 9×9 phase diagram heatmaps |

---

## Key Results at a Glance

```
Paper 1 (FLRW):    slow-roll preferred      E = 0.889 < E_deSitter
Paper 2 (1D):      kink wins 81/81          weight 99.96%
Paper 3 (2D):      selection transition     s★_c ∈ (0.20, 0.42)
Paper 6 (ensemble):CCI ≈ F_struct           Spearman r = 0.878
Paper 8 (scaling): CCI ∝ Ṡ/I^α            r_s = 0.734, α ≈ 2.5
```

---

## License

MIT License — free to use, modify, and share with attribution.
