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

### Paper 01 — General Framework
**A Structural Selection Principle for Dynamically Consistent Field Configurations**
*Maximum-Entropy Weighting on Physical Solution Spaces*

**Core idea:** δS[Φ]=0 defines possible configurations; a structural energy functional E[Φ] selects the physically preferred ones via a MaxEnt Boltzmann weight.

---

### Paper 02/03 — 1D and 2D Topological Selection
**Structural Selection in Nonlinear Field Theories:**
*Kink Selection in 1D and Activity-Dependent Domain-Wall Selection in 2D*

**Core results:**

| System | Result |
|--------|--------|
| **1D** | Kinks dominate 81/81 parameter settings, weight 99.96% |
| **2D** | Selection transition at s★_c ∈ (0.20, 0.42) |

---

### Paper 04 — CCI Framework
**A Structural Free Energy on the Solution Manifold**
*Coherence Diagnostics, Critical Coherence Index, and Regime Classification*

**Core idea:** CCI = destabilising / stabilising — a Reynolds-type number for field systems.

```
CCI < τ₁       →  Ordered   (solitons, kinks)
τ₁ ≤ CCI < τ₂  →  Critical  (intermittency)
CCI ≥ τ₂       →  Chaotic   (instability-dominated)
```

---

### Paper 05 — Physical Grounding of the CCI
**Physical Grounding of the Critical Coherence Index:**
*Entropy Production and Structural Information*

**Core result:** CCI correlates with Ṡ/I ratio:

| Quantity | Spearman r |
|----------|-----------|
| **Ṡ/I ratio** | **0.760** |
| Ṡ alone | 0.530 |
| I alone | 0.246 |

---

### Paper 06 — Continuous Ensemble Study
**Structural Free Energy and the CCI in a Continuous Ensemble of φ⁴ Field Configurations**
*An Extended Numerical Study with Mixed Initial Conditions*

**Core result:** CCI is primarily a proxy for structural free energy (160-trajectory ensemble):

| Quantity | Spearman r |
|----------|-----------|
| **F_struct** | **0.878** ← dominant |
| Ṡ/I ratio | 0.288 |
| I_nn alone | −0.341 |

> Note:
> The entropy–information ratio should be understood as an empirical
> scaling proxy, while the structural free energy provides a more
> global effective description of solution selection.

---

### Paper 07 — Entropy–Information Scaling (1D)
**Entropy–Information Scaling of the Critical Coherence Index in Nonlinear Field Dynamics**

**Core result:** CCI consistent with power-law scaling in 1D:

```
CCI ∝ Ṡ_cg⁺ / I_nn^α,   α ≈ 2.5
Spearman r_s = 0.734,  p ≈ 1.6×10⁻²¹
```

Phase diagrams over (σ_π, σ_noise) reveal a clear ordered → chaotic transition.
The exponent α is empirical; its theoretical origin is an open problem.

---

### Paper 08 — Landau-Type Theory
**Landau-Type Theory of Structural Selection in Nonlinear Field Systems**
*Entropy–Information Scaling and the Critical Coherence Index*

**Core idea:** CCI and structural free energy fit a Landau framework:

```
F_struct[O] = a(μ)·O² + b·O⁴ + γ(∇O)²

Order parameter:   O ~ I_nn  (structural information)
Control parameter: μ ~ CCI

μ < μ_c  →  Ordered   (high coherence)
μ = μ_c  →  Critical  (phase boundary)
μ > μ_c  →  Chaotic   (low coherence)
```

The Landau coefficients are not derived from first principles —
this is a theoretical proposal and interpretive framework.

---

### Paper 09 — Entropy–Information Scaling (2D)
**Entropy–Information Scaling of the CCI in Two-Dimensional φ⁴ Field Dynamics**
*Dimensional Dependence and Robustness of Structural Instability*

**Core result:** Scaling persists in 2D with shifted exponent:

| Quantity | 1D (Paper 07) | 2D (Paper 09) |
|----------|--------------|--------------|
| Spearman r_s | 0.734 | **0.870** |
| Best α | ≈ 2.5 | ≈ 3.0 (within tested range) |
| Dominant factor | Ṡ/I balance | **I_nn loss** |

Structural information (I_nn) becomes the dominant factor in 2D.
Log-linear regression: R² ≈ 0.74.

---

## Repository Structure

```
structural-selection-principle/
│
├── structural_selection_flrw_v4.py        ← Paper 01: FLRW benchmark
├── maat_structural_selection_study_v2.py  ← Paper 02: 1D kink sweep
├── maat_structural_selection_2d_sstar.py  ← Paper 03: 2D s★ sweep
├── cci_entropy_test.py                    ← Paper 05: CCI vs Ṡ/I
├── cci_continuous_ensemble.py             ← Paper 06: 160-run ensemble
├── cci_alpha_scaling.py                   ← Paper 07: α-scan + log-fit
├── cci_phase_diagram.py                   ← Paper 07: phase diagrams
├── cci_entropy_scaling_2d.py              ← Paper 09: 2D scaling test
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
| `structural_selection_flrw_v4.py` | 01 | FLRW 3-case MaxEnt benchmark |
| `maat_structural_selection_study_v2.py` | 02 | 1D kink dominance, 81 settings |
| `maat_structural_selection_2d_sstar.py` | 03 | 2D s★ phase diagram |
| `cci_entropy_test.py` | 05 | CCI vs Ṡ/I, discrete classes |
| `cci_continuous_ensemble.py` | 06 | CCI vs F_struct, 160 runs |
| `cci_alpha_scaling.py` | 07 | α-scan, log-fit, scatter plots |
| `cci_phase_diagram.py` | 07 | 9×9 phase diagram heatmaps |
| `cci_entropy_scaling_2d.py` | 09 | 2D α-scan + log-fit |


---

## Key Results at a Glance

```
Paper 01 (SSP):    slow-roll preferred        E = 0.889 < E_deSitter
Paper 02 (1D):     kink wins 81/81            weight 99.96%
Paper 03 (2D):     selection transition       s★_c ∈ (0.20, 0.42)
Paper 05 (ground): CCI vs Ṡ/I               Spearman r = 0.760
Paper 06 (ensemb): CCI ≈ F_struct            Spearman r = 0.878
Paper 07 (1D sc.): CCI ∝ Ṡ/I^α             r_s = 0.734, α ≈ 2.5
Paper 08 (Landau): O ~ I_nn, μ ~ CCI        theoretical proposal
Paper 09 (2D sc.): scaling persists in 2D    r_s = 0.870, α ≈ 3.0
```

---

## License

MIT License — free to use, modify, and share with attribution.
