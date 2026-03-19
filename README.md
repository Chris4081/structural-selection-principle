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

### Paper 10 — Dimensional Dependence of Scaling
**Dimensional Dependence of Entropy–Information Scaling in Nonlinear Field Systems**
*Evidence for a Dimension-Dependent Universality Class*

**Core result:** Systematic dimensional trend of the scaling exponent:

| Dimension | α_eff | Spearman r_s | Dominant factor |
|-----------|-------|-------------|-----------------|
| **1D** | ≈ 2.5 | 0.734 | Ṡ/I balance |
| **2D** | ≳ 3.0 (within tested range) | **0.870** | I_nn loss |

The effective exponent increases with dimension; the 2D value is a lower bound.
The Spearman correlation strengthens from 1D to 2D, indicating increased robustness.

> The dimensional dependence of α is a property of the projection from the
> structural free energy onto the entropy–information plane, rather than a
> fundamental property of the selection principle itself.

---


### Paper 11 — Entropy–Information Scaling (3D)
**Persistence of Entropy–Information Scaling in Three Dimensions**
*Breakdown of a Simple Exponent Trend in Nonlinear Field Systems*

**Core result:** Scaling persists in 3D but the exponent structure changes:

| Quantity | 1D | 2D | 3D |
|----------|-----|-----|-----|
| Spearman r_s | 0.734 | 0.870 | **0.850** |
| α_eff | ≈ 2.5 | ≳ 3.0 | plateau [2.8, 3.5] |
| R² | moderate | 0.74 | **0.37** |
| Dominant | Ṡ/I balance | I_nn loss | balanced |

The notion of a single effective exponent breaks down in 3D: instead of a sharp
optimum, a broad plateau of near-maximal correlation appears over α ∈ [2.8, 3.5].
The reduced R² indicates a degeneracy of effective scaling descriptions.

> The entropy–information scaling remains robust across dimensions,
> but transitions from a well-defined power-law exponent to a broad
> effective scaling regime in higher dimensions.

---

### Paper 12 — Cross-Dimensional Collapse Analysis
**Collapse Analysis of Entropy–Information Scaling Across Spatial Dimensions**
*Shared Structure and Dimensional Degeneracy in Nonlinear Field Systems*

**Core idea:** Direct visual comparison of the scaling relation across 1D, 2D, and 3D
via collapse plots, testing whether a unified cross-dimensional description exists.

**Core result:** Partial collapse observed — shared monotonic structure, no universal curve.

| View | Observation |
|------|-------------|
| Log-x | shared monotonic trend across 1D, 2D, 3D |
| Log-log | similar trend with increasing 3D dispersion |
| Normalized | partial overlap; 3D scatter broader |

> *Scaling survives, but the exponent does not.*
> The partial collapse reflects a shared structural principle
> expressed through dimension-dependent effective parameters.

**Script:** `Collapse-Plot.py`

**Requires the following CSV files** (how to generate them):

| CSV file | Generated by |
|----------|-------------|
| `cci_entropy_information_test.csv` | `python cci_alpha_scaling.py` |
| `2d_cci_entropy_information_test.csv` | `python cci_entropy_scaling_2d.py` |
| `3d_cci_entropy_information_test.csv` | `python cci_entropy_scaling_3d.py` |

---

### Paper 13 — Multi-Parameter Scaling Analysis
**Beyond Single-Exponent Scaling: Multi-Parameter Entropy–Information Relations**
*in Nonlinear Field Systems*

**Core idea:** Extends the single-exponent ansatz CCI ∝ Ṡ^α / I_nn^α to a
two-parameter model CCI ∝ Ṡ^a · I_nn^{-b}, testing whether the 3D exponent
plateau reflects a deeper multi-channel structure.

**Core result:** The entropy and information weights vary systematically with dimension:

| Dimension | a | b | R² |
|-----------|---|---|-----|
| **1D** | −1.60 | 0.33 | 0.76 |
| **2D** | −0.08 | 0.50 | 0.74 |
| **3D** | +0.30 | 0.30 | 0.37 |

The entropy contribution **changes sign** from 1D to 3D — a qualitative shift
in the role of dynamical activity across dimensions.

> The exponent description fails not because the scaling breaks, but because
> the system transitions to a higher-dimensional multi-parameter scaling structure.

**Script:** `multiparameter_fit.py`

**Requires the following CSV files** (see Paper 12 for how to generate them):

| CSV file | Generated by |
|----------|-------------|
| `cci_entropy_information_test.csv` | `python cci_alpha_scaling.py` |
| `2d_cci_entropy_information_test.csv` | `python cci_entropy_scaling_2d.py` |
| `3d_cci_entropy_information_test.csv` | `python cci_entropy_scaling_3d.py` |

---

### Paper 14 — Scaling Manifold Framework
**Scaling Laws as Projections: Evidence for a Dimension-Dependent Scaling Manifold**
*in Nonlinear Field Systems*

**Core idea:** The observed power-law scaling is not fundamental, but represents a
local projection of a higher-dimensional geometric object — the *scaling manifold*:

```
M_scale = { (Ṡ, I_nn, d, CCI) | CCI = F(Ṡ, I_nn, d) }
```

Effective exponents are local tangent directions of this manifold; exponent degeneracy
in 3D reflects its higher-dimensional geometry.

**Empirical evidence:**

| Evidence | Result |
|----------|--------|
| PCA of log-space observables | PC1 dominates in 1D/2D; PC2/PC3 gain in 3D |
| Joint PCA (1D+2D+3D) | shared structure with dimension-dependent spread |
| ξ_aniso quartile-split (3D) | plateau width Δα varies from 0.20 to 2.20 across quartiles |
| Continuous correlation | r_s(ξ_aniso, R_α) = −0.721, p = 0.0001 |

**Key result:**
> Scaling laws are coordinate-dependent projections of a scaling manifold.
> Universality may extend from scaling exponents to geometric structures.

**Scripts:** `manifold_geometry_plot.py`, `xi_aniso_full_test.py`

**Requires the following CSV files** (see Paper 12 for how to generate them):

| CSV file | Generated by |
|----------|-------------|
| `cci_entropy_information_test.csv` | `python cci_alpha_scaling.py` |
| `2d_cci_entropy_information_test.csv` | `python cci_entropy_scaling_2d.py` |
| `3d_cci_entropy_information_v2.csv` | `python cci_entropy_scaling_3d_v2.py` |

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
├── cci_entropy_scaling_2d.py              ← Paper 09: 2D α-scan + log-fit
├── figure_dimensional_trend.py            ← Paper 10: Figure 1 (α trend)
├── figure_scaling_comparison.py           ← Paper 10: Figure 2 (1D vs 2D)
├── cci_entropy_scaling_3d.py              ← Paper 11: 3D ensemble + figures
├── Collapse-Plot.py                       ← Paper 12: cross-dimensional collapse
├── multiparameter_fit.py                  ← Paper 13: multi-parameter scaling fit
├── cci_entropy_scaling_3d_v2.py           ← Paper 14: 3D with directional MI
├── manifold_geometry_plot.py              ← Paper 14: PCA + log-space figures
├── xi_aniso_full_test.py                  ← Paper 14: ξ_aniso quartile analysis
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
| `figure_dimensional_trend.py` | 10 | α(d) and r_s(d) trend figures |
| `figure_scaling_comparison.py` | 10 | 1D vs 2D scaling side-by-side |
| `cci_entropy_scaling_3d.py` | 11 | 3D ensemble, α-scan, dimension comparison |
| `Collapse-Plot.py` | 12 | cross-dimensional collapse (requires all 3 CSVs) |
| `multiparameter_fit.py` | 13 | two-parameter fit a,b across 1D–3D |
| `cci_entropy_scaling_3d_v2.py` | 14 | 3D sim with directional MI (I_x, I_y, I_z) |
| `manifold_geometry_plot.py` | 14 | PCA + log-space manifold figures |
| `xi_aniso_full_test.py` | 14 | ξ_aniso quartile-split + regression |

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
Paper 10 (dim.):   α increases with d        r_s: 0.734 → 0.870
Paper 11 (3D):     plateau in 3D              r_s = 0.850, α ∈ [2.8, 3.5]
Paper 12 (collapse):partial collapse 1D–3D     shared structure, no universal exponent
Paper 13 (multi-p): sign flip a: -1.6→+0.3       dimension-dependent scaling family
Paper 14 (manifold):scaling = projection          r_s(ξ_aniso,R_α) = -0.721, p=0.0001
```

---

## License

MIT License — free to use, modify, and share with attribution.
