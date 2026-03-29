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

### Paper 16 — Plateau Degeneracy Measure
**Plateau Degeneracy in Entropy–Information Scaling:**
*A Quantitative Measure for the Breakdown of Single-Exponent Descriptions in Nonlinear Field Systems*

**Core idea:** Formalises the breakdown of single-exponent scaling observed in 3D
by introducing a rigorous plateau degeneracy framework.

**Key definitions:**

| Measure | Formula | Meaning |
|---------|---------|---------|
| Plateau set | P_ε = {α \| r_s(α) ≥ r_s^max − ε_plat} | near-optimal exponents |
| Width | D_width = α_max^(ε) − α_min^(ε) | plateau size |
| Normalised | D_norm = D_width / (α_max − α_min) | in [0, 1] |
| Flatness | D_flat = Var(r_s) within P_ε | 0 = perfectly flat |
| Combined | D = D_width / (1 + D_flat) | joint measure |

**Exact results** (ε_plat = 0.01 · r_s^max, scan range α ∈ [0.5, 3.5], step h = 0.1):

| Dimension | α* | r_s^max | Plateau | D_norm | n |
|-----------|-----|---------|---------|--------|---|
| 1D | 3.5 | 0.748 | [2.9, 3.5] | 0.200 | 7 |
| 2D | 3.5 | 0.881 | [3.1, 3.5] | 0.133 | 5 |
| 3D | 2.8 | 0.849 | [2.6, 3.5] | 0.300 | 10 |

**Key finding:**
> D_norm(2D) < D_norm(1D) < D_norm(3D) — non-monotonic dimensional dependence.
> The breakdown of single-exponent scaling is not gradual but occurs via a
> dimensional transition between 2D and 3D.

**Script:** `plateau_degeneracy_exact.py`

**Requires the following CSV files:**

| CSV file | Used for |
|----------|----------|
| `cci_entropy_information_test.csv` | 1D r_s(α) scan (Paper 07) |
| `2d_cci_entropy_information_test.csv` | 2D r_s(α) scan (Paper 09) |
| `3d_cci_alpha_scan.csv` | 3D r_s(α) scan (Paper 11) |

---

### Paper 17 — Predictive Power of the CCI
**Predictive Power of the Critical Coherence Index:**
*From Structural Diagnostics to Regime Classification*

**Core idea:** Tests whether the CCI can predict the dynamical regime of a system
using a threshold-based classifier on a 120-sample 1D benchmark dataset.

**Core results:**

| Classifier | Accuracy | τ₁ | τ₂ |
|------------|----------|-----|-----|
| **CCI-only** | **1.000 (120/120)** | 0.2785 | 0.3527 |
| F_struct-only | 0.8417 (101/120) | 0.4204 | 0.4984 |

5-fold cross-validation (CCI-only): **Ā = 0.992 ± 0.019**
(folds: [1.000, 1.000, 1.000, 0.958, 1.000])

> The CCI alone is sufficient to fully separate ordered, critical, and chaotic
> regimes. The near-perfect cross-validated accuracy indicates robustness
> across data splits.

**Script:** `paper17_analysis.py`

**Requires:** `cci_entropy_information_test.csv` (1D benchmark dataset)

---

### Paper 18 — Dynamical Origin of Structural Selection
**A Minimal Dynamical Model for Structural Selection:**
*A Field Equation for the Critical Coherence Index*

**Core idea:** Promotes the Critical Coherence Index (CCI) from a diagnostic observable to a dynamical field C(x,t), governed by a minimal nonlinear relaxation equation with gradient-flow structure.

**Model equation:**

| Term | Expression | Meaning |
|------|------------|--------|
| Dynamics | ∂ₜC = D ∇²C + aC − bC³ | coherence-field evolution |
| Free energy | F[C] = ∫ [ (D/2)\|∇C\|² − (a/2)C² + (b/4)C⁴ ] dx | structural potential |
| Gradient flow | ∂ₜC = −δF/δC | energy minimisation |

**Stationary structure:**

| Regime | Condition | Solution |
|--------|----------|----------|
| Disordered | a < 0 | C = 0 (stable) |
| Critical point | a = 0 | pitchfork bifurcation |
| Ordered | a > 0 | C = ±√(a/b) |

**Numerical results** (a = b = 1, D = 0.8):

| Quantity | Value | Interpretation |
|----------|--------|----------------|
| Fixed points | ±1.0 | stable coherence states |
| ⟨C⟩ (1D) | -0.080 | domain balance |
| std(C) (1D) | 0.939 | saturation near ±1 |
| Energy drop | 0.057 → -27.076 | monotonic decay |

**Key findings:**
> The coherence-field equation reproduces spontaneous symmetry breaking,
> domain formation, and energy minimisation within a minimal framework.
> Structural selection emerges dynamically as gradient descent in a quartic free-energy landscape.

**Script:** `coherence_simulation.py`

**Features:**

| Component | Description |
|----------|-------------|
| 0D dynamics | relaxation to fixed points |
| 1D field | domain formation with periodic BC |
| Energy tracking | verifies gradient-flow behaviour |

---

### Paper 19 — Active Structural Control
**Active Structural Control in Nonlinear Field Systems:**
*Parameter-Dependent Steering of Coherence*

**Core idea:** Extends the CCI framework from a passive diagnostic into an active
control paradigm. A feedback control term is introduced into the coherence-field
dynamics, formalising the complete loop φ → C → U → φ. Structural coherence
becomes a controllable quantity rather than a purely emergent observable.

**Control scheme:**

| Step | Expression | Role |
|------|------------|------|
| Field → Coherence | C = G_σ[φ²] − (G_σ[φ])² | coarse-grained variance |
| Coherence → Control | U = −(C − C*) | feedback signal |
| Control → Field | φ̈ = ∇²φ − (φ³−φ) − γφ̇ + λU | controlled dynamics |

**Best-performing configuration** (γ = 0.02, λ = 0.24):

| Observable | Baseline | Controlled | Change |
|------------|----------|------------|--------|
| CCI | 0.0579 | 0.0540 | −6.7% |
| F_struct | 3.3433 | 3.1322 | −6.3% |
| I_nn | 0.1080 | 0.2746 | +154% |

**Two operational regimes:**

| Regime | Damping | Effect |
|--------|---------|--------|
| Structure refinement | γ ≳ 0.05 | F_struct ↓, I_nn ↑, CCI stable |
| Regime steering | γ ≲ 0.03 | CCI ↓, genuine regime transition |

> Structural selection can be formulated as a control problem.
> The feedback loop φ → C → U → φ enables target-driven steering
> of the system in coherence space, with qualitatively distinct
> effectiveness regimes depending on damping strength.

**Script:** `active_control_phi4.py`

**Requires:** 2D φ⁴ field simulation with RK4 integration (see repository)

---

### Paper 20 — Universal Stability Test
**Universal Stability of Structural Observables Across Dimensions:**
*A Comparative Test of CCI, Structural Free Energy, and Entropy–Information Scaling*

**Core idea:** Tests which structural observables remain dimensionally robust
across 1D, 2D, and 3D φ⁴ systems. Introduces relative drift D_rel = max(μ_d)/min(μ_d)
as a scale-invariant measure of dimensional stability.

**Core results:**

| Observable | 1D Mean | 2D Mean | 3D Mean | Max/Min Drift |
|------------|---------|---------|---------|---------------|
| **CCI** | 0.303 | 0.457 | 0.601 | **1.98** |
| **F_struct** | 0.430 | 0.641 | 0.824 | **1.92** |
| MI | 0.574 | 0.449 | 0.393 | 1.46 |
| Ṡ⁺ | 0.263 | 0.080 | 0.111 | 3.30 |
| **Ratio** | 6.81 | 71.64 | 13563 | **1992** |

**Key finding:**
> CCI and F_struct remain comparatively stable across dimensions (drift < 2),
> while the entropy–information ratio loses consistency as a dimensionally robust
> observable (drift ≈ 2000). The instability is not merely quantitative but
> structural: it reflects a projection-dependent construction that cannot define
> a consistent observable across dimensions.

**Hierarchical interpretation:**

| Observable | Role |
|------------|------|
| CCI | comparatively robust structural coordinate |
| F_struct | stable global effective functional |
| Ratio | structurally unstable projection-dependent proxy |

> Structural observables (CCI, F_struct) define intrinsic coordinates on the
> solution manifold. Scaling ratios correspond to dimension-dependent projections.

**Script:** `uni_stability_test.py`

**Requires the following CSV files** (how to generate them):

| CSV file | Generated by |
|----------|-------------|
| `cci_entropy_information_test.csv` | `python cci_alpha_scaling.py` |
| `2d_cci_entropy_information_test.csv` | `python cci_entropy_scaling_2d.py` |
| `3d_cci_entropy_information_test.csv` | `python cci_entropy_scaling_3d.py` |

---

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
├── plateau_degeneracy_exact.py            ← Paper 16: plateau degeneracy measures
├── paper17_analysis.py                    ← Paper 17: CCI regime classifier + CV
├── coherence_simulation.py                ← Paper 18: coherence-field dynamics
├── active_control_phi4.py                 ← Paper 19: active structural control
├── uni_stability_test.py                  ← Paper 20: cross-dimensional stability test
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
| `plateau_degeneracy_exact.py` | 16 | exact plateau degeneracy D_width, D_norm, D_flat, D |
| `paper17_analysis.py` | 17 | CCI threshold classifier, F_struct comparison, 5-fold CV |
| `coherence_simulation.py` | 18 | coherence-field dynamics (0D + 1D, domain formation, energy decay) |
| `active_control_phi4.py` | 19 | active structural control, parameter sweep, heatmaps |
| `uni_stability_test.py` | 20 | cross-dimensional stability test, drift analysis, figure |

---

## Key Results at a Glance

```
Paper 01 (SSP):      slow-roll preferred        E = 0.889 < E_deSitter
Paper 02 (1D):       kink wins 81/81            weight 99.96%
Paper 03 (2D):       selection transition        s★_c ∈ (0.20, 0.42)
Paper 05 (ground):   CCI vs Ṡ/I                Spearman r = 0.760
Paper 06 (ensemb):   CCI ≈ F_struct             Spearman r = 0.878
Paper 07 (1D sc.):   CCI ∝ Ṡ/I^α              r_s = 0.734, α ≈ 2.5
Paper 08 (Landau):   O ~ I_nn, μ ~ CCI          theoretical proposal
Paper 09 (2D sc.):   scaling persists in 2D     r_s = 0.870, α ≈ 3.0
Paper 10 (dim.):     α increases with d         r_s: 0.734 → 0.870
Paper 11 (3D):       plateau in 3D               r_s = 0.850, α ∈ [2.8, 3.5]
Paper 12 (collapse): partial collapse 1D–3D      shared structure, no universal exponent
Paper 13 (multi-p):  sign flip a: -1.6→+0.3     dimension-dependent scaling family
Paper 14 (manifold): scaling = projection         r_s(ξ_aniso,R_α) = -0.721, p=0.0001
Paper 16 (degeneracy):D_norm non-monotonic        D(2D)<D(1D)<D(3D), transition at 2D→3D
Paper 17 (predict.): CCI perfect classifier      Acc=1.00, CV Ā=0.992±0.019
Paper 18 (dynamics): coherence field emerges     symmetry breaking, domains, F ↓ monotonic
Paper 19 (control):  active coherence steering   CCI ↓6.7%, I_nn ↑154%, two control regimes
Paper 20 (robust.):  CCI/F_struct stable           drift<2 vs ratio drift≈2000, projection vs intrinsic
```

---

www.maat-research.com

## License

MIT License — free to use, modify, and share with attribution.
