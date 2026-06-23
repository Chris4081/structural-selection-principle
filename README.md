# Structural Selection Principle

Companion code and papers for an ongoing research series on
structural selection, coherence, and chaos in nonlinear field theories.

> *This work is part of an ongoing research programme on structural
> selection principles for dynamically consistent field configurations.*

**Author:** Christof Krieg — Independent Researcher
[Academia.edu](https://independent.academia.edu/KriegChristof) ·
[GitHub: Chris4081](https://github.com/Chris4081)

---

## Start Here — Core Formalism Primer

If you are new to the series, start with **Paper 00**:

**Core Formalism of Structural Selection**  
*A Reader's Primer for the MAAT / Structural Selection Paper Series*

This primer collects the common notation, master functional, support
geometry, covariance-response closure, emergent robustness closure, CCI
diagnostics, and cosmological projection layer used throughout the later
papers. It is intended as the conceptual and mathematical entry point before
reading the numbered research papers.

| File | Role |
|------|------|
| `documentation/00_Core_Formalism_of_Structural_Selection.pdf` | compiled reader primer |
| `documentation/00_Core_Formalism_of_Structural_Selection.tex` | LaTeX source |

---

## Foundational Support Note

**Foundations of MAAT Structural Supports**  
*Why H, B, S, V, Log-Support Selection, and Emergent Robustness*

This note addresses the foundational question left open by the MaxEnt and
defect-field papers: why the canonical MAAT support language uses four
primitive sectors `(H, B, S, V)`, why defects are mapped to bounded supports
`Gamma_a = 1/(1+d_a)`, why the selection functional is logarithmic in sector
supports, and why robustness is treated as an emergent closure rather than a
fifth primitive field.

Core conclusion:

> The paper does not prove MAAT as a microscopic fundamental theory. It gives
> an operational axiomatic basis: `H,B,S,V` are a minimal diagnostic basis for
> consistency, balance, controlled activity, and connectivity; the log-support
> form follows from multiplicative viability and MaxEnt; and `R_rob` is a
> bottleneck closure over the primitive supports.

| File | Role |
|------|------|
| `documentation/MAAT_Foundations_of_Structural_Supports.pdf` | compiled foundational support note |
| `documentation/MAAT_Foundations_of_Structural_Supports.tex` | LaTeX source |

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

### Paper 22 — Numerical Validation of Structural Selection
**Numerical Validation of Structural Selection:**
*Evidence Beyond Energy-Based Ranking in a Minimal φ⁴ Model*

**Core idea:** First static lattice validation of the structural-selection
functional in a minimal 1+1 dimensional φ⁴ model. Tests whether
`F_struct` does more than reproduce ordinary energy ranking by
benchmarking vacua, the unstable homogeneous saddle, the discrete kink, and a
sector-matched distortion ensemble.

**Core results:**

| Benchmark | Result | Interpretation |
|----------|--------|----------------|
| Vacuum discrimination | **Δ_vac = 2.583852** | true vacua cleanly separated from unstable saddle |
| Kink discrimination | **p_K = 1.000** | discrete kink outranks all tested distortions |
| Ranking test (20%) | **0.010363 < 0.012128** | structural ranking outperforms energy ranking at the 20% percentile |

**Key finding:**
> Structural selection becomes numerically testable in a minimal field-theory
> setting. The structural layer is not merely a relabeling of on-shellness or
> static energy: it distinguishes stable from unstable exact solutions and adds
> measurable ranking power within a nontrivial topological sector.

**Scripts:**

| Script | Role |
|--------|------|
| `structural_selection_phi4_protocol.py` | runs the static lattice benchmark and generates `structural_selection_phi4_results.json` |
| `structural_selection_phi4_plots.py` | reads `structural_selection_phi4_results.json` and generates the validation figures |

**To reproduce the results:**

```bash
python3 structural_selection_phi4_protocol.py --json-output structural_selection_phi4_results.json --pretty-json
python3 structural_selection_phi4_plots.py --input structural_selection_phi4_results.json --output-dir structural_selection_phi4_plots
```

This produces:

| Output | Generated by |
|--------|--------------|
| `structural_selection_phi4_results.json` | `structural_selection_phi4_protocol.py` |
| `structural_selection_phi4_plots/` | `structural_selection_phi4_plots.py` |

**Documentation PDF:** `documentation/22_Numerical_Validation_of_Structural_Selection__Evidence_Beyond_Energy_Based_Ranking_in_a_Minimal_phi4_Model.pdf`

---

### Paper 26 — Structural Selection of Effective Constants
**Structural Selection of Effective Constants:**
*From MAAT Basins to MaxEnt-Weighted RG Bridge Tests*

**Core idea:** Extends structural selection from nonlinear field configurations
to effective constants. The paper tests whether MAAT-type structural scores
define low-defect basins in constant space and whether Standard-Model-like
RG flow preserves cross-sector structural constraints. The v12/v13 extension
calibrates the sector weights `lambda_a` through a maximum-entropy procedure.

**Scientific status:** This is a phenomenological benchmark, not a
first-principles derivation of the constants, not a precision Standard Model
fit, and not a solution of the cosmological-constant problem.

**Core results:**

| Benchmark | Result | Interpretation |
|----------|--------|----------------|
| Natural-constants basin | observed proxy and MAAT optimum lie in same low-defect basin | basin-level compatibility, not exact prediction |
| SM RG bridge | `F_bridge = 0.210446 < F_obs = 0.231789` | selected bridge point lies in same structural region as observed SM proxy |
| v11 holdout: α_em | predicted/observed = 1.581 | direct electromagnetic term can be removed while cross-sector constraints remain informative |
| v11 holdout: sin²θ_W | predicted/observed = 1.543 | weak-sector proxy remains order-of-magnitude constrained |
| v11 holdout: λ_H | predicted/observed = 0.170 | Higgs quartic sector is the main limitation of the current toy bridge |
| v12 MaxEnt weights | `R > V ≈ S > B > H` | equal sector weighting is replaced by calibrated effective weights |
| v13 MaxEnt bridge | stability = 0.9936 | MaxEnt-weighted run preserves the low-defect constants basin |

**Key finding:**
> Structural selection identifies a compatible basin of effective constants
> rather than a unique point. The v11 holdout test provides a first
> falsifiable cross-sector prediction protocol, while v12/v13 reduce the
> earlier equal-weight assumption by MaxEnt-calibrating effective sector
> weights.

**Reference data sources:**

| Source | Used for |
|--------|----------|
| [CODATA/NIST fundamental constants](https://www.nist.gov/programs-projects/codata-values-fundamental-physical-constants) | low-energy constants and dimensionless comparison values |
| [Particle Data Group 2024](https://pdg.lbl.gov/index-2024.html) | Standard-Model proxy values near the electroweak scale |

**Data attribution and license note:** PDG 2024 content is licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) except where
otherwise noted. NIST web information is generally public information unless
marked otherwise, but NIST Standard Reference Data may have separate copyright
and licensing conditions. The benchmark uses only reference comparison values;
users should cite CODATA/NIST and PDG when reusing or discussing the numerical
inputs.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/natural_constants_selection/` | v1--v13 natural-constants basin, robustness, ablation, gradient-flow, and MaxEnt-weight tests |
| `experiments/standard_model_bridge/` | one-loop SM-like RG bridge, four-panel summary figure, and v11 holdout test |

**Documentation PDF:** `documentation/26_Structural_Selection_of_Effective_Constants.pdf`

---

### Paper 27 — Boundary-Aware Calibration of MAAT Structural Weights
**Boundary-Aware Calibration of MAAT Structural Weights:**
*A Reproducible Benchmark for Constraint-Dominated Structural Selection*

**Core idea:** Calibrates MAAT structural weights on a fused dataset containing
SAT hardness instances, unconstrained MAAT-Core states, and explicit
constraint-boundary regimes. The key question is whether boundary information
changes the inferred selection hierarchy.

**Core results:**

| Quantity | Result |
|----------|--------|
| Fused samples | 3400 |
| Sources | SAT hardness atlas, MAAT-Core, MAAT-Core boundary |
| Lambda hierarchy | `R > V > H > S ~= B` |
| Dominant share | `R = 0.3917` |
| Fit loss | `0.0029543310` |

**Key finding:**
> Robustness / Respect becomes the dominant structural selector once explicit
> boundary, margin, and violation information are included.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/boundary_aware_lambda_calibration/` | fused defect dataset, closed lambda fit, plots, and result JSON |

**Documentation PDF:** `documentation/27_Boundary_Aware_Calibration_of_MAAT_Structural_Weights.pdf`

---

### Paper 28 — Cosmological Critical Coherence Index
**Cosmological Critical Coherence Index:**
*A Structural-Stress Observable for Cosmic Evolution*

**Core idea:** Defines a dimensionless cosmological CCI diagnostic combining
expansion stress, redshift activity, linear growth coherence, and structural
imbalance. The goal is to visualise the balance between expansion history and
structure-growth history, not to replace standard cosmological inference.

**Core results:**

| Quantity | Result |
|----------|--------|
| Reference cosmology | Planck-2018 flat LCDM |
| Chronometer points | 31 |
| Model range | `0 <= z <= 10` |
| Normalised CCI at `z=1` | `~6.66` |
| Normalised CCI at `z=2` | `~22.6` |
| Normalised CCI at `z=10` | `~803.8` |
| Chronometer residual RMS | `~1.01` |

**Key finding:**
> Cosmic evolution can be represented not only as expansion history, but as a
> structural-stress history comparing expansion/activity against growth
> coherence.

**Operational convention:** This first CCI-cosmology projection uses unit
weights, fixes `kappa = 1`, and does not fit regime cutoffs. Full MAAT
connectivity (`V`) and robustness (`R`) sectors are specified in the paper as
operational targets for future multi-probe work, not as measured sectors in
the chronometer-only projection.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/cosmological_cci/` | Cosmological CCI script, chronometer input table, generated CSVs, and plots |

**Data attribution and license note:** Planck-2018 parameter values and Cosmic
Chronometer measurements are external scientific data and should be cited to
the original publications/collaborations when reused or discussed. The CSV
tables and figures in this repository are derived analysis artifacts generated
for reproducibility of the CCI diagnostic. No endorsement by the Planck
Collaboration or the chronometer-data authors is implied.

**Documentation PDF:** `documentation/28_Cosmological_Critical_Coherence_Index.pdf`

---

### Paper 29 — Cosmological CCI with Growth Connectivity and Robustness
**Cosmological Critical Coherence Index with Growth Connectivity and Robustness Margins:**
*A Multi-Sector Structural-Stress Benchmark Using H(z) and f sigma_8(z)*

**Core idea:** Extends the chronometer-only cosmological CCI into a
multi-sector diagnostic by making the previously open connectivity (`V`) and
robustness (`R`) sectors operational. Cosmic Chronometer H(z) residuals define
expansion consistency, BOSS DR12 f sigma_8 residuals define growth
connectivity, and joint expansion-growth consistency defines robustness.

**Core results:**

| Quantity | Result |
|----------|--------|
| Chronometer points | 31 |
| BOSS DR12 f sigma_8 points | 3 |
| H(z) pull RMS | `0.693` |
| f sigma_8 pull RMS | `0.697` |
| Transition proxy | `z_c ~= 1.114` |
| Normalised v0.3 CCI at `z=1` | `~6.60` |
| Normalised v0.3 CCI at `z=2` | `~15.42` |
| Lambda hierarchy | `S > H > R > V` |

**Key finding:**
> The cosmological CCI can be extended from an expansion-growth scalar into a
> multi-sector observable with measured growth connectivity, robustness
> margins, companion MaxEnt weights, and a reproducible transition marker.

**Scientific status:** This is a diagnostic benchmark, not precision
cosmological inference. It does not fit cosmological parameters, replace
LCDM, or claim evidence for modified gravity.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/cosmological_cci_v03/` | v0.3 cosmological CCI script, H(z) and f sigma_8 input tables, generated CSV/JSON outputs, and plots |

**Data attribution and license note:** Paper 29 uses Planck-2018 reference
parameters, Cosmic Chronometer H(z) values, and BOSS DR12 consensus
f sigma_8 measurements from the cited literature. Repository CSV/PNG files are
derived reproducibility artifacts only. No endorsement by the Planck
Collaboration, SDSS/BOSS Collaboration, or the chronometer-data authors is
implied.

**Documentation PDF:** `documentation/29_Cosmological_CCI_Growth_Connectivity_and_Robustness.pdf`

---

### Paper 30 — Response-Based Derivation of MAAT Structural Weights
**Response-Based Derivation of MAAT Structural Weights:**
*From Covariance Geometry to Selection Pressure*

**Core idea:** Provides the missing response-theoretic interpretation of
`lambda_a`.  Instead of treating the weights as arbitrary fitted constants,
the addendum derives them as linear-response coefficients of the empirical
defect ensemble:

```text
lambda = (Cov_mu0[d] + eta tr(C)/A I)^(-1)(<d>_mu0 - <d>_target)
```

Here `mu0` is the empirical reference measure over the fused defect ensemble,
`Cov_mu0[d]` is the defect covariance matrix, and `<d>_target` is the selected
target sub-ensemble.

**Core results:**

| Target ensemble | Main hierarchy | Interpretation |
|----------------|----------------|----------------|
| Low-defect 20% | `R > B > S > H > V` | generic low-defect selection gives dominant but not exclusive robustness |
| Safe + core-safe | `R >> H > B > V` | explicit safety targets are overwhelmingly robustness-driven |
| Safe boundary only | `R > S > B > V > H` | boundary-safe states reproduce a balanced but R-dominant hierarchy |
| Not violated | `R > S` | excluding violations activates mainly robustness and activity |

**Key finding:**
> `R` dominance emerges as a covariance-response effect of safety/boundary
> target selection.  The weights are still effective and ensemble-dependent,
> but they are no longer merely hand-chosen fit parameters.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/lambda_response_closure/` | response-theoretic lambda derivation, result JSON, and comparison plots |

**Documentation PDF:** `documentation/30_Response_Theoretic_Closure_of_MAAT_Structural_Weights.pdf`

---

### Paper 31 — Dynamic Structural Selection
**Dynamic Structural Selection:**
*From Response Fields to Observable Cosmology in the MAAT Framework*

**Core idea:** Consolidates the MAAT versioned development from static
structural ranking into a dynamical, local, gravitational, and observable
effective-theory pipeline:

```text
d_a -> C_ab -> lambda_a(t) -> lambda_a(x,t) -> T_MAAT -> FLRW -> observables
```

The paper combines response-field dynamics, local selection-pressure fields,
effective gravitational coupling, a scalar worked example, an FLRW stability
scan, and an observable-projection layer.

**Core results:**

| Layer | Result | Interpretation |
|-------|--------|----------------|
| v0.5 response flow | final tracking error `9.91e-7` | lambda follows covariance-response fixed point |
| v0.6 local fields | final residual ratio `0.1648` | local fields reduce equation residuals with energy unchanged |
| v0.9 FLRW scan | `666/784` stable, `618/784` background-safe | broad stable region in the toy scalar scan |
| v0.10 observables | max `|Delta H/H| = 0.0399` | small structured expansion deviation |
| v0.10 growth proxy | max `|Delta f sigma_8/f sigma_8| = 0.0336` | first approximate growth-sector signal |

**Key finding:**
> Structural selection can be represented as a reproducible effective-theory
> ladder from defect diagnostics to observable cosmological signatures, while
> remaining explicitly toy-level and not a completed cosmological model.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/maat_dynamic_fields_v05_v09/` | v0.5 response flow, v0.6 local field benchmark, v0.9 FLRW stability scan |
| `experiments/maat_observable_predictions_v10/` | observable projection used in Paper 31 |

**Documentation PDF:** `documentation/31_Dynamic_Structural_Selection_From_Response_Fields_to_Observable_Cosmology.pdf`

---

### Paper 32 — First H(z) Data Comparison
**First H(z) Data Comparison of MAAT Structural-Selection Cosmology:**
*Cosmic Chronometer Diagnostic and Chi-Square Scan*

**Core idea:** Moves the MAAT cosmology stack from internal observable
projection to a first direct comparison with observational expansion data.
The paper uses Cosmic Chronometer H(z) measurements as a diagnostic
goodness-of-fit test for the fixed v0.10 branch and for a two-parameter
MAAT scalar scan.

**Core results:**

| Quantity | Result |
|----------|--------|
| Chronometer points | `31` |
| LCDM reference | `chi2 = 14.8759`, `chi2/N = 0.4799` |
| Fixed v0.10 MAAT | `chi2 = 15.3661`, `chi2/N = 0.4957` |
| Best scanned MAAT branch | `chi2 = 15.3674`, `Delta chi2 = +0.4915` vs LCDM |
| Best scan reduced chi-square | `chi2_nu = 0.5299` for `N-k = 29` |
| Stable scan points | `581/616` |
| Best branch maximum `Omega_MAAT` | `0.02668` |
| Best branch maximum `|Delta H/H|` | `0.02347` |

**Key finding:**
> In this first H(z) diagnostic, MAAT remains stable, subdominant, and close
> to LCDM, but it is not favoured over LCDM by the Cosmic Chronometer chi-square.
> The value of the result is therefore falsifiability and pipeline closure, not
> a detection claim.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/maat_hz_chi2_paper32/` | fixed v0.10 chronometer comparison and two-parameter H(z) chi-square scan |

**Documentation PDF:** `documentation/32_First_Hz_Data_Comparison_of_MAAT_Structural_Selection_Cosmology.pdf`

---

### Paper 33 — CCI as Projection Observable
**The Critical Coherence Index as a Projection Observable:**
*Breadth--Depth Compression and Transition Robustness in a Minimal Cosmological Benchmark*

**Core idea:** Clarifies what the CCI represents in the cosmological setting:
not an energy density and not a fitted cosmological parameter, but a projection
observable measuring stress between expansion-driven breadth and accessible
structural depth.

**Core results:**

| Quantity | Result |
|----------|--------|
| v0.3 reference background | Planck-like flat LCDM |
| v0.3 transition estimate | `z_tr = 0.84985` |
| v0.3 `C_norm(z_tr)` | `12.0466` |
| v0.3 final `C_norm(z=3)` | `474.184` |
| v0.4 sensitivity points | `6930` |
| v0.4 mean transition redshift | `0.84491` |
| v0.4 median transition redshift | `0.88889` |
| v0.4 fraction in `0.5 <= z_tr <= 1.1` | `0.55657` |

**Key finding:**
> The CCI can be operationally interpreted as a stabilized projection stress:
> it grows when expansion breadth outpaces accessible coherent depth. The
> resulting transition marker is reproducible across a broad proxy-parameter
> grid, but it is not a measured cosmological transition or a data-fit result.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/maat_cci_projection_paper33/` | v0.2 stabilization predecessor, v0.3 LCDM-growth projection benchmark, and v0.4 sensitivity scan |

**Documentation PDF:** `documentation/33_CCI_as_Projection_Observable.pdf`

---

### Paper 34 — Projection Observable vs Growth and Expansion Data
**MAAT Projection Observable vs Growth and Expansion Data:**
*Transition Marker and Parameter Sensitivity in v0.11*

**Core idea:** Tests whether the MAAT projection observable `C_proj` collapses
onto standard growth data. It does not: the large mismatch against
`f sigma_8` is interpreted as evidence that `C_proj` is a projection-level
diagnostic, not a matter-growth observable.

**Core results:**

| Quantity | Result |
|----------|--------|
| v0.11 baseline transition marker | `z_tr = 1.04356` |
| `C_norm(z_tr)` | `22.64463` |
| `R_proj(z_tr)` | `1.11151e-03` |
| `C_norm(z=3)` | `647.49295` |
| LCDM `f sigma_8` chi-square/dof | `1.04813` |
| Rescaled `C_proj` vs `f sigma_8` chi-square/dof | `42.76835` |
| LCDM `H(z)` chi-square/dof | `0.47944` |
| Gamma-alpha scan points | `625` |
| Scan median transition | `1.18715` |
| Scan transition range | `[0.28073, 1.58353]` |

**Key finding:**
> Paper 34 demonstrates that the MAAT projection observable does not reduce
> to known growth observables, establishing it as a distinct layer of
> cosmological information. The result is diagnostic, not a cosmological fit
> or evidence claim.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/maat_projection_growth_paper34/` | v0.11 projection residual, growth-data mismatch diagnostic, H(z) reference, and gamma-alpha transition scan |

**Documentation PDF:** `documentation/34_MAAT_Projection_Observable_vs_Growth_and_Expansion_Data.pdf`

---

### Paper 35 — Linear Growth Embedding
**Linear Growth Embedding of MAAT Structural Selection:**
*Perturbation Stability, Positivity, and Sub-Percent Growth Signatures*

**Core idea:** Moves beyond the direct `C_proj` versus `f sigma_8` template
comparison of Paper 34 by embedding the representative MAAT branch into a
minimal linear growth equation. The goal is to test whether MAAT produces
small and stable growth-sector effects when treated as a forward model.

**Core results:**

| Quantity | Result |
|----------|--------|
| `Omega_MAAT,0` | `0.00331` |
| `w_MAAT` | `-0.801` |
| Max `|Delta D/D|` for `z < 2` | `0.0591 %` |
| Max `|Delta f sigma_8/f sigma_8|` for `z < 2` | `0.4570 %` |
| NEC margin | `+0.199 rho` |
| Positivity stress tests | `3/3 positive` |
| Fixed-point mode stability | all tested `gamma_k > 0` |

**Key finding:**
> A minimal MAAT growth embedding produces sub-percent growth-sector
> deviations while preserving linear stability, positivity, and acceptable
> energy-condition behaviour for the tested branch. This is a first
> perturbative toy embedding, not a precision cosmological constraint.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/maat_growth_perturbation_paper35/` | corrected LCDM reference, MAAT growth equation, selection-pressure perturbations, positivity stress test, and energy-condition check |

**Documentation PDF:** `documentation/35_Linear_Growth_Embedding_of_MAAT_Structural_Selection.pdf`

---

### Paper 36 — The Maturation of MAAT
**The Maturation of MAAT:**
*From Projection Observables (v0.11) to Structural Selection (v1.0) to Emergent Robustness Closure (v1.2.1)*

**Core idea:** Reconstructs the conceptual evolution of the MAAT framework
from projection diagnostics to structural selection and finally to the
v1.2.1 robustness closure. The paper resolves the ambiguity of the old
independent `R` sector by replacing it with derived support-level quantities:
`R_resp` for structural respect and `R_rob` for emergent robustness.

**Core results:**

| Quantity | Result |
|----------|--------|
| Primary support sectors | `H, B, S, V` |
| Old v1.0 input | independent `R` |
| v1.2.1 respect closure | `R_resp = (H B V)^(1/3)` |
| v1.2.1 robustness closure | `R_rob = min(R_resp, (H B S V)^(1/4))` |
| Stability convention | `Stability = R_rob` |
| Closure coefficient | `lambda_cl`, covariance/response-derived |
| Primary-field count | reduced from `5` to `4` |

**Key finding:**
> MAAT v1.2.1 does not add another sector. It reduces the framework by
> reconstructing robustness as a derived closure quantity on normalised
> supports. The result is a more closed and interpretable selection
> framework with the same baseline observable pipeline.

**Scientific status:** This is a conceptual closure and synthesis paper.
It does not introduce a new numerical benchmark. Its role is to make the
notation, support-level interpretation, and R-closure logic internally
consistent across Papers 30--35 and the Paper 37 v1.2.1 benchmark.

**Documentation PDF:** `documentation/36_The_Maturation_of_MAAT.pdf`

---

### Paper 37 — MAAT v1.2.1 Observable Proxy Predictions
**MAAT v1.2.1 Observable Proxy Predictions and Stability Landscape:**
*Structural Fields, Projection Layer, CCI Diagnostics, and Parameter Sensitivity*

**Core idea:** Provides the numerical companion benchmark for Paper 36.
The paper applies the v1.2.1 closure convention with four primary support
sectors `H, B, S, V`, derived structural respect `R_resp=(H B V)^(1/3)`,
and emergent robustness `R_rob=min(R_resp,(H B S V)^(1/4))`. It then tests
the baseline observable proxy and scans a two-parameter projection/damping
landscape.

**Core results:**

| Quantity | Result |
|----------|--------|
| Trajectory points | `382` |
| Redshift range | `0 <= z <= 2.33` |
| Baseline max `|Delta H/H|` | `1.70%` |
| Baseline max `Omega_MAAT` | `3.31%` |
| Baseline max `|Delta f sigma_8 / f sigma_8|` | `1.83%` |
| Mean structural respect | `<R_resp> = 0.818` |
| Mean emergent robustness | `<R_rob> = 0.805` |
| Mean CCI v1.2.1 | `0.208` |
| Stability scan | `1089/1089` internally stable points |
| Best scan point | `eta = 0.000`, `gamma = 0.050` |
| SAT companion correlations | `rho(V, runtime)=0.647`, `rho(R_rob, runtime)=0.484` |

**Key finding:**
> The v1.2.1 closure does not destabilise the observable proxy pipeline.
> Emergent robustness remains finite across the scanned parameter space,
> and the baseline is not an isolated fine-tuned point within the tested
> two-parameter landscape.

**Scientific status:** This is a toy/proxy benchmark, not an observational
cosmology fit. Growth deviations are simplified proxy quantities, not
Boltzmann-code predictions.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/maat_v121_observables_stability_paper37/` | Paper 37 observable proxy run, stability landscape scan, SAT companion validation, output CSV/JSON files, and figures |

**Documentation PDF:** `documentation/37_MAAT_v121_Observable_Proxy_Predictions_and_Stability_Landscape.pdf`

---

### Paper 38 — MAAT v1.2.1 Robustness Closure in Linear Growth
**MAAT v1.2.1 Robustness Closure in Linear Growth:**
*Sub-Percent Growth Signatures, Emergent Robustness, and Selection-Field Stability*

**Core idea:** Applies the v1.2.1 closure convention to a minimal
linear-growth benchmark. The paper updates the Paper 35 growth pipeline with
the Paper 36/Paper 37 closure logic:
`R_resp=(H B V)^(1/3)`, `R_rob=min(R_resp,(H B S V)^(1/4))`,
and `Stability=R_rob`. It also checks selection-field perturbations,
positivity of the lambda-relaxation equation, fixed-point decay, and basic
energy-condition diagnostics.

**Core results:**

| Quantity | Result |
|----------|--------|
| Max `|Delta D/D|` for `z < 2` | `0.0591%` |
| Max `|Delta f sigma_8 / f sigma_8|` for `z < 2` | `0.4570%` |
| Mean `R_resp` | `0.9543` |
| Mean `R_rob` | `0.9313` |
| Mean `CCI_min` | `0.3114` |
| Mean `CCI_diag` | `0.2348` |
| Selection-field positivity | passed all tested scenarios |
| Fixed-point perturbations | stable for all tested modes |

**Key finding:**
> The v1.2.1 robustness closure can be inserted into a growth-sensitive
> benchmark without producing large observable deviations or obvious
> selection-field instabilities. The result is an internal consistency check,
> not an observational detection.

**Scientific status:** This is an effective-theory growth benchmark, not a
Boltzmann-code calculation and not a cosmological likelihood.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/maat_paper38_v121_robustness_closure/` | Paper 38 linear-growth closure benchmark, selection-field perturbation test, output CSV/JSON files, and figures |

**Documentation PDF:** `documentation/38_MAAT_v121_Robustness_Closure_in_Linear_Growth.pdf`

---

### Paper 39 — MAAT v1.2.1 Observable Growth Signature Proxy
**MAAT v1.2.1 Observable Growth Signature Proxy:**
*Projection-Modulated f sigma_8, Emergent Robustness, and a Boundary-Limited Diagnostic Scan*

**Core idea:** Tests whether the MAAT projection observable can be inserted
as a small explicit signature template in the growth observable `f sigma_8(z)`
while preserving the v1.2.1 closure convention:
`R_resp=(H B V)^(1/3)`, `R_rob=min(R_resp,(H B S V)^(1/4))`,
and `Stability=R_rob`.

**Core results:**

| Quantity | Result |
|----------|--------|
| Growth comparison points | `13` |
| Projection transition estimate | `z_tr ~= 0.6508` |
| Epsilon scan range | `[-0.01, 0.01]` |
| Best epsilon | `-0.0100` |
| LCDM chi2 | `12.4373` |
| MAAT proxy chi2 | `12.3772` |
| Delta chi2 | `-0.0601` |
| Max `|Delta f sigma_8 / f sigma_8|` | `0.9891%` |
| Mean `|Delta f sigma_8 / f sigma_8|` | `0.5725%` |
| Mean `R_resp` | `0.7247` |
| Mean `R_rob` | `0.6673` |
| Mean `CCI_min` | `0.2662` |
| Mean `CCI_diag` | `0.2036` |

**Key finding:**
> A bounded MAAT projection template can modulate `f sigma_8(z)` at the
> sub-percent level while preserving the v1.2.1 robustness closure. The best
> epsilon value lies at the scan boundary, so the result is a diagnostic
> compatibility check, not a measured parameter or evidence for modified
> growth.

**Scientific status:** This is an observable-signature proxy, not a full
Boltzmann-code calculation and not a precision cosmological likelihood.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/maat_paper39_observable_growth_signature/` | Paper 39 projection-modulated `f sigma_8` proxy, epsilon chi2 scan, v1.2.1 closure diagnostics, output CSV/JSON files, and figures |

**Documentation PDF:** `documentation/39_MAAT_v121_Observable_Growth_Signature_Proxy.pdf`

---

### Paper 40 — MAAT v1.2.1 Structural Signature Test in Growth Data
**MAAT v1.2.1 Structural Signature Test in Growth Data:**
*Projection CCI, Residual Magnitudes, and Exploratory Null Tests*

**Core idea:** Tests whether MAAT v1.2.1 structural diagnostics align with
the residual structure of a compact `f sigma_8(z)` growth comparison set
relative to a Planck-normalised LCDM baseline. Unlike Paper 39, this paper
does not ask whether MAAT improves the fit. It asks whether structural
diagnostics track where the reference model shows larger residual stress.

**Core results:**

| Quantity | Result |
|----------|--------|
| Growth comparison points | `13` |
| Projection transition estimate | `z_tr ~= 0.6508` |
| Spearman `R_proj` vs signed residual | `0.6319`, `p = 0.0228` |
| Spearman `V` vs signed residual | `-0.6319`, `p = 0.0228` |
| Spearman `CCI_diag` vs `|residual_sigma|` | `0.5934`, `p = 0.0338` |
| Random-field null for `CCI_diag` | `p = 0.0363` |
| Redshift-shuffle null for `CCI_diag` | `p = 0.0359` |
| Spearman `B` vs `|residual_sigma|` | `-0.7912`, `p = 0.0017` |

**Key finding:**
> The diagnostic CCI tracks residual magnitude more clearly than residual
> direction in this compact growth-data benchmark. Because the balance support
> `B` is residual-sensitive, the result is a semi-supervised structural
> consistency test, not a blind prediction or detection claim.

**Scientific status:** This is an exploratory residual-structure diagnostic,
not a Boltzmann-code calculation and not a full cosmological likelihood.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/maat_paper40_structural_signature_test/` | Paper 40 CCI/residual structural signature test, permutation null tests, output CSV/JSON files, and figures |

**Documentation PDF:** `documentation/40_MAAT_v121_Structural_Signature_Test_in_Growth_Data.pdf`

---

### Paper 41 — MAAT v1.2.1 Variable Closure and Measurement Definitions
**MAAT v1.2.1 Variable Closure and Measurement Definitions:**
*Formal Definition of Structural Supports, Activity, Response-Derived Closure,
and Projection Observable Derivation*

**Core idea:** Provides the formal variable-definition layer for the
v1.2.1 framework used in Papers 33--40. It defines the four primary support
sectors `H`, `B`, `S`, `V` from non-negative defects, keeps robustness
emergent through `R_resp=(H B V)^(1/3)` and
`R_rob=min(R_resp,(H B S V)^(1/4))`, and makes clear that `lambda_cl` is a
derived closure coefficient rather than a fifth primary sector weight.

**Core definitions:**

| Quantity | Definition / role |
|----------|-------------------|
| Support scores | `Gamma_a = 1 / (1 + d_a)` |
| Primary sectors | `H`, `B`, `S`, `V` only |
| Implicit respect | `R_resp = (H B V)^(1/3)` |
| Emergent robustness | `R_rob = min(R_resp,(H B S V)^(1/4))` |
| Stability | `Stability = R_rob` |
| Closure coefficient | `lambda_cl = sum_a w_a lambda_a`, sensitivity-projected from primary weights |
| CCI diagnostic | `CCI_diag = S / (H + B + V + R_rob + epsilon)` |

**Key finding:**
> MAAT v1.2.1 is operationally closed at the variable-definition level:
> no independent fifth robustness sector is introduced, `R_resp` is not
> double-counted in the CCI denominator, and projection parameters are fixed
> by the stated response-map convention rather than fitted as independent
> observables.

**Scientific status:** This is a formal closure and measurement-definition
note. It introduces no new numerical data set, likelihood analysis, or
observational fit.

**Reproducibility:** Paper 41 has no standalone code folder. It documents the
variable conventions used by the existing reproducibility folders:
`experiments/lambda_response_closure/`,
`experiments/maat_cci_projection_paper33/`,
`experiments/maat_paper39_observable_growth_signature/`, and
`experiments/maat_paper40_structural_signature_test/`.

**Documentation PDF:** `documentation/41_MAAT_v121_Variable_Closure_and_Measurement_Definitions.pdf`

---

### Paper 42 — Blind Projection Test
**Blind Projection Test:**
*Response-Derived Projection, No Projection-Shape Tuning, and Residual Structure Diagnostics*

**Core idea:** Removes the projection-shape tuning used in earlier residual
diagnostics. The observable parameters are derived from response-closed
MAAT v1.2.1 sector weights:

```text
defects -> covariance C -> lambda -> response shares pi_a
        -> gamma_lambda, Bstar_lambda, alpha_lambda
        -> C_proj_lambda(z)
```

The test uses a compact `f sigma_8(z)` growth comparison set and asks whether
the response-derived projection observable and CCI diagnostics retain a
residual-structure alignment without fitting `epsilon`, `gamma`, `Bstar`, or
`alpha`. The small denominator regulator is fixed numerically and is not
treated as a physical parameter.

**Core results:**

| Quantity | Result |
|----------|--------|
| Growth comparison points | `13` |
| Response shares | `pi_H=0.0960`, `pi_B=0.1257`, `pi_S=0.3879`, `pi_V=0.3903` |
| Derived `gamma_lambda` | `1.3305` |
| Derived `Bstar_lambda` | `4.5091` |
| Derived `alpha_lambda` | `1.9937` |
| Transition marker | `z_tr = 0.6986` |
| Spearman `CCI_diag` vs `|residual_sigma|` | `0.4286`, permutation `p = 0.1456` |
| Spearman `C_proj_lambda` vs `|residual_sigma|` | `0.4286`, permutation `p = 0.1432` |

**Key finding:**
> The residual-structure direction survives blind projection-shape closure,
> but the signal is reduced to a positive, non-significant exploratory trend.
> Paper 42 is therefore a robustness stress test, not a detection claim.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/maat_paper42_blind_projection_test/` | Paper 42 blind projection test, response-derived parameters, residual/CCI correlations, permutation null tests, output CSV/JSON files, and figures |

**Data attribution and license note:** The Planck-normalised reference
parameters and compact `f sigma_8` comparison points are external scientific
data and should be cited to the original publications/collaborations. The
repository CSV/PNG files are derived reproducibility artifacts only. No
endorsement by the Planck Collaboration, survey collaborations, or original
data authors is implied.

**Documentation PDF:** `documentation/42_Blind_Projection_Test.pdf`

---

### Paper 43 — Linear Perturbations and Structure Growth
**Linear Perturbations and Structure Growth:**
*From Projection Signatures to Effective Growth Dynamics in Structural-Selection Cosmology*

**Core idea:** Promotes the MAAT projection layer from a direct growth
signature template to an explicit effective linear matter-growth equation.
The bounded projection kernel `C_hat_proj(z)` modifies the sub-horizon growth
source through

```text
mu(z) = G_eff/G = 1 + eta * C_hat_proj(z)
```

with `C_hat_proj` normalised over `0 <= z <= 3`. This makes `eta` the maximum
fractional growth-coupling deformation in the benchmark interval.

**Core results:**

| Quantity | Result |
|----------|--------|
| Growth comparison points | `13` |
| Eta scan | `eta in [0, 0.08]`, `41` grid points |
| Stable branches | `41/41` |
| Best eta | `0.0000` |
| LCDM chi2 | `12.2021` |
| Representative eta | `0.0200` |
| Max `|Delta D/D|` for `z <= 3` at eta `0.02` | `0.7656%` |
| Max `|Delta f sigma_8 / f sigma_8|` for `z <= 3` at eta `0.02` | `0.5441%` |
| Max `|mu - 1|` for `z <= 3` at eta `0.02` | `1.9966%` |
| Representative eta chi2 | `12.5737` |

**Key finding:**
> The compact growth comparison does not prefer a nonzero MAAT coupling:
> the best value is the LCDM limit `eta=0`. However, nonzero perturbative
> branches remain stable and produce only small, structured deviations. Paper
> 43 therefore establishes a controlled bridge from projection observables to
> explicit linear growth dynamics.

**Scientific status:** This is an effective sub-horizon linear-growth
benchmark, not a Boltzmann-code calculation, not a CMB or weak-lensing
likelihood, and not evidence for modified gravity.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/maat_linear_perturbations_paper43/` | Paper 43 linear perturbation benchmark, eta scan, compact `f sigma_8` comparison, output CSV/JSON files, and figures |

**Data attribution and license note:** The Planck-normalised reference
parameters and compact `f sigma_8` comparison points are external scientific
data and should be cited to the original publications/collaborations. The
repository CSV/PNG files are derived reproducibility artifacts only. No
endorsement by the Planck Collaboration, survey collaborations, or original
data authors is implied.

**Documentation PDF:** `documentation/43_Linear_Perturbations_and_Structure_Growth.pdf`

---

### Paper 44 — Active Structural Significance
**Active Structural Significance:**
*A Toy Extension of the MAAT v1.2.1 Robustness Closure*  
MAAT Structural Selection Series --- Paper 44 / Diagnostic Note

**Core idea:** Extends the v1.2.1 interpretation without returning to a fifth
primitive `R` sector. The paper treats `R` as an emergent hierarchy over the
four primary support sectors `H, B, S, V`:

```text
R_resp = (H B V)^(1/3)
R_rob  = min(R_resp, (H B S_eff V)^(1/4))
R_sig  = R_resp^(1-alpha) S_eff^alpha
```

Here `S_eff` is controlled activity in an optimal window, not raw activity:

```text
S_eff(A) = exp[-0.5 ((A - A_star) / sigma_A)^2]
```

**Core results:**

| Quantity | Result |
|----------|--------|
| Ensemble size | `6000` |
| Activity optimum | `A_star = 1.0` |
| Activity width | `sigma_A = 0.28` |
| Significance exponent | `alpha = 0.45` |
| Max `R_sig` location | `A = 0.9977` |
| Max `R_sig` | `0.9918` |
| Dormant coherent `R_resp` / `R_sig` | `0.6361 / 0.0980` |
| Stable-star `R_resp` / `R_sig` | `0.9850 / 0.9917` |
| Chaotic-burst `R_resp` / `R_sig` | `0.4451 / 0.0806` |
| High `R_resp` but low `R_sig` fraction | `13.07%` |

**Key finding:**
> MAAT v1.2.1 should be read as a robustness closure, not a complete theory of
> structural significance. Active significance peaks only when passive
> adherence is combined with controlled activity in the optimal window.

**Scientific status:** This is a diagnostic note and toy conceptual
simulation, not a physical stellar model, not a cosmological fit, and not a
new primitive-field proposal. The use of `S_eff` in the robustness boundary is
part of this toy extension and does not retroactively modify the v1.2.1 papers.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/active_respect_significance/` | Paper 44 active significance toy simulation, archetype comparison, ensemble CSV/JSON outputs, and figures |

**Documentation PDF:** `documentation/44_Active_Structural_Significance.pdf`

---

### Paper 45 — Effective Metric Response in MAAT Structural Cosmology
**Effective Metric Response in MAAT Structural Cosmology:**
*Gravitational Slip, Weyl-Potential Proxies, and Growth-Lensing Consistency*

**Core idea:** Extends the Paper-43 linear-growth benchmark from a
growth-only response channel to a minimal metric-response benchmark. The
bounded MAAT projection kernel `C_hat_proj(z)` now sources both the effective
Newtonian growth coupling and a gravitational-slip channel:

```text
mu(z)       = G_eff/G = 1 + eta_g * C_hat_proj(z)
eta_slip(z)= Phi/Psi = 1 + beta_s * C_hat_proj(z)
Sigma(z)   = mu(z) * [1 + eta_slip(z)] / 2
```

Here `Sigma` parameterizes the Weyl/lensing response channel. A nonzero
`eta_slip` is equivalent, at the effective perturbation level, to introducing
a metric anisotropic-stress channel.

Two diagnostic consistency proxies are reported:

```text
Weyl_proxy(z) = Sigma(z) * D_MAAT(z) / D_LCDM(z)
EG_proxy(z)   = Sigma(z) * f_LCDM(z) / f_MAAT(z)
```

**Core results:**

| Quantity | Result |
|----------|--------|
| `eta_g` scan | `[0.00, 0.08]`, `41` points |
| `beta_s` scan | `[-0.06, 0.06]`, `61` points |
| Total metric-response branches | `2501` |
| Stable / positive branches | `2501 / 2501` |
| Growth-only best `eta_g` | `0.0000` |
| Growth-only `beta_s` | unconstrained without lensing data |
| Representative branch | `eta_g = 0.02`, `beta_s = 0.03` |
| Representative max `|mu - 1|` | `1.9966%` |
| Representative max `|eta_slip - 1|` | `2.9949%` |
| Representative max `|Sigma - 1|` | `3.5240%` |
| Representative Weyl-proxy deviation | `2.7314%` |
| Representative `E_G`-proxy deviation | `2.3019%` |

**Key finding:**
> MAAT projection structure can be coupled to matter growth and metric
> potentials as two distinct bounded response channels. Growth-only data still
> prefer the LCDM limit for `eta_g` and cannot constrain `beta_s`, because
> slip changes the metric potentials rather than the matter-growth source term;
> lensing or metric-potential data are required for the next step.

**Scientific status:** This is an effective metric-response benchmark, not a
Boltzmann-code implementation, not a weak-lensing likelihood, not a CMB
anisotropy calculation, and not evidence for modified gravity.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/maat_metric_response_paper45/` | Paper 45 metric-response benchmark, gravitational-slip scan, Weyl/`E_G` proxy outputs, CSV/JSON files, and figures |

**Data attribution and license note:** The Planck-normalised reference
parameters and compact `f sigma_8` comparison points are external scientific
data and should be cited to the original publications/collaborations. The
repository CSV/PNG files are derived reproducibility artifacts only. No
endorsement by the Planck Collaboration, survey collaborations, or original
data authors is implied.

**Documentation PDF:** `documentation/45_Effective_Metric_Response_in_MAAT_Structural_Cosmology.pdf`

---

### Paper 46 — Environmental Screening in MAAT Structural Cosmology
**Environmental Screening in MAAT Structural Cosmology:**
*Emergent Suppression of Metric Response in High-Coherence Regions*

**Core idea:** Extends Paper 45 by making the bounded MAAT projection response
environment-dependent. The same projection kernel that sources growth and
metric response is multiplied by an environmental suppression factor:

```text
C_env(z, Delta, Sigma_env) = C_hat_proj(z) * S_env(Delta, Sigma_env)

S_env = [1 + alpha_rho Delta_+^n + alpha_sigma Sigma_env^m]^(-1)
Delta_+ = max(Delta - 1, 0)
```

The screened response channels are:

```text
mu_env(z)        = 1 + eta_g  C_env(z)
eta_slip_env(z)  = 1 + beta_s C_env(z)
Sigma_lens_env   = mu_env(z) * [1 + eta_slip_env(z)] / 2
```

Because `0 <= C_hat_proj <= 1` and `0 <= S_env <= 1`, the environmental
response is bounded by construction. The GR-recovery limit is explicit:

```text
S_env -> 0  =>  mu_env, eta_slip_env, Sigma_lens_env -> 1
```

**Core results:**

| Quantity | Result |
|----------|--------|
| Representative response | `eta_g = 0.02`, `beta_s = 0.03` |
| Screening parameters | `alpha_rho = 0.15`, `alpha_sigma = 1.0`, `n = 0.75`, `m = 2.0` |
| Synthetic environments | `void`, `sheet`, `field`, `filament`, `cluster`, `local_dense` |
| Stable environments | `6 / 6` |
| `S_env(void)` | `0.9975` |
| `S_env(cluster)` | `0.1555` |
| `S_env(local_dense)` | `0.0002107` |
| Void max `|Sigma_lens - 1|` | `3.5151%` |
| Cluster max `|Sigma_lens - 1|` | `0.5441%` |
| Local-dense max `|Sigma_lens - 1|` | `0.000736%` |
| Screening transition at `Sigma_env=0.2` | `S_env ~= 0.5` at `Delta ~= 12.35` |

**Key finding:**
> MAAT metric response can remain percent-level in void-like regions while
> being strongly suppressed in dense or highly organized environments. The
> construction provides an effective structural screening layer, not a
> microscopic chameleon, Vainshtein, or symmetron derivation.

**Scientific status:** This is an effective environmental-screening benchmark,
not a local-gravity test, not an N-body or halo-model analysis, and not a
first-principles derivation of a screening mechanism.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/maat_environmental_screening_paper46/` | Paper 46 environmental screening benchmark, environment archetypes, phase-space grid, CSV/JSON outputs, and figures |

**Data attribution and license note:** The Planck-normalised reference
parameters and compact `f sigma_8` comparison points are external scientific
data and should be cited to the original publications/collaborations. The
repository CSV/PNG files are derived reproducibility artifacts only. No
endorsement by the Planck Collaboration, survey collaborations, or original
data authors is implied.

**Documentation PDF:** `documentation/46_Environmental_Screening_in_MAAT_Structural_Cosmology.pdf`

---

### Paper 47 — From Environmental Screening to Observable Void Signatures
**From Environmental Screening to Observable Void Signatures:**
*A Halo-Environment Interpretation of MAAT Metric Response*

**Core idea:** Addresses the main limitation of Paper 46 by translating the
environmental screening axis into large-scale-structure language. Instead of
using only an abstract `Sigma_env`, Paper 47 uses overdensity
`Delta = rho / rho_bar` and a bounded tidal/shear proxy:

```text
Sigma_env(q) = q^2 / (q^2 + q_star^2)
```

The screened kernel becomes:

```text
C_env(z, Delta, q) = C_hat_proj(z) * S_env(Delta, q)

S_env(Delta, q)
  = [1 + alpha_rho Delta_+^n + alpha_sigma Sigma_env(q)^m]^(-1)

Delta_+ = max(Delta - 1, 0)
```

The observable target is the lensing-to-growth response ratio:

```text
R_LG(z) = Sigma_lens_env(z) / mu_env(z)
        = [1 + eta_slip_env(z)] / 2
```

**Core results:**

| Quantity | Result |
|----------|--------|
| Synthetic populations | `void`, `field`, `filament`, `cluster`, `local_dense` |
| Samples per population | `6000` |
| Random seed | `47` |
| `mean S_env(void)` | `0.998821` |
| `mean S_env(cluster)` | `0.160656` |
| `mean S_env(local_dense)` | `0.000213` |
| `mean max R_LG - 1` in voids | `1.498231%` |
| `mean max R_LG - 1` in clusters | `0.240984%` |
| `mean max R_LG - 1` in local-dense proxy | `0.000320%` |
| Void--cluster `R_LG` contrast | `1.257248 percentage points` |
| Void--local `R_LG` contrast | `1.497912 percentage points` |
| Minimum stable fraction | `1.0` |
| Screening sensitivity cases | `11` |
| Positive void--cluster contrast under sensitivity | `11 / 11` |
| Void--cluster contrast range under sensitivity | `0.571843` to `1.400886` percentage points |

**Key finding:**
> If the MAAT metric response is environmentally screened, the residual
> lensing-to-growth response should be largest in void-like environments,
> strongly suppressed in clusters, and essentially absent in high-density
> local-test proxies.

**Scientific status:** This is a phenomenological bridge benchmark. It is not
a full halo model, not a weak-lensing likelihood, not a Boltzmann-code
calculation, and not a microscopic derivation of a screening mechanism.

**Scripts and reproducibility:**

Repository URL:

```text
https://github.com/Chris4081/structural-selection-principle/tree/main/experiments/maat_void_environment_paper47
```

| Folder | Role |
|--------|------|
| `experiments/maat_void_environment_paper47/` | Paper 47 void/cluster environment benchmark, Monte-Carlo samples, tidal phase-space grid, CSV/JSON outputs, and figures |

**Data attribution and license note:** The benchmark uses Planck-normalised
reference cosmological parameters as external literature values and synthetic
environment populations generated by the script. No external survey catalogue
is redistributed. Repository CSV/PNG files are derived reproducibility
artifacts only. No endorsement by the Planck Collaboration, survey
collaborations, or cited authors is implied.

**Documentation PDF:** `documentation/47_From_Environmental_Screening_to_Observable_Void_Signatures.pdf`

---

### Paper 48 — Maximum-Entropy Derivation of Universal Structural Selection
**A Maximum-Entropy Derivation of Universal Structural Selection:**
*From Constraint Defects to Exponential Structural Measures*

**Core idea:** Provides the formal closure behind the structural-selection
measure. If admissible configurations are selected under finite information
about structural defect means, Jaynes maximum entropy implies an exponential
selection measure:

```text
P[X] = Z^{-1} exp[-F_struct[X]]

F_struct[X] = sum_a lambda_a d_a[X]
```

The MAAT support form used in later papers is recovered as a logarithmic
multiplicative-support representation:

```text
F_MAAT[X] = - sum_a lambda_a log(epsilon + Gamma_a[X])
Gamma_a[X] = 1 / (1 + d_a[X])
```

**Core result:**

| Quantity | Result |
|----------|--------|
| Formal status | Conditional MaxEnt theorem |
| Required inputs | admissible space `X_adm`, reference measure `mu0`, defect observables `d_a` |
| Derived measure | `P[X] proportional to exp[-F_struct[X]]` |
| Linear form | Gibbs--Jaynes structural cost |
| Support form | bounded multiplicative viability representation |
| Weight interpretation | Lagrange multipliers / covariance-response pressures |
| No new numerical solver | formal derivation; empirical tests live in existing experiment folders |

**Key finding:**
> The exponential MAAT/structural-selection form is not an arbitrary scoring
> ansatz. It is the least-biased probability measure generated by fixed
> structural defect constraints.

**Scientific status:** This paper proves the form of the measure conditional
on the chosen defects, reference measure, and constraints. It does not prove
that the MAAT defect basis is unique or fundamental. Universality still
requires recurrence across domains and falsification against simpler
alternatives such as energy-only, entropy-only, stability-only, random,
anthropic, or complexity-only rankings.

**Reproducibility:** This is a formal derivation with no standalone code. The
empirical evidence ladder is reproduced by the existing benchmark folders in
`experiments/`, including cosmological CCI, growth/projection, string-landscape,
SAT/constraint, and environmental-response tests.

**Documentation PDF:** `documentation/48_Maximum_Entropy_Derivation_of_Universal_Structural_Selection.pdf`

---

### Paper 49 — Universal Structural Selection Benchmarks
**Universal Structural Selection Benchmarks:**
*Testing MAAT Against Energy, Stability, Random, and Domain-Specific Baselines*

**Core idea:** Turns the universality claim following Paper 48 into an
explicit falsification matrix. The paper audits existing experiments across
fields, string landscapes, cosmology, SAT/constraint systems, AI/safety
boundary calibration, SM-like constants, and quantum measurement.

Each domain is checked against competitor classes:

```text
energy-only
entropy/activity-only
stability-only
random or shuffled nulls
domain-specific baselines
cross-domain / holdout transfer
```

**Core results:**

| Quantity | Result |
|----------|--------|
| Domains audited | `7` |
| Complete domain evidence | `2` |
| Partial domain evidence | `4` |
| Missing domain evidence | `1` |
| Complete matrix cells | `8` |
| Partial matrix cells | `13` |
| Missing matrix cells | `17` |
| Benchmark-readiness index | `0.3816` |
| Strongest current successes | field fixed-energy and string/path energy-only baselines |
| Main open stress tests | SAT and quantum measurement |
| Immediate roadmap | SAT Structural Hardness II, Quantum Pointer-State Selection, Cross-Domain Transfer |

**Key finding:**
> Energy-only baselines are already defeated in field and string benchmarks,
> but the universal MAAT claim is not yet earned. Entropy/activity-only,
> stability-only, shuffled-null, SAT, and quantum tests remain the decisive
> next layer.
> In fixed-energy fields the structural score changes by nearly `8x` at
> sub-permille energy variation, while in the KKLT path graph the top-20
> structural and energy transitions have `0/20` overlap.

**Scientific status:** This is a meta-benchmark registry and audit paper. It
does not introduce a new solver and does not claim universal proof. It turns
the post-MaxEnt universality claim into a reproducible benchmark programme.
The low readiness index is interpreted as an audit warning, not as a success
metric. SAT is treated explicitly as an Achilles test: the current SAT result
is poor and motivates new graph-local defect sectors and future cross-solver
validation.

**Scripts and reproducibility:**

Repository URL:

```text
https://github.com/Chris4081/structural-selection-principle/tree/main/experiments/universal_structural_selection_benchmarks_paper49
```

| Folder | Role |
|--------|------|
| `experiments/universal_structural_selection_benchmarks_paper49/` | Paper 49 meta-benchmark aggregator, evidence registry, readiness matrix, JSON summary, and figures |

**Data attribution and license note:** The script aggregates derived outputs
already present in the repository. External scientific data used by the
underlying experiments remain attributed in their respective experiment folders
and papers. No new external dataset is redistributed here.

**Documentation PDF:** `documentation/49_Universal_Structural_Selection_Benchmarks.pdf`

---

### Paper 50 — SAT Structural Hardness II
**SAT Structural Hardness II:**
*Graph-Local Defects and Structural Hardness in Random 3-SAT*

**Core idea:** Directly addresses the SAT weakness identified in Paper 49.
The older SAT validation used coarse global fields and performed poorly.
Paper 50 generates a deterministic synthetic random-3-SAT ensemble and tests
new graph-local defects against density-only, standard graph, stability-only,
scalar MAAT, MAAT graph-local, and shuffled-null baselines.

**New structural diagnostics:**

```text
constraint-graph expansion
clause-variable interaction cycles
backdoor-density proxy
local covariance conditioning
unit-propagation depth
contradiction rate
literal and degree balance
```

**Core results:**

| Quantity | Result |
|----------|--------|
| Synthetic random-3-SAT instances | `520` |
| Variable range | `24` to `36` |
| Clause-density range | `3.17 <= alpha <= 5.22` |
| SAT fraction | `0.6346` |
| Timeout fraction | `0.0000` |
| Density-only CV `R2` | `0.1886` |
| Standard graph CV `R2` | `0.2469` |
| Stability-only CV `R2` | `-0.0145` |
| Scalar `F_MAAT` CV `R2` | `-0.0145` |
| MAAT graph-local CV `R2` | `0.1469` |
| MAAT graph-local + density CV `R2` | `0.2347` |
| MAAT+density gain over density-only | `+0.0461` |
| MAAT+density gap to standard graph | `-0.0122` |
| Shuffled-null 95% Spearman | `0.0623` |
| MAAT graph-local Spearman | `0.3914` |

**Key finding:**
> SAT hardness is graph-local and phase-aware. The new MAAT graph-local
> diagnostics clearly beat shuffled-null structure and improve over
> density-only once the SAT phase-density channel is retained, but they do not
> yet outperform a simple standard graph baseline.
> The key quantitative gain is `Delta R2 = +0.0461` over density-only,
> indicating nontrivial structural information beyond phase-transition
> proximity alone.

**Scientific status:** This is a partial success and useful failure. It does
not solve SAT hardness, prove `P != NP`, or establish MAAT universality.
It narrows the SAT problem: scalar `F_MAAT` and stability-only scores fail,
while local covariance conditioning and graph-local diagnostics carry real
but still insufficient hardness signal.
The graph-local diagnostics are solver-independent inputs, but the target is
a deterministic DPLL node-count proxy; solver-independent hardness remains a
future cross-solver validation task.

**Scripts and reproducibility:**

Repository URL:

```text
https://github.com/Chris4081/structural-selection-principle/tree/main/experiments/sat_structural_hardness_paper50
```

| Folder | Role |
|--------|------|
| `experiments/sat_structural_hardness_paper50/` | Paper 50 synthetic SAT benchmark, graph-local features, DPLL target, model comparisons, JSON/CSV outputs, and figures |

**Data attribution and license note:** The benchmark generates synthetic
random-3-SAT formulas with a fixed random seed. No external SAT benchmark
dataset is redistributed. Repository CSV/PNG files are derived
reproducibility artifacts generated by the script.

**Documentation PDF:** `documentation/50_SAT_Structural_Hardness_II.pdf`

---

### Paper 50b — SAT Structural Hardness IIb
**SAT Structural Hardness IIb:**
*Local Frustration Fields and the Failure of Scalar Compression*

**Core idea:** Companion test to Paper 50. Paper 50 showed that scalar
`F_MAAT` and global robustness are weak for SAT hardness. Paper 50b tests the
sharper thesis that SAT hardness is damaged by aggressive global scalar
compression and is better represented by local frustration fields, rare-event
tails, and hotspot geometry.

**Scientific status:** This is not the Paper 51 roadmap item announced in
Paper 49. It is a Paper-50 companion experiment. It does not solve SAT
hardness, prove `P != NP`, or provide a modern CDCL benchmark. The target is
still a deterministic DPLL node-count proxy.

**Core results:**

| Quantity | Result |
|----------|--------|
| Synthetic SAT instances | `550` |
| Formula families | random, planted, modular |
| Variable range | `24` to `38` |
| Clause-density range | `3.26 <= alpha <= 5.18` |
| SAT fraction | `0.7018` |
| Timeout fraction | `0.0000` |
| RF density-only CV `R2` | `0.0362` |
| RF density + formula type CV `R2` | `0.1834` |
| RF standard graph CV `R2` | `0.2959` |
| RF global MAAT scalar/supports CV `R2` | `0.1375` |
| RF single `F_MAAT^multi` scalar CV `R2` | `-0.0806` |
| RF `F_MAAT^multi` components CV `R2` | `-0.0224` |
| RF `F_MAAT^multi` components + density/type CV `R2` | `0.2833` |
| RF local frustration + density/type CV `R2` | `0.2946` |
| RF raw graph + frustration CV `R2` | `0.3152` |
| Local frustration + density/type gain over global MAAT | `+0.1571` |
| Multi components + density/type gain over global MAAT | `+0.1458` |
| Raw graph + frustration gain over standard graph | `+0.0194` |

**Key finding:**
> SAT hardness is not well described by global robustness or scalar
> `F_MAAT`. The useful signal lives in graph-local, phase-aware, multi-scale
> frustration geometry. Local frustration fields only become useful when the
> density/formula-family channel is retained. The Paper-50c test confirms that
> `F_MAAT^multi` should not be re-compressed into one scalar; its components
> carry the useful signal when context is retained.

**Scripts and reproducibility:**

Repository URL:

```text
https://github.com/Chris4081/structural-selection-principle/tree/main/experiments/sat_frustration_fields_paper50b
```

| Folder | Role |
|--------|------|
| `experiments/sat_frustration_fields_paper50b/` | Paper 50b SAT frustration-field compression test, synthetic instances, model ladder, feature importance, JSON/CSV outputs, and figures |

**Data attribution and license note:** The benchmark generates synthetic SAT
formulas with a fixed random seed. No external SAT benchmark dataset is
redistributed. Repository CSV/PNG files are derived reproducibility artifacts
generated by the script.

**Documentation PDF:** `documentation/50b_SAT_Frustration_Field_Compression_Test.pdf`

---

### Paper 50c — SAT Structural Hardness IIc
**SAT Structural Hardness IIc:**
*From Global Robustness to Multi-Scale Defect-Field Structural Selection*

**Core idea:** Theoretical companion note to Papers 50 and 50b. It does not
replace the old MAAT formula, but reinterprets it as a low-frequency macro
projection of a richer defect-field functional. SAT motivates the update
because scalar `F_MAAT` and `R_rob` lose rare local frustration, tails,
hotspots, and cluster structure.

**Updated master form:**

```text
F_MAAT^multi = F_mean + F_tail + F_cluster + F_scale
```

with:

```text
F_mean    = sum_a lambda_a^(0) <d_a(x)>_G
F_tail    = sum_a lambda_a^tail q_p(d_a)
F_cluster = sum_a lambda_a^cl C_max^(a)/|G|
F_scale   = sum_{a,k} lambda_{a,k} d_a^(k)
```

**Key finding:**
> Global MAAT is the macro projection; full MAAT is the defect-field
> functional. The old scalar formula remains useful for smooth or
> self-averaging systems, while SAT-like domains require local, tail, cluster,
> and scale-resolved defect structure.
> The companion code test shows the same lesson numerically: a single
> `F_MAAT^multi` scalar is weak, but its components plus density/type reach
> `R2 = 0.2833`, far above global MAAT compression.

**Scientific status:** Effective-theory update, not a microscopic derivation
and not a new benchmark. It is the theoretical synthesis of the Paper 50/50b
SAT failure analysis.

**Documentation PDF:** `documentation/50c_Multi_Scale_Defect_Field_Structural_Selection.pdf`

---

### Paper 51 — Quantum Pointer-State Selection
**Quantum Pointer-State Selection:**
*A MAAT v1.4 Defect-Field Benchmark in a Non-Commuting Lindblad System*

**Core idea:** Tests the MAAT v1.4 defect-field thesis in a minimal open
quantum system. A qubit evolves under non-commuting Lindblad channels:
`z`-dephasing, `x`-dephasing, and amplitude damping. The target is initial
state survival fidelity after Lindblad evolution. All features are computed
from the initial state and generator parameters, not from the evolved final
state.

**Core results:**

| Model | CV `R2` |
|-------|--------:|
| MAAT scalar/supports | `0.5383 ± 0.0151` |
| MAAT v1.4 field features | `0.7424 ± 0.0057` |
| Rates only | `0.0542 ± 0.0547` |
| Basis geometry only | `-0.1134 ± 0.0331` |
| Hamiltonian only | `-0.0406 ± 0.0241` |
| Domain combo | `0.1362 ± 0.0291` |
| v1.4 gain over scalar MAAT | `+0.2041` |
| v1.4 gain over domain combo | `+0.6062` |

**Key finding:**
> Pointer-state robustness is not captured by rates alone. Scalar MAAT retains
> moderate signal, but the v1.4 field-feature extension substantially improves
> prediction by combining environment alignment, channel conflict, damping
> pressure, coherence exposure, and robustness closure.

**Scientific status:** Toy benchmark only. This is not a derivation of quantum
measurement, not a Born-rule derivation, not an experimental fit, and not a
full decoherence theory.

**Scripts and reproducibility:**

Repository URL:

```text
https://github.com/Chris4081/structural-selection-principle/tree/main/experiments/quantum_pointer_state_selection_paper51
```

| Folder | Role |
|--------|------|
| `experiments/quantum_pointer_state_selection_paper51/` | Paper 51 quantum pointer-state toy benchmark, generated ensemble, model comparison, feature importance, JSON/CSV outputs, and figures |

**Data attribution and license note:** The benchmark generates synthetic
quantum toy-model data with a fixed random seed. No external experimental
dataset is redistributed. Repository CSV/PNG files are derived reproducibility
artifacts generated by the script.

**Documentation PDF:** `documentation/51_Quantum_Pointer_State_Selection.pdf`

---

### Paper 52 — Cross-Domain Transfer of Structural Selection
**Cross-Domain Transfer of Structural Selection:**
*A Frozen-Weight Two-Domain Pilot Study*

**Core idea:** Tests whether a MAAT defect architecture calibrated in one
domain transfers to a second internal benchmark domain without target-domain
retuning. The common architecture uses only the shared supports
`H,B,S,V,R_rob`. Source-domain weights are learned once, frozen, and then
evaluated in the target domain with a predeclared score direction and
shuffled-defect nulls.

**Core results:**

| Transfer | Frozen Spearman | Equal/scalar Spearman | Interpretation |
|----------|----------------:|----------------------:|----------------|
| SAT -> Quantum | `0.6034` | `0.3771` | strong asymmetric transfer above shuffled-defect null |
| Quantum -> SAT | `-0.0408` | `-0.2683` | weak/failed positive transfer, but less bad than scalar baselines |

**Key finding:**
> Transfer is partial and asymmetric. SAT-trained connectivity transfers
> strongly to quantum pointer robustness, while quantum-trained structure does
> not solve SAT hardness. The dominant transferred signal is `V`
> connectivity, not a balanced multi-sector architecture. The conservative
> interpretation is that the two internal benchmarks share a
> connectivity-sensitive structural mode; this is not yet independent evidence
> for universal cross-domain transfer.
> The result supports the v1.4/v1.4.1 reading that universality, if present,
> is more likely to live in partially transferable structural modes than in one
> compressed scalar score.

**Scientific status:** Two-domain internal transfer benchmark, not a
universality proof and not independent external validation. The protocol is
intentionally adversarial: frozen transfer is allowed to fail. The result is
deliberately reported as mixed: one direction succeeds strongly, the reverse
direction remains a failure mode for global scalar transfer. The SAT side uses
DPLL node count as a transparent controlled proxy, not modern CDCL conflict
statistics.

**Scripts and reproducibility:**

Repository URL:

```text
https://github.com/Chris4081/structural-selection-principle/tree/main/experiments/cross_domain_transfer_paper52
```

| Folder | Role |
|--------|------|
| `experiments/cross_domain_transfer_paper52/` | Paper 52 frozen cross-domain transfer runner, shuffled-null validation, CSV/JSON outputs, and figures |

**Data attribution and license note:** This experiment reuses repository-internal
synthetic/derived artifacts from Paper 50b and Paper 51. No new external
dataset is introduced. Repository CSV/PNG files are reproducibility artifacts
generated by the script.

**Documentation PDF:** `documentation/52_Cross_Domain_Transfer_Structural_Selection.pdf`

---

### Paper 53 — Structural Mode Transfer Robustness
**Robustness of Structural Mode Transfer in MAAT v1.4:**
*Alternative Frozen Fits, Regularisation Tests, and Multi-Domain Null Validation*

**Core idea:** Tests whether the Paper 52 `V`/connectivity result is stable
under alternative source-only fit rules, regularisation schemes, and additional
internal benchmark domains. This is an internal robustness audit, not an
external cross-domain validation. The protocol keeps the frozen-transfer rule:
source weights are learned once, frozen, and transferred without target-domain
retuning or sign flips.

**Domains and fit rules:**

| Category | Included |
|----------|----------|
| Domains | SAT-Frustration, Quantum-Pointer, Active-Significance, Societal-CCI |
| Fit rules | equal, NNLS, positive ridge, positive Lasso, response-top20, V-only |
| Nulls | shuffled-defect nulls for cross-domain transfers |

**Core results:**

| Result | Value |
|--------|------:|
| SAT -> Quantum, NNLS | `rho = 0.6034`, `lambda_V = 1.0000` |
| SAT -> Quantum, ridge_pos | `rho = 0.5981`, `lambda_V = 0.8394` |
| SAT -> Quantum, lasso_pos | `rho = 0.6034`, `lambda_V = 1.0000` |
| SAT -> Quantum, response_top20 | `rho = 0.5975`, `lambda_V = 0.8264` |
| Mean cross-domain Spearman, equal | `0.5156` |
| Overall V-dominant frozen architectures | `~58%` |

**Key finding:**
> The Paper 52 SAT -> Quantum connectivity result is not merely an NNLS
> artifact: it survives ridge, Lasso, covariance-response closure, and a V-only
> baseline. But `V` is not universal. Across all internal domains, broad equal
> support averaging performs best on support-composite toy targets, and SAT
> remains resistant to simple incoming transfer.

The equal-weight result is an important negative control rather than a cosmetic
detail: in the present partly support-composite domain set, simple broad
averaging remains competitive with more sophisticated source-specific fitting.
Thus Paper 53 supports conditional structural portability, not universal
transfer.

**Scientific status:** Robustness/ablation paper, not external validation. The
added Active-Significance and Societal-CCI domains are internal toy domains
whose targets are partly support-composite, so aggregate cross-domain averages
must be interpreted cautiously. The Societal-CCI ensemble is deliberately
synthetic and is not a real social dataset, political ranking, or independent
external validation target.

**Scripts and reproducibility:**

Repository URL:

```text
https://github.com/Chris4081/structural-selection-principle/tree/main/experiments/structural_mode_transfer_paper53
```

| Folder | Role |
|--------|------|
| `experiments/structural_mode_transfer_paper53/` | Paper 53 multi-fit/multi-domain transfer robustness runner, shuffled-null validation, CSV/JSON outputs, and figures |

**Data attribution and license note:** This experiment reuses repository-internal
synthetic/derived artifacts from earlier papers. No new external dataset is
introduced. Repository CSV/PNG files are reproducibility artifacts generated by
the script.

**Documentation PDF:** `documentation/53_Structural_Mode_Transfer_Robustness.pdf`

---

### Paper 54 — External Validation of Structural Modes
**External Validation of Structural Modes:**  
*MAAT Defect Features on Public Datasets Accessed via scikit-learn*

**Core idea:** Breaks the internal-artifact loop of Papers 52--53 by testing
MAAT structural supports on public external machine-learning datasets loaded
via `scikit-learn`. Sample-level hardness is predeclared as the fraction
of repeated cross-validation runs in which a sample is misclassified by a
fixed Random Forest classifier.

**External datasets:**

| Dataset | Samples | Features | Classes |
|---------|--------:|---------:|--------:|
| `breast_cancer` | 569 | 30 | 2 |
| `wine` | 178 | 13 | 3 |
| `digits` | 1797 | 64 | 10 |

**Core results:**

| Model / score | Mean Spearman | Mean AUROC |
|---------------|--------------:|-----------:|
| external geometry feature set | `0.3889` | `0.6685` |
| combined geometry feature set | `0.3824` | `0.6666` |
| MAAT supports | `0.3241` | `0.6393` |
| MAAT v1.4 field features | `0.2795` | `0.6191` |
| label disagreement scalar | `0.4299` | `0.6719` |
| MAAT robustness-loss scalar | `0.3301` | `0.6708` |
| MAAT scalar cost | `0.3226` | `0.6706` |

**Key finding:**
> MAAT scalar cost and robustness loss carry positive external
> sample-hardness signal on all three public datasets, but simple external
> geometry and local label-disagreement baselines remain stronger in this
> first public-data benchmark.

**Scientific status:** First external validation benchmark, not a proof of
universal transfer. The MAAT v1.4 field features beat shuffled-defect nulls on
breast cancer and digits, but not on wine. The result is therefore mixed and
useful: it establishes measurable external signal while preventing internal
consistency from being mistaken for external universality.

**Scripts and reproducibility:**

Repository URL:

```text
https://github.com/Chris4081/structural-selection-principle/tree/main/experiments/external_ml_hardness_paper54
```

| Folder | Role |
|--------|------|
| `experiments/external_ml_hardness_paper54/` | Paper 54 external ML sample-hardness benchmark, scikit-learn dataset loader, CV hardness target, CSV/JSON outputs, and figures |

**Data attribution and license note:** The experiment accesses public datasets
via `scikit-learn`. No raw external dataset files are redistributed in the
experiment folder. Generated CSV/JSON/PNG files are derived reproducibility
artifacts. Cite `scikit-learn` and the original dataset sources as appropriate.
No endorsement by `scikit-learn` maintainers, UCI dataset curators, or original
dataset authors is implied.

**Documentation PDF:** `documentation/54_External_Validation_Structural_Modes.pdf`

---

### Paper 55 — CDCL Hardness Prediction on Standard SAT Families
**CDCL Hardness Prediction on Standard SAT Families:**  
*A Benchmark for Graph-Local and Defect-Field Features*

**Core idea:** Replaces the earlier transparent DPLL SAT-hardness target with
real CDCL solver statistics from PySAT implementations of Glucose3 and
MiniSat22. The benchmark uses deterministically generated standard SAT
families rather than previous repository-internal artifacts: random 3-SAT,
planted 3-SAT, modular 3-SAT, 3-XOR CNF, graph coloring, and pigeonhole
formulas.

**Target:**

```text
y_CDCL = log(1 + mean_conflicts + 0.10 mean_decisions + 0.01 mean_propagations)
```

**Core results:**

| Feature set | RF R2 | RF Spearman |
|-------------|------:|------------:|
| Standard graph | `0.2742` | `0.4254` |
| Multi-scale defect field | `0.2621` | `0.3984` |
| Graph-local structural | `0.2573` | `0.4170` |
| Scalar structural | `0.1062` | `0.3400` |
| Shuffled defect null | `-0.0609` | `-0.0070` |
| Density-only | `-0.0894` | `0.2503` |

**Key finding:**
> Multi-scale defect-field features carry real CDCL-hardness signal and
> strongly beat shuffled-defect nulls and scalar compression, but they do not
> beat a strong standard graph baseline. The decisive negative result is
> leave-family-out transfer: all feature sets fail to extrapolate cleanly across
> canonical SAT families.

**Scientific status:** This is a standard-family CDCL benchmark, not a SAT
Competition or industrial SAT benchmark. It strengthens the local
defect-field direction while weakening broad universality claims.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/sat_cdcl_external_paper55/` | canonical SAT-family generator, PySAT/CDCL solver statistics, feature extraction, shuffled-null tests, model comparisons, figures |

**Data attribution and license note:** No external CNF dataset is redistributed.
Instances are generated from standard SAT benchmark families by the script.
Solver statistics are produced locally through PySAT-compatible Glucose3 and
MiniSat22 solvers. Generated CSV/JSON/PNG files are derived reproducibility
artifacts.

**Documentation PDF:** `documentation/55_External_SAT_CDCL_Hardness_Validation.pdf`

---

### Paper 56 — From Scalar Scores to Defect-Field Geometry
**From Scalar Scores to Defect-Field Geometry:**  
*Formal Definitions for Multi-Scale Structural Selection*

**Core idea:** Consolidates the operator language needed by MAAT v1.4.1:
local defect fields, multi-scale coarse-graining, defect-field dynamics,
observation projections, and no-retuning universality tests. The paper is a
formal definition layer rather than a new numerical benchmark.

**Key definitions:**

```text
d_a(x) = D_a(Phi|N(x), C|N(x), partial N(x), mu|N(x)) >= 0
d_a^(ell) = Pi_ell[d_a]
partial_logell d_a^(ell) = beta_a[d^(ell)]
partial_t d_a = F_a[d] + xi_a
Pi_obs: Omega_full -> Omega_obs
U(theta) = E_{D_test notin D_train}[1{Delta I > 0}]
```

**v1.4.1 anchor:** The paper explicitly uses the v1.4.1 lesson that
structural sectors must be operationally separated rather than merely renamed
correlated summaries. In the quantum two-level anchor case:

```text
B = 1 - |p0 - p1|
S = |rho_01| / sqrt(rho_00 rho_11)
```

This keeps balance and activity distinct: balanced-inactive,
balanced-coherent, imbalanced-incoherent, and imbalanced-coherent states are
not collapsed into the same structural mode.

**Scientific status:** This is a consolidation and definition paper. It does
not claim a microscopic derivation, a complete RG theorem, or universal
cross-domain transfer. It defines the formal objects needed to make those
claims testable.

**Documentation PDF:** `documentation/56_Defect_Field_Geometry_Projection_and_Universality.pdf`

---

### Paper 57 — Where Does SAT Hardness Live?
**Where Does SAT Hardness Live? Mean, Tail, Cluster, and Scale Projections of Defect-Field Geometry**

**Core idea:** Tests the formal projection language of Paper 56 on the Paper-55
SAT/CDCL benchmark. Instead of asking whether one global structural scalar
predicts hardness, the experiment decomposes the defect field into mean, tail,
cluster, and scale projections and asks which component carries the CDCL
hardness signal.

**Core results:**

| Feature set / projection | 5-fold RF R2 | Spearman |
|--------------------------|-------------:|---------:|
| Standard graph baseline | `0.2169` | `0.4234` |
| All defect projections | `0.1996` | `0.3751` |
| Macro scalar | `0.1898` | `0.4193` |
| Tail + scale | `0.1725` | `0.3793` |
| Scale projection | `0.1381` | `0.3304` |
| Tail projection | `0.0937` | `0.2982` |
| Mean projection | `0.0575` | `0.2944` |
| Shuffled all-projection null | `-0.0666` | `-0.0134` |

**Projection importance inside the all-defect model:**

| Projection | Grouped R2 drop |
|------------|----------------:|
| Tail | `0.4044` |
| Scale | `0.3414` |
| Mean | `0.1567` |
| Cluster | `0.0229` |

**Key finding:**
> SAT/CDCL hardness signal is not primarily localized in mean scalar
> compression. In this benchmark, tail and scale projections carry most of the
> internal defect-field signal, while leave-family-out transfer remains weak.

**Scientific status:** This is a projection-decomposition benchmark, not a
solution of SAT hardness and not evidence for universal transfer. The standard
graph baseline remains slightly stronger in ordinary cross-validation. The
positive result is signal localization; the negative result is continued
family-transfer failure.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/sat_projection_decomposition_paper57/` | Paper-57 projection ablation, Paper-55 CSV reader, mean/tail/cluster/scale tests, shuffled-null output, figures |

**Data attribution and license note:** No external CNF files are redistributed.
The experiment uses derived reproducibility artifacts generated by Paper 55
from deterministic canonical SAT-family generators and PySAT-compatible solver
statistics.

**Documentation PDF:** `documentation/57_Where_Does_SAT_Hardness_Live.pdf`

---

### Paper 58 — Stationarity Balance Beats Population Balance
**Stationarity Balance Beats Population Balance:**  
*A Defect-Field Benchmark on Open Quantum Pointer States*

**Core idea:** Tests structural defect-field diagnostics on a fully classical
one-qubit Lindblad simulation. No quantum computer, Qiskit runtime, cloud
hardware, or paid backend access is required. The benchmark asks whether
stationarity-sensitive balance predicts robust pointer-like states better than
raw population balance.

**Simulation setup:**

| Component | Description |
|----------|-------------|
| System | one-qubit density matrix |
| Dynamics | Lindblad master equation integrated by RK4 |
| Channel families | z-dephasing, x-dephasing, amplitude damping, thermal relaxation, mixed z/x, depolarizing |
| Instances | `3600` sampled trajectories reused from a `5200`-trajectory generated ensemble |
| Target | `0.50 fidelity_retention + 0.35 probability_stability + 0.15 purity_final` |

**Key definitions:**

```text
B_pop  = 1 - |p0 - p1|
B_stat = 1 / (1 + 2 |dp/dt|_0 + 0.35 |d Tr(rho^2)/dt|_0)
```

**Core results:**

| Feature set | 5-fold RF R2 | Spearman |
|-------------|-------------:|---------:|
| Standard quantum baseline | `0.9011` | `0.9241` |
| Defect-field stationary | `0.8527` | `0.8856` |
| Scalar stationarity | `0.6457` | `0.7719` |
| Scalar population balance | `0.4561` | `0.6221` |
| Rates only | `0.1641` | `0.4646` |
| State geometry only | `0.0008` | `0.2945` |
| Shuffled defect null | `-0.0261` | `-0.0003` |

**Leave-channel-family-out results:**

| Feature set | LFO R2 | LFO Spearman |
|-------------|------:|-------------:|
| Defect-field stationary | `0.7075` | `0.8286` |
| Standard quantum baseline | `0.6637` | `0.8105` |
| Scalar stationarity | `0.4928` | `0.7284` |
| Scalar population balance | `0.3108` | `0.5604` |

**Key finding:**
> Pointer robustness is better captured by stationarity-sensitive balance than
> by raw population equality. The scalar stationarity model improves over the
> scalar population-balance model by `Delta R2 = +0.1896`, and the full
> stationarity-sensitive defect-field model performs best under
> leave-channel-family-out transfer.

**Scientific status:** This is a reproducible classical open-system
simulation. It is not a proof of quantum measurement, not a collapse theory,
and not a hardware experiment. The result supports a narrower refinement:
quantum balance should track stationarity under the generator, not mere
population equality.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/open_quantum_pointer_paper58/` | Paper-58 Lindblad simulation, pointer-robustness target, B_pop/B_stat comparison, model outputs, figures |

**Data attribution and license note:** No external quantum dataset is
redistributed. All CSV/JSON/PNG files are derived reproducibility artifacts
generated locally by the script.

**Documentation PDF:** `documentation/58_Open_Quantum_Pointer_Stability.pdf`

---

### Extra Framework Paper — Generator Stationarity and Dynamical Defect-Field Geometry
**Generator Stationarity and Dynamical Defect-Field Geometry:**  
*Refining Structural Selection for Evolving Systems (MAAT v1.5)*

**Core idea:** Consolidates the next MAAT version after v1.4/v1.4.1.
The central update is that structural supports should be defined relative to
the generator that preserves, damps, transports, or destabilises them. Static
support definitions are therefore not enough; they must pass a
generator-stationarity sanity check.

**Key rule:**

```text
A support coordinate should measure what the dynamics preserve,
not merely what looks balanced in a static snapshot.
```

**v1.5 master functional:**

```text
F_MAAT^(1.5)
= F_mean + F_tail + F_cluster + F_scale + F_stat
```

where `F_stat` measures stationarity defects relative to a domain generator
`G`, such as a Lindbladian, Hamiltonian, field-equation residual, solver
transition rule, coarse-graining map, or nonlinear flow.

**Core result from Paper 58 motivating v1.5:**

```text
B_pop  = 1 - |p0 - p1|
B_stat = 1 / (1 + 2 |dp/dt|_0 + 0.35 |d Tr(rho^2)/dt|_0)
```

Stationarity balance outperforms population balance for open quantum
pointer robustness:

| Comparison | Result |
|------------|-------:|
| Scalar stationarity R2 | `0.6457` |
| Scalar population balance R2 | `0.4561` |
| Improvement | `Delta R2 = +0.1896` |

**Key finding:**
> MAAT v1.5 turns the defect-field programme from a static structural language
> into a generator-aware structural language. Balance should track
> stationarity under the relevant dynamics, not mere static equality.

**Scientific status:** This is a formal framework note, not a new empirical
benchmark and not a microscopic derivation of the defect sectors. It tightens
the operational definition of supports by requiring generator-aware
validation.

**Reproducibility anchor:**

| Folder | Role |
|--------|------|
| `experiments/open_quantum_pointer_paper58/` | Paper-58 Lindblad benchmark motivating the v1.5 stationarity rule |

**Documentation PDF:** `documentation/MAAT_v15_Generator_Stationarity_and_Dynamical_Defect_Fields.pdf`

---

### Extra Framework + Experiment Paper — Structural Search Geometry
**Structural Search Geometry:**  
*A MAAT v1.6 Framework Note with a Toy Search Benchmark*

**Core idea:** Reframes MAAT from state ranking alone into a search-geometry
layer over complex candidate spaces. Defect supports become
priority-modulating coordinates for proof search, SAT search, mathematical discovery, AI
representation search, physical configuration search, and scientific
hypothesis search.

**v1.6 design rule:**

```text
Use structural supports to modulate verification-driven search.
Do not let structural support replace the target, proof obligation,
or validation criterion.
```

**Search functional:**

```text
F_search(x)
= lambda_D d_ver(x)
 + lambda_E E(x)
 + lambda_M[-sum_a log(epsilon + Gamma_a(x))]
 + lambda_R(1 - R_rob(x))
 - lambda_V V(x)
```

where `d_ver` is the distance to the verification target, `E` is a domain
objective or energy proxy, and `V` rewards bridge-like or connective search
positions.

**Toy benchmark:** Synthetic graph-search spaces with four families:
`smooth`, `trap`, `bridge`, and `trap_bridge`.

**Core results:**

| Algorithm | Success rate | Median expansions |
|-----------|-------------:|------------------:|
| distance-only | `1.000` | `59.0` |
| connectivity-only | `1.000` | `59.0` |
| MAAT v1.6 anchored | `1.000` | `92.0` |
| energy-only | `1.000` | `106.0` |
| MAAT v1.6 heavy | `1.000` | `234.0` |
| random frontier | `0.741` | `702.5` |

**Key finding:** Anchored structural search improves over energy-only search in
trap-like landscapes (`84.5` vs `150.5` median expansions in `trap`; `92.0`
vs `183.5` in `trap_bridge`; roughly 44--50% fewer median expansions) but does
not beat simple distance/connectivity baselines on clean grid tasks. This is
expected because the target is visible, Euclidean distance is cheap, and the
bridges are simple graph bottlenecks. Heavy structural weighting is harmful.
This turns v1.6 into a constrained search-geometry programme rather than a
broad claim about universal discovery.

**Scientific status:** Toy benchmark and framework note. It is not an
algorithmic optimality proof and not evidence that MAAT outperforms classical
search heuristics in general. The next meaningful tests should use richer
spaces where local evaluation is costly or noisy and bridges are not visible
as simple geometric bottlenecks: theorem-search graphs, CDCL/SAT state spaces,
graph-construction problems, model-configuration spaces, or optimisation tasks
with expensive validation.

**Reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/structural_search_geometry_v16/` | MAAT v1.6 structural search toy benchmark, CSV/JSON outputs, and figures |

Run:

```bash
cd experiments/structural_search_geometry_v16
python3 structural_search_geometry_v16.py
```

**Data attribution and license note:** No external dataset is redistributed.
All CSV/JSON/PNG files are derived reproducibility artifacts generated from
synthetic search-space instances.

**Documentation PDF:** `documentation/MAAT_v16_Structural_Search_Geometry.pdf`

---

### Extra Framework Specification — MAAT v1.6-T Triage Instantiation
**MAAT v1.6-T: The Triage Instantiation:**  
*Structural Search Geometry as a Pre-Inference Layer*

**Core idea:** Defines the policy/triage face of MAAT v1.6.  The v1.6 measure
face ranks paths through structural search spaces; the v1.6-T triage face
ranks candidate continuations before full inference, proof, simulation, or
solving is committed.

**Search arena:**

```text
(N, C, cost, goal)
```

where `N` is the search-state space, `C(n)` the candidate-continuation set,
`cost` the expansion cost, and `goal` the domain validation criterion.

**Unification equation:**

```text
pi(c | n) proportional exp[-beta * h_MAAT^(h)(c)]
```

At `h -> infinity`, the finite-horizon heuristic approaches a
quasipotential-like cost-to-go.  At `h = 0`, it becomes greedy structural
triage.

**Goal-blindness guardrail:**

```text
F_MAAT measures coherence, not progress.
```

Therefore v1.6-T permits two sanctioned usage modes:

| Mode | Role |
|------|------|
| Mode F | feasibility prior; structural support filters or down-weights broken candidates |
| Mode A | combined acquisition rule: `A(c)=prog(c)-tau*h_MAAT(c)` |

**Predeclared obligations:**

| Obligation | Test |
|------------|------|
| T1 | insert Mode A into CDCL and measure compute-to-solution regret against VSIDS/LRB |
| T2 | compare Mode F, Mode A, progress-only, and MAAT-only ablations |
| T3 | frozen-weight triage transfer across domains |
| T4 | verify finite-horizon `h_MAAT` approaches quasipotential orderings in small arenas |

**Scientific status:** Framework specification, not evidence.  Paper 62
implements the first T1 benchmark and shows that global `h=0` Mode A is mixed:
family-sensitive and not yet globally compute-positive.

**Documentation PDF:** `documentation/MAAT_v16_T_Triage_Specification.pdf`

---

### Paper 59 — MAAT v1.6 Guided SAT Search
**MAAT v1.6 Guided SAT Search:**  
*Structural Search Geometry for Transparent DPLL Branching*

**Core idea:** Moves the SAT programme from passive hardness prediction to
active search steering. A transparent DPLL solver with unit propagation uses
MAAT v1.6 structural-search coordinates as a branching heuristic and compares
them against classical baselines.

**Scope:** This is not a proof of `P != NP`, not a CDCL solver competition, and
not a replacement for modern solvers such as Glucose, MiniSat, CaDiCaL, or
Kissat. It is a controlled first test of whether structural search geometry
can provide active branching signal.

**Tested heuristics:**

```text
first, random, degree, Jeroslow-Wang, MOMS,
MAAT v1.6 branch, MAAT v1.6 heavy
```

**Core aggregate results:**

| Heuristic | Solved rate | Median decisions | Median conflicts | Median log cost |
|-----------|------------:|-----------------:|-----------------:|----------------:|
| MOMS | `1.000` | `9.0` | `2.0` | `2.6283` |
| Jeroslow-Wang | `1.000` | `11.0` | `2.0` | `2.7788` |
| MAAT v1.6 branch | `1.000` | `12.0` | `2.0` | `2.9178` |
| MAAT v1.6 heavy | `1.000` | `12.0` | `2.0` | `2.9232` |
| Degree | `1.000` | `14.0` | `3.0` | `3.0587` |
| First | `0.993` | `17.0` | `6.0` | `3.2696` |
| Random | `1.000` | `18.0` | `5.0` | `3.3087` |

**Key finding:**
> MAAT v1.6 provides a usable structural branching signal: it beats degree,
> first-unassigned, and random branching in aggregate. In its present form,
> however, it does not yet compete with the strongest SAT-specific heuristics
> tested here, MOMS and Jeroslow-Wang.

**Important negative control:** The heavy MAAT variant is slightly worse than
the anchored MAAT brancher (`2.9232` vs. `2.9178` median log cost). This
supports the v1.6 lesson that structural search geometry must remain anchored
to immediate verification pressure; increasing structural weight alone does not
automatically improve search.

**Graphical representation:** Paper 59 now also includes SAT graph
visualisations analogous in spirit to dense mathematical graph drawings:
variable co-occurrence graphs, clause-variable bipartite graphs, and a
branching-component plot. These figures act as root-state feature maps:
they show visible structural pressure points and connect them to verification
gain, short-clause pressure, connectivity, and robustness closure. The numerical
results then show the limit of that signal: useful, but still weaker than
short-clause-specialised SAT heuristics.

**Scientific status:** First active SAT-search benchmark in the MAAT series.
The solver is DPLL with unit propagation only; there is no clause learning,
restart policy, VSIDS, LRB, CHB, or phase saving.
The generated benchmark is light-to-moderate rather than adversarially hard:
most nontrivial heuristics solve almost all instances, so the experiment
primarily measures branching-cost ordering rather than hard timeout separation.
The improvement over degree, first-unassigned, and random is therefore a
moderate positive signal, not a scalability claim.

**Reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/sat_v16_guided_search_paper59/` | Paper-59 DPLL branching benchmark, generated SAT families, heuristic comparisons, figures, CSV/JSON outputs |

Run:

```bash
cd experiments/sat_v16_guided_search_paper59
python3 sat_v16_guided_search_paper59.py
python3 sat_v16_structure_visualization.py
```

**Data attribution and license note:** No external CNF dataset is redistributed.
Instances are generated locally from standard synthetic SAT-family generators.
CSV/JSON/PNG outputs are derived reproducibility artifacts.

**Documentation PDF:** `documentation/59_MAAT_v16_Guided_SAT_Search.pdf`

---

### Paper 60 — Fluid Coherence and Blow-up Diagnostics
**Fluid Coherence and Blow-up Diagnostics:**  
*A MAAT Early-Warning Benchmark for Burgers and Navier--Stokes Toy Flows*

**Core idea:** Tests whether MAAT-style structural supports can act as
early-warning diagnostics for fluid steepening and turbulence-like activity.
The paper is explicitly not a Navier--Stokes regularity proof. It uses
viscous Burgers as a finite-gradient warning proxy and unforced 2D
Navier--Stokes as a regular control case.

**Operational supports:**

```text
H      equation-residual consistency
B      energy/enstrophy balance consistency
S_eff  controlled activity pressure
V      low-mode spectral coherence
R_rob  v1.2.1 robustness closure
W_MAAT robustness-loss warning with spectral-tail pressure
```

**Core results:**

| Quantity | Result |
|----------|-------:|
| Burgers inviscid shock time | `1.3334` |
| Burgers first 70% warning time | `1.1600` |
| Burgers warning lead time | `0.1734` |
| Spearman warning vs max gradient | `0.8413` |
| Pre-shock warning vs inverse time-to-shock | `0.9952` |
| 2D Navier-Stokes false high-warning fraction | `0.0000` |
| 2D Navier-Stokes minimum robustness | `0.8956` |

**Key finding:**
> MAAT fluid supports can track loss of structural coherence in a Burgers
> shock-warning setting while remaining bounded in a regular 2D
> Navier--Stokes control. This supports structural early-warning diagnostics,
> not a proof of fluid regularity.

**Interpretation note:** The Burgers correlation is partly expected because the
warning includes activity pressure and spectral-tail growth. The stronger
control result is the bounded no-false-alarm behaviour in the regular 2D
Navier--Stokes run.

**Scientific status:** Toy diagnostic benchmark. The Burgers test provides a
known finite-gradient warning target; the 2D Navier--Stokes run is a regular
control, not a 3D blow-up test. A serious next step would use forced turbulence,
vortex-stretching proxies, adaptive simulations, or public DNS datasets.

**Reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/fluid_blowup_diagnostics_paper60/` | Paper-60 Burgers and 2D Navier--Stokes diagnostic benchmark, CSV/JSON outputs, and plots |

Run:

```bash
cd experiments/fluid_blowup_diagnostics_paper60
python3 fluid_blowup_diagnostics_paper60.py
```

**Data attribution and license note:** No external fluid dataset is
redistributed. All fields are generated locally from synthetic initial
conditions. CSV/JSON/PNG outputs are derived reproducibility artifacts.

**Documentation PDF:** `documentation/60_MAAT_Fluid_Coherence_Blowup_Diagnostics.pdf`

---

### Paper 61 — SPARC Pilot Test of Structural Dark-Matter Residuals
**SPARC Pilot Test of Structural Dark-Matter Residuals:**  
*A MAAT v1.6 Real-Data Bridge from Rotation-Curve Residuals to NFW/MOND Baselines*

**Core idea:** Converts the previous synthetic dark-matter residual note into
a real SPARC pilot test.  It prepares the Zenodo-hosted SPARC `Rotmod_LTG`
rotation-curve files, computes baryonic decomposition and MAAT residual
supports, and compares the structural residual against NFW-like and MOND/RAR
baseline mismatch channels plus shuffled-galaxy nulls.

**Scientific status:** First real-data pilot, not a precision halo fit and not
a dark-matter solution.  The NFW channel is lightweight, stellar
mass-to-light factors are fixed, and the result must be checked against
published SPARC halo catalogues.
This is a feasibility and signal test: it shows applicability to real data and
non-random structure, not superiority over professional halo or MOND/RAR
analyses.

**Core diagnostic:**

```text
D_struct(r) = f_res(r) * R_rob(r) * V(r)
```

with

```text
v_res^2(r) = max(v_obs^2(r) - v_bar^2(r), 0)
```

where `v_bar` is computed from gas, disk, and bulge components using
mass-to-light factors.

**External baselines:**

| Baseline | Role |
|----------|------|
| NFW-like residual fit | lightweight halo-shape comparison channel |
| MOND/RAR proxy | acceleration-relation comparison channel |
| shuffled-galaxy null | tests whether galaxy-level pairing carries signal |
| future shuffled-radius null | tests whether radial structure carries signal |

**SPARC pilot results:**

| Quantity | Result |
|----------|-------:|
| SPARC galaxies | `175` |
| radius rows | `3391` |
| mean residual fraction | `0.5592` |
| mean `D_struct` | `0.1652` |
| mean `R_rob` | `0.5814` |
| mean `V` | `0.4661` |
| Spearman `D_struct` vs NFW-like mismatch | `-0.6378` |
| Spearman `D_struct` vs RAR mismatch | `-0.3256` |
| mean shuffled-null `|rho|` for NFW channel | `0.0608` |
| mean shuffled-null `|rho|` for RAR channel | `0.0593` |

**Key finding:**
> Higher structural residual support anticorrelates with both NFW-like and
> RAR mismatch in the SPARC pilot.  The real galaxy pairing is much stronger
> than shuffled-galaxy nulls, indicating that the diagnostic carries
> nontrivial SPARC-level signal.

Here smaller mismatch means lower NFW-like residual-shape error or lower RAR
velocity-residual score under simplified comparison channels.

**Reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/sparc_structural_residual_pilot_paper61/` | Paper-61 SPARC structural residual pilot, SPARC data preparation, CSV/JSON outputs, and plots |

Prepare SPARC data:

```bash
cd experiments/sparc_structural_residual_pilot_paper61
python3 prepare_sparc_rotmod.py
```

Run SPARC pilot:

```bash
cd experiments/sparc_structural_residual_pilot_paper61
python3 sparc_structural_residual_pilot.py --input data/sparc_rotation_curves.csv --output outputs_sparc_real
```

**Data attribution and license note:** SPARC data were obtained from the
Zenodo record `10.5281/zenodo.16284118`, which lists CC-BY-4.0.  Cite Lelli,
McGaugh, and Schombert (2016), cite the Zenodo DOI when using the Zenodo
record, respect the CC-BY-4.0 attribution requirement, include the VizieR/CDS
acknowledgment if using VizieR, and clearly distinguish original SPARC
measurements from derived MAAT artifacts.

**Scope and compliance checklist:** No endorsement by the SPARC authors,
Zenodo, VizieR, CDS, or original data providers is claimed or implied.  The
reported numbers are pilot correlations, not a dark-matter model fit.  The
NFW-like and MOND/RAR-like baselines are comparison diagnostics only, not full
professional halo or modified-dynamics fits.  Fixed
`Upsilon_disk=0.5` and `Upsilon_bulge=0.7` are reproducible pilot
simplifications, not fitted astrophysical parameters.

**Documentation PDF:** `documentation/61_SPARC_Pilot_Test_Structural_Dark_Matter_Residuals.pdf`

---

### Paper 62 — Mode-A CDCL Branching
**Mode-A CDCL Branching Pilot:**  
*Compute-to-Solution Regret against VSIDS in a Transparent Conflict-Learning Solver*

**Core idea:** Paper 59 tested MAAT v1.6 as an active DPLL branching
heuristic. The T1 pilot moves one layer closer to modern SAT practice by
inserting the Mode-A acquisition rule into a small transparent CDCL solver:

```text
A(c) = prog(c) - tau * h_MAAT(c)
```

The benchmark measures paired compute-to-solution regret against VSIDS under
fixed seeds, fixed beta/tau/horizon, shared propagation, and shared conflict
learning.

**Core results:**

| Quantity | Result |
|----------|--------|
| CNF instances | 99 |
| Policies | VSIDS, progress-only, MAAT-only, Mode A, random |
| Conflict budget | 3500 |
| Mode-A median regret vs VSIDS | 0.000 |
| Mode-A mean regret vs VSIDS | +7.4619 |
| Mode-A wins/losses/ties vs VSIDS | 44 / 47 / 8 |
| Main positive signal | Mode A improves median regret on modular SAT |
| Main negative signal | Mode A does not beat progress-only overall |

**Interpretation:** Paper 62 is the first policy-level test of the v1.6 search
geometry layer. The result is mixed and useful: Mode A is an active branching
signal, but the structural penalty does not yet earn compute globally. The
family-level pattern supports a sector-/scale-dependent view of structural
search rather than a universal scalar-heuristic claim.

**Reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/sat_cdcl_mode_a_t1/` | transparent CDCL solver, fixed Mode-A parameters, generated SAT families, paired regret CSV/JSON outputs, and figures |

Run:

```bash
cd experiments/sat_cdcl_mode_a_t1
python3 sat_cdcl_mode_a_t1.py
```

**Data attribution and license note:** No external CNF dataset is redistributed.
Instances are generated locally from standard synthetic SAT-family generators.
CSV/JSON/PNG outputs are derived reproducibility artifacts.

**Documentation PDF:** `documentation/62_Mode_A_CDCL_Branching_Pilot.pdf`

---

### Paper 63 — Structure-Gated Mode-A CDCL Branching
**Structure-Gated Mode-A CDCL Branching:**  
*Instance-Level Triage of Structural Search in MAAT v1.6-T*

**Core idea:** Paper 62 showed that global `h=0` Mode A is family-sensitive:
it can help in structure-rich cases but also creates concentrated regret in
families where structural steering is misleading. Paper 63 tests the next
v1.6-T hypothesis: gate the structural penalty by an instance-level diagnostic
before applying it.

**Gated acquisition:**

```text
A_g(c) = prog(c) - g(F) * tau * h_MAAT(c)
```

where `g(F)` is computed once at the root from clause-variable graph
diagnostics: co-occurrence clustering, V-mode spread, degree-tail pressure,
and a symmetry penalty.  The gate uses no family labels.

**Core results:**

| Quantity | Result |
|----------|-------:|
| CNF instances | `99` |
| Policies | VSIDS, progress-only, MAAT-only, Mode A, gated Mode A, random |
| Gated Mode-A mean regret vs VSIDS | `+7.1023` |
| Global Mode-A mean regret vs VSIDS | `+7.6727` |
| Gated vs global Mode-A mean regret | `-0.5703` |
| Gated vs global Mode-A median regret | `-0.0080` |
| Gated vs global wins/losses/ties | `50 / 45 / 4` |
| Pigeonhole median gate | `0.0000` |
| Random 3-SAT median gate | `0.5810` |

**Key finding:**
> Structure-gating reduces mean regret relative to global Mode A and correctly
> suppresses the structural penalty on pigeonhole formulas.  However, the
> current hand-built gate remains too permissive on random and XOR-like
> families and still does not beat progress-only globally.

**Scientific status:** Partial repair, not a victory.  Paper 63 supports
structure-gated triage as the right direction, but shows that the current gate
is too coarse.  The next step should test learned or predeclared multi-feature
gates against stronger CDCL baselines and harder external SAT instances.

**Reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/sat_cdcl_structure_gated_mode_a_paper63/` | Paper-63 transparent CDCL benchmark, structure-gate diagnostics, paired regret CSV/JSON outputs, and figures |

Run:

```bash
cd experiments/sat_cdcl_structure_gated_mode_a_paper63
python3 sat_cdcl_structure_gated_mode_a_paper63.py
```

**Data attribution and license note:** No external CNF dataset is redistributed.
Instances are generated locally from standard synthetic SAT-family generators.
CSV/JSON/PNG outputs are derived reproducibility artifacts.

**Documentation PDF:** `documentation/63_Structure_Gated_Mode_A_CDCL_Branching.pdf`

---

### Paper 64 — Correlated Stationarity for Two-Qubit Entanglement Robustness
**Correlated Stationarity for Two-Qubit Entanglement Robustness:**  
*A Two-Qubit Pointer-Selection Benchmark*

**Core idea:** Extends the one-qubit pointer-state programme of Papers 51 and
58 to two qubits. The benchmark asks which entangled states survive open-system
Lindblad evolution, and whether raw local population balance should be replaced
by correlated generator stationarity when the protected structure is
entanglement rather than a one-qubit pointer axis.

**Simulation setup:**

| Component | Description |
|----------|-------------|
| System | two-qubit density matrix |
| Dynamics | Lindblad master equation integrated by RK4 |
| Initial-state families | product-z, Bell-phi, Bell-psi, partial-phi, partial-psi, Haar-entangled |
| Channel families | local z-dephasing, collective z-dephasing, pair ZZ dephasing, amplitude damping, thermal mixed, depolarizing |
| Instances | `1800` sampled trajectories |
| Target | final concurrence + concurrence retention + mutual-information retention + final purity |

**Key definitions:**

```text
B_pop = (1 - |<Z1>|) (1 - |<Z2>|)

B_corr-stat =
1 / (1 + (0.65 D_local + 1.15 D_corr
           + 0.90 D_concurrence + 0.25 D_purity) / 2.2)
```

**Core results:**

| Feature set | RF R2 | RF Spearman | Leave-channel R2 |
|-------------|------:|------------:|-----------------:|
| defect-field stationary | `0.9160` | `0.9477` | `0.6033` |
| short-time generator baseline | `0.8060` | `0.8773` | `0.7800` |
| scalar correlated stationarity | `0.7831` | `0.8154` | `0.7286` |
| einselection baseline | `0.7574` | `0.8672` | `-0.0110` |
| scalar population | `0.7232` | `0.7899` | `0.6179` |
| rates only | `0.4871` | `0.8520` | `0.0357` |
| shuffled defect null | `-0.0493` | `-0.0048` | `-0.2176` |

**Key finding:**
> Two-qubit entanglement robustness is better aligned with correlated
> generator stationarity than with raw local population balance.  The full
> stationarity-sensitive defect field performs best in ordinary
> cross-validation, but the compact short-time generator baseline remains best
> under leave-channel-family-out transfer.

**Scientific status:** Classical simulation benchmark only. This is not a
quantum-hardware experiment, not a collapse model, not a derivation of quantum
measurement, and not a replacement for decoherence/einselection theory.

**Reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/two_qubit_pointer_paper64/` | Paper-64 two-qubit Lindblad simulation, entanglement-robustness target, model comparisons, CSV/JSON outputs, and figures |

Run:

```bash
cd experiments/two_qubit_pointer_paper64
python3 two_qubit_pointer_paper64.py
```

**Data attribution and license note:** No external quantum dataset is
redistributed. All trajectories are generated locally from the script using a
fixed random seed. CSV/JSON/PNG outputs are derived reproducibility artifacts.

**Documentation PDF:** `documentation/64_Two_Qubit_Pointer_Selection.pdf`

---

### Paper 65 — SPARC II
**SPARC II:**  
*Shuffled-Radius Nulls, Halo-Catalogue Cross-Check Interface, and Morphology Splits*

**Core idea:** Extends the Paper-61 SPARC pilot with a stronger
within-galaxy shuffled-radius null, baryonic morphology/proxy splits, and an
optional Li-et-al.-style published halo-catalogue cross-check interface.
The goal is not to fit a dark-matter model, but to test whether the
structural residual keeps radius-resolved information under a null that
preserves each galaxy's one-point residual distribution while destroying
radial pairing.

**Core results:**

| Quantity | Result |
|----------|-------:|
| SPARC galaxies | `175` |
| Radius rows | `3391` |
| Radius-shuffle nulls | `500` |
| pooled rho(`D_struct`, NFW-like residual) | `0.3930` |
| NFW-like shuffled-radius `|rho|_95` | `0.0355` |
| NFW-like exceeds shuffled null? | `yes` |
| pooled rho(`D_struct`, RAR residual) | `0.0123` |
| RAR shuffled-radius `|rho|_95` | `0.0336` |
| RAR exceeds shuffled null? | `no` |
| median galaxywise rho(`D_struct`, NFW-like residual) | `0.3000` |
| published halo-catalogue status | `not run; no local catalogue supplied` |

**Key finding:**
> The structural residual retains radius-resolved SPARC information in the
> NFW-like comparison channel under a shuffled-radius null, but the same does
> not hold for the RAR channel. The signal is therefore channel- and
> population-dependent, not a dark-matter solution.

**Scientific status:** Robustness follow-up to Paper 61. The optional
published-halo-catalogue interface is implemented, but the present repository
run reports `not_run_no_local_catalogue` because no local Li-et-al. halo-fit
table is supplied. The morphology splits are baryonic proxy classes derived
from available SPARC components, not official Hubble types.

**Reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/sparc_ii_paper65/` | Paper-65 shuffled-radius nulls, proxy morphology splits, optional halo-catalogue cross-check interface, generated CSV/JSON outputs, and figures |

Run:

```bash
cd experiments/sparc_ii_paper65
python3 sparc_ii_paper65.py
```

Optional published-halo-catalogue interface:

```bash
python3 sparc_ii_paper65.py --halo-catalog path/to/li_halo_catalog.csv
```

**Data attribution and license note:** Paper 65 reuses the SPARC-derived
Paper-61 inputs. SPARC data are external scientific data released via Zenodo
under CC-BY-4.0 and should be cited to the SPARC source and Lelli et al. when
reused. If a published halo-fit catalogue is supplied locally, cite the
original catalogue source, e.g. Li et al. for SPARC halo fits. The exploratory
`Q_bar` comparison uses a collaborator-supplied derived HGD-GSR descriptor from
Ali Alhawarat; it is not a SPARC measurement and not a MAAT-derived quantity.
Unless broader redistribution terms are separately supplied by the HGD-GSR
author, treat the `Q_bar` table as a neutral comparison input for this pilot and
cite/acknowledge Ali Alhawarat and HGD-GSR when discussing or reusing it.
Repository CSV/JSON/PNG files are derived reproducibility artifacts only. No
endorsement by the SPARC team, Lelli et al., Li et al., HGD-GSR, Ali Alhawarat,
or any catalogue authors is implied.

**Documentation PDF:** `documentation/65_SPARC_II_Shuffled_Radius_Nulls.pdf`

---

### Collaboration Pilot — SPARC MAAT x HGD-GSR Cross-Framework Test
**SPARC Cross-Framework Structural Signal Test:**  
*Comparing MAAT `D_struct` with externally supplied HGD-GSR `Q_bar` values*

**Core idea:** Provides a collaboration-ready pilot for testing whether the
MAAT structural residual signal from Paper 65 and the HGD-GSR `Q_bar` signal
are related, independent, complementary, or redundant on the same SPARC galaxy
basis.

The repository includes a collaborator-supplied `Q_bar` table for the neutral
pilot comparison. If that CSV is removed or replaced, the script can still write
a 175-galaxy template and a data contract. With `data/hgd_gsr_qbar.csv`
present, the pipeline performs the join, correlation tests, bootstrap
intervals, and cross-validated predictive comparisons against the Paper-65
mismatch targets.

**Reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/sparc_hgd_gsr_cross_framework_pilot/` | Collaboration-ready MAAT x HGD-GSR SPARC pilot; writes Q-bar template if external data are absent |

Run:

```bash
cd experiments/sparc_hgd_gsr_cross_framework_pilot
python3 sparc_hgd_gsr_cross_framework_pilot.py
```

Expected external input:

```text
data/hgd_gsr_qbar.csv
columns: galaxy,Q_bar
```

**Scientific status:** Exploratory cross-framework pilot. The present run uses
real collaborator-supplied HGD-GSR `Q_bar` values for all 175 SPARC galaxies.
It is not a confirmatory validation of either framework and should be read as a
neutral relation/complementarity test.

**Data attribution and license note:** MAAT/SPARC-derived inputs follow the
Paper-65 SPARC attribution and CC-BY-4.0 notes. HGD-GSR `Q_bar` values are
external collaborator-supplied derived descriptors from Ali Alhawarat's HGD-GSR
framework. They are included only for this neutral cross-framework pilot unless
broader redistribution terms are separately supplied by the HGD-GSR author.
Cite/acknowledge Ali Alhawarat and HGD-GSR when discussing or reusing the
`Q_bar` comparison. No endorsement by SPARC, HGD-GSR, Ali Alhawarat, Zenodo, or
any original data provider is implied.

---

### Extra Phenomenological Paper — BKM-Aware MAAT Diagnostics for Navier-Stokes
**BKM-Aware MAAT Diagnostics for Navier--Stokes:**  
*A Phenomenological Structural-Action Test on Taylor--Green Vortices*

**Core idea:** Operationalises the symbolic expression

```text
ToE_MAAT = integral [ (H+B+S+V+R) * Z * cascade_resistance ]
                    / [ DeltaE + DeltaQ + D0 ] dt
```

as a reproducible 3D incompressible Navier--Stokes diagnostic. The symbolic
`cascade_resistance` term is implemented as cascade/tail resistance coupled to
vorticity-stretching pressure, while `D0=11` is a fixed dimensionless baseline
penalty.

**Scientific status:** Phenomenological diagnostic note, not a proof of
Navier--Stokes regularity.
The purpose is to test whether MAAT v1.6 structural coordinates can act as
BKM-aware warning coordinates along a Taylor-Green vortex trajectory.

**Core results:**

| Scenario | min `R_rob` | max warning | final `ToE_MAAT` action | warning vs stretching |
|----------|------------:|------------:|------------------------:|----------------------:|
| moderate viscosity | `0.9029` | `0.1677` | `0.9567` | `0.9912` |
| low viscosity / high stress | `0.7578` | `0.5692` | `0.6786` | `1.0000` |

**Key finding:**
> The MAAT warning does not simply duplicate the vorticity infinity norm. In
> this toy 3D setting it aligns much more strongly with vortex-stretching
> pressure, while the integrated `ToE_MAAT` action is larger for the more
> coherent moderate-viscosity trajectory.

**Structural value added:** The diagnostic points to a genuinely
three-dimensional stress channel instead of merely renaming a standard
vorticity monitor.

**Interpretation guardrail:** The observed correlations should not be
interpreted as evidence for singularity prediction. They indicate only that the
structural diagnostic appears sensitive to the dynamically relevant stretching
channel in these controlled toy trajectories.

**Scope note:** Only two Taylor-Green trajectories are tested. The high
correlations are proof-of-concept evidence for the diagnostic channel, not broad
turbulence evidence.

**Reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/navier_stokes_maat_regularization/` | 3D Taylor-Green Navier--Stokes solver, MAAT/BKM diagnostics, symbolic action outputs, and plots |

Run:

```bash
cd experiments/navier_stokes_maat_regularization
python3 navier_stokes_maat_regularization.py
```

**Data attribution and license note:** No external fluid dataset is
redistributed. Velocity fields are generated locally from Taylor-Green initial
conditions. CSV/JSON/PNG files are derived reproducibility artifacts.

**Documentation PDF:** `documentation/Phenomenological_MAAT_Navier_Stokes_Regularity_Functional.pdf`

---

### Extra Phenomenological Paper — Dark Matter Residuals as Structural Coherence Diagnostics
**Dark Matter Residuals as Structural Coherence Diagnostics in MAAT v1.6:**  
*A Phenomenological Interpretation of Galactic Rotation Residuals*

**Core idea:** Interprets the galactic dark-matter signal as a structural
coherence residual, not as a claim against particle dark matter.  In MAAT v1.6,
the effective residual

```text
v_dark^2(r) = max(v_obs^2(r) - v_baryon^2(r), 0)
```

marks where the visible baryonic state fails to close the observed
gravitational trajectory.  The paper asks how this residual can be read as a
structural support channel in gravitational search geometry.

**Scientific status:** Phenomenological diagnostic note, not a galaxy fit, not
a dark-matter particle model, and not a replacement for Lambda-CDM or cold dark
matter.  The result is compatible with particle dark matter: a particle halo
would be one possible microscopic carrier of the residual channel.

**Important caveat:** The benchmark is intentionally weak and synthetic.  It
does not add a new empirical dark-matter result beyond the familiar
rotation-curve residual.  Its value is to make the residual explicit as a
MAAT v1.6 structural diagnostic and to define what a nontrivial real-data test
would require.

**Core diagnostic:**

```text
D_struct(r) = f_res(r) * R_rob(r) * V(r)
```

where `f_res` is the inferred dark residual fraction, `R_rob` is the MAAT
robustness closure, and `V` measures connectedness between baryonic structure
and residual-density shape.

**Core results:**

| Quantity | Result |
|----------|-------:|
| max observed velocity | `204.9972 km/s` |
| outer mean residual fraction | `0.8802` |
| outer mean robustness | `0.7806` |
| transition peak structural dark support | `0.3486` |
| peak-support radius | `20.3782 kpc` |
| mean `H` | `1.0000` |
| mean `B` | `0.9060` |
| mean `S_eff` | `0.7506` |
| mean `V` | `0.5919` |
| mean `R_rob` | `0.7611` |

**Key finding:**
> Dark matter can be described, at the phenomenological level, as the
> gravitationally inferred structural residual required for galactic dynamical
> coherence.  This statement is structural, not microscopic.

**Interpretation guardrail:** The paper does not claim that dark matter is
``only structure''.  It says that any successful microscopic dark-sector
theory must carry the residual gravitational structure while preserving
coherence, balance, connectedness, and robustness.

**Required next tests:** Real rotation-curve samples, comparison against halo
and MOND-like baselines, lensing residuals, cluster systems, CMB constraints,
and large-scale-structure growth.  Without those tests this remains a
conceptual translation layer, not a competing dark-matter model.

**SPARC/NFW/MOND upgrade protocol:** The paper now specifies a concrete
real-data path: use SPARC rotation curves and baryonic mass models, compare the
MAAT structural residual against NFW halo fits, MOND/RAR baselines, and
published SPARC multi-halo catalogues, and predeclare shuffled-radius and
shuffled-galaxy null tests.

**Future data attribution:** This repository version does not redistribute
SPARC data.  If SPARC is added later, cite Lelli, McGaugh, and Schombert
(2016), the Zenodo DOI `10.5281/zenodo.16284118` when using the Zenodo record,
respect its CC-BY-4.0 attribution requirement, and include the standard
VizieR/CDS acknowledgment if the VizieR route is used.  Derived MAAT CSV/PNG
artifacts must be clearly separated from original SPARC measurements.

**Reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/dark_matter_structural_coherence/` | Synthetic rotation-curve toy model, inferred dark residual, MAAT support diagnostics, CSV/JSON outputs, and plots |

Run:

```bash
cd experiments/dark_matter_structural_coherence
python3 dark_matter_structural_coherence.py
```

**Data attribution and license note:** No external astronomical dataset is
redistributed.  The rotation curve, baryonic profile, residual profile, CSV
tables, JSON summary, and figures are synthetic reproducibility artifacts
generated locally by the script.

**Documentation PDF:** `documentation/Phenomenological_MAAT_Dark_Matter_Structural_Coherence.pdf`

---

### Extra Phenomenological Paper — Global Coherence versus Local Projection Defects
**Global Coherence versus Local Projection Defects:**  
*Rare B-Decay Anomaly Diagnostics in MAAT v1.6*

**Core idea:** Treats rare electroweak-penguin B-meson tensions as a
structural anomaly-classification problem. The paper does **not** claim new
physics, does **not** fit Wilson coefficients, and does **not** redistribute
LHCb/CMS/Belle/HEPData measurements. It proposes a MAAT v1.6 pre-fit
diagnostic layer for distinguishing coherent global residual modes from
localized hadronic or projection-like residuals.

**Diagnostic question:**

```text
Is the anomaly geometry closer to a coherent global coefficient-space defect
or to a local q^2-dependent long-distance hadronic projection defect?
```

Here a projection defect means a charm-loop / charming-penguin-like
long-distance hadronic residual projected onto a small set of measured angular
observables and q² bins, not a new fundamental field.

**Toy benchmark results:**

| Scenario | Globality | Locality |
|----------|-----------|----------|
| Coherent Wilson-shift toy pattern | `0.9161` | `0.1586` |
| Localized charm-like projection toy pattern | `0.0139` | `0.6894` |
| Mixed coherent plus local toy pattern | `0.6594` | `0.4178` |

**Bin-shuffled null control:** The toy code additionally preserves residual
amplitudes while shuffling q²-bin coherence. The shuffled null raises the mean
diagnostic cost by `+0.1823` for the coherent Wilson-shift toy pattern,
`+2.0578` for the localized charm-like pattern, and `+1.8338` for the mixed
pattern.

**Scientific status:** Synthetic diagnostic protocol only. The benchmark uses
generated toy residual patterns and is not a particle-physics data analysis.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/rare_b_decay_anomaly_diagnostics/` | Synthetic q^2 residual-pattern benchmark, MAAT v1.6 diagnostic supports, CSV/JSON outputs, and plots |

Run:

```bash
cd experiments/rare_b_decay_anomaly_diagnostics
python3 rare_b_decay_anomaly_diagnostics.py
```

**Data attribution and license note:** No external experimental data are
redistributed. No CERN, LHCb, Nature, Physical Review Letters, or HEPData
figures, tables, images, or long text passages are reused. Any future use of
LHCb, CMS, Belle, Belle II, or HEPData tables, covariance matrices, plots, or
figures must cite the relevant collaboration papers and comply with the
corresponding data and figure licenses. No endorsement by CERN, LHCb, CMS,
Belle, Belle II, HEPData, Physical Review Letters, Nature, or the cited authors
is implied.

**Documentation PDF:** `documentation/Phenomenological_MAAT_Rare_B_Decay_Anomaly_Diagnostics.pdf`

---

### Extra Methodology Paper — Ranking Mathematical Candidates
**Ranking Mathematical Candidates:**  
*A Structural Selection Layer for Discovery and Counterexample Search*

**Core idea:** Uses the 2026 OpenAI unit-distance breakthrough as motivation
for a general MAAT-style discovery layer. It proposes that MAAT can help rank
mathematical candidates before proof by scoring logical coherence, constraint
balance, constructive novelty, cross-domain connectivity, and verification
robustness. The paper does **not** reproduce or verify the OpenAI proof.

**Discovery functional:**

```text
F_MAAT^disc[X]
= - sum_{a in {H,B,S,V}} lambda_a log(epsilon + Gamma_a[X])
  - lambda_cl log(epsilon + R_rob[X])
```

where lower score means higher structural support for further proof effort,
not mathematical truth.

**Key workflow:**

```text
generate broadly
-> select structurally
-> prove rigorously
-> verify independently
```

**Key finding:**
> MAAT does not prove the theorem; it ranks which mathematical structures are
> worth trying to prove.

**Bridge-detection update:** The paper expands the connectivity sector into
four bridge quantities: reach, preservation, yield, and translation cost. This
turns the OpenAI-style lesson into an operational question: which external
structure preserves the conjecture's constraints while changing the growth
mechanism?

**Toy example:** The paper includes a minimal extremal-graph benchmark
template. Candidate graph families are ranked by forbidden-condition
violations, edge-target balance, construction degrees of freedom,
cross-domain bridge support, and perturbation robustness, then compared
against known-family and random-search baselines. A small-n implementation
could enumerate labelled or isomorphism-reduced candidates for `n <= 10`,
enumerate restricted construction families up to `n <= 12`, sample larger
cases up to about `n <= 20`, and use subgraph counting, single-edge
perturbation tests, random edge-sampling, graph-feature distance from known
families, and edge-baseline improvements.

**Scientific status:** Conceptual methodology paper. It is a search-prior and
robustness-audit proposal for AI-assisted mathematics, not a theorem prover,
not a verification system, and not a claim to reproduce the OpenAI
unit-distance construction.

**Documentation PDF:** `documentation/Structural_Selection_for_Mathematical_Discovery.pdf`

---

### Extra Mathematical Diagnostic Paper — Riemann Zeta Critical-Line Diagnostic
**A Structural Selection Diagnostic on the Riemann Zeta Critical Line**  
*A purely diagnostic, non-proof toy experiment with controls, ablations, and pair-correlation tests*

**Core idea:** Applies a simple heuristic structural diagnostic to known
nontrivial Riemann zeta zeros. The test asks whether known zero ordinates have
lower diagnostic cost on the critical line than under small off-line shifts,
and how much of that behaviour is built into the scoring rule itself. This is
a mathematical toy benchmark, not evidence for or a proof of the Riemann
Hypothesis. The pair-correlation component compares normalised zeta-zero
spacings with a GUE spacing proxy and Poisson/uniform controls.

**Diagnostic terms:**

| Term | Definition |
|------|------------|
| `H` | `1 / (1 + |zeta(sigma + i t)|)` |
| `B` | `1 / (1 + alpha |sigma - 1/2|)` |
| `R` | symmetric real-direction zeta-response robustness |

**Core results:**

| Set | Config | Best mean delta | Mean F at best |
|-----|--------|----------------:|---------------:|
| known zeros | HBR | `0.000` | `0.00818` |
| known zeros | HR, no explicit B | `0.000` | `0.00818` |
| known zeros | BR, no H | `0.000` | `0.00818` |
| jittered zeros | HBR | `0.000` | `0.16544` |
| random critical-line points | HBR | `0.000` | `0.86054` |
| random critical-line points | HR, no explicit B | `0.050` | `0.83426` |

Pair-correlation diagnostic:

| Series | Spacings | L1 to GUE | V_pair | F_pair |
|--------|---------:|----------:|-------:|-------:|
| jittered zeta spacings | 199 | `0.2059` | `0.9844` | `0.0157` |
| zeta normalised spacings | 199 | `0.2938` | `0.9797` | `0.0205` |
| uniform-points control | 199 | `0.8103` | `0.9573` | `0.0437` |
| poisson control | 199 | `0.8225` | `0.9559` | `0.0451` |

**Key finding:**
> Known zeta-zero locations exhibit a low-cost minimum on the critical line
> under the proposed diagnostic score. However, the score partly
> builds in critical-line preference through the balance term and directly
> rewards known zeros through the zeta-null term. The result is therefore a
> diagnostic benchmark, not an RH result.

**Scientific status:** Mathematical toy diagnostic with controls and
ablations. The full HBR and no-balance HR diagnostics select `delta=0` for
known zeros across the tested parameter sweep, while random critical-line
controls reveal how much of the critical-line preference is score-induced.
The pair-correlation diagnostic is conceptually closer to spectral-universality
tests than to direct zero-location tests, but it remains finite-sample and uses
a GUE/Wigner-surmise proxy rather than a theorem.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/zeta_critical_line_selection/` | zeta critical-line diagnostic script, zero/control tables, CSV/JSON outputs, and figures |

**Data attribution and license note:** No external dataset file is
redistributed. The script computes zeta zeros using `mpmath.zetazero` and
evaluates the zeta function directly. Generated CSV/JSON/PNG files are derived
reproducibility artifacts.

**Documentation PDF:** `documentation/Structural_Selection_on_the_Critical_Line.pdf`

---

### Extra Experiment Paper — Societal Critical Coherence Index
**Societal Critical Coherence Index:**
*A Toy Framework for Structural Stress and Constructive Transformation*

**Core idea:** Transfers the CCI/MAAT diagnostic logic to a deliberately
synthetic social-system toy model. The goal is not to rank real people,
parties, countries, organisations, or communities. The paper tests whether a
structural index can separate raw activity, passive coherence, robustness,
and constructive transformation.

**Master definitions:**

```text
S_eff(A)    = exp[-0.5 ((A_raw - A_star) / sigma_A)^2]
R_resp,soc  = (H B V)^(1/3)
R_rob,soc   = min(R_resp,soc, (H B S_eff V)^(1/4))
R_sig,soc   = R_resp,soc^(1-alpha) S_eff^alpha
CCI_soc     = A_raw / (H + B + V + R_rob,soc + epsilon)
ASI_soc     = R_sig,soc / (1 + CCI_soc)
```

**Core results:**

| Quantity | Result |
|----------|--------|
| Synthetic archetypes | `6` |
| Random ensemble size | `5000` |
| Fixed seed | `44` |
| Activity optimum | `A_star = 1.0` |
| Activity width | `sigma_A = 0.30` |
| Significance exponent | `alpha = 0.45` |
| Best `ASI_soc` archetype | `creative_democratic_renewal` |
| Best `ASI_soc` | `0.6940` |
| Highest `R_sig,soc` archetype | `creative_democratic_renewal` |
| Highest `R_sig,soc` | `0.9055` |
| Highest `CCI_soc` archetype | `polarized_mobilization` |
| Highest `CCI_soc` | `1.0177` |
| Ensemble mean `CCI_soc` | `0.4530` |
| Ensemble mean `ASI_soc` | `0.2813` |

**Key finding:**
> Low structural stress can still be stagnant, and high social activity can
> still be polarising or fragmenting. Constructive transformation appears only
> when controlled activity is supported by coherence, balance, and
> connectedness.

**Ethical status:** This is a conceptual diagnostics toy framework. It is not
an instrument for political scoring, population ranking, surveillance,
persuasion, or normative judgement of real communities.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/societal_cci/` | Synthetic societal CCI script, archetype table, random ensemble, JSON/CSV outputs, and figures |

**Documentation PDF:** `documentation/Societal_Critical_Coherence_Index.pdf`

---


### Extra Phenomenological Paper — Structural Selection in SO(10)-Motivated Unified Field Theories
**Structural Selection in SO(10)-Motivated Unified Field Theories:**
*A Phenomenological MAAT Layer for Gauge and Yukawa Benchmarks*

**Core idea:** Treats SO(10)-motivated grand unification as the dynamical
base and MAAT v1.2.1 as a structural-selection layer over gauge and Yukawa
parameter regions. The paper does not construct a complete SO(10) model,
prove precision unification, or derive the Standard Model spectrum from first
principles. It tests whether structurally coherent GUT-scale regions can be
ranked reproducibly.

**Core results:**

| Quantity | Result |
|----------|--------|
| Gauge benchmark | one-loop SM running, no threshold corrections |
| Best `alpha_GUT` | `0.022676468338` |
| Best `M_GUT` | `6.639198e15 GeV` |
| Gauge `SM chi2` | `100.789861` |
| Gauge `R_rob` | `0.60301490` |
| Yukawa benchmark | SO(10)-motivated `b-tau` boundary |
| Yukawa `M_GUT` | `1.858006e16 GeV` |
| `Delta_b` | `0.0506451697` |
| `chi2_yukawa` | `0.00024341` |
| Yukawa `R_rob` | `0.99916912` |

**Key finding:**
> The one-loop gauge benchmark locates a plausible GUT-scale region but fails
> precision unification, while the third-generation Yukawa benchmark shows
> that an SO(10)-motivated `b-tau` boundary can be compatible with low-energy
> effective masses using a modest bottom-threshold correction. The result is a
> structural compatibility benchmark, not a first-principles derivation.

**Scientific status:** This is an extra phenomenological bridge paper. It is
not a complete SO(10) construction and not a precision GUT fit.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/maat_so10_structural_selection/` | SO(10)-motivated gauge and Yukawa benchmarks, JSON summaries, and selection-landscape figure |

**Data attribution and license note:** The electroweak-scale inputs use
standard reference values for `M_Z`, gauge couplings, and third-generation
fermion masses cited in the paper to the Particle Data Group. Repository
JSON/PNG/PDF artifacts are derived analysis outputs only. No endorsement by
the Particle Data Group or any external collaboration is implied.

**Documentation PDF:** `documentation/Structural_Selection_in_SO10_Unified_Field_Theories.pdf`

---

### Extra Phenomenological Paper — MAAT--SO(10) Extended Benchmark Framework
**MAAT--SO(10) Extended Benchmark Framework:**  
*From Phenomenological Bridge to Response-Closed Multi-Scale GUT Selection*

**Core idea:** Extends the first SO(10)-motivated MAAT bridge paper into a
versioned benchmark programme. The paper does not claim a complete SO(10)
construction, precision unification, or a derivation of the Standard Model
spectrum. Instead, it converts the limitations of the first bridge into a
testable roadmap: two-loop and threshold-aware gauge diagnostics, expanded
Yukawa and neutrino tests, response-closed structural weights, v1.2.1
emergent robustness, multi-scale defect-field selection, RG
generator-stationarity, CCI-style regime diagnostics, and explicit
falsification baselines.

**Starting evidence from the first SO(10) bridge:**

| Sector | Result | Interpretation |
|--------|--------|----------------|
| Gauge one-loop | `M_GUT ~= 6.64e15 GeV`, `chi2 ~= 100.79`, `R_rob ~= 0.603` | plausible scale, no precision unification |
| Yukawa `b-tau` | `M_GUT ~= 1.86e16 GeV`, `Delta_b ~= 0.0506`, `R_rob ~= 0.999` | coherent third-generation compatibility valley |

**Key finding:**
> The SO(10) bridge should not be inflated into a result claim. Its strongest
> next form is a falsifiable benchmark ladder: MAAT must show added predictive
> information beyond chi-square-only, threshold-only, complexity-only,
> random-weight, and shuffled-defect baselines.

**Scientific status:** Extended benchmark specification. It is a roadmap and
formalisation paper, not a new numerical fit and not a complete GUT model.

**Reproducibility:** Uses the existing SO(10) experiment outputs:

| Folder | Role |
|--------|------|
| `experiments/maat_so10_structural_selection/` | Existing gauge and Yukawa benchmark scripts, JSON summaries, and Yukawa selection-landscape figure |

**Documentation PDF:** `documentation/MAAT_SO10_Extended_Benchmark_Framework.pdf`

---

### Extra Phenomenological Paper — Structural Selection in the String Landscape
**Structural Selection in the String Landscape:**
*A MAAT-Based Phenomenological Framework for Vacuum Ranking*

**Core idea:** Treats string theory as the dynamical base and MAAT as a
structural ranking layer over candidate string backgrounds. The paper does not
replace string theory, derive a unique string measure, construct a Standard
Model vacuum, or solve the landscape problem. It proposes a reproducible
multi-sector ranking architecture over admissible or approximately admissible
backgrounds.

**Master formula:**

```
F_MAAT^landscape[X] =
- sum_a lambda_a log(epsilon + Gamma_a[X]),
Gamma_a[X] = 1 / (1 + d_a[X])
```

**Benchmark evidence:** Existing `string_landscape_selection/` scripts test
10D-inspired tadpole closure, KKLT-style bridges, period-controlled flux
backgrounds, and backreaction / Standard-Model-sector proxies. Across these
benchmarks, structural ranking differs from energy-only ordering and can reduce
tadpole, stability, or phenomenological obstruction in the selected subsets.

**Scripts and reproducibility:**

| Folder | Role |
|--------|------|
| `experiments/string_landscape_selection/` | string-landscape toy and bridge benchmarks, result JSON files, and generated plots |

**Documentation PDF:** `documentation/Structural_Selection_in_the_String_Landscape_MAAT_Framework.pdf`

---

### Extra Phenomenological Paper — MAAT String Selection II
**MAAT String Selection II:**
*Structural Path Selection in a Reduced KKLT Landscape*

**Core idea:** Extends the string-landscape vacuum-ranking benchmark from
static endpoint selection to transition-path selection. Candidate KKLT bridge
vacua become nodes in a reduced landscape graph, and near-neighbour
transitions become directed edges with a structural path action.

**Path action:**

```text
A_ij =
  B_ij
  + lambda_F max(F_MAAT[B_j] - F_MAAT[B_i], 0)
  + lambda_c C_ij
  + lambda_Q Q_ij
  + lambda_R R_loss,ij

P(B_i -> B_j) ∝ Gamma_ij exp[-A_ij],
Gamma_ij = exp[-B_ij]
```

**Core results:**

| Quantity | Result |
|----------|--------|
| Minima nodes | `77` |
| Directed near-neighbour edges | `462` |
| AdS minima | `45` |
| dS minima | `32` |
| Spearman `abs_energy` vs `F_bridge` among minima | `0.1271` |
| Spearman `Delta_Q_10` vs `F_bridge` among minima | `0.9486` |
| Fraction structurally downhill edges | `0.539` |
| Top-20 structural/energy edge overlap | `0 / 20` |
| Best-energy to best-structural path | `60 -> 22 -> 30` |
| Structural path cost | `0.907889` |
| Energy-action cost to same target | `1.719058` |

**Key finding:**
> MAAT extends static string vacuum ranking into structural landscape flow:
> preferred endpoints and preferred transition paths need not be determined by
> energy alone.

**Scientific status:** This is a phenomenological transition-graph benchmark.
It is not a Coleman-De Luccia tunnelling calculation, not a microscopic string
decay-rate derivation, and not a complete landscape dynamics model.

**Scripts and reproducibility:**

Repository URL:

```text
https://github.com/Chris4081/structural-selection-principle/tree/main/experiments/string_landscape_path_selection
```

| Folder | Role |
|--------|------|
| `experiments/string_landscape_path_selection/` | MAAT String Selection II path graph, transition edges, selected paths, CSV/JSON outputs, and figures |

**Data attribution and license note:** This experiment reuses derived
synthetic KKLT-bridge candidates generated by the repository's
`string_landscape_selection` benchmark scripts. No external survey data or
proprietary dataset is redistributed. The benchmark is not endorsed by the
authors of the original string-theory literature cited in the paper.

**Documentation PDF:** `documentation/MAAT_String_Selection_II_Structural_Path_Selection.pdf`

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
├── structural_selection_phi4_protocol.py  ← Paper 22: static φ⁴ validation benchmark
├── structural_selection_phi4_plots.py     ← Paper 22: figure generation from benchmark JSON
│
├── experiments/
│   ├── natural_constants_selection/       ← Paper 26: natural-constants v1--v13 benchmark
│   ├── standard_model_bridge/             ← Paper 26: SM-like RG bridge and v11 holdout
│   ├── boundary_aware_lambda_calibration/ ← Paper 27: fused boundary-aware λ calibration
│   ├── lambda_response_closure/           ← addendum: response-theoretic λ closure
│   ├── cosmological_cci/                  ← Paper 28: cosmological CCI observable
│   ├── cosmological_cci_v03/              ← Paper 29: growth connectivity + robustness CCI
│   ├── maat_dynamic_fields_v05_v09/       ← Paper 31: response/local/gravity/FLRW tests
│   ├── maat_observable_predictions_v10/   ← Paper 31: observable projection layer
│   ├── maat_hz_chi2_paper32/              ← Paper 32: H(z) chronometer comparison and χ² scan
│   ├── maat_cci_projection_paper33/       ← Paper 33: CCI projection observable + sensitivity scan
│   ├── maat_projection_growth_paper34/    ← Paper 34: projection observable vs growth-data diagnostic
│   ├── maat_growth_perturbation_paper35/  ← Paper 35: linear growth embedding + perturbation stability
│   ├── maat_v121_observables_stability_paper37/ ← Paper 37: v1.2.1 observable proxy + stability landscape
│   ├── maat_paper38_v121_robustness_closure/ ← Paper 38: v1.2.1 closure in linear growth
│   ├── maat_paper39_observable_growth_signature/ ← Paper 39: v1.2.1 observable growth signature proxy
│   ├── maat_paper40_structural_signature_test/ ← Paper 40: v1.2.1 CCI residual signature test
│   ├── maat_paper42_blind_projection_test/ ← Paper 42: response-derived blind projection test
│   ├── maat_so10_structural_selection/   ← extra paper: SO(10)-motivated gauge and Yukawa benchmarks
│   └── string_landscape_selection/        ← extra paper: MAAT string-landscape ranking
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
| `structural_selection_phi4_protocol.py` | 22 | static φ⁴ benchmark, vacua vs saddle, kink vs distortions, writes result JSON |
| `structural_selection_phi4_plots.py` | 22 | generates Paper 22 validation plots from `structural_selection_phi4_results.json` |
| `experiments/natural_constants_selection/naturkonstante_v7_landscape_heatmap.py` | 26 | natural-constants basin visualisation |
| `experiments/natural_constants_selection/naturkonstante_v8_gradient_flow.py` | 26 | gradient-flow attractor test in constant space |
| `experiments/natural_constants_selection/naturkonstante_v12_maxent_lambda.py` | 26 | maximum-entropy calibration of sector weights |
| `experiments/natural_constants_selection/naturkonstante_v13_maxent_sm_bridge.py` | 26 | MaxEnt-weighted constants bridge using calibrated sector weights |
| `experiments/standard_model_bridge/standard_model_rg_maat_bridge.py` | 26 | one-loop SM-like RG bridge from UV parameters to IR effective observables |
| `experiments/standard_model_bridge/standard_model_rg_maat_summary_figure.py` | 26 | generates the four-panel Paper 26 summary figure |
| `experiments/standard_model_bridge/standard_model_rg_maat_v11_holdout.py` | 26 | direct-term holdout benchmark for cross-sector predictivity |
| `experiments/boundary_aware_lambda_calibration/fit_closed_maat_lambda_v1.py` | 27 | closed boundary-aware MAAT lambda calibration over fused defect data |
| `experiments/boundary_aware_lambda_calibration/plot_closed_maat_lambda_v2.py` | 27 | generates Paper 27 lambda and defect-comparison figures |
| `experiments/lambda_response_closure/lambda_response_closure.py` | addendum | derives effective lambda weights from defect covariance and target-response geometry |
| `experiments/cosmological_cci/maat_cci_cosmology_v02.py` | 28 | generates the cosmological CCI model grid, chronometer projection, and plots |
| `experiments/cosmological_cci_v03/maat_cci_cosmology_v03_growth.py` | 29 | adds f sigma_8 growth connectivity, robustness margins, MaxEnt companion weights, and curvature transition proxy |
| `experiments/maat_dynamic_fields_v05_v09/v05_dynamic_lambda_flow/lambda_dynamic_flow_v05.py` | 31 | covariance-driven response-field dynamics |
| `experiments/maat_dynamic_fields_v05_v09/v06_local_selection_fields/local_selection_phi4_v06.py` | 31 | local selection-pressure fields in a perturbed 1D phi4 benchmark |
| `experiments/maat_dynamic_fields_v05_v09/v09_flrw_stability_scan/maat_flrw_stability_scan_v09.py` | 31 | toy scalar FLRW stability scan |
| `experiments/maat_observable_predictions_v10/maat_observable_predictions_v10.py` | 31 | observable projection: H(z), Delta H/H, w(z), Omega_MAAT, and f sigma_8 proxy |
| `experiments/maat_hz_chi2_paper32/maat_hz_data_comparison_v01.py` | 32 | fixed v0.10 MAAT comparison against Cosmic Chronometer H(z) data |
| `experiments/maat_hz_chi2_paper32/maat_hz_chi2_fit_v01.py` | 32 | two-parameter MAAT H(z) chi-square scan |
| `experiments/maat_v121_observables_stability_paper37/paper37_observables_emergent_robustness.py` | 37 | v1.2.1 baseline observable proxy with derived respect and emergent robustness |
| `experiments/maat_v121_observables_stability_paper37/paper37_stability_landscape.py` | 37 | two-parameter v1.2.1 proxy stability landscape scan |
| `experiments/maat_v121_observables_stability_paper37/sat_validation/maat_v121_sat_validation.py` | 37 | companion SAT correlation validation for the empirical discussion |
| `experiments/maat_paper38_v121_robustness_closure/maat_paper38_v121_robustness_closure.py` | 38 | v1.2.1 robustness closure in a linear-growth benchmark plus selection-field perturbation tests |
| `experiments/universal_structural_selection_benchmarks_paper49/universal_structural_selection_benchmarks_paper49.py` | 49 | aggregates existing benchmark outputs into a universality evidence registry and competitor-readiness matrix |
| `experiments/sat_structural_hardness_paper50/sat_structural_hardness_paper50.py` | 50 | synthetic random-3-SAT graph-local hardness benchmark with density, graph, MAAT, and shuffled-null baselines |
| `experiments/sat_cdcl_mode_a_t1/sat_cdcl_mode_a_t1.py` | 62 | Mode-A v1.6-T CDCL branching pilot with paired compute-to-solution regret against VSIDS |
| `experiments/sat_cdcl_structure_gated_mode_a_paper63/sat_cdcl_structure_gated_mode_a_paper63.py` | 63 | structure-gated Mode-A CDCL branching benchmark with instance-level gate diagnostics |
| `experiments/two_qubit_pointer_paper64/two_qubit_pointer_paper64.py` | 64 | two-qubit Lindblad pointer-selection benchmark for entanglement robustness under correlated generator stationarity |
| `experiments/sparc_ii_paper65/sparc_ii_paper65.py` | 65 | SPARC-II shuffled-radius nulls, morphology/proxy splits, and optional published halo-catalogue cross-check |
| `experiments/sparc_hgd_gsr_cross_framework_pilot/sparc_hgd_gsr_cross_framework_pilot.py` | pilot | collaboration-ready SPARC MAAT x HGD-GSR cross-framework test using external Q_bar input |

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
Paper 22 (validate): struct. beats energy rank     Δ_vac=2.583852, p_K=1.000, gain at 20% percentile
Paper 26 (constants): basin-level SM compatibility, MaxEnt weights R>V≈S>B>H, v11/v13 test predictivity
Paper 27 (boundary): R dominates closed λ fit       R share=0.3917, λ_R=8.078, N=3400 fused samples
Lambda closure:       λ from Cov[d] response         safe target: R share≈0.833, low-defect target: R share≈0.297
Paper 28 (cosmo CCI): structural-stress history     CCI_norm(z=10)≈803.8, chronometer RMS≈1.01
Paper 29 (cosmo v03): V/R measured + λ + z_c        z_c≈1.114, λ:S>H>R>V, CCI_v03(z=2)≈15.42
Paper 31 (dynamic):  diagnostics -> observables      track err=9.91e-7, residual=0.1648, max |DeltaH/H|=0.0399
Paper 32 (H(z)):     first chronometer chi2 test      LCDM χ²=14.8759, MAAT χ²=15.3661, χ²ν(best)=0.5299
Paper 33 (CCI proj): breadth-depth projection stress z_tr=0.8499, scan median=0.8889, N=6930
Paper 34 (proj/growth): distinct info layer           z_tr=1.0436, C_proj χ²ν=42.77 vs fσ8, scan N=625
Paper 35 (growth):  forward perturbation embedding    max |ΔD/D|=0.0591%, max |Δfσ8/fσ8|=0.4570%
Paper 36 (closure): R becomes derived robustness       R_resp=(HBV)^(1/3), R_rob=min(R_resp,(HBSV)^(1/4))
Paper 37 (v1.2.1): observable proxy + stability scan   max |ΔH/H|=1.70%, max |Δfσ8/fσ8|=1.83%, stable=1089/1089
Paper 38 (growth closure): v1.2.1 in linear growth      <R_rob>=0.9313, max |ΔD/D|=0.0591%, λ positivity passed
Paper 39 (growth signature): projection-modulated fσ8   ε_best=-0.0100, Δχ²=-0.0601, max |Δfσ8/fσ8|=0.9891%
Paper 40 (residual signature): CCI vs residual stress    ρS(CCI_diag, |rσ|)=0.5934, p=0.0338, null p≈0.036
Paper 41 (variable closure): definitions + measurement map H,B,S,V primary; R_rob emergent; no standalone code
Paper 42 (blind projection): response-derived projection  ρS(CCI_diag, |rσ|)=0.4286, p=0.1456, no shape tuning
Paper 48 (MaxEnt):   structural measure theorem           P[X]∝exp[-F_struct], support form = bounded viability
Paper 49 (benchmarks): universality audit matrix           complete/partial/missing domains=2/4/1, readiness=0.3816
Paper 50 (SAT II): graph-local SAT hardness benchmark      MAAT+density R²=0.2347 vs density 0.1886, graph 0.2469
Paper 62 (CDCL):   Mode-A active branching vs VSIDS         median regret=0.000, mean regret=+7.4619, wins/losses/ties=44/47/8
Paper 63 (gated):  structure-gated Mode-A vs global Mode-A  mean regret=-0.5703, median=-0.0080, wins/losses/ties=50/45/4
Paper 64 (2-qubit): correlated stationarity improves pointer-entanglement prediction  scalar R²: 0.7232→0.7831, field R²=0.9160
Paper 65 (SPARC II): shuffled-radius null supports NFW-like radial signal rho=0.3930 vs |rho|95=0.0355; RAR rho=0.0123 not significant
SO(10) extra: gauge one-loop + Yukawa bridge           M_GUT≈1.86e16 GeV, Δb≈0.0506, Yukawa R_rob≈0.999
SO(10) extended: response-closed multi-scale GUT benchmark ladder specified; falsification baselines declared
```

---

www.maat-research.com

## License
Code:
MIT License — free to use, modify, and share with attribution.

External data:
External scientific datasets used or redistributed in experiment folders retain
their original licences and attribution requirements. In particular, SPARC
rotation-curve data are attributed to Lelli, McGaugh, and Schombert and to the
Zenodo record `10.5281/zenodo.16284118`, listed under CC-BY-4.0. Repository
CSV/JSON/PNG files derived from external data are reproducibility artifacts and
do not imply endorsement by the original data providers.
