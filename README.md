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
SO(10) extra: gauge one-loop + Yukawa bridge           M_GUT≈1.86e16 GeV, Δb≈0.0506, Yukawa R_rob≈0.999
```

---

www.maat-research.com

## License
Code:
MIT License — free to use, modify, and share with attribution.
