# MAAT--SO(10) Structural Selection Benchmark

This folder contains the reproducibility artifacts for the extra
phenomenological paper:

**Structural Selection in SO(10)-Motivated Unified Field Theories:**
*A Phenomenological MAAT Layer for Gauge and Yukawa Benchmarks*

## Scientific Status

This is a phenomenological bridge benchmark, not a complete SO(10) model.
It does not construct a full Higgs sector, prove precision unification,
derive the Standard Model spectrum from first principles, or solve vacuum
selection. It tests whether MAAT v1.2.1 structural selection can rank
SO(10)-motivated gauge and Yukawa parameter regions reproducibly.

## Scripts

Run from this folder:

```bash
python3 maat_so10_gauge_2d_fit.py
python3 maat_so10_yukawa_coupled_rg_v04.py
```

The scripts generate:

- `outputs/gauge_2d_summary.json`
- `outputs/yukawa_v04_summary.json`
- `outputs/maat_selection_landscape.png`

## Main Results

Gauge-sector one-loop diagnostic:

- `alpha_GUT = 0.022676468338`
- `M_GUT = 6.639198e15 GeV`
- `SM chi2 = 100.789861`
- `R_rob = 0.60301490`

Interpretation: the one-loop gauge benchmark locates a plausible GUT-scale
region but does not achieve precision unification. Threshold corrections,
intermediate scales, or model-specific spectra are required.

Yukawa-sector benchmark:

- `M_GUT = 1.858006e16 GeV`
- `y_t(MGUT) = 0.4136766812`
- `y_btau(MGUT) = 0.0077389774`
- `Delta_b = 0.0506451697`
- `chi2_yukawa = 0.00024341`
- `R_rob = 0.99916912`

Interpretation: an SO(10)-motivated high-scale `b-tau` boundary condition can
be made compatible with low-energy effective third-generation masses using a
modest fitted bottom-threshold correction.

## Data Attribution and Licence Note

The electroweak-scale inputs use standard reference values for `M_Z`, gauge
couplings, and third-generation fermion masses, cited in the paper to the
Particle Data Group review. These are external scientific reference data.
Repository JSON/PNG/PDF files are derived reproducibility artifacts only.
No endorsement by the Particle Data Group or any external collaboration is
implied.

