# MAAT Growth Perturbation Benchmark (Paper 35)

This folder contains the reproduction script and generated artifacts for:

**Paper 35 -- Linear Growth Embedding of MAAT Structural Selection**  
*Perturbation Stability, Positivity, and Sub-Percent Growth Signatures*

## Purpose

Paper 35 moves beyond the direct projection-template comparison of Paper 34.
Instead of comparing `C_proj` directly to `f sigma_8`, it embeds a small MAAT
source term into a linear growth equation and tests whether the resulting
growth-sector modifications remain stable and small.

## Main Command

Run from this folder:

```bash
python3 maat_paper35_reproduction.py
```

The script creates the `paper35_outputs/` folder and writes all CSV, JSON, and
PNG artifacts used in the paper.

## Key Results

| Quantity | Result |
|---|---:|
| `Omega_MAAT,0` | `0.00331` |
| `w_MAAT` | `-0.801` |
| Max `|Delta D/D|` for z < 2 | `0.0591 %` |
| Max `|Delta f sigma_8/f sigma_8|` for z < 2 | `0.4570 %` |
| NEC margin | `+0.199 rho` |
| Positivity stress tests | `3/3 positive` |
| Fixed-point mode stability | all tested `gamma_k > 0` |

## Outputs

| File | Role |
|---|---|
| `paper35_summary.json` | Machine-readable summary of all key numerical results. |
| `paper35_growth_curves.csv` | Growth factor and `f sigma_8` curves for LCDM and MAAT. |
| `paper35_dlambda_evolution.csv` | Fourier-space `delta lambda_a(k,t)` evolution. |
| `fig_paper35_growth.png` | Growth-factor and `f sigma_8` comparison. |
| `fig_paper35_perturbation.png` | Selection-pressure perturbation evolution and decay rates. |
| `fig_paper35_positivity.png` | Positivity and fixed-point stability plots. |

## Data and Attribution Note

The script includes a small representative `f sigma_8` table for visual
orientation only. The paper does not perform a likelihood fit. When discussing
observational growth data, cite the original survey publications and
compilations. The repository contains derived analysis artifacts only; no
endorsement by the original collaborations is implied.

## Scientific Status

This is a first perturbative toy embedding, not a full relativistic
perturbation solver and not a precision cosmological constraint. Its purpose is
to show that the representative MAAT branch can be inserted into a linear
growth equation without producing large or unstable effects.
