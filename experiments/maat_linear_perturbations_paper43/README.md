# Paper 43 — Linear Perturbations and Structure Growth

This folder contains the first perturbation-level benchmark for MAAT
structural-selection cosmology.

Paper 43 moves beyond growth-index proxies of the form
`f(z) ~= Omega_m(z)^gamma` by solving an explicit linear matter-growth
equation:

```text
D_NN + [2 + d ln H / dN] D_N - 3/2 Omega_m(a) mu(a) D = 0
```

with `N = ln(a)` and

```text
mu(z) = G_eff/G = 1 + eta * C_hat_proj(z).
```

The raw Paper-42 projection observable grows strongly with redshift, so the
perturbative coupling uses a bounded projection kernel:

```text
C_hat_proj(z) = [C_proj(z) - C_proj(0)] / [C_proj(zmax) - C_proj(0)]
```

over the tested interval `0 <= z <= 3`. This makes `eta` the maximum fractional
change in the effective growth coupling over the benchmark interval.

## Run

```bash
cd experiments/maat_linear_perturbations_paper43
python3 maat_linear_growth_solver_v01.py
```

## Outputs

The script writes into `outputs_paper43/`:

- `paper43_summary.json`
- `paper43_growth_curves.csv`
- `paper43_eta_scan.csv`
- `paper43_fsigma8_comparison.csv`
- `fig1_growth_perturbation_summary.png`
- `fig2_eta_scan.png`

## Main Results

| Quantity | Result |
| --- | --- |
| Eta scan | `eta in [0, 0.08]`, `41` points |
| Stable branches | `41/41` |
| Best eta against compact `f sigma_8` table | `0.0000` |
| LCDM chi-square | `12.2021` |
| Representative eta | `0.0200` |
| Max `abs(Delta D/D)` for `z <= 3` at eta `0.02` | `0.7656 %` |
| Max `abs(Delta f sigma_8/f sigma_8)` for `z <= 3` at eta `0.02` | `0.5441 %` |
| Max `abs(mu-1)` for `z <= 3` at eta `0.02` | `1.9966 %` |
| Representative eta chi-square | `12.5737` |

The benchmark is therefore stable and perturbatively controlled. It is not
favoured over LCDM by this compact comparison set; its value is the explicit
growth-equation bridge from structural projection to perturbation dynamics.

## Scientific Status

This is a controlled effective linear-growth benchmark. It is not a
Boltzmann-code perturbation solver, not a CMB or weak-lensing likelihood, and
not a precision cosmological fit. Its purpose is to test whether a
projection-derived structural coupling can be inserted into a standard linear
growth equation without producing large or unstable effects.

## Data Attribution and License Notes

The Planck-normalised reference parameters and compact `f sigma_8` comparison
points are external scientific data and should be cited to the original
publications/collaborations when reused or discussed. The CSV tables and
figures generated here are derived reproducibility artifacts only.

No endorsement by the Planck Collaboration, survey collaborations, or original
data authors is implied. Repository code is released under the repository
license.
