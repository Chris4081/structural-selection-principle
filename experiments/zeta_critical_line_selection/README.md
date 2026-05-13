# Extra Paper — Riemann Zeta Critical-Line Diagnostic

**A Structural Selection Diagnostic on the Riemann Zeta Critical Line**  
*A purely diagnostic, non-proof toy experiment with controls, ablations, and pair-correlation tests*

This experiment is a mathematical diagnostic toy benchmark.
It is not a proof of, or evidence for, the Riemann Hypothesis.

The test asks whether known nontrivial Riemann zeta zeros exhibit lower
diagnostic cost on the critical line than under small off-line shifts, and how
much of that behaviour is built into the scoring rule itself.

The pair-correlation component adds a less pointwise test: the normalised
nearest-neighbour spacing distribution of the first 200 zeta zeros is compared
with a GUE spacing proxy and Poisson/uniform controls.

## Diagnostic Terms

For `s = sigma + i t`, the default diagnostic uses:

```text
H(sigma,t) = 1 / (1 + |zeta(sigma + i t)|)
B(sigma)   = 1 / (1 + alpha |sigma - 1/2|)
R(sigma,t) = 1 / (1 + |zeta(sigma + dsigma + i t)
                    - zeta(sigma - dsigma + i t)|)
```

The diagnostic cost is:

```text
F_diag = -log(eps + H) - log(eps + B) - log(eps + R)
```

Lower values mean lower diagnostic stress under the proposed score.

## Controls

The script includes:

- known zeta zeros,
- jittered-zero controls,
- random critical-line points,
- no-B ablation (`HR`),
- no-H ablation (`BR`),
- parameter sweeps over `alpha` and `dsigma`.

## Headline Results

Default run:

```text
n_zeros  = 40
n_random = 40
alpha    = 10.0
dsigma   = 0.002
```

| Set | Config | Best mean delta | Mean F at best |
|---|---|---:|---:|
| known_zero | HBR | 0.000 | 0.00818 |
| known_zero | HR | 0.000 | 0.00818 |
| known_zero | BR | 0.000 | 0.00818 |
| jittered_zero | HBR | 0.000 | 0.16544 |
| random_critical_line | HBR | 0.000 | 0.86054 |
| random_critical_line | HR | 0.050 | 0.83426 |

Parameter-sweep result for known zeros:

| Config | Fraction with delta=0 best | Mean abs(best delta) |
|---|---:|---:|
| HBR | 1.0 | 0.000 |
| HR | 1.0 | 0.000 |
| BR | 0.8 | 0.010 |

Pair-correlation diagnostic result:

| Series | Spacings | L1 to GUE | V_pair | F_pair |
|---|---:|---:|---:|---:|
| jittered zeta spacings | 199 | 0.2059 | 0.9844 | 0.0157 |
| zeta normalised spacings | 199 | 0.2938 | 0.9797 | 0.0205 |
| uniform-points control | 199 | 0.8103 | 0.9573 | 0.0437 |
| poisson control | 199 | 0.8225 | 0.9559 | 0.0451 |

## Interpretation

Known zeta zeros show a strong low-cost minimum on the critical line under the
full HBR diagnostic.
However, this is not RH evidence because:

- `H` directly rewards small `|zeta|`,
- `B` explicitly rewards `sigma = 1/2`,
- random critical-line points also inherit the `B/R` critical-line minimum.

The legitimate conclusion is therefore:

```text
Known zeta-zero locations exhibit a low-cost minimum on the critical line, but
the present diagnostic partly builds in critical-line preference.
```

The pair-correlation diagnostic is less tautological and conceptually closer to
spectral-universality tests: it shows that the zeta spacing geometry is closer
to the GUE proxy than Poisson/uniform controls. However, jittered zeta spacings
remain similarly strong, so the diagnostic should be interpreted as a
spacing-geometry signal, not exact zero-location evidence.

## Reproduction

Run from the repository root:

```bash
python3 experiments/zeta_critical_line_selection/zeta_critical_line_selection_v02.py
```

Outputs are written to:

```text
experiments/zeta_critical_line_selection/outputs/
```

## Data and License Note

No external dataset file is redistributed.  The script computes the first
Riemann zeta zeros using `mpmath.zetazero` and evaluates the zeta function
directly. Generated CSV/JSON/PNG files are derived reproducibility artifacts.
