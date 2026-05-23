# Paper 60 - Fluid Coherence and Blow-up Diagnostics

This experiment tests whether MAAT-style structural supports can act as
early-warning diagnostics for fluid steepening and turbulence-like activity.

It is deliberately conservative:

- not a proof of Navier-Stokes regularity
- not a 3D turbulence simulation
- not a replacement for PDE analysis
- not a claim that MAAT predicts singularities in physical fluids

## Test Systems

The script runs two controlled toy systems:

- `viscous Burgers`: a 1D shock-warning proxy with known inviscid shock time
- `2D Navier-Stokes`: an unforced incompressible vorticity simulation used as a
  regularity/control case

## Operational Supports

The diagnostic uses:

- `H`: equation-residual consistency
- `B`: energy/enstrophy balance consistency
- `S_eff`: controlled activity pressure
- `V`: low-mode spectral coherence
- `R_rob`: v1.2.1 robustness closure
- `MAAT_warning`: robustness loss weighted by activity and spectral tail pressure

## Run

```bash
python3 fluid_blowup_diagnostics_paper60.py
```

## Core Results

| Quantity | Result |
|---|---:|
| Burgers inviscid shock time | `1.3334` |
| Burgers first 70% warning time | `1.1600` |
| Burgers warning lead time | `0.1734` |
| Spearman warning vs max gradient | `0.8413` |
| Pre-shock warning vs inverse time-to-shock | `0.9952` |
| 2D NS false high-warning fraction (`warning > 1`) | `0.0000` |
| 2D NS minimum robustness | `0.8956` |

## Main Finding

The MAAT warning functional tracks Burgers gradient steepening strongly and
crosses a high-warning threshold before the inviscid shock time. This is partly
expected because the diagnostic explicitly includes activity pressure and
spectral-tail growth. The stronger control result is that, in the regular 2D
Navier-Stokes case, the warning remains bounded despite increasing
palinstrophy.

This supports a narrow interpretation:

> MAAT fluid supports may be useful as structural early-warning diagnostics,
> not as a proof mechanism for Navier-Stokes regularity.

## Outputs

The script writes:

- `outputs_paper60/paper60_burgers_diagnostics.csv`
- `outputs_paper60/paper60_ns2d_diagnostics.csv`
- `outputs_paper60/paper60_summary.json`
- `outputs_paper60/fig1_burgers_warning_timeseries.png`
- `outputs_paper60/fig2_ns2d_regular_control.png`
- `outputs_paper60/fig3_warning_scatter.png`

## Data and License Note

No external fluid dataset is redistributed. All fields are generated locally by
the script from synthetic initial conditions. CSV/JSON/PNG outputs are derived
reproducibility artifacts.
