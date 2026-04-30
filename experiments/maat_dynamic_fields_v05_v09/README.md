# MAAT Dynamic Selection Fields v0.5-v0.9

This folder contains the reproducibility package for the MAAT dynamic
selection-field sequence developed in Papers 31-36.

The package is intentionally scoped as an effective-theory benchmark. It does
not claim a completed fundamental theory. Its purpose is to make the transition
from static structural selection weights to dynamic, local, and gravitationally
coupled selection fields reproducible.

## Structure

```text
maat_dynamic_fields_v05_v09/
├── v05_dynamic_lambda_flow/
├── v06_local_selection_fields/
└── v09_flrw_stability_scan/
```

Papers v0.7 and v0.8 are analytic bridge papers and therefore do not have a
separate numerical test script in this folder. Their PDF and LaTeX sources are
stored in `documentation/` as Papers 33 and 34.

## Reproducibility

Run the scripts from their respective subfolders:

```bash
cd experiments/maat_dynamic_fields_v05_v09/v05_dynamic_lambda_flow
python3 lambda_dynamic_flow_v05.py
```

```bash
cd experiments/maat_dynamic_fields_v05_v09/v06_local_selection_fields
python3 local_selection_phi4_v06.py
```

```bash
cd experiments/maat_dynamic_fields_v05_v09/v09_flrw_stability_scan
python3 maat_flrw_stability_scan_v09.py
```

Each script regenerates the corresponding CSV/JSON output files and figures in
its local `outputs/` directory.

## Paper Map

 Version | Role |
|-------|---------|
| v0.5 | Dynamic covariance-response flow for global structural weights |
| v0.6 | Local selection-pressure fields in a 1D nonlinear field benchmark |
| v0.7 | Effective gravitational coupling through `T_MAAT` |
| v0.8 | Worked scalar-field example and Lorentzian stability branch |
| v0.9 | FLRW stability scan for the kinetic selection branch |
| v0.5-v0.9 | Synthesis paper connecting the full chain |

## Key Results

- v0.5: selection weights are treated as dynamical response variables rather
  than fixed constants.
- v0.6: local MAAT fields reduce post-perturbation structural residuals in a
  perturbed 1D phi4 kink benchmark while leaving energy nearly unchanged.
- v0.7: local selection fields can be coupled to gravity through an effective
  stress-energy tensor derived by metric variation.
- v0.8: a scalar kinetic worked example shows that the Lorentzian sign branch
  must be fixed by no-ghost, sound-speed, and energy-positivity conditions.
- v0.9: the FLRW toy scan finds 666/784 dynamically stable trajectories and
  618/784 background-safe trajectories.

## Notes

The v0.9 scan is a stability and consistency benchmark. It is not a fit to
cosmological data and should not be interpreted as evidence for a completed
dark-energy model.
