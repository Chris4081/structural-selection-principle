# Phenomenological String-Landscape Selection Experiments

This directory contains the reproducibility bundle for the paper

**A Phenomenological Structural Selection Measure for String Backgrounds**.

The experiments are exploratory benchmark tests for extending the structural
selection framework to string-background ranking. They are not claimed to be a
complete string-theoretic derivation, a unique vacuum measure, a Standard Model
compactification, or a final Theory of Everything. Their purpose is narrower:
to test whether layered structural ranking carries information beyond
energy-only ordering in progressively more structured toy and bridge ensembles.

## Layout

```text
experiments/
├── phenomenological_structural_selection_string_measure.pdf
└── string_landscape_selection/
    ├── structural_selection_10d_tadpole_toy.py
    ├── structural_selection_iib_kklt_scan.py
    ├── structural_selection_iib_kklt_bridge.py
    ├── structural_selection_iib_period_kklt_bridge.py
    ├── structural_selection_iib_exact_period_kklt_bridge.py
    ├── structural_selection_iib_backreaction_sm_bridge.py
    ├── *_results.json
    └── *_plots/
```

## Requirements

The scripts require Python 3 and the following packages:

```bash
python3 -m pip install numpy scipy matplotlib mpmath
```

No network access is required once the Python packages are installed.

## Reproducing the Results

Run the scripts from the `string_landscape_selection/` directory:

```bash
cd experiments/string_landscape_selection

python3 structural_selection_10d_tadpole_toy.py
python3 structural_selection_iib_kklt_bridge.py
python3 structural_selection_iib_period_kklt_bridge.py
python3 structural_selection_iib_exact_period_kklt_bridge.py
python3 structural_selection_iib_backreaction_sm_bridge.py
```

Each script writes a JSON result file and a plot directory in the same folder.
The default seeds are fixed inside the scripts, so the benchmark outputs are
deterministic for the same Python/package versions.

## Script Overview

| Script | Purpose | Main outputs |
| --- | --- | --- |
| `structural_selection_10d_tadpole_toy.py` | Type-IIB-inspired flux toy ensemble with an explicit D3 tadpole proxy. | `structural_selection_10d_tadpole_toy_results.json`, `structural_selection_10d_tadpole_toy_plots/` |
| `structural_selection_iib_kklt_scan.py` | Shared reduced KKLT potential and structural diagnostic utilities. | Helper module used by later scripts. |
| `structural_selection_iib_kklt_bridge.py` | Flux-to-KKLT bridge ensemble with stationarity and stability tests. | `structural_selection_iib_kklt_bridge_results.json`, `structural_selection_iib_kklt_bridge_plots/` |
| `structural_selection_iib_period_kklt_bridge.py` | Large-complex-structure period proxy connected to KKLT ranking. | `structural_selection_iib_period_kklt_bridge_results.json`, `structural_selection_iib_period_kklt_bridge_plots/` |
| `structural_selection_iib_exact_period_kklt_bridge.py` | Mirror-quintic Picard-Fuchs period benchmark connected to KKLT ranking. | `structural_selection_iib_exact_period_kklt_bridge_results.json`, `structural_selection_iib_exact_period_kklt_bridge_plots/` |
| `structural_selection_iib_backreaction_sm_bridge.py` | Adds phenomenological backreaction and Standard-Model-sector proxy layers. | `structural_selection_iib_backreaction_sm_bridge_results.json`, `structural_selection_iib_backreaction_sm_bridge_plots/` |

## Scientific Status

These scripts implement toy and bridge models. In particular:

- The tadpole and backreaction terms are operational proxies, not full 10D
  supergravity solutions.
- The Standard-Model-sector layer is a phenomenological diagnostic proxy, not a
  constructed chiral compactification.
- The Picard-Fuchs benchmark uses a controlled mirror-quintic period model, but
  does not constitute a full global compactification scan.
- The structural measure is intended as a testable ranking architecture, not as
  a final microscopic string measure.

This status is intentional. The experiments provide a compact reproducibility
loop:

```text
model assumptions -> generated ensemble -> structural ranking -> JSON results -> plots
```

## Citation / Paper Reference

If using or discussing these experiments, cite the accompanying paper:

Christof Krieg,
*A Phenomenological Structural Selection Measure for String Backgrounds*,
2026.

