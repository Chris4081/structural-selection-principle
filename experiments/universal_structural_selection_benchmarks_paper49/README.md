# Universal Structural Selection Benchmarks — Paper 49

This folder contains the reproducibility bundle for Paper 49:

```text
Universal Structural Selection Benchmarks:
Testing MAAT Against Energy, Stability, Random, and Domain-Specific Baselines
```

## Status

This is a meta-benchmark registry. It does not introduce a new solver. Instead,
it aggregates results from existing experiments and audits whether the MAAT /
Structural Selection framework has been tested against simple alternatives
across independent domains.

The goal is deliberately conservative:

```text
Paper 48 derives the exponential structural-selection measure.
Paper 49 asks where that measure has actually beaten simpler competitors.
```

## Reproduce

Run from this folder:

```bash
python3 universal_structural_selection_benchmarks_paper49.py
```

The script reads existing outputs from the repository and writes:

| Output | Role |
| --- | --- |
| `outputs_paper49/paper49_domain_evidence_registry.csv` | Domain-by-domain evidence table |
| `outputs_paper49/paper49_benchmark_readiness_matrix.csv` | Competitor-test readiness matrix |
| `outputs_paper49/paper49_summary.json` | Machine-readable summary |
| `outputs_paper49/fig1_benchmark_readiness_matrix.png` | Heatmap of complete / partial / missing tests |
| `outputs_paper49/fig2_domain_evidence_status.png` | Domain-level evidence status chart |

## Competitor Classes

Paper 49 audits whether each domain has been compared against:

- energy-only ranking,
- entropy/activity-only ranking,
- stability-only ranking,
- random or shuffled null rankings,
- domain-specific simple baselines,
- cross-domain or holdout transfer tests.

## Initial Findings

The current registry is intentionally mixed:

- field and string benchmarks already contain explicit energy-only failures;
- cosmology has null/permutation protocols but remains diagnostic;
- boundary/safety calibration shows robustness dominance, but not yet external
  AI-system generalization;
- SAT validation is currently a stress-test/open problem rather than a success;
- quantum measurement / pointer-state selection is missing and should become a
  future universality benchmark.

## Data and License Note

This script aggregates derived outputs already present in the repository.
External scientific data used by the underlying experiments remain attributed
in their respective experiment folders and papers. No new external dataset is
redistributed here. The generated CSV/JSON/PNG files are derived
reproducibility artifacts only.
