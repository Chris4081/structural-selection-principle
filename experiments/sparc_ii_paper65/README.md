# Paper 65 -- SPARC II

This folder contains the Paper-65 SPARC-II robustness benchmark and the
exploratory external `Q_bar` comparison.

## Main run

```bash
python3 sparc_ii_paper65.py
```

Outputs are written to:

```text
outputs_paper65/
```

The main run performs the shuffled-radius null test, proxy morphology splits,
and optional halo-catalogue cross-check interface.

## External Q_bar comparison

```bash
python3 sparc_qbar_external_comparison.py
```

Outputs are written to:

```text
outputs_qbar_external/
```

The input table is:

```text
external_qbar_alhawarat.csv
```

Locked definition supplied for the pilot:

```text
Q_bar = median(V_bar, inner 20% radial points) /
        median(V_bar, outer 20% radial points)
```

`Q_bar` is treated as an external HGD-GSR descriptor, not as a MAAT-derived
quantity. The comparison is exploratory because it was added after the original
Paper-65 shuffled-radius null protocol.

The Q-bar comparison now includes the requested independence checks:

- direct correlation between `Q_bar` and `D_struct`;
- partial rank correlation controlling for an approximate baryonic mass scale
  and galaxy size;
- five-fold cross-validated NFW/RAR mismatch prediction using `D_struct` only,
  `Q_bar` only, and `D_struct + Q_bar`.

The baryonic mass control is the enclosed baryonic mass-scale proxy
`V_bar(R_max)^2 R_max / G`, derived from the available SPARC rotation rows. It
is not a full photometric stellar-mass estimate.

## Data attribution and license note

SPARC-derived inputs follow the Paper-61/Paper-65 attribution notes. SPARC
rotation-curve data are external scientific data attributed to Lelli, McGaugh,
and Schombert and are available through the Zenodo-hosted SPARC record under
CC-BY-4.0.

If a published halo-fit catalogue is supplied locally, cite the original
catalogue source and respect its distribution terms.

The `Q_bar` table is a collaborator-supplied derived HGD-GSR descriptor from
Ali Alhawarat. It is not an original SPARC measurement and not a MAAT-derived
quantity. Unless broader redistribution terms are separately supplied by the
HGD-GSR author, treat the table as a neutral comparison input for this pilot and
cite/acknowledge Ali Alhawarat and HGD-GSR when discussing or reusing the
comparison.

No endorsement by SPARC, Lelli et al., Li et al., Zenodo, VizieR, CDS, HGD-GSR,
Ali Alhawarat, or any original data provider is implied.
