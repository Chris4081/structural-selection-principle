# HGD-GSR Q_bar Data Contract

Place a collaborator-provided CSV at `data/hgd_gsr_qbar.csv` with columns:

```text
galaxy,Q_bar
NGC2403,0.123
...
```

Accepted galaxy aliases: galaxy, Galaxy, name, Name, SPARC, galaxy_name.
Accepted Q_bar aliases: Q_bar, q_bar, Qbar, qbar, hgd_qbar, HGD_Q_bar.

The pipeline joins by exact galaxy name after stripping whitespace.
No HGD-GSR data are redistributed by this repository unless explicitly
provided under clear permission/licence terms by the data owner.
