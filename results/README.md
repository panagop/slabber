# Generated results

The PyNite `SSSS` pilot writes `pynite_ssss.json`. The execution log is
`pynite_ssss.log`. See the [runbook](../docs/REPRODUCIBILITY.md) for the
exact command and the [project status](../docs/PROJECT_STATUS.md) for results.

The nine-case comparison writes `four_edge_benchmark.json` with each
mesh and `four_edge_benchmark.csv` with finest-mesh comparisons. The
execution log is `four_edge_benchmark.log`.

The Hahn three-edge matrix writes `three_edge_benchmark.json` with
all refinements and `three_edge_benchmark.csv` with the finest comparisons.
The execution log is `three_edge_benchmark.log`. The focused runs at
ratio 0.50 use the `three_edge_0p50_refinement` filename stem.

The benchmark files named above are included so the app can load these
comparisons on a fresh checkout. Write new solver outputs here. Other
generated JSON, CSV and log files are ignored by Git until explicitly
curated. The inherited finite-difference baseline is in
`czerny_dataset/reference_ssss_fd*.csv`.
