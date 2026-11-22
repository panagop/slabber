# Codex project handoff

## Context in one paragraph

We are digitizing **Czerny's four-edge rectangular slab coefficient tables**
and **Hahn's six three-edge tables** into long-form CSV, computing
dimensions/lookup using `czerny.py`, and independently reconstructing
equivalent plate responses. Four-edge original table grid: 9 configurations,
21 ratio nodes each, 1,701 values. Secondary Hahn set: 6 cases, 456
values with a `K=p*lx*ly` moment scale. Four-edge Czerny uses `p*lx^2`.
The SciPy FD and PyNite `Rect` comparison covers all nine four-edge cases
at three aspect ratios. The Rayleigh–Ritz and PyNite comparison covers
all six Hahn cases at three printed ratios. All 31 tests passed. The
independent audit of the full source tables remains open.
[Project status](PROJECT_STATUS.md)

## Next benchmark

Review the outlying [four-edge comparisons](PROJECT_STATUS.md#nine-four-edge-cases)
and [Hahn comparisons](PROJECT_STATUS.md#six-hahn-three-edge-cases).
The [runbook](REPRODUCIBILITY.md) records models, refinements and commands.

## Three-edge and solver work

1. Investigate the numerical behavior of the fixed/free junction moments
   and the largest discrepancies at printed Hahn stations.
2. Extend comparisons to additional printed knots as needed.
3. Compare another element or solver following the numerical roadmap.

## Do not do automatically

Avoid broad refactors, uncontrolled source CSV edits, "correcting" historical
values to match FEM, assuming pointwise shear equals support reactions,
or inferring engineering-grade validation from unit tests.
