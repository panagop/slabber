# Proposed numerical coefficient-generation procedure

## Goal

Create independent *numerical estimates* of the specific historical
coefficients, not just a generic FEM output field. Preserve a clear
distinction between original printed values and numerical estimates.

## Stage 0 — baseline (done for SSSS only)

- Input: homogeneous isotropic Kirchhoff–Love slab, `nu=0`, uniformly
  distributed full-panel load `p`, constant `h`, relevant support codes.
- Solve simply supported (`SSSS`) with independent SciPy finite differences;
  store 3-ratio refinement study and 21-ratio comparison.
- The FD reference evaluates `My_max` over all interior grid nodes at
  ratios 1, 1.5 and 2.

## Stage 1 — PyNite Rect SSSS (complete)

The executed pilot used the pinned PyNiteFEA 3.2.1 `Rect` API and solved
all three ratios on short-span meshes 7, 11, 15, 21 and 31. The comparison
uses center `Mx`, center deflection and two positive `My` peaks. The
whole-panel peak matches the scope of the FD reference. The peak on
`x=lx/2` matches the station shown in the Czerny sketch on printed page
281. See the [reproducibility runbook](REPRODUCIBILITY.md) and
[pilot results](PROJECT_STATUS.md#pynite-rect-ssss-pilot).

## Stage 2 — nine four-edge support cases

The square `FFFF` check was extended to all nine four-edge cases at the
three pilot ratios. The FD reference adds fixed-edge ghost conditions.
Source comparisons use printed stations, including both fixed-edge minima
in `FF_adj`. Review the discrepancies listed in the
[matrix results](PROJECT_STATUS.md#nine-four-edge-cases) before extending
the comparison to the rest of the printed aspect-ratio grid.

## Stage 3 — three edges and one free edge

The PyNite and independent Kirchhoff Rayleigh–Ritz comparison covers
all six cases at three common printed knots. The free boundary is `y=0`.
The H17 benchmark distinguishes the integrated edge reactions from
the signed corner forces. Review the published `My` differences, the
mesh-sensitive fixed/free corner moments and the assumptions listed in
the [three-edge results](PROJECT_STATUS.md#six-hahn-three-edge-cases).

## Stage 4 — alternative solvers, one at a time

- Compare PyNite `Rect` against its `Quad`/DKMQ if appropriate.
- Then consider xdfem2D (plate FE) and scikit-fem (Kirchhoff plate
  example) as independent numerical reconstructions.
- Preserve each solver's own output and provenance; do not overwrite
  source-table CSVs with numerical estimates.

## Data/QA format for future numerical coefficients

Store new long-form numerical estimates separately, e.g.
`results/numerical_coefficients_pynite.csv` with columns for:

```text
solver,solver_version,case_id,ratio_ly_lx,quantity,
sample_location,mesh_nx,mesh_ny,nu,load,E,thickness,
moment_or_deflection,normalized_coefficient,source_coefficient,
error_percent,converged,run_id,notes
```

Include raw JSON output, machine OS/Python version and command-line
arguments. Do not mix moment denominators with physical multipliers.
Explicitly label coefficients that lack a match in the published table.
