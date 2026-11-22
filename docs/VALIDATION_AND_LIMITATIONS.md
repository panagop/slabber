# Validation checklist, open risks and engineering limits

## Historical source verification

- [ ] Independently re-read all 1,701 four-edge and 456 three-edge cells
  from licensed/available printed scans, including unusual outliers.
- [ ] Confirm all edge-condition sketches and x/y orientation from print.
- [ ] Confirm exact point/line of every maximum and support coefficient.
- [ ] Identify whether a printed denominator is based on `p*lx^2`,
  `p*lx`, total `K`, or some alternative quantity.
- [ ] Verify original three-edge Czerny chapter instead of calling Hahn
  an exact substitute; record every disagreement as a separate issue.
- [ ] Record two-person audit metadata if used for production engineering.

## Numerical cross-validation

- [x] Run `SSSS` independent finite-difference benchmark on 21 ratios.
- [x] Execute the `SSSS` PyNite benchmark on three ratios and five meshes
      using PyNiteFEA 3.2.1.
- [x] Confirm the `SSSS` local axes, applied load and deflection signs,
      force balance and physical units. Check peak moments and mesh changes
      against Czerny and the independent finite-difference reference.
- [x] Execute PyNite and an independent finite-difference solution for all
      nine four-edge cases at ratios 1.0, 1.5 and 2.0. Compare fixed-edge
      hogging at the printed midpoints or edge minima.
- [x] Execute PyNite and an independent Kirchhoff Rayleigh–Ritz solution
      for all six three-edge cases at ratios 0.50, 1.00 and 1.40. Check
      the H17 force partition and the released interior free-edge nodes.
- [ ] Check force equilibrium and spatial moment extrema for every case.
- [ ] Check convergence of **each** requested source coefficient, not
  merely global displacement; pay special attention to support edges.
- [ ] Separate transverse shear, effective shear, line reactions and
  corner concentrated forces.
- [ ] Check singular/discontinuous corner moments and reaction artifacts.
- [ ] Compare FEM with independent Kirchhoff/scikit-fem solution before
  asserting high accuracy for fixed/free edge cases.

The `SSSS` pilot includes bilateral vertical restraints at all four
corners. The negative nodal corner reactions are the forces needed to
hold those corners down in this idealized model. The JSON records them
separately from the open-edge support resultants and element shear. A
model with corner lift-off needs a different support condition. See the
[pilot results](PROJECT_STATUS.md#pynite-rect-ssss-pilot).

The expanded matrix uses the same corner restraint. It screens agreement
at three aspect ratios and reports outlying source and numerical
comparisons as discrepancies. The reported nodal force resultants and
corner forces are in kN. Source transverse shear and line reaction
density comparisons remain open.
[Matrix results](PROJECT_STATUS.md#nine-four-edge-cases)

The three-edge benchmark assumes `nu=0`. The compared Hahn pages lack
an explicit Poisson ratio. H17 line resultants and concentrated corner
forces retain distinct definitions. Fixed/free junction moments can
change substantially on the last mesh or polynomial refinement. Treat
those values as corner estimates with measured sensitivity.
[Three-edge results](PROJECT_STATUS.md#six-hahn-three-edge-cases)

## Structural engineering limits

Historical plate tables assume ideal elastic homogeneous isotropic plates,
idealized restraints and uniform full-area loading (in these extracts).
They do NOT automatically account for:

- Cracking, creep, shrinkage, concrete nonlinear stiffness, or creep
  deflections; steel reinforcement nonlinear contribution;
- Real joint/beam flexural stiffness, partial edge restraint, uplift,
  slab-column interactions, and construction-stage support movement;
- Reinforced-concrete design code ULS/SLS checks, detailing,
  reinforcement anchorage, durability, or seismic effects;
- Model error caused by thick plates, shear deformation or finite edge
  support width when thin Kirchhoff theory is assumed.

**Do not use this workspace directly as a structural design tool.**
It is an audit/research program for reproducing legacy coefficients and
exploring numerical plate models; verify against trusted references and
professional engineering practice before any real application.
