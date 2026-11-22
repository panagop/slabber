# Numerical reconstruction pilot — Czerny SSSS (November 2026)

## What is actually executed?

- **Executed in this package environment:** `finite_difference_reference.py` (SciPy/Numpy finite-difference solver), for the four-edge simply-supported Kirchhoff–Love plate (`SSSS`) with Poisson ratio `nu=0`, 4 m short span, uniform `p=10 kN/m²`, `E=33e6 kN/m²`, `h=0.18 m`. Meshes of 8,16,32,64 elements along the short span at aspect ratios 1, 1.5, 2. The 21 historical aspect ratios 1.00,1.05,...,2.00 were also independently computed on meshes of 40 divisions along the short span. Raw results are in `reference_ssss_fd.csv` and `reference_ssss_fd_all_ratios.csv`.
- **Prepared and syntax-checked, but *not run*:** `pynite_czerny_benchmark.py`. This environment has no PyNiteFEA and outbound package download was unsuccessful. NO PyNite numerical results are claimed. Install PyNiteFEA in a connected local environment to run the file and compare its output with the executed FD results.
- Neither xdfem2D nor scikit-fem has been run. They are separate possible next-phase solvers, intentionally deferred.

## Why the `Rect` plate in PyNite?

PyNite's documentation says its `Rect` element uses a 12-term polynomial for out-of-plane bending and is suitable for thin-plate bending. Its `Quad` element uses a DKMQ formulation and includes shear deformation; for thin slabs like Czerny's it is useful to start with the `Rect` approximation. Rect elements allow direct queries of `moment(x,y)`, `shear(x,y)` and nodal reactions. See: https://pynite.readthedocs.io/en/stable/plate.html and https://pynite.readthedocs.io/en/stable/FEModel3D.html .

## Original coefficient normalization

For a **four-edge** Czerny slab with `a=lx`, uniform loading `p`:

```
Mx_mid = p*a²/kx                 => kx = p*a²/abs(Mx_mid)
My_max = p*a²/ky                => ky = p*a²/abs(My_max)
center_deflection = c*p*a⁴/(E*h³)  => c = |w_mid|*E*h³/(p*a⁴)
```

For **three-edge** Hahn tables with the convention `K=p*lx*ly`, bending denominators are `k=K/abs(M)`; they must not be mixed with the four-edge `p*lx²/k` convention.

The `My_max` value is the greatest positive bending moment anywhere on the panel, NOT necessarily the middle value. For ratio 2, the numerical field shows this maximum away from the center. This distinction is essential: use spatial samples and source sketch locations, not just a center evaluation.

## Reproducible independent finite-difference benchmark

The finite-difference solver discretizes `D*(∂xx+∂yy)^2 w = p` (linear elastic Kirchhoff plate) with the simply supported Navier conditions `w=0` and `∇²w=0` on the four edges, where `D=E*h³/[12(1-ν²)]`. It uses two sparse Poisson solves. Bending moments are recovered by second finite differences:

```
Mx=-D*(w_xx + nu*w_yy),  My=-D*(w_yy + nu*w_xx)
```

Current executable implementation checks only `nu=0`. It is *not* a finite-element method or PyNite and does not model any fixed/free edges. It is an independent benchmark with truncation errors.

With 64 short-span divisions:

| b/a | Czerny kx | FD kx | Czerny ky_max | FD ky_max | FD deflection multiplier c |
|--:|--:|--:|--:|--:|--:|
| 1.00 | 27.2 | 27.1528 | 27.2 | 27.1528 | 0.04875 |
| 1.50 | 13.7 | 13.7470 | 34.7 | 34.6345 | 0.09269 |
| 2.00 | 10.4 | 10.3682 | 40.3 | 40.1527 | 0.12155 |

The 21-point independent comparison at 40 short-span divisions had maximum *absolute* deviations (relative to printed Czerny coefficients) of **0.973% for kx**, **0.515% for ky_max**, and **0.654% for deflection multiplier c**. This validates the basic source interpretation and demonstrates coefficient reconstruction for SSSS; it does not establish numerical results for the other 14 configurations.

## Run locally

Prerequisites: Python >=3.11 and, for the independent benchmark, `numpy`, `scipy`; for the prepared PyNite script, `PyNiteFEA==3.2.1`:

```bash
python -m pip install numpy scipy "PyNiteFEA==3.2.1"
cd czerny_dataset
python finite_difference_reference.py
python run_fd_all_ssss.py
python pynite_czerny_benchmark.py --case SSSS --ratio 1.0 --meshes 7 11 15
python pynite_czerny_benchmark.py --case SSSS --ratio 1.5 --meshes 7 11 15
python pynite_czerny_benchmark.py --case SSSS --ratio 2.0 --meshes 7 11 15
```

The PyNite script additionally reads `coverage.csv` and can build the other 14 support configurations, including the six three-edge-supported/free-edge cases. **These other conditions and PyNite execution have not yet been independently validated**. If you run those cases, check each source edge sketch and convergence carefully.

## Mesh / boundary setup in the PyNite script

- Plate in global XY plane, rectangular elements wound bottom-left, bottom-right, top-right, top-left. Uniform downward pressure applied in local Z as `-p`.
- Simply supported edge: constrain transverse `DZ`, leave out-of-plane rotations free.
- Fixed/clamped edge: constrain `DZ`, `RX`, `RY` (this extra restriction does not conflict with a classical smooth fixed straight edge).
- Free edge: leave `DZ`, `RX`, `RY` free.
- Only for this *pure-bending, no membrane-load* pilot, constrain `DX`, `DY`, `RZ` at **every** node to remove otherwise unrestrained membrane/drilling modes. This artificial membrane condition must not be used as a model for true in-plane slab behavior.
- Element center moment samples are collected across the full plate to locate approximate extrema; center deflection is a four-node displacement average. Both results require mesh refinement; center and edge moments should not be assumed mesh-independent.
- Reaction balance `sum(RxnFZ)` against total applied load `p*a*b` provides an important physical sanity check, but not an independent guarantee of correct stiffness or edge moments.
- Printed moment denominator signs are not automatically comparable to the program's signed local components; code compares absolute values for the SSSS pilot.

## Open issues before automated generation of ALL Czerny values

1. For fixed-edge cases, identify exactly the location of printed negative moment maxima and sign convention; check if some values are averaged across support lines.
2. For free-edge cases, distinguish center / free-edge moments and torsion from global extrema. Edge/corner singularities and concentrated reactions may not converge as simple maxima.
3. Recover appropriate `Qx`,`Qy` at published stations separately from the **Kirchhoff effective boundary shear** (reaction includes twisting-moment derivatives). Nodal reactions summed along support lines do not directly equal pointwise shear.
4. Mesh refinement must be checked per output and support case. Compare `Rect` and `Quad` at matching meshes only after the first `Rect` run works.
5. Treat the historical CSVs as research transcriptions, not independently certified engineering values.

## Source links

- Czerny scanned four-edge tables (1996 excerpt): https://www.klancek.si/sites/default/files/datoteke/files/tabele.pdf
- Hahn three-edge tables: https://fgg-web.fgg.uni-lj.si/KMLK/Sebastjan/masivni%20objekti/Hahnove%20tabele.pdf
- PyNite documented plate elements/API: https://pynite.readthedocs.io/en/stable/plate.html
- PyNite FEModel3D API: https://pynite.readthedocs.io/en/stable/FEModel3D.html
- PyNite PyPI: https://pypi.org/project/PyNiteFEA/3.2.1/
