# Runbook — uv, Python, VS Code

## Prerequisites and lockfile

- Slabber defaults to **Python 3.12**, pinned in `.python-version`.
- Install `uv`: https://docs.astral.sh/uv/getting-started/installation/.
- Recommended IDE: VS Code with Python/Pylance and Codex extension.
- Run from the **repository root**, never from an unrelated working dir.

```bash
uv python install 3.12
uv sync
uv run python --version
uv run python -m unittest discover -s czerny_dataset -p 'test_*.py' -v
uv run python scripts/verify_data.py
```

`uv sync` uses `uv.lock` to install NumPy, SciPy and PyNiteFEA 3.2.1
into `.venv/`. If installation fails, check the Python architecture,
SSL configuration and network connection.

## Historical table examples

```bash
uv run python czerny_dataset/czerny.py --case SSSS --lx 4 --ly 6 --load 10 --E 33000000 --thickness .18
uv run python czerny_dataset/czerny.py --case H17_SSS --lx 2.5 --ly 1.5 --load 9.5 --E 33000000 --thickness .18
uv run python czerny_dataset/czerny.py --case H22_ALL_FIXED --lx 2.5 --ly 1.5 --load 7.5 --E 33000000 --thickness .18
```

Prior example JSONs are in `examples/`. `lookup()` returns physical
factor + printed provenance; `calculate()` returns dimensioned results.
Note that H17/H22 use `K=p*lx*ly`, not `p*lx**2`.

## FD numerical validation — repeat a quick comparison

```bash
uv run python czerny_dataset/finite_difference_reference.py \
  --ratios 1 1.5 2 --refinements 8 16 32 64 \
  --output results/fd_recheck.csv
```

`czerny_dataset/reference_ssss_fd.csv` and
`reference_ssss_fd_all_ratios.csv` preserve the executed historical pilot.
To regenerate *all* original ratios (writes its original result path):

```bash
uv run python czerny_dataset/run_fd_all_ssss.py
```

Before running a data-rebuilding script (`_build_dataset.py` or
`_build_three_edge.py`), make a copy or use Git; these overwrite CSVs.

## PyNite Rect `SSSS` pilot

```bash
uv run --locked python -m scripts.run_pynite_ssss --meshes 7 11 15 21 31 --output results/pynite_ssss.json
```

This command runs ratios `ly/lx = 1, 1.5, 2` using `lx = 4 m`,
uniform downward pressure `p = 10 kN/m²`, `E = 33e6 kN/m²`,
`h = 0.18 m` and `nu = 0`. The script chooses an odd long-span mesh near
`nx*ly/lx`, so the geometric center is inside one element. The JSON stores
runtime versions, per-mesh dimensions, responses, comparisons and criteria.

PyNite's `Rect` elements are wound from the lower-left corner toward `+X`
and then `+Y`. The script records the local axes expressed in global
coordinates and applies pressure along negative Z.
All four edges have `DZ=0` with `RX` and `RY` free. The `DX`, `DY` and `RZ`
constraints stabilize the unloaded in-plane DOFs of this pure-bending pilot.
The model uses `add_plate`, `add_plate_surface_pressure`, `analyze_linear`,
element `moment`, `shear` and `d`, and nodal `RxnFZ` with PyNiteFEA 3.2.1.
[PyNite plate API](https://pynite.readthedocs.io/en/stable/plate.html)

The center deflection is evaluated from the Rect displacement polynomial
using the four nodal displacements and bending rotations. The unsmoothed
`Mx` and `My` fields are bilinear within each element. Checking all four
corners of every element gives their panel-wide maxima, including both
element limits at shared boundaries. The center-strip endpoints give the
positive `My` peak along `x=lx/2`, the mid-section in Czerny's sketch.
The inherited finite-difference CSV provides the 64 short-span division
reference and its maximum over interior grid nodes.
[Czerny notation, page 277](https://www.klancek.si/sites/default/files/datoteke/files/tabele.pdf#page=1)
· [Czerny table and sketch, page 281](https://www.klancek.si/sites/default/files/datoteke/files/tabele.pdf#page=5)

For this sign convention sagging moments are positive and `w` is negative
under the downward load. The four-edge coefficient denominators are
`p*lx²/M` and the deflection multiplier is `-w*E*h³/(p*lx⁴)`. Comparisons
include signed percent differences for both the physical response and its
coefficient. The signed vertical residual is `sum(RxnFZ) - p*lx*ly`,
divided by `p*lx*ly` for the percentage. Each corner appears once in
that sum. Open-edge resultants are signed sums of their non-corner nodal
forces. Element shear is reported in kN/m and reactions in kN.

The pilot requires an absolute relative reaction residual of `1e-6%` or
less on every mesh, absolute finest-mesh errors of `2%` or less against
both references, and an absolute change of `1%` or less on the last
refinement for the compared responses and their coefficients.
[Observed results](PROJECT_STATUS.md#pynite-rect-ssss-pilot)

## Nine four-edge cases

Run from the repository root.

```bash
uv run --locked python -m scripts.run_four_edge_benchmark --output results/four_edge_benchmark.json
```

The extension uses the pilot's material, load and short span. It runs
the nine cases in `coverage.csv` at exact printed ratios `1.0`, `1.5` and
`2.0`. The PyNite short-span counts are `7, 11, 15, 21, 31` and the FD
counts are `32, 64, 128`. All four borders restrain vertical displacement.
A fixed border also restrains both nodal bending rotations. The corners
remain vertically restrained. The JSON records edge ordering, the two
grids at each level, source provenance, dimensional values, sample
locations and nodal reactions. The companion CSV gives the finest-grid
comparisons and the changes on the last refinement in each solver.

The Streamlit four-edge page offers every printed ratio from 1.00 to
2.00 in steps of 0.05. For a selection without saved results, it runs
PyNite at `nx=21,31` and FD at `nx=64,128`, then caches the comparison in
the app process. The FD long-span count is the nearest even
number to `nx*ly/lx`, so the center remains a grid node. The existing
default ratios keep their original grids.

To reproduce one on-demand selection and write a separate report:

```bash
uv run --locked python -m scripts.run_four_edge_benchmark --cases FFFF --ratios 1.05 --fem-meshes 21 31 --fd-meshes 64 128 --output results/four_edge_ffff_1p05.json
```

Use `--all-ratios` in place of `--ratios` for the full 21-ratio grid.

The FD reference solves `D Δ²w = p` with `nu = 0`. The simply supported
edges use odd ghost reflection for zero normal bending moment. The fixed
edges use even ghost reflection for zero normal slope. Fixed-edge moments
use the resulting boundary normal curvature. The `SSSS` regression
compares this implementation with the inherited two-Poisson solution.

Each solver reports the center moments and deflection and the positive
panel-wide `Mx` and `My` maxima. Source `my_max` is sampled along `x=lx/2`
and source `my_mid` at the center. The fixed-edge normal moment is sampled
at the midpoint or searched for the most negative value along the edge,
according to the printed row. For `FF_adj` the printed p. 286 rows
`m_xermin` and `m_yermin` map to edge minima. The CSV keys `mx_fixed_mid`
and `my_fixed_mid` are kept in the result metadata for traceability.
[Czerny tables and sketches, pp. 281–289](https://www.klancek.si/sites/default/files/datoteke/files/tabele.pdf#page=5)

At corners shared by two simply supported edges, the table comparison
uses the magnitudes of the twisting moment and corner force. PyNite
supplies signed nodal vertical corner reactions. The FD corner force is
derived from the corner twisting moment using the Kirchhoff relation
`|R| = 2|Mxy|`. Line support resultants in the JSON are signed sums over
the non-corner support nodes. Pointwise transverse shear and line reaction
densities have separate source rows and are outside this comparison.

Results can be screened using the pilot's `2%` reference difference and
`1%` last-refinement thresholds. The maximum error considers the physical
quantity and its normalized coefficient. See the
[matrix results](PROJECT_STATUS.md#nine-four-edge-cases) for exceptions.

## Six Hahn three-edge cases

```bash
uv run --locked python -m scripts.run_three_edge_benchmark --output results/three_edge_benchmark.json
```

The default ratios `0.50`, `1.00` and `1.40` are printed in all six
tables. The runner accepts other exact, in-range printed knots for the
selected cases. The benchmark uses `lx=4 m`, `p=10 kN/m²`,
`E=33e6 kN/m²`, `h=0.18 m` and `nu=0`. The benchmark adopts `nu=0`
explicitly. The printed pages used for comparison lack a Poisson-ratio
value. It uses PyNite `x`-axis divisions `7, 11, 15, 21, 31` and
Rayleigh–Ritz polynomial degrees `8, 14, 22` per coordinate.

PyNite places no bending or transverse displacement restraint on the
interior nodes of `y=0`. Each corner still follows the adjoining side
support. The independent Ritz solution minimizes Kirchhoff–Love bending
energy using tensor-products of Legendre polynomials. It imposes zero
displacement at supported edges and zero normal slope at fixed edges.
Zero normal moment and zero effective shear on the open free edge arise
from the variational problem. Ritz extrema use a 241 × 241 sampling grid.
PyNite scans the bilinear moment field at element corners, including
each element's limit on shared boundaries.
[Kirchhoff–Love variational boundary conditions](https://docu.ngsolve.org/ngs24/SaS/plates_derivation.html)

The source comparison uses `K=p*lx*ly` for moment denominators and
`K*lx²/(omega*E*h³)` for H17 free-edge midpoint deflection. For H18 the
profile on `x=lx/2` provides the positive `My` peak and the first zero
crossing measured from the rear, which defines `beta`. For H19 and H20,
the table notes switch `My` from the center to a peak near the rear.
The benchmark samples a peak over the rear half for ratios above 0.8
and 0.6 respectively. The CSV gives both solver station coordinates.
[Hahn tables and station sketches, pp. 373–379](https://fgg-web.fgg.uni-lj.si/KMLK/Sebastjan/masivni%20objekti/Hahnove%20tabele.pdf#page=13)

H17's `Kx` is the force along one side and `Ky` is the rear force. The
PyNite line resultants sum only non-corner support-node reactions. Its
corner reactions are retained with their signs and can include part of
the distributed edge force. The Ritz H17 line resultants integrate the
effective shear along each supported edge. The Ritz corner forces come
from the signed `2*Mxy` twisting jumps. The JSON reports the Ritz force
residual from this recovery separately from the PyNite nodal balance.
The H17 worked example on p. 374 gives the force partition used here.
[Hahn worked example, p. 374](https://fgg-web.fgg.uni-lj.si/KMLK/Sebastjan/masivni%20objekti/Hahnove%20tabele.pdf#page=14)

The focused refinement of the short-panel cases is reproducible with
this command.

```bash
uv run --locked python -m scripts.run_three_edge_benchmark --cases H18_REAR_FIXED H19_ONE_SIDE_FIXED H20_BOTH_SIDES_FIXED H21_ADJ_FIXED H22_ALL_FIXED --ratios 0.5 --meshes 31 45 --degrees 22 32 --output results/three_edge_0p50_refinement.json
```

[Three-edge results](PROJECT_STATUS.md#six-hahn-three-edge-cases)

## VS Code

```bash
code .
```

Select the interpreter in `.venv`. The `.vscode/tasks.json` includes
data integrity checks, unit tests and an FD run. The launch configs
cover table lookup and the PyNite pilot. The Codex extension may be installed from
the VS Code extension marketplace; tell Codex to read `AGENTS.md` and
`docs/CODEX_HANDOFF.md`.

## Troubleshooting

- If the `Pynite` import fails, run `uv sync` and check that your IDE uses
  the interpreter in `.venv/`.
- A missing NumPy/SciPy means `uv sync` was not completed.
- Error about absent CSV: keep scripts under `czerny_dataset/`, because
  original code resolves data paths relative to its own file.
- If a PyNite run fails, do **not** claim earlier FD coefficients came
  from FEM. First verify library signatures and coordinate conventions.
