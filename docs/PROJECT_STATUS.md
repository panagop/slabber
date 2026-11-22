# Slabber project status

Snapshot: 2026-11-22. PyNite benchmarks cover all nine four-edge and all
six Hahn three-edge cases at three aspect ratios per case. The independent
references are finite differences and Rayleigh–Ritz, respectively.

## Inventory

| Item | Count | What the count means |
|---|---:|---|
| Four-edge support configurations | 9 | Fixed and simply-supported combinations from Czerny 1996 |
| Four-edge aspect ratios | 21 per case | 1.00–2.00 inclusive, increments 0.05 |
| Four-edge numeric entries | 1,701 | Long-form CSV lines excluding header, not 1,701 individually checked cells |
| Three-edge support configurations | 6 | Hahn Tafeln 17–22, one free edge each |
| Three-edge numeric entries | 456 | All uniform-load rows transcribed, irregular source grid |
| Total entries | 2,157 | 1,701 + 456 |
| Inherited tests | 19 | Unit tests, not a full source audit |
| PyNite pilot tests | 4 | Signed force balance, normalization, displacement recovery and an actual solve |
| Four-edge extension tests | 4 | FD agreement with the inherited solution, rotation symmetry, fixed-edge stations and a square `FFFF` solve |
| Three-edge extension tests | 4 | Ritz cross-check, Hahn scaling, H17 force partition, free-edge restraints and station selection |
| Numerical SSSS comparison ratios | 21 | Finite-difference benchmark, at 40 divisions of short span |
| Expanded benchmark configurations | 27 | Nine four-edge cases at ratios 1.0, 1.5 and 2.0 |
| Hahn benchmark configurations | 18 | Six three-edge cases at ratios 0.50, 1.00 and 1.40 |

## What has actually been done

- **Python environment**: `uv sync` resolved Slabber's `uv.lock` and
  installed PyNiteFEA 3.2.1 as a main dependency under CPython 3.12.12.
  The PyNite import, all 19 inherited tests and `scripts/verify_data.py`
  passed in this environment.
- **Four-edge CSV**: source page references, formulas and table samples
  reconstructed from scanned *Beton-Kalender* 1996 pp. 281–289.
- **Three-edge CSV**: secondary-source compilation transcribed from Hahn
  1970, Tafeln 17–22, pp. 373–379. This is **not** Czerny's complete
  original three-edge table chapter, though Hahn's Tafel 22 cites Czerny.
- **Python lookup**: `lookup` interpolates physical reciprocal factors by
  default, `calculate` applies the original normalization and units.
- **Numerical SSSS pilot executed** with NumPy/SciPy finite differences.
  Its computed max absolute discrepancy across the historical 21-ratio grid
  (40 short-span divisions) was about 0.973% for `kx`, 0.515% for
  `ky_max`, and 0.654% for the deflection multiplier. These are
  discretized **FD** results, not FEM.
- **PyNite `Rect` pilot executed** at three aspect ratios and five mesh
  levels. See below for the quantitative results and [the runbook](REPRODUCIBILITY.md)
  for the model and exact command.
- xdfem2D and scikit-fem remain deferred.
- The Streamlit explorer displays the saved four-edge and three-edge
  comparisons and computes PyNite centerline moments for the selected case
  and aspect ratio. Launch it with `uv run --locked streamlit run app.py`.
- The four-edge selector covers all 21 printed ratios from 1.00 to 2.00.
  Selections without a saved result use on-demand PyNite and finite-difference
  runs on the two finest mesh levels in the benchmark. Coarse-mesh smoke
  checks passed for all 189 case-ratio combinations. The app's two finest
  mesh levels were also exercised for `FFFF` at 1.05 and `FFF_y` at 1.95.
  All 33 project tests passed.
- The full-slab Streamlit views show signed PyNite `Mx`, `My`, `Mxy` and
  vertical deflection at element centers. The optional `rendering` extra
  enables PyNite's nodally averaged moment contours and deformed-shape PNGs.
  Verification covered both pages with Streamlit 1.67.0 and native PNGs
  with PyVista 0.49.2 and VTK 9.7.2. All 31 numerical tests passed.

## PyNite Rect SSSS pilot

The run used PyNiteFEA 3.2.1, CPython 3.12.12, NumPy 2.5.5 and SciPy
1.18.1 on Windows. It used 7, 11, 15, 21 and 31 short-span divisions at
each ratio. All 15 meshes balanced their signed vertical loads to within
`9.71e-10%`. The short-span model parameters were `lx = 4 m`, uniform
downward `p = 10 kN/m²`, `E = 33e6 kN/m²`, `h = 0.18 m` and `nu = 0`.

| `ly/lx` | Finest `nx × ny` | Center `Mx` (kN·m/m) | Peak `My` (kN·m/m) | Center deflection down (mm) | Max error vs Czerny | Max error vs FD | Max change from `nx=21` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.00 | 31 × 31 | 5.89024 | 5.89026 | 0.64898 | 0.183% | 0.087% | 0.099% |
| 1.50 | 31 × 47 | 11.63784 | 4.61822 | 1.23353 | 0.352% | 0.047% | 0.081% |
| 2.00 | 31 × 63 | 15.42979 | 3.98652 | 1.61723 | 0.411% | 0.044% | 0.133% |

The maximum errors include both the physical response and its normalized
coefficient. The `My` peaks on the source mid-section and across the whole
PyNite panel agree to the precision shown. At ratio 2.00 the finest
mid-section peak is at `y ≈ 1.397 m` or its symmetric location
`y ≈ 6.603 m`, with `ly = 8 m`. Some unsmoothed element maxima are attained
along a segment, and the JSON records one maximizing location.

All 23 tests and the CSV integrity check passed. All 18 inherited file
hashes matched. See the [raw mesh results](../results/pynite_ssss.json)
and [solver log](../results/pynite_ssss.log).

## Nine four-edge cases

The expanded benchmark used the same loading and material parameters as
the pilot and the same five PyNite mesh levels, plus independent FD grids
with 32, 64 and 128 short-span divisions. It compares center and peak
moments, fixed-edge hogging, center deflection and eligible corner values.
All 135 PyNite solves balanced their signed vertical reactions within
`9.71e-10%` of the applied load. The 81 FD solves used the same geometry
and boundary codes. The finest-grid CSV has 135 source comparisons and
78 additional solver-to-solver comparisons.

| Case | Max PyNite vs Czerny | Max FD vs Czerny | Max PyNite vs FD at source stations |
|---|---:|---:|---:|
| `SSSS` | 0.922% | 0.392% | 0.932% |
| `F_y0` | 1.073% | 0.728% | 1.111% |
| `F_x0` | 1.230% | 0.580% | 1.111% |
| `FF_y` | 0.596% | 0.411% | 0.195% |
| `FF_x` | 1.015% | 1.108% | 0.178% |
| `FF_adj` | 2.226% | 2.259% | 0.164% |
| `FFF_x` | 0.963% | 0.809% | 0.352% |
| `FFF_y` | 25.734% | 25.574% | 0.397% |
| `FFFF` | 2.700% | 2.833% | 0.397% |

Each maximum takes the greater absolute difference in physical value or
normalized coefficient across all three ratios. Both solvers were within
the `2%` source screen on 132 of 135 comparisons. The three exceptions
are mid-section `My` peaks at ratio 2.0 in `FF_adj`, `FFF_y` and `FFFF`.
For `FFF_y`, printed page 288 gives denominator `97.0`, or `1.6495 kN·m/m`.
PyNite gives `2.0740` and FD gives `2.0713 kN·m/m` near `x=2.00 m`,
`y=7.10 m`. [Czerny table and section sketch, p. 288](https://www.klancek.si/sites/default/files/datoteke/files/tabele.pdf#page=12)

The largest PyNite-to-FD difference across all 213 comparisons was
`1.111%`. The last FD refinement changed each reported quantity by less
than `0.325%`. Four PyNite corner-force comparisons exceeded the `1%`
last-refinement screen, with a maximum of `1.186%`. These were `SSSS`,
`F_x0` and `F_y0` at ratio 1.0, and `F_x0` at ratio 2.0.

All 27 tests and the CSV integrity check passed. The 18 inherited files
matched their recorded hashes.

[Finest-grid comparison CSV](../results/four_edge_benchmark.csv) ·
[Raw JSON](../results/four_edge_benchmark.json) ·
[Execution log](../results/four_edge_benchmark.log) ·
[Runbook](REPRODUCIBILITY.md#nine-four-edge-cases)

## Six Hahn three-edge cases

The main run used 90 PyNiteFEA 3.2.1 `Rect` solves and 54 independent
Kirchhoff Rayleigh–Ritz solves. It covers printed ratios 0.50, 1.00 and
1.40 and uses `nu=0`. All 90 PyNite solves balanced the signed vertical
forces within `1.21e-9%`. The open free-edge nodes had zero vertical
reaction. The H17 Ritz forces recovered from effective shear and signed
corner jumps balanced within `0.0093%` on the three finest models.

| Case | Max PyNite vs Hahn | Max Ritz vs Hahn | Max PyNite vs Ritz at source stations |
|---|---:|---:|---:|
| `H17_SSS` | 33.05% | 33.01% | 3.10% |
| `H18_REAR_FIXED` | 10.55% | 10.58% | 1.70% |
| `H19_ONE_SIDE_FIXED` | 21.49% | 25.23% | 5.01% |
| `H20_BOTH_SIDES_FIXED` | 11.69% | 16.46% | 4.59% |
| `H21_ADJ_FIXED` | 74.58% | 74.41% | 6.47% |
| `H22_ALL_FIXED` | 52.76% | 52.69% | 6.18% |

These are absolute differences in physical quantities across the tested
ratios. The CSV also shows errors in normalized coefficients. Both solvers
were within `2%` of Hahn on 51 of 108 source rows when both measures were
screened. The two solvers agreed within `2%` on 101 of 108 source rows.

The largest source difference was center `My` for `H21_ADJ_FIXED` at
ratio 0.50. Printed p. 378 gives `K/80 = 1.000 kN·m/m`. PyNite gave
`0.2542` and Ritz gave `0.2559 kN·m/m`. A Ritz spot check at `nu=0.2`
gave `0.2455 kN·m/m`. The source value is retained in the CSV.
[Hahn Tafel 21, p. 378](https://fgg-web.fgg.uni-lj.si/KMLK/Sebastjan/masivni%20objekti/Hahnove%20tabele.pdf#page=18)

The main refinement screen flags the `Mx` values at fixed/free corners,
especially at ratio 0.50. PyNite also shows sensitivity in H17's open-edge
resultant and free-corner twisting estimate. Among the diagnostics without
a source value, `H18_REAR_FIXED` center `My` at ratio 0.50 differs by
`0.00291 kN·m/m` on the main finest grids, which is `32.44%` of the
near-zero Ritz value.

A focused run at ratio 0.50 refined H18–H22 from 31 to 45 PyNite `x`
divisions and from Ritz degree 22 to 32. The fixed/free corner moments
in H19–H22 still changed by `1.65–2.42%` in PyNite and `1.83–2.84%`
in Ritz. The two finer solutions differed by `4.76–6.85%` at those
corners. The H18 center-moment difference fell to `0.00181 kN·m/m`.
[Focused refinement CSV](../results/three_edge_0p50_refinement.csv) ·
[Focused raw JSON](../results/three_edge_0p50_refinement.json)

All 31 tests and the CSV integrity check passed. All 18 inherited file
hashes matched.
[Finest-grid comparison CSV](../results/three_edge_benchmark.csv) ·
[Raw JSON](../results/three_edge_benchmark.json) ·
[Execution log](../results/three_edge_benchmark.log) ·
[Runbook](REPRODUCIBILITY.md#six-hahn-three-edge-cases)

## Evidence levels / remaining uncertainty

1. **Bibliographic**: citations/URLs to primary and secondary scans.
2. **Transcribed**: historical CSV rows entered from sources; multiple
   checks but not a second person's independent source-by-source audit.
3. **Software tested**: 19 inherited and 12 added regression tests
   passed under the current lockfile.
4. **Numerically executed**: PyNite and FD for the nine four-edge cases,
   and PyNite and Ritz for the six Hahn cases at the selected ratios.
5. **Cross-checked**: see the matrix comparison for the remaining source
   discrepancies and fine-grid sensitivity.

No result qualifies as verified for reinforced-concrete structural design.

## Key folders

- `czerny_dataset/` — frozen-ish inherited CSVs, scripts and tests.
- `docs/original/` — untouched copies of prior research notes.
- `examples/` — saved earlier sample outputs.
- `results/` — results generated on the user's machine.
