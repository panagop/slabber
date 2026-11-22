# Slabber

Slabber is a **research-only** Python project for reproducing Czerny and
Hahn's published rectangular-plate coefficient tables and comparing them
against numerical solutions.

**Not a structural design package.** The CSV numbers are manually transcribed
from historical scans and have **not** undergone an independent cell-by-cell
audit. PyNite and independent numerical comparisons cover all nine
four-edge support cases and all six three-edge cases at three aspect
ratios per case. The four-edge reference uses finite differences, and the
three-edge reference uses Rayleigh–Ritz.
Do not use results for reinforced-concrete design without independent
verification and appropriate code checks.

## Quick start — VS Code / Codex

Install [uv](https://docs.astral.sh/uv/getting-started/installation/),
[VS Code](https://code.visualstudio.com/), and Python/Codex extensions as
desired. Unzip this folder, open a terminal in its root and run:

```bash
uv python install 3.12
uv sync
uv run python -m unittest discover -s czerny_dataset -p 'test_*.py' -v
uv run python scripts/verify_data.py
code .
```

Slabber defaults to Python 3.12 via `.python-version`. `uv sync` creates
`.venv/` and installs the project dependencies, including `PyNiteFEA[pyvista]==3.2.1`.

### Launch the slab explorer

```bash
uv run --locked streamlit run app.py
```

The sidebar switches between four-edge and three-edge slabs. Select a
support case and aspect ratio to see the support plan, signed PyNite
centerline moment diagrams, and tables for physical values, normalized
coefficients and mesh refinement. The four-edge page offers all 21
printed ratios from 1.00 through 2.00 in steps of 0.05. It uses saved
results when available and calculates any other selection on demand.
The three-edge page uses the saved Hahn benchmark ratios. If its report
is absent, the page shows its regeneration command and draws the PyNite
diagrams.

Both pages also show full-slab maps of `Mx`, `My`, `Mxy` and vertical
deflection. The interactive maps show signed values at element centers,
with deflection in millimetres, positive downward. The **PyNite render** tab
uses PyVista/VTK to draw nodally averaged moment contours and the exaggerated
deformed shape. On Streamlit Community Cloud, the root `packages.txt` installs
the Debian libraries for off-screen rendering via EGL, with OSMesa as a
software fallback.

### First table calculation

```bash
uv run python czerny_dataset/czerny.py --case SSSS \
  --lx 4 --ly 6 --load 10 --E 33000000 --thickness 0.18
```

### Independently executed numerical benchmark (repeat it locally)

```bash
uv run python czerny_dataset/finite_difference_reference.py \
  --ratios 1 1.5 2 --refinements 8 16 32 64 \
  --output results/fd_recheck.csv
```

### Reproduce the PyNite `SSSS` pilot

```bash
uv run --locked python -m scripts.run_pynite_ssss --meshes 7 11 15 21 31 --output results/pynite_ssss.json
```

The pilot reports signed reaction balance, moments, deflections and mesh
convergence against Czerny and the independent finite-difference reference.
[Benchmark results](docs/PROJECT_STATUS.md#pynite-rect-ssss-pilot)

### Compare all nine four-edge cases

```bash
uv run --locked python -m scripts.run_four_edge_benchmark --output results/four_edge_benchmark.json
```

The run produces JSON for every mesh and a CSV for the finest-grid
comparisons. Pass `--ratios 1.05 1.10` to select printed ratios or
`--all-ratios` to run the full grid. [Four-edge results](docs/PROJECT_STATUS.md#nine-four-edge-cases)

### Compare all six Hahn three-edge cases

```bash
uv run --locked python -m scripts.run_three_edge_benchmark --output results/three_edge_benchmark.json
```

The benchmark compares PyNite with an independent Kirchhoff Rayleigh–Ritz
solution and the printed Hahn coefficients. [Three-edge results](docs/PROJECT_STATUS.md#six-hahn-three-edge-cases)

## Current contents

| Artifact | Scope | Status |
|---|---|---|
| `czerny_dataset/czerny_coefficients.csv` | Four-edge Czerny, nine cases, 21 ratios, 1,701 entries | Transcribed, **not cell-by-cell audited** |
| `czerny_dataset/czerny_three_edge_coefficients.csv` | Hahn, six three-edge cases, 456 entries | Transcribed secondary source; not full original Czerny |
| `czerny_dataset/czerny.py` | Lookup, interpolation, and dimensional conversion | 33 project tests pass, including 19 inherited checks |
| `czerny_dataset/finite_difference_reference.py` | Kirchhoff SSSS numerical solution | **Executed**; comparison CSV supplied |
| `scripts/run_pynite_ssss.py` | Reproducible PyNite `Rect` comparison | **Executed** at three aspect ratios and five mesh levels |
| `scripts/run_four_edge_benchmark.py` | PyNite, FD and Czerny comparisons | Nine four-edge cases at ratios 1.0, 1.5 and 2.0 |
| `scripts/run_three_edge_benchmark.py` | PyNite, Rayleigh–Ritz and Hahn comparisons | Six three-edge cases at ratios 0.50, 1.00 and 1.40 |
| `czerny_dataset/pynite_czerny_benchmark.py` | Inherited multi-case PyNite harness | Preserved for provenance |
| `czerny_dataset/_build*.py`, `_full_grid.py` | Reconstruction of source CSV transcriptions | Retained for provenance |
| `examples/*.json` | Example values from prior calculations | Retained reference outputs |

## Documentation index

- [Project status and evidence levels](docs/PROJECT_STATUS.md)
- [Sources and original references](docs/SOURCES.md)
- [Dataset schema, notation, physical units](docs/DATA_DICTIONARY.md)
- [Research chronology and decisions](docs/RESEARCH_HISTORY.md)
- [How to reproduce the tests and benchmarks](docs/REPRODUCIBILITY.md)
- [Detailed numerical methodology and next stages](docs/NUMERICAL_ROADMAP.md)
- [Technical limitations and validation checklist](docs/VALIDATION_AND_LIMITATIONS.md)
- [Codex handoff and next benchmark](docs/CODEX_HANDOFF.md)
- [Unmodified earlier notes](docs/original/README.md),
  [three-edge source notes](docs/original/THREE_EDGE_SOURCES.md), and
  [numerical validation report](docs/original/NUMERICAL_VALIDATION.md)

[`AGENTS.md`](AGENTS.md) contains project rules for Codex or other coding
assistants. The `.vscode/` folder provides launch configurations and tasks.

## Table conventions (do not mix them)

- Four-edge Czerny coefficients use `M = p*lx^2/k`, where `lx <= ly`.
  Shear uses `p*lx/k`; deflection uses a printed multiplier times
  `p*lx^4/(E*h^3)`.
- Three-edge **Hahn** coefficients generally use `K = p*lx*ly` and
  `M = K/k`; `y=0` is **free**. Ratios may be less than one.
- `F` in a boundary **case code** means **fixed**, whereas the free edge is
  explicitly named `FREE`. Boundary order: `x=0, y=0, x=lx, y=ly`.
- Printed rows have different physical locations. A spatial maximum is not
  necessarily at the center; e.g. `My_max` for long simply supported slabs.
- Missing tabulated quantities are unknown, **not zero**. Neither table
  family includes every possible combination of shears and reactions.

## External references

Original scans are **linked, not bundled**, to preserve provenance and
avoid redistributing potentially copyrighted works:

- [Czerny, Beton-Kalender 1996, four-edge table excerpt](https://www.klancek.si/sites/default/files/datoteke/files/tabele.pdf)
- [Hahn, Tafeln 17–22, three supported edges](https://fgg-web.fgg.uni-lj.si/KMLK/Sebastjan/masivni%20objekti/Hahnove%20tabele.pdf)
- [PyNite documentation](https://pynite.readthedocs.io/en/stable/plate.html)
- [uv documentation](https://docs.astral.sh/uv/guides/projects/)
