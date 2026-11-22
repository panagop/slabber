# Czerny / Hahn rectangular plate coefficient library

**Research-only release, November 2026.** CSV transcriptions plus a dependency-free
Python 3.9+ lookup and unit-scaled evaluator for two collections:

| Source | CSV | Cases | Original numeric entries | Aspect ratio |
|---|---|---:|---:|---|
| Czerny, *Beton-Kalender* 1996, four-edge pp. 281–289 | `czerny_coefficients.csv` | 9 | 1,701 | `ly/lx = 1.00,1.05,...,2.00` |
| Hahn, 1970, Tafeln 17–22, three-edge pp. 373–379 | `czerny_three_edge_coefficients.csv` | 6 | 456 | `ly/lx = 0.25,...,1.50` or `2.00`, irregular spacing |
| **Total** | | **15** | **2,157** | |

The three-edge numerical values come from Hahn's **secondary tables**, not from a
verified complete transcription of Czerny's own 1996 Chapter 4. Tafel 22
explicitly attributes its moment numbers to Czerny. All three-edge cases
have one free edge and are normalized differently from the four-edge cases.

## Sources

- Czerny, F., *Tafeln für Rechteckplatten*, *Beton-Kalender* (1996), pp. 277–289
  freely accessible excerpt: https://www.klancek.si/sites/default/files/datoteke/files/tabele.pdf
- Hahn, J., *Durchlaufträger, Rahmen, Platten und Balken auf elastischer
  Bettung* (1970), Tafeln 17–22: https://fgg-web.fgg.uni-lj.si/KMLK/Sebastjan/masivni%20objekti/Hahnove%20tabele.pdf
- Details and limitations of Hahn transcription: [`THREE_EDGE_SOURCES.md`](THREE_EDGE_SOURCES.md)
- Coverage and source page per case: `coverage.csv`

## How the numbers work

For all **four-edge** cases, the x span `lx` is the shorter plate side,
`ly >= lx`; `nu=0`, uniform full load `p`:

- Bending/twisting moments and corner point forces: `p*lx²/k`, signed as tabulated.
- Transverse shear and distributed support reactions: `p*lx/k`.
- Elastic middle deflection: `c*p*lx⁴/(E*h³)` (the printed `c` is a multiplier).

For all **three-edge** cases, `lx` is horizontal, `ly` is vertical,
`y=0` is the free edge; `lx` need not be shorter:

- Set `K = p*lx*ly`, the **total panel load**.
- All tabulated bending/twisting moments: `M = ±K/k` (moment per unit width).
- Tafel 17 support reaction resultants: `Kx = nu_x*K`, `Ky = nu_y*K`.
- Tafel 17 free-edge deflection: `w_r = K*lx²/(omega_r*E*h³)`.
- Tafel 18 negative-bending region reach: `beta*ly`.

The tables include different **quantities**, not a full output suite per
case. In particular, no transverse shear or distributed reaction data was
transcribed for Hahn Tafeln 18–22 because those rows are not provided in the
published `Lastfall 1` tables. Do not interpret missing quantities as zero.
The `H19`/`H20` row named `my_center_or_max` switches interpretation from
center moment to maximum near the rear edge at larger aspect ratios; see
source notes.

`F` in case IDs and edge codes means fixed/clamped; `FREE` unrestrained; `S`
simply supported. Edges are ordered `(x=0, y=0, x=lx, y=ly)` in `coverage.csv`:

| Three-edge case ID | Left | Bottom | Right | Rear (top) |
|---|---|---|---|---|
| `H17_SSS` | S | FREE | S | S |
| `H18_REAR_FIXED` | S | FREE | S | F |
| `H19_ONE_SIDE_FIXED` | F | FREE | S | S |
| `H20_BOTH_SIDES_FIXED` | F | FREE | F | S |
| `H21_ADJ_FIXED` | F | FREE | S | F |
| `H22_ALL_FIXED` | F | FREE | F | F |

The original four-edge IDs (`SSSS`, `F_y0`, `F_x0`, `FF_y`, `FF_x`, `FF_adj`,
`FFF_x`, `FFF_y`, `FFFF`) and their orientations are unchanged from the
previous package. Mirror images of asymmetric boundary cases require users
to rotate the plate/interpret x-y moment location accordingly.

## Python example

Standard library only. From the extracted directory:

```bash
python czerny.py --case SSSS --lx 4 --ly 6 --load 10 --E 33000000 --thickness .18
python czerny.py --case H17_SSS --lx 2.5 --ly 1.5 --load 9.5 --E 33000000 --thickness .18
python czerny.py --case H22_ALL_FIXED --lx 2.5 --ly 1.5 --load 7.5 --E 33000000 --thickness .18
python -m unittest discover -v -p 'test_*.py'
```

```python
from czerny import lookup, calculate, load_all_tables

four_edge_factor, info = lookup('SSSS', 1.5, 'mx_mid')  # 1 / 13.7
hahn_factor, info = lookup('H17_SSS', 0.6, 'mx_free_mid')  # 1 / 9.2
print(info['scale_basis'])  # 'K', with K = p*lx*ly

r = calculate('H17_SSS', lx=2.5, ly=1.5, load=9.5,
              E=33e6, thickness=0.18)
print(r['responses']['mx_free_mid'])  # ~3.87 kNm/m
print(r['responses']['reaction_x_total'])  # ~11.05 kN per side
```

`lookup` returns the *physical multiplier* and provenance metadata, not the
printed denominator. `calculate` scales to dimensional outputs and retains
source/sign/normalization per result. Cases are selected by ID, not inferred
from an unordered collection of edge flags.

**Interpolation:** exact published nodes are never marked interpolated. For
any intermediate ratio, default method linearly interpolates the physical
coefficient `1/k` (or direct multiplier); use `interpolate='denominator'` to
interpolate printed denominators instead. This method has **not** been
independently certified by the source authors. Values beyond each original
table range raise `ValueError`.

## Reproducibility

- Four-edge CSV: `python _build_dataset.py`, using `_full_grid.py`.
- Three-edge CSV: `python _build_three_edge.py`.
- Tests: `python -m unittest discover -v -p 'test_*.py'`.

These are transcription and scaling tests, **not** a second independent
validation of every printed coefficient. Some table entries are notably
non-monotone; none has been silently smoothed. Check original page scans and
source sketches before evaluating any real plate.

## Engineering limitations

NOT engineering/design-certified. Assumes elastic plate behavior for the
source cases and does not check reinforced-concrete crack stages, corner
uplift, real edge restraint, linear/nonlinear support flexibility, code
ULS/SLS, creep, deflection limits, reinforcement, or seismic actions.
Signs marked `±` describe magnitudes only; actual direction and symmetry
require reviewing the source sketch. Neither code nor CSV gives all spatial
moment/shear values across the slab. This material is for research and
independent cross-checking, not an automatic concrete-slab design.


## Numerical reconstruction pilot (new)

See [`NUMERICAL_VALIDATION.md`](NUMERICAL_VALIDATION.md) for an independently
**executed** SciPy finite-difference reconstruction of all 21 four-edge SSSS
aspect ratios, comparison CSVs, and a PyNiteFEA 3.2.x `Rect` plate benchmark
script prepared for local execution. The PyNite script has been syntax-checked
but **has not been run** because PyNite could not be installed in this
execution environment. Do not mistake the finite-difference results for PyNite.

```bash
# after installing numpy/scipy locally
python run_fd_all_ssss.py

# after installing PyNiteFEA==3.2.1 locally
python pynite_czerny_benchmark.py --case SSSS --ratio 1.5 --meshes 7 11 15
```
