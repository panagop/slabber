# Datasets, quantities and normalization

## Historical source files

`czerny_dataset/czerny_coefficients.csv` (four-edge) — 1,701 rows,
columns: `case_id`, `ratio_ly_lx`, `quantity`, `original_symbol`,
`measure_kind`, `sign`, `table_value`, `value_kind`, `source_page`, `source_url`.

`czerny_dataset/czerny_three_edge_coefficients.csv` (three-edge) — 456
rows, including `scale_basis`, `source_table`, `source_pdf_page` and `notes`.

`czerny_dataset/coverage.csv` — the 15 cases, their boundary codes, source
pages, ratio ranges, parameter counts and coverage status.

`czerny_dataset/reference_ssss_fd.csv` — 3 ratios × 4 refinement levels;
`reference_ssss_fd_all_ratios.csv` — all 21 ratios with a 40-division
short-side mesh. These are **executed finite-difference results**, not
Czerny original coefficients and not PyNite output.

## Physical and coordinate convention

- `lx`, `ly`, `h`, `w`: metres; `E`: kN/m²; `p` or `load`: kN/m².
- Moment per unit width: kNm/m; transverse shear per width: kN/m;
  total support reactions/corner forces: kN.
- Four-edge: `lx <= ly`; aspect ratio `ly/lx` spans 1.0–2.0.
- Three-edge: `ly/lx` may be below 1; **y=0 is free**; boundary order
  `x=0,y=0,x=lx,y=ly`. No rotation/orientation is implicit.
- `S`: simply supported; `F`: fixed/clamped; `FREE`: free edge. `F` in
  **case IDs** is not free.
- Historical values are commonly printed as denominator `k`, e.g.
  `M = p*lx^2/k`; the code's `lookup()` returns a **signed multiplier**
  (approximately `1/k`) and metadata containing the printed `k`.

## Source-specific formulas

**Four-edge Czerny:**

```text
Mx, My, Mxy, corner forces  = signed p*lx^2/k   [see original quantity]
Qx, Qy, line reactions       = signed p*lx/k
w_mid                        = c*p*lx^4/(E*h^3)
```

**Three-edge Hahn:**

```text
K                         = p*lx*ly
M                         = signed K/k
H17 support resultants    = nu_x*K, nu_y*K
H17 free-edge deflection  = K*lx^2/(omega_r*E*h^3)
H18 negative-region reach = beta*ly
```

Conventions are saved in `value_kind`/`scale_basis`. `sign=±` denotes a
**magnitude only**: source diagrams are required to establish actual
direction. A missing row is **not zero**; avoid constructing a fake full
output suite from partial Hahn tables.

## Interpolation and the audit boundary

Exact published ratios are reported as `interpolated=False`.
Intermediate ratios are linearly interpolated by default in the *physical
response factor* (the reciprocal printed denominator), with
`interpolate='denominator'` as an explicit alternative. Never extrapolate
beyond the printed ratio range. An interpolated value is **not** an
original published value.

Original source rows have been visually transcribed and regression tested,
but have **not** received independent cell-by-cell transcription review.
A source audit should record every scanned page/row/column and a second
reader's verification. See `docs/VALIDATION_AND_LIMITATIONS.md`.
