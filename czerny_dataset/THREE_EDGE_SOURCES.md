# Three-edge source and verification notes

## Dataset status

**Six three-supported-edge cases transcribed from a secondary source.** The file
`czerny_three_edge_coefficients.csv` contains 456 printed numeric entries for
**uniform full-area loading** from J. Hahn, *Durchlaufträger, Rahmen, Platten
und Balken auf elastischer Bettung* (1970), Tafeln 17–22. These are not the
full original Czerny 1996 tables; Hahn says Tafel 22's moment coefficients
are after Czerny, *Beton-Kalender* (1961).

Source scan (University of Ljubljana):
https://fgg-web.fgg.uni-lj.si/KMLK/Sebastjan/masivni%20objekti/Hahnove%20tabele.pdf

| Case ID | Table | Printed pages | PDF page index, 0-based | Fourth edge | Side/rear edge conditions | Ratios `ly/lx` | Numeric entries |
|---|---:|---:|---:|---|---|---|---:|
| `H17_SSS` | 17 | 373–374 | 13 | free at y=0 | S left, S right, S rear | 0.25–1.50 | 112 |
| `H18_REAR_FIXED` | 18 | 375 | 14 | free at y=0 | S left, S right, F rear | 0.25–1.50 | 84 |
| `H19_ONE_SIDE_FIXED` | 19 | 376 | 15 | free at y=0 | F left, S right, S rear | 0.25–2.00 | 55 |
| `H20_BOTH_SIDES_FIXED` | 20 | 377 | 16 | free at y=0 | F left, F right, S rear | 0.25–2.00 | 55 |
| `H21_ADJ_FIXED` | 21 | 378 | 17 | free at y=0 | F left, S right, F rear | 0.25–2.00 | 66 |
| `H22_ALL_FIXED` | 22 | 379 | 18 | free at y=0 | F left, F right, F rear | 0.25–1.50 | 84 |

The orientation is explicit: `(x=0, y=0, x=lx, y=ly)`, with `y=0` free;
`F` always means fixed/clamped, never free. `ly/lx` is not forced to be >=1.
Mirroring left and right gives equivalent restraints but changes where to
interpret edge-specific results.

## Exactly what the three-edge CSV captures

Only `Lastfall 1` (uniform full-area load). Hahn also has triangular loads,
edge loads, and applied edge moments, all excluded to avoid mixing incompatible
loading bases. All numeric rows of the *uniform-load* table sections have
been entered, including reaction/deflection parameters available for Tafel 17
and the dimensionless beta in Tafel 18. Other tables do not provide a
full shear-force or support-reaction table; **a missing result is not zero**.

For Tafel 17 and Tafeln 18–22 (except direct multipliers), Hahn defines the
total distributed load `K = p*lx*ly` and `M = K/k`. These coefficients **cannot
be directly pasted** into the four-edge Czerny convention `M = p*lx^2/k`.
The Python module keeps the datasets separate and automatically uses the
correct scaling convention. `±` is a magnitude only; use source diagrams to
assign actual directions and corner signs.

Specific nonmoment entries:

- Tafel 17 `nu_x`, `nu_y`: multipliers `K_x = nu_x*K`, `K_y = nu_y*K`
  for support-resultant forces (see distribution in original drawing).
- Tafel 17 `omega_r`: reciprocal deflection factor `w_r = K*lx²/(omega_r * E*h³)`
  from source formula; check the original reference-point convention.
- Tafel 17 corner force magnitudes: `R1 = 2 * m_xy1` (free corner),
  `R2 = 2 * m_xy2` (rear supported corner). The code returns moments only,
  from which the point forces can be derived with sign conventions.
- Tafel 18 `beta`: dimensionless reach of negative `My` along the side;
  the program also reports `beta*ly` in metres.
- Tafeln 19–20: the row labelled `m_y` corresponds to an off-center `My,max`
  near the rear edge for `eps >~ 0.8` (T19) or `eps >~ 0.6` (T20), as noted
  below the original tables.

Original sign and k conventions are stored in the CSV in `sign`, `value_kind`,
`scale_basis`, `original_symbol` and `notes` columns. The program does *not*
calculate shear-force diagrams or distinguish positive and negative corner
force directions beyond the printed annotations.

## Completeness and QA limitations

- **456 source-table entries transcribed**, using the visually readable PDF
  page images, with multiple anchor-value checks and Hahn's worked example.
- The 15 Python regression tests cover data shape, representative printed
  nodes, load scaling, signs, interpolation bounds, and the Tafel 17 worked
  example. They **do not** constitute independent cell-by-cell verification.
- The printed Tafel 21 `m_ym` row is particularly non-monotone. Values were
  kept as seen, including a sharp increase at `eps=0.25`, rather than
  smoothing or fitting them. Check directly against a print copy before use.
- The Hahn three-edge grids are **0.1 or 0.2-spaced and nonuniform**, not
  the 0.05 Czerny four-edge grid. Values between printed nodes returned by
  the program are *interpolated estimates*, never labelled printed entries.
- No three-edge coefficients have been fabricated for the ratios outside
  each table's source range or for missing outputs.

## Where the original Czerny tables may be found

- F. Czerny, *Tafeln für Rechteckplatten*, *Beton-Kalender* (1996), Chapter 4
  starts around printed p. 309. Publicly indexed reproduction appears to show
  the three-edge cases, but scanned download/complete access has not been
  verified:
  https://www.studocu.com/de-at/document/universitat-fur-bodenkultur-wien/konstruktion-flachentragwerke/czerny-plattentafel/97365742
- F. Czerny, *Tafeln für vierseitige und dreiseitig gelagerte Rechteckplatten*
  (1978), Ernst & Sohn, library record:
  https://books.google.com/books/about/Tafeln_f%C3%BCr_vierseitige_und_dreiseitig_g.html?id=lYj-PgAACAAJ
- The available original Czerny PDF fragment, printed pp. 277–289, only
  includes nine **four-edge** cases:
  https://www.klancek.si/sites/default/files/datoteke/files/tabele.pdf

**Engineering caution:** Use only for studying or checking elastic plate
coefficients; not validated for structural design. RC cracking, restraint
stiffness, uplift, torsional edge effects, detailing and code checks are not
modelled. Independently review all numbers and sign conventions.
