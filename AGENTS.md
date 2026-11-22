# Codex instructions for Slabber

Read `README.md` and `docs/PROJECT_STATUS.md` first. **This is a research
workspace**, not a certified structural-analysis library.

## Non-negotiable rules

1. Preserve the original CSVs and the inherited scripts. Before changing a
   source coefficient, compare it visually with its **printed page**, record
   page/table/cell, and add a regression test. Prefer append-only annotations
   or a new versioned derived CSV instead of silently editing source values.
2. Distinguish *source-transcribed*, *interpolated*, *independent finite-
   difference*, and *PyNite FEM* outputs. Do not describe an unexecuted
   script as a passed numerical validation.
3. Do not conflate Czerny's four-edge `p*lx^2` moment normalization with
   Hahn's three-edge `K=p*lx*ly` normalization.
4. Distinguish `Qx/Qy`, edge reaction distributions, integrated support
   resultants, Kirchhoff effective shear, and corner point forces.
5. Check signs, axes, units, physical reference locations, edge ordering,
   free-corner behavior, and mesh convergence before comparing coefficients.
6. Use `uv` for dependency/environment commands. Save new generated outputs
   in `results/` and give every accepted benchmark reproducible parameters.
7. No extrapolation past original aspect-ratio ranges without explicitly
   implementing a new numerical model and labelling the output as such.
8. Do not claim code compliance or provide reinforcement design as a direct
   inference from elastic historical coefficients.

## Existing checks

```bash
uv sync
uv run python -m unittest discover -s czerny_dataset -p 'test_*.py' -v
uv run python scripts/verify_data.py
```

`uv sync` installs `PyNiteFEA==3.2.1` with the other project dependencies.
The PyNite and SciPy FD four-edge comparison has been executed at three
aspect ratios for every four-edge case. The PyNite and Rayleigh–Ritz
comparison covers the six Hahn three-edge cases at three printed ratios.
See `docs/PROJECT_STATUS.md` and `docs/REPRODUCIBILITY.md`.

## Next numerical milestone

Review the remaining Czerny and Hahn discrepancies identified in
`docs/PROJECT_STATUS.md`. Follow `docs/NUMERICAL_ROADMAP.md` for
subsequent validation work.

## Documentation expectations

- Update `docs/PROJECT_STATUS.md` after new verified achievements.
- Update `docs/NUMERICAL_ROADMAP.md` and `docs/VALIDATION_AND_LIMITATIONS.md`
  when methods or risk assessments change.
- Record package versions, scripts, cases, units, mesh sizes, and test results.
- Preserve source-level traceability; link to the underlying scan where possible.
