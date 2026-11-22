# Research timeline / decisions

This document summarizes the work performed before the uv handoff, not a
claim that all linked material has been independently verified.

1. **Identify primary Czerny**: a 1955 *Bautechnik-Archiv* publication,
   subsequent *Beton-Kalender* tabulations, and a 1978 book record.
   Accessible 1996 excerpt covers the four-edge pages 281–289.
2. **Transcribe nine four-edge cases**: started as a partial dataset,
   completed 21 original aspect-ratio nodes (1.00, 1.05, ..., 2.00) for
   all nine cases: **1,701 printed numeric entries**. Denominators and
   original meaning retained in long-form CSV.
3. **Find three-edge source**: located 1970 Hahn tables 17–22 (six
   boundary configurations). They are secondary-source tables; Tafel 22
   references Czerny 1961. **456 numeric uniform-load entries** added;
   no invented 0.05 increments on Hahn's irregular grid.
4. **Write lookup and scaling module**: `czerny.py`, including provenance
   metadata, `lookup`, `calculate`, bounded interpolation, and two
   incompatible-but-explicit source normalizations.
5. **Select numerical approach**: PyNiteFEA chosen for first thin-plate
   `Rect` FEM check; xdfem2D and scikit-fem explicitly deferred.
6. **Execute independent reference**: SciPy finite-difference solver
   implemented and run for `SSSS`. Saved 3-ratio refinement study and
   21-ratio comparison. Moment **extremum location** proved important:
   `My_max` can occur away from plate center for large aspect ratios.
7. **Prepare PyNite harness**: model builder and command line created,
   source codes mapped to restraints; script syntax-checked only.
   Neither PyNite installation nor FEM simulation was completed in the
   former environment; first actual run is the recommended next step.
8. **Create uv/VS Code/Codex handoff**: preserve the raw research artifacts
   and 19 inherited passing tests, add project metadata, docs, task
   launchers and agent instructions, without mutating historical rows.

## Open questions for further research

- Recover the **original** complete three-edge *Czerny* 1996 chapter;
  reconcile overlapping rows against Hahn's secondary-source tables.
- Verify every transcribed cell with a second read of the original scans.
- Confirm `Rect` element sign/orientation, shear-output conventions,
  reaction balance, mesh convergence and the PyNite 3.2.1 API.
- Recover support-line moment/shear reference locations from each printed
  plate sketch, especially for fixed/free edges and corners.
- Decide when scikit-fem (Kirchhoff elements) or xdfem2D should be
  independently trialled, after the first validated PyNite pilot.
