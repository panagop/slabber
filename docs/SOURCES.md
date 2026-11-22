# Historical sources and software references

Accessed/rechecked as of 2026-11-22 where public search records were
available. Links may later require access or return different scans.

## Original Czerny publication and subsequent editions

- **F. Czerny (1955)**, *Tafeln für gleichmäßig vollbelastete Rechteckplatten*,
  *Bautechnik-Archiv*, Heft 11, Wilhelm Ernst & Sohn, Berlin.
  Bibliographic corroboration in older engineering references and secondary
  [biographical bibliography](https://de.wikipedia.org/wiki/Friedrich_Czerny).
  **The 1955 original is NOT bundled or transcribed in this project.**
- **F. Czerny, Beton-Kalender (1996)**, *Tafeln für Rechteckplatten*.
  The [public excerpt](https://www.klancek.si/sites/default/files/datoteke/files/tabele.pdf)
  contains printed pages 277–289, notably nine four-edge tabulations
  on pages **281–289**. The full three-edge chapter is **not in this file**.
  Dataset: `czerny_dataset/czerny_coefficients.csv`.
- **F. Czerny (1978)**, *Tafeln für vierseitige und dreiseitig gelagerte
  Rechteckplatten*, Ernst & Sohn.
  [Google Books bibliographic record](https://books.google.com/books/about/Tafeln_f%C3%BCr_vierseitige_und_dreiseitig_g.html?id=lYj-PgAACAAJ).
  A possible route to a definitive primary three-edge source; **not
  available as a verified complete transcription here**.
- A searchable [1996 chapter reproduction](https://www.studocu.com/de-at/document/universitat-fur-bodenkultur-wien/konstruktion-flachentragwerke/czerny-plattentafel/97365742)
  reportedly includes original three-edge pages from approximately p. 309;
  download access, completeness and OCR correctness remain **unverified**.

## Secondary source for three supported edges

- **J. Hahn (1970)**, *Durchlaufträger, Rahmen, Platten und Balken auf
  elastischer Bettung*, Werner Verlag, Tafeln 17–22, printed pp. 373–379:
  [University of Ljubljana-hosted scan](https://fgg-web.fgg.uni-lj.si/KMLK/Sebastjan/masivni%20objekti/Hahnove%20tabele.pdf).
  These six cases have a free edge. Tafel 22 attributes its bending
  coefficient values to Czerny (1961 *Beton-Kalender*). The code preserves
  **Hahn** as source for all these values, without mislabelling them original
  1996 Czerny coefficients. Dataset:
  `czerny_dataset/czerny_three_edge_coefficients.csv`.

## Solver / workflow documentation

- [PyNite plate element documentation](https://pynite.readthedocs.io/en/stable/plate.html).
- [PyNite FEModel3D API](https://pynite.readthedocs.io/en/stable/FEModel3D.html).
- [PyNiteFEA 3.2.1 distribution](https://pypi.org/project/PyNiteFEA/3.2.1/).
- [uv: working on projects](https://docs.astral.sh/uv/guides/projects/).
- [uv: dependency groups / optional dependencies](https://docs.astral.sh/uv/concepts/dependencies/).
- [Codex IDE workflow](https://learn.chatgpt.com/docs/codex/ide).
- [Agent instruction files](https://learn.chatgpt.com/docs/agent-configuration/agents-md).
- Alternative future solvers:
  [xdfem2D](https://github.com/pcachim/xdfem2D),
  [scikit-fem examples](https://scikit-fem.readthedocs.io/en/stable/listofexamples.html).

## Copyright / availability policy

PDFs have not been repackaged or redistributed. This project preserves
citations, reconstructed research data, source-specific page references,
and the Python work. If a source goes offline, use a licensed library
copy rather than silently changing the provenance of the dataset.
