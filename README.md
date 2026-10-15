# The p53 pulse encoder counts by topological necessity

Companion code, results and manuscript source for the paper:

> H. Li. *The p53 pulse encoder counts by topological necessity: theory, evidence, and the molecular origin of the 5.5-hour clock.* Companion manuscript (2026), split from the umbrella project [cell-communication-audit](https://github.com/lihuanlz/cell-communication-audit) (archived at Zenodo, doi:10.5281/zenodo.23008824).

## Layout

- `manuscript/` — LaTeX source, compiled PDF, four figure PNGs, and one figure-generation script per figure (`make_fig1.py` … `make_fig4.py`; fig3/fig4 read the audited code89/code92 result JSONs under `results/molecular_layer/`).
- `code/` — analysis scripts (ours): Fisher stratification (code13), cascade self-check (code50), Fisher audit series (code51r1–r4), discrimination map (code52), full-cascade audit (code81), mechanism audit (code82), digital-limit and cascade-model studies (code83–88), molecular-constant decomposition and cross-validation (code89–93), plus the cell-line-1 bifurcation-and-discrimination suite cited in Methods.
- `results/` — derived outputs (ours): JSON result files, verdict figures, molecular-layer tables.
- `data_external/` — third-party data used; see `data_external/PROVENANCE.md`. Large source archives are not re-hosted; the file lists exact sources.

## Running

Scripts are plain Python 3 (numpy/scipy/matplotlib). Each `codeNN_*.py` is self-contained against the tables in `results/` and `data_external/`; figure scripts write the manuscript PNGs.

## Ownership

All code, derived results and the manuscript are the author's work product, released under MIT (code) with the manuscript text under CC BY 4.0. Third-party datasets are redistributed only where small and license-clear, and always with provenance; everything else is referenced, never copied. See `data_external/PROVENANCE.md`.

## Companion papers

- *Amplitude and midpoint statistics of the same chemotaxis cells demand incompatible receptor models* → [chemotaxis-receptor-layer](https://github.com/lihuanlz/chemotaxis-receptor-layer)
- *The quantitative currency of GPCR bias signalling fails its own algebraic identities* → [gpcr-field-audit](https://github.com/lihuanlz/gpcr-field-audit)
- *The price of absolute scale: degeneracy groups, surviving statistics, and the reference economics of cellular signalling* → [scale-free-framework](https://github.com/lihuanlz/scale-free-framework)
