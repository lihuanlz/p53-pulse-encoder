# Third-party data provenance

This paper's analysis re-uses published single-cell p53 records. **No large third-party archive is re-hosted in this repository.** The pipeline consumes extracted tables that are part of our own audit trail (`results/`), traceable to the following sources:

| Source | Used for | Access |
|---|---|---|
| Batchelor, Loewer, Mock & Lahav, *Mol. Syst. Biol.* 7, 488 (2011) | p53 pulse amplitudes, dose-invariance of amplitude | Publisher supplementary; our extraction is cross-checked against the tables in `results/` |
| Stewart-Ornstein & Lahav, *Sci. Signal.* 10, eaah6671 (2017) (PMC5504473) | 12-cell-line period invariance (~5 h) | Publisher supplementary |
| Geva-Zatorsky et al., *Mol. Syst. Biol.* 2, 33 (2006) | Pulse counting across doses | Publisher supplementary |
| Loewer et al., *Cell* 142, 89–100 (2010) (PMC3003696) | Spontaneous-pulse baseline and dispersion-ratio reference (code85/87/88) | Publisher site |
| I. Belgacem, arXiv:2605.23722 (2026) | Closed-form cross-validation target | arXiv |
| Wang, Liu & Zhou, *IET Syst. Biol.* 13, 180–185 (2019); Monk, *Curr. Biol.* 13, 1409–1413 (2003); Batchelor et al., *Mol. Cell* 30, 277–289 (2008); Mönke et al., *Sci. Rep.* 7, 46571 (2017) | Molecular-layer constants: transcription/translation delays, half-lives (constant registry v02) | Publisher sites |

If you re-run the extraction layer, download the source supplementary files from the publishers above; the derived numbers checked against the manuscript are all under `results/`.
