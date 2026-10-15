# Third-party data provenance

This paper's analysis re-uses published single-cell p53 records. **No large third-party archive is re-hosted in this repository.** The pipeline consumes extracted tables that are part of our own audit trail (`results/`), traceable to the following sources:

| Source | Used for | Access |
|---|---|---|
| Batchelor, Loewer, Mock & Lahav, *Mol. Syst. Biol.* 7, 488 (2011) | p53 pulse amplitudes, dose-invariance of amplitude | Publisher supplementary; our extraction is cross-checked against the tables in `results/` |
| Stewart-Ornstein & Lahav, *Sci. Signal.* 10, eaah6671 (2017) (PMC5504473) | 12-cell-line period invariance (~5 h) | Publisher supplementary |
| Geva-Zatorsky et al., *Mol. Syst. Biol.* 2, 33 (2006) | Pulse counting across doses | Publisher supplementary |
| Loewer et al., *Science* 336, 1442–1444 (2010) | Spontaneous pulses | Publisher supplementary |
| Belgacem & Ismail, arXiv:2605.23722 (2026) | Closed-form cross-validation target | arXiv |
| Gillies et al., *Cell Syst.* 11, 445–460 (2020) / Finzel dissertation (2016) | Mdm2/Wip1 layer constants | Publisher / institutional repository |

If you re-run the extraction layer, download the source supplementary files from the publishers above; the derived numbers checked against the manuscript are all under `results/`.
