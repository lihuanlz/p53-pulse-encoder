# Pipeline: the paper in order

The script numbers are **historical identifiers**. They come from the original
cell-line numbering of the parent project and are cited verbatim in the Methods
(`codes 13, 50--52, 81, 82, 89--93`). Treat them as stable citation keys, not as
a sequence: renumbering would break the manuscript, the SI pipeline map and the
archive alias table at once.

This file gives the order in which the paper presents the work.

## At a glance

| Paper position | Scripts | Lands in |
|---|---|---|
| §2 separation theorem, Lemma 1--3 | `code/cellline1_bifurcation_suite/code2`, `code3` | Fig. 1a--d |
| §3 the published record decides the class | `code/code52` + `data_external/batchelor2011/digitized_Fig1GH.json` | Eq. (3), Fig. 1e, Table 1 |
| §4 channel census | `code/code51r1`--`r4`, `code/code81` | Fig. 2, Eqs. (4)--(6) |
| §5 molecular layer | `code/code82`, `code/code89`--`93` | Figs. 3--4, Table 3 |
| §6 the decisive experiment | — | no computation |

## Step by step

| # | Paper position | Script | Reads | Writes |
|---|---|---|---|---|
| 1 | Lemma 1--2 numerics, Fig. 1a--d | `cellline1_bifurcation_suite/code2` | — | `TrackB_bifurcation_validation.png` (beside the script) |
| 2 | discrimination grid | `cellline1_bifurcation_suite/code3` | — | figure beside the script |
| 3 | Eq. (3): the two alpha estimates | `code52` | `data_external/batchelor2011/digitized_Fig1GH.json` | `results/code52_TrackC2_results.json`, `results/code52_TrackC2_alpha_nu_discrimination.png` |
| 4 | Fisher spectrum, Eq. (6) | `code51r3`, `code51r4` | — | `results/code51r3_p53_Fisher_audit_results.json`, `results/code51r4A_noise_calibration_robustness.json`, `results/code51r4B_Fisher_M4000.json` |
| 5 | channel shares, Eq. (4) | `code51r3` | — | (same JSON; `fisher_h2.share`) |
| 6 | section 4.5 digital limit | `code51r1`, `code51r2` | — | `results/code51_p53_Fisher_audit_results.json`, `results/code51r2_...json` |
| 7 | cascade information identity | `code81` | — | `results/code81_cascade_audit_results.json` |
| 8 | mechanism-level recheck (the 45.1%/54.9% split) | `code82` | — | `results/code82_mechanism_audit_results.json` |
| 9 | period and phase shares | `code89` | — | `results/molecular_layer/code89_molecular_decomposition_results.json` |
| 10 | 14 static links, onset decomposition | `code90` | — | `results/molecular_layer/code90_cascade_audit_results.json` |
| 11 | FWHM decomposition | `code91` | — | `results/molecular_layer/code91_pulse_width_results.json` |
| 12 | closed-form cross-validation | `code92` | — | `results/molecular_layer/code92_closedform_crossvalidation_results.json` |
| 13 | TCS--logistic correspondence | `code93` | — | `results/molecular_layer/code93_tcs_logistic_correspondence_results.json` |

Supporting analyses that back individual threshold choices but generate no
figure: `code50` (cascade framework self-check), `code83`--`code88`
(digital-limit variational argument, heterogeneity recheck, foci trigger and
termination, clock architecture).

## Running

```bash
python reproduce.py              # every step, in the order above
python reproduce.py section4     # one group
python reproduce.py --list       # what the groups are
```

Environment: Python 3.10+, `numpy`, `scipy`, `matplotlib`. No network access is
needed; every input is inside the repository.

## Two things worth knowing

- **`code51r4` runs in modes.** `python code51r4...py A` (the default) does the
  robustness sweep in seconds; `... B lo hi` recomputes the `M=4000` spectrum in
  segments with a resumable cache, and is the only slow step; `... FIG` redraws
  the four-panel figure from the two JSONs. `reproduce.py` runs A and FIG; run B
  only if you want to regenerate `code51r4B` (the shipped one is archived).
- **The archived result files are the output of these exact scripts.** Re-running
  them reproduces every number; a handful of human-readable strings in the JSON
  differ from the shipped ones where the prose was edited after generation.
- **`code82` is a Monte Carlo audit.** The mechanism-level split that the paper
  quotes is stable across runs (`share_time` 0.451137 vs 0.451139), but
  near-singular quantities move: the Fisher condition number of the 23-parameter
  model varies by roughly 10% between runs. The paper quotes the shares, not
  that condition number.
