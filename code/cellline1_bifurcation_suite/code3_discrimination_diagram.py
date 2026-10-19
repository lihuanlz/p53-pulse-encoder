# -*- coding: utf-8 -*-
"""
code 3: discrimination diagram (pulse observables vs circuit generic class)
==========================================
draw the 'three-signature discrimination theorem' of manuscript section 4 as an operational diagram:
  x-axis alpha = dlnA/dlnD (amplitude-dose sensitivity, log axis)
  y-axis nu = pulse-number increment per decade of dose (count sensitivity)

points fall in two classes:
  - model points (from code 2 / Track B computations) - define the two generic-class regions
  - experimental points (semi-quantitative, from published statistics) - test the theorem's predictions

derivation of the experimental point values (same as manuscript Methods):
  p53/γ-IR (Lahav 2004): 0.1–10 Gy (100×) amplitude 'dose-independent',
      taking a 15% amplitude-variation cap -> |alpha| <= ln1.15/ln100 ~ 0.03;
      count ~1.5 -> 5.5 (2 decades) -> nu ~ 2
  p53/NCS (Mönke 2017): 25→400 ng/ml (16×) mean amplitude independent -> |alpha| <~ 0.05; count increases -> nu ~ 1.8
  p53/UV (Batchelor 2011): amplitude proportional to dose (graded) -> alpha ~ 1; single pulse -> nu ~ 0
      - built-in control: same protein, same cell, experimentally judged 'non-excitable', signature lands in the analog region

Run: python3 code3_discrimination_diagram.py   (a few seconds, outputs the discrimination-map PNG)
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
REPO = Path(__file__).resolve().parent.parent.parent

fig, ax = plt.subplots(figsize=(8.5, 7))

# ---------- regions of the two generic classes ----------
ax.axvspan(0.3, 15, alpha=0.09, color='red')     # NF analog region
ax.axvspan(5e-4, 0.15, alpha=0.09, color='blue') # EXC digital region
ax.text(2.2, 3.4, 'NF family (analog)\n$\\alpha\\gtrsim0.5$\n(subcritical: +hysteresis)',
        color='darkred', fontsize=10, ha='center')
ax.text(0.008, 3.4, 'EXC family (digital)\n$\\alpha\\lesssim0.15$', color='darkblue', fontsize=10, ha='center')
ax.axvline(0.22, color='gray', ls='--', lw=1)

# ---------- model points (Track B results, nu in pulses/decade) ----------
ax.errorbar([8.76], [0.02], fmt='s', color='red', ms=9, label='NF model, near onset (alpha->inf)')
ax.errorbar([0.7], [0.05], fmt='s', color='red', ms=9, mfc='none', label='NF model, far from onset')
ax.errorbar([0.014], [4.05], fmt='o', color='blue', ms=9, label='FHN slow drive (pure digital)')
ax.errorbar([0.26], [2.0], fmt='^', color='green', ms=9, label='Morris-Lecar fast drive (leaky)')

# ---------- experimental points (semi-quantitative) ----------
ax.errorbar([0.03], [2.0], xerr=[0.02], yerr=[0.5], fmt='*', color='purple',
            ms=19, capsize=4, label='p53 / γ-IR (Lahav 2004) : 0.1–10 Gy')
ax.errorbar([0.05], [1.8], xerr=[0.04], yerr=[0.5], fmt='*', color='magenta',
            ms=17, capsize=4, label='p53 / NCS (Mönke 2017) : 25–400 ng/ml')
ax.errorbar([1.0], [0.05], xerr=[0.3], yerr=[0.05], fmt='D', color='darkorange',
            ms=10, capsize=4, label='p53 / UV (Batchelor 2011): graded, non-excitable')
ax.annotate('built-in control: same protein, same cell\ntopology<->signature strictly covary', (0.8, 1.35), fontsize=8.5, color='darkorange')

ax.set_xscale('log'); ax.set_xlim(5e-4, 15); ax.set_ylim(-0.3, 4.6)
ax.set_xlabel('$\\alpha = \\partial\\ln A/\\partial\\ln D$(amplitude-dose sensitivity)')
ax.set_ylabel('$\\nu$ (pulses per decade of dose) ')
ax.set_title('discrimination map: experimental points (semi-quantitative; see SI for precise digitization)', fontsize=11)
ax.legend(fontsize=8.5, loc='center right')
plt.tight_layout()
plt.savefig(str(REPO / "results" / "TrackC_discrimination_map_experimental_points.png"), dpi=200)
print('saved TrackC_discrimination_map_experimental_points.png')
print('\nreading: purple/magenta stars fall deep in the EXC region (digital), orange diamonds in the analog region (control) - no quadrant misplacement.')
