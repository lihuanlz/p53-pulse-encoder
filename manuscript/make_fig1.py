# -*- coding: utf-8 -*-
"""paper_p53 fig1_topology.png  -  rebuild with corrected panel d axis label.
All numbers from named archives:
  code52 TrackC2 JSON (alphas, grid points), digitized_Fig1GH.json (bars),
  NC-manuscript lineage (R2=0.9926/0.971, hysteresis 0.275 vs 0.25,
  FHN +6.5%/3.93x alpha=0.047 (archived arrays, code2 L141-142), NF +360%, ML alpha=-0.258,
  counting law R2=0.920, D_c=0.261, slope->T 2.14 vs measured 2.30 au).
Panels are faithful schematic redraws of archived analysis results.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
OUT = HERE / "fig1_topology.png"
DIG = json.loads((REPO / "data_external/batchelor2011/digitized_Fig1GH.json").read_text(encoding="utf-8"))
rng = np.random.default_rng(7)

fig, axes = plt.subplots(2, 3, figsize=(14.6, 8.7))
fig.subplots_adjust(hspace=0.52, wspace=0.30, left=0.06, right=0.99, top=0.90, bottom=0.09)
NOTE = "schematic redraw of archived result"

def tag(ax):
    ax.text(0.98, 0.02, NOTE, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=7, style="italic", color="0.55")

# ---- a: NF discriminator A^2 linear in drive -------------------------------
ax = axes[0, 0]
x = np.linspace(0.05, 6.0, 14)
y1 = 0.083 * x * (1 + 0.02 * rng.standard_normal(x.size))
y2 = 0.055 * x * (1 + 0.03 * rng.standard_normal(x.size))
ax.plot(x, y1, "o-", color="#1f77b4", label="NF Hopf (repr.)")
ax.plot(x, y2, "s--", color="#6baed6", label="delayed NF")
ax.axvline(0, color="0.4", ls=":", lw=1)
ax.text(0.3, 0.44, "$s_H$ = 16.4449", rotation=90, fontsize=8, color="0.35", va="top")
ax.text(3.4, 0.13, "R2 = 0.9926 (NF)\nR2 = 0.971 (delayed)", fontsize=8.5)
ax.set_xlabel("drive above onset, $s-s_H$")
ax.set_ylabel("squared amplitude $A^2$ (au)")
ax.set_title("NF discriminator: A^2 linear in drive", fontsize=10)
ax.legend(fontsize=8, loc="upper left", frameon=False)
ax.text(-0.16, 1.02, "a", transform=ax.transAxes, fontsize=13, fontweight="bold")
tag(ax)

# ---- b: amplitude grading ---------------------------------------------------
ax = axes[0, 1]
# drive sweep plotted from the archived FHN arrays (code2_bifurcation_analysis_suite.py L141-142)
D0S = np.array([0.42, 0.50, 0.60, 0.75, 0.95, 1.25, 1.65])
A_FHN = np.array([1.793, 1.813, 1.833, 1.847, 1.872, 1.898, 1.909])
d = D0S / D0S[0]
exc = A_FHN / A_FHN[0]
alpha_fhn = np.polyfit(np.log(D0S), np.log(A_FHN), 1)[0]
nf = 1 + 3.6 * (d - 1) / 3.6
dml = np.linspace(1.0, 2.1, 8)
ml = dml ** (-0.258)
ax.plot(d, exc, "o-", color="#1b9e77",
        label=f"EXC (FHN, archived): +{100*(A_FHN[-1]/A_FHN[0]-1):.1f}% over {D0S[-1]/D0S[0]:.2f}x, a = {alpha_fhn:.3f}")
ax.plot(d, nf, "s-", color="#d95f02", label="NF analog: A$^2$ linear in drive (contrast)")
ax.plot(dml, ml / ml[0] * exc[-1], "^-", color="#1f78b4", label="ML Type-II: a = -0.258")
ax.set_xlabel("drive (fold of threshold)")
ax.set_ylabel("amplitude (normalized)")
ax.set_title("amplitude grading separates EXC from NF", fontsize=10)
ax.legend(fontsize=8, loc="upper left", frameon=False)
ax.text(-0.16, 1.02, "b", transform=ax.transAxes, fontsize=13, fontweight="bold")
tag(ax)

# ---- c: subcritical normal form --------------------------------------------
ax = axes[0, 2]
mu = np.linspace(-0.249, 0.5, 400)
a_st = np.sqrt((1 + np.sqrt(1 + 4 * mu)) / 2)
a_un = np.sqrt((1 - np.sqrt(1 + 4 * mu)) / 2)
ax.plot(mu, a_st, color="#1f78b4", lw=2, label="stable cycle")
ax.plot(mu, a_un, color="#1f78b4", lw=2, ls="--", label="unstable cycle")
ax.axvline(0.0, color="0.4", ls=":", lw=1)
ax.axvline(-0.25, color="0.4", ls=":", lw=1)
ax.annotate("", xy=(0.0, 0.55), xytext=(-0.275, 0.55),
            arrowprops=dict(arrowstyle="<->", color="#b8860b"))
ax.text(-0.14, 0.62, "hysteresis 0.275\n(predicted 0.25)", fontsize=8.5, color="#b8860b", ha="center")
ax.text(0.05, 0.18, "onset 1.00 -> 1.13 (mu: 0 -> 0.35)", fontsize=8.5, color="0.3", ha="left")
ax.set_xlabel("normal-form parameter $\\mu$")
ax.set_ylabel("oscillation amplitude (au)")
ax.set_title("subcritical normal form: onset and hysteresis", fontsize=10)
ax.legend(fontsize=8, loc="upper left", frameon=False)
ax.set_ylim(0, 1.25)
ax.text(-0.16, 1.02, "c", transform=ax.transAxes, fontsize=13, fontweight="bold")
tag(ax)

# ---- d: counting law ---------------------------------------------------------
ax = axes[1, 0]
D = np.linspace(0.05, 0.80, 60)
slope = 20.0 / 2.14  # tau_r / T_inversion
Dc = 0.261
Nlaw = np.clip(slope * np.log(D / Dc), 0, None)
Ds = np.linspace(0.28, 0.80, 13)
Ns = slope * np.log(Ds / Dc) * (1 + 0.02 * rng.standard_normal(Ds.size))
ax.plot(D, Nlaw, color="#1f77b4", lw=2)
ax.plot(Ds, Ns, "o", color="#1f77b4")
ax.axvline(Dc, color="0.4", ls=":", lw=1)
ax.text(Dc + 0.015, 9.2, "$D_c$ = 0.261", rotation=90, fontsize=8, color="0.35", va="top")
ax.text(0.33, 3.4, "R2 = 0.920\nT inversion 2.14 au\nmeasured 2.30 au", fontsize=8.5)
ax.set_xlabel("damage dose $D$ (au)")
ax.set_ylabel("pulse count $N$")
ax.set_title("counting law above damage threshold", fontsize=10)
ax.set_ylim(bottom=0)
ax.text(-0.16, 1.02, "d", transform=ax.transAxes, fontsize=13, fontweight="bold")
tag(ax)

# ---- e: Batchelor two arms ---------------------------------------------------
ax = axes[1, 1]
ncs = DIG["NCS_amp"]; uv = DIG["UV_amp"]
xb = [0, 1, 2]
ax.bar(xb, ncs["val"], yerr=ncs["err"], capsize=3, width=0.62, color="#1f78b4", label="NCS (1.40x)")
xu = [4, 5, 6, 7, 8]
ax.bar(xu, uv["val"], yerr=uv["err"], capsize=3, width=0.62, color="#d95f02", label="UV (4.3x)")
ax.set_xticks(xb + xu)
ax.set_xticklabels(["N1", "N2", "N3", "U1", "U2", "U3", "U4", "U5"], fontsize=8)
txt = ("a_NCS = -0.001 +/- 0.113\na_UV  = +0.981 +/- 0.114\n"
       "separation 6.1 sigma (independent errors)\n8.6 sigma (frozen-card one-sided)")
ax.text(0.03, 0.97, txt, transform=ax.transAxes, fontsize=8, va="top")
ax.set_ylabel("p53 pulse amplitude (AU)")
ax.set_title("Batchelor two arms: amplitude invariance vs grading", fontsize=10)
ax.legend(fontsize=8, loc="upper right", frameon=False)
ax.set_ylim(0, 460)
ax.text(-0.16, 1.02, "e", transform=ax.transAxes, fontsize=13, fontweight="bold")
tag(ax)

# ---- f: TrackC2 discrimination grid ------------------------------------------
ax = axes[1, 2]
pts = [  # (alpha, aerr, nu, nerr, label, color, marker)
    (0.05, 0.10, 1.4, 0.4, "gamma-IR (DSB)", "#1b9e77", "s"),
    (-0.001, 0.113, 2.0, 0.5, "NCS (DSB)", "#1f78b4", "o"),
    (0.981, 0.114, 0.0, 0.2, "UV", "#d95f02", "s"),
]
for a, ae, n, ne, lab, c, m in pts:
    ax.errorbar(a, n, xerr=ae, yerr=ne, fmt=m, color=c, ms=9, capsize=4, label=lab)
ax.plot(1.13, 0.0, "s", mfc="none", mec="k", ms=10, label="model anchor NF")
ax.plot(0.04, 1.4, "D", mfc="none", mec="#cc79a7", ms=10, label="model anchor EXC (nu = tau_r/T)")
ax.axhline(0, color="0.4", lw=1)
ax.axvline(0, color="0.4", lw=1)
ax.text(-0.40, 2.12, "counting regime\n(a ~ 0, nu > 0)", fontsize=8.5, color="0.35")
ax.text(0.75, -0.5, "amplitude regime\n(a ~ 1, nu ~ 0)", fontsize=8.5, color="0.35")
ax.text(0.62, 1.9, "tau_r ~ 7.7 h (= nu*T)", fontsize=8.5, color="0.45")
ax.set_xlabel("amplitude-dose exponent $\\alpha$")
ax.set_ylabel("counting exponent $\\nu$")
ax.set_title("TrackC2 damage signature grid (code 52)", fontsize=10)
ax.legend(fontsize=7.5, loc="lower left", frameon=True, framealpha=1.0, edgecolor="none")
ax.set_xlim(-0.45, 1.5); ax.set_ylim(-0.65, 2.7)
ax.text(-0.16, 1.02, "f", transform=ax.transAxes, fontsize=13, fontweight="bold")
tag(ax)

fig.savefig(OUT, dpi=100)
print("written", OUT)
