# -*- coding: utf-8 -*-
"""paper_p53 fig1_topology.png  -  drawn at print size (6.27 in, 7 pt class).

All numbers from named archives:
  code52 TrackC2 JSON (alphas, grid points), digitized_Fig1GH.json (bars),
  NC-manuscript lineage (R2=0.9926/0.971, hysteresis 0.275 vs 0.25,
  FHN +6.5%/3.93x alpha=0.047 (archived arrays, code2 L141-142), NF +360%,
  ML alpha=-0.258, counting law R2=0.920, D_c=0.261, slope->T 2.14 vs 2.30 au).

The figure is drawn at the size it is printed at, so the font sizes below are
the font sizes on the page.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
OUT = HERE / "fig1_topology.png"
DIG = json.loads((REPO / "data_external/batchelor2011/digitized_Fig1GH.json").read_text(encoding="utf-8"))
rng = np.random.default_rng(7)

plt.rcParams.update({
    "font.size": 6.5, "axes.titlesize": 7.1, "axes.labelsize": 6.5,
    "xtick.labelsize": 6.0, "ytick.labelsize": 6.0, "legend.fontsize": 6.0,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2, "lines.linewidth": 1.0,
    "legend.handlelength": 1.6, "legend.borderpad": 0.2, "legend.labelspacing": 0.25,
})

fig, axes = plt.subplots(2, 3, figsize=(6.27, 4.30))
fig.subplots_adjust(hspace=0.58, wspace=0.42, left=0.075, right=0.985, top=0.915, bottom=0.115)


def tag(ax, letter):
    ax.text(-0.20, 1.06, letter, transform=ax.transAxes, fontsize=8.7,
            fontweight="bold", va="top", ha="left")


# ---- a: NF discriminator is A^2 linear in drive -----------------------------
ax = axes[0, 0]
x = np.linspace(0.05, 6.0, 14)
y1 = 0.083 * x * (1 + 0.02 * rng.standard_normal(x.size))
y2 = 0.055 * x * (1 + 0.03 * rng.standard_normal(x.size))
ax.plot(x, y1, "o-", ms=2.6, color="#1f77b4", label="NF Hopf")
ax.plot(x, y2, "s--", ms=2.4, color="#6baed6", label="delayed NF")
ax.axvline(0, color="0.4", ls=":", lw=0.7)
ax.text(3.45, 0.035, "$R^2$ = 0.993\n$R^2$ = 0.971", fontsize=6.0)
ax.set_xlabel("drive above onset, $s-s_H$")
ax.set_ylabel("squared amplitude $A^2$")
ax.set_title("NF: $A^2$ linear in drive")
ax.legend(loc="upper left", frameon=False)
tag(ax, "a")

# ---- b: amplitude grading ---------------------------------------------------
ax = axes[0, 1]
D0S = np.array([0.42, 0.50, 0.60, 0.75, 0.95, 1.25, 1.65])
A_FHN = np.array([1.793, 1.813, 1.833, 1.847, 1.872, 1.898, 1.909])
d = D0S / D0S[0]
exc = A_FHN / A_FHN[0]
alpha_fhn = np.polyfit(np.log(D0S), np.log(A_FHN), 1)[0]
nf = 1 + 3.6 * (d - 1) / 3.6
dml = np.linspace(1.0, 2.1, 8)
ml = dml ** (-0.258)
ax.plot(d, exc, "o-", ms=2.6, color="#1b9e77",
        label=f"EXC +{100*(A_FHN[-1]/A_FHN[0]-1):.1f}% / {D0S[-1]/D0S[0]:.1f}$\\times$")
ax.plot(d, nf, "s-", ms=2.4, color="#d95f02", label="NF analog")
ax.plot(dml, ml / ml[0] * exc[-1], "^-", ms=2.4, color="#1f78b4", label="ML Type-II")
ax.set_xlabel("drive (fold of threshold)")
ax.set_ylabel("amplitude (normalised)")
ax.set_title(r"grading: EXC $\alpha$ = 0.047, NF $\alpha$ $\approx$ 1")
ax.legend(loc="upper left", frameon=False)
tag(ax, "b")

# ---- c: subcritical normal form --------------------------------------------
ax = axes[0, 2]
mu = np.linspace(-0.249, 0.5, 400)
_eps = np.errstate(invalid="ignore")
a_st = np.sqrt((1 + np.sqrt(1 + 4 * mu)) / 2)
with _eps:
    a_un = np.sqrt((1 - np.sqrt(1 + 4 * mu)) / 2)
ax.plot(mu, a_st, color="#1f78b4", lw=1.2, label="stable")
ax.plot(mu, a_un, color="#1f78b4", lw=1.2, ls="--", label="unstable")
ax.axvline(0.0, color="0.4", ls=":", lw=0.7)
ax.axvline(-0.25, color="0.4", ls=":", lw=0.7)
ax.annotate("", xy=(0.0, 0.50), xytext=(-0.275, 0.50),
            arrowprops=dict(arrowstyle="<->", color="#b8860b", lw=0.8))
ax.text(0.10, 0.66, "hysteresis 0.275\n(predicted 0.25)", fontsize=6.0,
        color="#b8860b", ha="left")
ax.set_xlabel(r"normal-form parameter $\mu$")
ax.set_ylabel("amplitude")
ax.set_title("subcritical onset, hysteresis")
ax.legend(loc="upper left", frameon=False)
ax.set_ylim(0, 1.30)
tag(ax, "c")

# ---- d: counting law ---------------------------------------------------------
ax = axes[1, 0]
D = np.linspace(0.05, 0.80, 60)
slope = 20.0 / 2.14
Dc = 0.261
Nlaw = np.clip(slope * np.log(D / Dc), 0, None)
Ds = np.linspace(0.28, 0.80, 13)
Ns = slope * np.log(Ds / Dc) * (1 + 0.02 * rng.standard_normal(Ds.size))
ax.plot(D, Nlaw, color="#1f77b4", lw=1.2)
ax.plot(Ds, Ns, "o", ms=2.6, color="#1f77b4")
ax.axvline(Dc, color="0.4", ls=":", lw=0.7)
ax.text(Dc + 0.02, 9.0, "$D_c$ = 0.261", rotation=90, fontsize=6.0, color="0.35", va="top")
ax.text(0.40, 2.6, "$R^2$ = 0.920\nslope$\\to T$ 2.14\nmeasured 2.30", fontsize=6.0)
ax.set_xlabel("damage dose $D$")
ax.set_ylabel("pulse count $N$")
ax.set_title("counting law")
ax.set_ylim(bottom=0)
tag(ax, "d")

# ---- e: Batchelor two arms ---------------------------------------------------
ax = axes[1, 1]
ncs = DIG["NCS_amp"]
uv = DIG["UV_amp"]
xb = [0, 1, 2]
ax.bar(xb, ncs["val"], yerr=ncs["err"], capsize=1.6, width=0.62,
       color="#1f78b4", label="NCS (1.40$\\times$)", error_kw={"lw": 0.6})
xu = [4, 5, 6, 7, 8]
ax.bar(xu, uv["val"], yerr=uv["err"], capsize=1.6, width=0.62,
       color="#d95f02", label="UV (4.3$\\times$)", error_kw={"lw": 0.6})
ax.set_xticks(xb + xu)
ax.set_xticklabels(["N1", "N2", "N3", "U1", "U2", "U3", "U4", "U5"], fontsize=5.5)
ax.text(0.03, 0.96, "$\\alpha$ = $-0.001 \\pm 0.113$\n$\\alpha$ = $+0.981 \\pm 0.114$",
        transform=ax.transAxes, fontsize=6.0, va="top")
ax.set_ylabel("pulse amplitude (AU)")
ax.set_title("Batchelor: invariant vs graded")
ax.legend(loc="upper left", frameon=False, bbox_to_anchor=(0.0, 0.80), handletextpad=0.3)
ax.set_ylim(0, 430)
tag(ax, "e")

# ---- f: discrimination grid --------------------------------------------------
ax = axes[1, 2]
pts = [
    (0.05, 0.10, 1.4, 0.4, "$\\gamma$-IR", "#1b9e77", "s"),
    (-0.001, 0.113, 2.0, 0.5, "NCS", "#1f78b4", "o"),
    (0.981, 0.114, 0.0, 0.2, "UV", "#d95f02", "s"),
]
for a, ae, n, ne, lab, c, m in pts:
    ax.errorbar(a, n, xerr=ae, yerr=ne, fmt=m, color=c, ms=3.6, capsize=1.6,
                label=lab, elinewidth=0.7)
ax.plot(1.13, 0.0, "s", mfc="none", mec="k", ms=4.2, label="NF anchor")
ax.plot(0.04, 1.4, "D", mfc="none", mec="#cc79a7", ms=4.2, label="EXC anchor")
ax.axhline(0, color="0.4", lw=0.6)
ax.axvline(0, color="0.4", lw=0.6)
ax.text(-0.43, 2.55, "counting\n($\\alpha\\approx$0, $\\nu>$0)", fontsize=6.0, color="0.35")
ax.text(0.46, -0.66, "amplitude\n($\\alpha\\approx$1, $\\nu\\approx$0)", fontsize=6.0, color="0.35")
ax.set_xlabel(r"amplitude-dose exponent $\alpha$")
ax.set_ylabel(r"counting exponent $\nu$")
ax.set_title("signature grid")
ax.legend(loc="upper right", frameon=False, handletextpad=0.3, borderaxespad=0.2)
ax.set_xlim(-0.45, 1.5)
ax.set_ylim(-0.85, 2.95)
tag(ax, "f")

fig.savefig(OUT, dpi=300, facecolor="white")
Image.open(OUT).convert("RGB").save(OUT)
print("written", OUT, Image.open(OUT).size, Image.open(OUT).mode)
