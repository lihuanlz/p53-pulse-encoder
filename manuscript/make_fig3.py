# -*- coding: utf-8 -*-
"""Manuscript Fig. 3 (paper_p53): emergent loop constants decomposed into
molecular constants, drawn at print size (6.27 in, 7 pt class).

Figure layer of code89; numbers read from
results/molecular_layer/code89_molecular_decomposition_results.json so the
figure cannot drift from the audited pipeline output.
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "fig3_molecular.png")


def find_json(name):
    d = HERE
    for _ in range(6):
        cand = os.path.join(d, "results", "molecular_layer", name)
        if os.path.exists(cand):
            return cand
        d = os.path.dirname(d)
    raise SystemExit("cannot locate results/molecular_layer/%s" % name)


with open(find_json("code89_molecular_decomposition_results.json"), encoding="utf-8") as f:
    res = json.load(f)

plt.rcParams.update({
    "font.size": 6, "axes.titlesize": 6.5, "axes.labelsize": 6,
    "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.5,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.handlelength": 1.5, "legend.borderpad": 0.2, "legend.labelspacing": 0.25,
})

K3_pred = res["T1_K3_mdm2_lag"]["predicted_h"]
K4_pred = res["T2_K4_wip1_delay"]["predicted_h"]
S = res["T3_K1_period"]["S"]
T_pred = res["T3_K1_period"]["T_pred_h"]

fig, axes = plt.subplots(2, 2, figsize=(6.27, 4.30))
fig.subplots_adjust(hspace=0.52, wspace=0.30, left=0.085, right=0.975, top=0.905, bottom=0.095)


def tag(ax, letter):
    ax.text(-0.18, 1.05, letter, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="top", ha="left")


# ---- a: K3 decomposition -----------------------------------------------------
ax = axes[0, 0]
comp = res["T1_K3_mdm2_lag"]["components_h"]
labels = ["transcription\n0.5 h", "mdm2 mRNA\nphase", "Mdm2 protein\nphase"]
vals = [comp["transcription"], comp["mdm2_mRNA_phase"], comp["Mdm2_protein_phase"]]
colors = ["#4C72B0", "#55A868", "#2F6B44"]
bottom = 0
for lb, v, c in zip(labels, vals, colors):
    ax.bar([0], [v], bottom=bottom, color=c, width=0.42, label=lb)
    ax.text(0, bottom + v / 2, f"{v:.2f}", ha="center", va="center",
            fontsize=5.5, color="white")
    bottom += v
ax.axhspan(1.5, 2.5, color="#C44E52", alpha=0.15)
ax.axhline(2.0, color="#C44E52", lw=0.8)
ax.text(0.26, 2.05, "measured\n2.0 $\\pm$ 0.5 h", fontsize=5.2, color="#C44E52")
ax.set_xlim(-0.45, 0.62)
ax.set_ylim(0, 3.0)
ax.set_xticks([])
ax.set_ylabel("lag (h)")
ax.set_title(f"K3 (Mdm2 lag): {K3_pred:.2f} h vs 2.0")
ax.legend(loc="lower right", frameon=False, handletextpad=0.3, fontsize=5)
tag(ax, "a")

# ---- b: K4 decomposition -----------------------------------------------------
ax = axes[0, 1]
comp = res["T2_K4_wip1_delay"]["components_h"]
vals = [comp["transcription"], comp["wip1_mRNA_phase"]]
bottom = 0
for lb, v, c in zip(["transcription\n0.5 h", "wip1 mRNA\nphase"], vals,
                    ["#4C72B0", "#8172B3"]):
    ax.bar([0], [v], bottom=bottom, color=c, width=0.42, label=lb)
    ax.text(0, bottom + v / 2, f"{v:.2f}", ha="center", va="center",
            fontsize=5.5, color="white")
    bottom += v
ax.axhline(1.25, color="#C44E52", lw=0.8)
ax.axhspan(1.0, 1.5, color="#C44E52", alpha=0.15)
ax.text(0.26, 1.30, "Wip1-arm delay\ntarget 1.25 h", fontsize=5.2, color="#C44E52")
ax.set_xlim(-0.45, 0.62)
ax.set_ylim(0, 2.2)
ax.set_xticks([])
ax.set_ylabel("delay (h)")
ax.set_title(f"K4 (Wip1 arm): {K4_pred:.2f} h vs 1.25")
ax.legend(loc="lower right", frameon=False, handletextpad=0.3, fontsize=5)
tag(ax, "b")

# ---- c: period phase shares --------------------------------------------------
ax = axes[1, 0]
sh = res["T3_K1_period"]["per_species_phase_share"]
order_s = sorted(sh.items(), key=lambda kv: -kv[1])
names = [k for k, _ in order_s]
vals = [100 * v for _, v in order_s]
ax.barh(range(len(vals)), vals, color="#4C72B0", height=0.68)
ax.set_yticks(range(len(vals)))
ax.set_yticklabels(names, fontsize=6)
ax.invert_yaxis()
for i, v in enumerate(vals):
    ax.text(v + 0.8, i, f"{v:.0f}", va="center", fontsize=5.5)
ax.set_xlim(0, 36)
ax.set_xlabel("share of loop phase sum (%)")
ax.set_title(f"K1 period: Jacobian {T_pred:.2f} h vs 5.48")
tag(ax, "c")

# ---- d: closure status of the Wip1 arm ---------------------------------------
ax = axes[1, 1]
rows = [
    ("$T_W$  Wip1 translation", "universal constant, 4 aa/s", "#55A868"),
    ("$d_W$  Wip1 degradation", "0.7-2.3 /h, distribution", "#DD8452"),
    ("$T_w$  wip1 transcription", "not measured, gap", "#C44E52"),
    ("$P$  Wip1 dephosphorylation", "in vitro only, gap", "#C44E52"),
]
y = np.arange(len(rows))[::-1]
for yi, (nm, st, c) in zip(y, rows):
    ax.scatter([0.03], [yi], s=26, marker="s", color=c, clip_on=False)
    ax.text(0.09, yi + 0.12, nm, fontsize=6, va="center")
    ax.text(0.09, yi - 0.26, st, fontsize=5.2, va="center", color="#444444")
ax.set_xlim(0, 1)
ax.set_ylim(-0.7, len(rows) - 0.25)
ax.set_xticks([])
ax.set_yticks([])
for sp in ("left", "bottom"):
    ax.spines[sp].set_visible(False)
ax.set_title("Wip1-arm closure: 1 constant, 2 gaps")
tag(ax, "d")

fig.savefig(OUT, dpi=300, facecolor="white")
Image.open(OUT).convert("RGB").save(OUT)
print("wrote", OUT, Image.open(OUT).size, Image.open(OUT).mode)
