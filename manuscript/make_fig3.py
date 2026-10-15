# -*- coding: utf-8 -*-
"""Manuscript Fig. 3 (paper_p53): emergent loop constants decomposed into
molecular constants. Figure layer of code89 (molecular decomposition);
numbers read from results/molecular_layer/code89_molecular_decomposition_
results.json so the figure cannot drift from the audited pipeline output.
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def find_json(name):
    """Search the script directory and its parents for results/molecular_layer."""
    d = HERE
    for _ in range(6):
        cand = os.path.join(d, "results", "molecular_layer", name)
        if os.path.exists(cand):
            return cand
        d = os.path.dirname(d)
    raise SystemExit(
        "cannot locate results/molecular_layer/%s; run from inside the "
        "p53-pulse-encoder repository" % name)


with open(find_json("code89_molecular_decomposition_results.json"),
          encoding="utf-8") as f:
    res = json.load(f)

plt.rcParams.update({"svg.fonttype": "none", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})

K3_pred = res["T1_K3_mdm2_lag"]["predicted_h"]
K4_pred = res["T2_K4_wip1_delay"]["predicted_h"]
S = res["T3_K1_period"]["S"]
T_pred = res["T3_K1_period"]["T_pred_h"]

fig, axes = plt.subplots(2, 2, figsize=(10, 7.2))

ax = axes[0, 0]
comp = res["T1_K3_mdm2_lag"]["components_h"]
labels = ["transcription\n(0.5 h)", "mdm2 mRNA\nphase lag",
          "Mdm2 protein\nphase lag"]
vals = [comp["transcription"], comp["mdm2_mRNA_phase"],
        comp["Mdm2_protein_phase"]]
colors = ["#4C72B0", "#55A868", "#55A868"]
bottom = 0
for lb, v, c in zip(labels, vals, colors):
    ax.bar([0], [v], bottom=bottom, color=c, width=0.45, label=lb)
    ax.text(0, bottom + v / 2, f"{v:.2f} h", ha="center", va="center",
            fontsize=8, color="white")
    bottom += v
ax.axhspan(1.5, 2.5, color="#C44E52", alpha=0.18)
ax.axhline(2.0, color="#C44E52", lw=1.2)
ax.text(0.28, 2.02, "measured 2.0 ± 0.5 h\n(Geva-Zatorsky 2006)",
        fontsize=8, color="#C44E52")
ax.set_xlim(-0.5, 0.6)
ax.set_ylim(0, 3.0)
ax.set_xticks([])
ax.set_ylabel("lag (h)")
ax.set_title("a  K3 decomposition: p53-to-Mdm2 lag\n"
             f"sum = {K3_pred:.2f} h  vs  measured 2.0 h  [closes]")
ax.legend(fontsize=7, loc="upper left", frameon=False)

ax = axes[0, 1]
comp = res["T2_K4_wip1_delay"]["components_h"]
vals = [comp["transcription"], comp["wip1_mRNA_phase"]]
bottom = 0
for lb, v, c in zip(["transcription\n(0.5 h)", "wip1 mRNA\nphase lag"], vals,
                    ["#4C72B0", "#8172B3"]):
    ax.bar([0], [v], bottom=bottom, color=c, width=0.45, label=lb)
    ax.text(0, bottom + v / 2, f"{v:.2f} h", ha="center", va="center",
            fontsize=8, color="white")
    bottom += v
ax.axhline(1.25, color="#C44E52", lw=1.2)
ax.axhspan(1.0, 1.5, color="#C44E52", alpha=0.18)
ax.text(0.28, 1.27, "Batchelor τi = 1.25 h\n(Wang 2019: 1.25 h)",
        fontsize=8, color="#C44E52")
ax.set_xlim(-0.5, 0.6)
ax.set_ylim(0, 2.2)
ax.set_xticks([])
ax.set_ylabel("delay (h)")
ax.set_title("b  K4 decomposition: p53-to-Wip1 delay\n"
             f"sum = {K4_pred:.2f} h  vs  target 1.25 h  [closes]")
ax.legend(fontsize=7, loc="upper left", frameon=False)

ax = axes[1, 0]
sh = res["T3_K1_period"]["per_species_phase_share"]
order_s = sorted(sh.items(), key=lambda kv: -kv[1])
names = [k for k, _ in order_s]
vals = [v for _, v in order_s]
ax.barh(range(len(vals)), vals, color="#4C72B0")
ax.set_yticks(range(len(vals)))
ax.set_yticklabels(names, fontsize=8)
ax.invert_yaxis()
for i, v in enumerate(vals):
    ax.text(v + 0.01, i, f"{100*v:.0f}%", va="center", fontsize=8)
ax.set_xlabel("share of loop phase sum")
ax.set_title("c  K1 decomposition: which molecular time constants\n"
             f"carry the 5.5 h period (Jacobian at S={S:.3f}, "
             f"T_pred = {T_pred:.2f} h vs 5.48 h)")

ax = axes[1, 1]
rows = [
    ("TW  Wip1 translation", "constant (universal 4 aa/s)", 1.0),
    ("dW  Wip1 degradation", "spread 0.7-2.3 /h -> distribution", 0.5),
    ("Tw  wip1 transcription rate", "not measured -> gap", 0.0),
    ("P   Wip1 dephosphorylation of ATM*", "in-vitro only -> gap", 0.0),
]
y = np.arange(len(rows))[::-1]
status_color = {1.0: "#55A868", 0.5: "#DD8452", 0.0: "#C44E52"}
for yi, (nm, st, sc) in zip(y, rows):
    ax.barh(yi, 1.0, color="#EEEEEE", height=0.62)
    ax.barh(yi, sc, color=status_color[sc], height=0.62)
    ax.text(0.02, yi, nm, fontsize=8, va="center", fontweight="bold")
    ax.text(0.02, yi - 0.31, st, fontsize=7, va="center", color="#333333")
ax.set_xlim(0, 1)
ax.set_ylim(-0.8, len(rows) - 0.2)
ax.set_xticks([])
ax.set_yticks([])
ax.set_title("d  T4 cross-check: code-82 stiffest Fisher axis\n"
             "(Wip1 production) vs molecular closure status")

fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig3_molecular.png"), dpi=300,
            bbox_inches="tight")
print("wrote fig3_molecular.png")
