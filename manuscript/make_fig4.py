# -*- coding: utf-8 -*-
"""Manuscript Fig. 4 (paper_p53): Belgacem closed-form Hopf period
cross-validated against our molecular constants, drawn at print size.

The closed form is recomputed inline (cheap); the DDE integration points are
read from results/molecular_layer/code92_closedform_crossvalidation_results.json
(240 h integrations are not re-run for a figure).
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "fig4_closedform.png")


def find_json(name):
    d = HERE
    for _ in range(6):
        cand = os.path.join(d, "results", "molecular_layer", name)
        if os.path.exists(cand):
            return cand
        d = os.path.dirname(d)
    raise SystemExit("cannot locate results/molecular_layer/%s" % name)


with open(find_json("code92_closedform_crossvalidation_results.json"), encoding="utf-8") as f:
    res = json.load(f)

plt.rcParams.update({
    "font.size": 6, "axes.titlesize": 6.5, "axes.labelsize": 6,
    "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.2,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2, "lines.linewidth": 1.0,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.handlelength": 1.5, "legend.borderpad": 0.2, "legend.labelspacing": 0.25,
})

LN2 = np.log(2.0)


def omega_c(g1, g2, AB):
    p = (-(g1**2 + g2**2) + np.sqrt((g1**2 - g2**2)**2 + 4 * AB**2)) / 2
    return np.sqrt(p) if p > 0 else np.nan


def tau_c(g1, g2, AB):
    w = omega_c(g1, g2, AB)
    return np.arctan2(w * (g1 + g2), w**2 - g1 * g2) / w


def T_of(g1, g2, AB):
    return 2 * np.pi / omega_c(g1, g2, AB)


def AB_from_tauc(g1, g2, tau_target):
    f = lambda AB: tau_c(g1, g2, AB) - tau_target
    return brentq(f, g1 * g2 * (1 + 1e-6), 1e6, xtol=1e-14, rtol=1e-14)


def T_from_measured_tau(g1, g2, tau_measured):
    return float(T_of(g1, g2, AB_from_tauc(g1, g2, tau_measured)))


dde_rows = res["P2_dde_integration"]["rows"]
T1 = T_from_measured_tau(LN2, LN2, 1.0)
T2 = T_from_measured_tau(LN2, LN2, 2.0)
T3 = T_from_measured_tau(LN2, LN2, 1.70)
T4 = T_from_measured_tau(LN2 / 0.25, 2.0, 2.0)

fig, axes = plt.subplots(1, 2, figsize=(6.27, 2.75))
fig.subplots_adjust(wspace=0.34, left=0.085, right=0.975, top=0.86, bottom=0.19)


def tag(ax, letter):
    ax.text(-0.16, 1.04, letter, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="top", ha="left")


# ---- a: closed form vs DDE ---------------------------------------------------
ax = axes[0]
taus = np.linspace(1.001, 3.0, 200)
Ts = [T_of(LN2, LN2, AB_from_tauc(LN2, LN2, tv)) for tv in taus]
ax.plot(taus, Ts, color="#4C72B0", lw=1.2, label="closed form")
xs = [r["tau_h"] for r in dde_rows if r["period_sim_h"]]
ys = [r["period_sim_h"] for r in dde_rows if r["period_sim_h"]]
ax.plot(xs, ys, "o", ms=3.2, color="#C44E52", label="our DDE")
ax.axhline(5.5, color="#333333", lw=0.7, ls="--", label="measured 5.5 h")
ax.axvline(2.0, color="#DD8452", lw=0.7, ls=":")
ax.text(2.06, 9.6, "K3 lag 2.0 h", fontsize=5.2, color="#B5641E")
ax.set_xlabel(r"loop delay $\tau$ (h)")
ax.set_ylabel("period (h)")
ax.set_ylim(5.2, 11.2)
ax.set_title("closed form vs our DDE")
ax.legend(loc="upper left", frameon=False, handletextpad=0.3)
tag(ax, "a")

# ---- b: four constant sets ----------------------------------------------------
ax = axes[1]
cases = ["Belgacem\n$\\tau$=1.0", "K3 meas.\n$\\tau$=2.0",
         "K3 decomp.\n$\\tau$=1.70", "stress\nregime"]
Ts4 = [T1, T2, T3, T4]
colors = ["#888888", "#4C72B0", "#55A868", "#DD8452"]
ax.bar(range(4), Ts4, color=colors, width=0.58)
for i, v in enumerate(Ts4):
    ax.text(i, v + 0.15, f"{v:.2f}", ha="center", fontsize=5.5)
ax.axhspan(5.2, 5.8, color="#C44E52", alpha=0.15)
ax.axhline(5.5, color="#C44E52", lw=0.7, ls="--", label="measured 5.5 h")
ax.set_xticks(range(4))
ax.set_xticklabels(cases, fontsize=5.2)
ax.set_ylabel("period (h)")
ax.set_ylim(0, 10.2)
ax.set_title("four constant sets")
ax.legend(loc="upper right", frameon=False, handletextpad=0.3)
tag(ax, "b")

fig.savefig(OUT, dpi=300, facecolor="white")
Image.open(OUT).convert("RGB").save(OUT)
print("wrote", OUT, Image.open(OUT).size, "T =", [round(t, 3) for t in Ts4])
