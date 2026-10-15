# -*- coding: utf-8 -*-
"""Manuscript Fig. 4 (paper_p53): Belgacem closed-form Hopf period
cross-validated against our molecular constants. Figure layer of code92.
The closed form is recomputed inline (cheap); the DDE integration points
are read from results/molecular_layer/code92_closedform_crossvalidation_
results.json (240 h integrations are not re-run for a figure).
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))


def find_json(name):
    d = HERE
    for _ in range(6):
        cand = os.path.join(d, "results", "molecular_layer", name)
        if os.path.exists(cand):
            return cand
        d = os.path.dirname(d)
    raise SystemExit(
        "cannot locate results/molecular_layer/%s; run from inside the "
        "p53-pulse-encoder repository" % name)


with open(find_json("code92_closedform_crossvalidation_results.json"),
          encoding="utf-8") as f:
    res = json.load(f)

plt.rcParams.update({"svg.fonttype": "none", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})

LN2 = np.log(2.0)


# Closed-form formulas (Belgacem & Ismail, arXiv:2605.23722, Section 4)
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
    AB = AB_from_tauc(g1, g2, tau_measured)
    return float(T_of(g1, g2, AB))


dde_rows = res["P2_dde_integration"]["rows"]
T1 = T_from_measured_tau(LN2, LN2, 1.0)          # Belgacem calibration
T2 = T_from_measured_tau(LN2, LN2, 2.0)          # K3 measured lag
T3 = T_from_measured_tau(LN2, LN2, 1.70)         # K3 decomposed (code89)
T4 = T_from_measured_tau(LN2 / 0.25, 2.0, 2.0)   # stress regime

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))

ax = axes[0]
taus = np.linspace(1.001, 3.0, 200)
Ts = [T_of(LN2, LN2, AB_from_tauc(LN2, LN2, tv)) for tv in taus]
ax.plot(taus, Ts, color="#4C72B0", lw=1.8,
        label="Belgacem closed form T(tau) at onset (g=ln2)")
xs = [r["tau_h"] for r in dde_rows if r["period_sim_h"]]
ys = [r["period_sim_h"] for r in dde_rows if r["period_sim_h"]]
ax.plot(xs, ys, "o", color="#C44E52", ms=7,
        label="our Python DDE integration (his model)")
ax.axhline(5.5, color="#333333", lw=1.0, ls="--")
ax.text(1.02, 5.62, "measured 5.5 h (MCF7)  |  our Monke simulation 5.48 h",
        fontsize=8, color="#333333")
ax.axhline(5.48, color="#55A868", lw=0.9, ls=":")
ax.axvline(2.0, color="#DD8452", lw=0.9, ls=":")
ax.text(2.03, 10.0, "K3 measured lag 2.0 h", fontsize=8, color="#B5641E")
ax.set_xlabel("loop delay tau (h)")
ax.set_ylabel("period (h)")
ax.set_title("a  Closed-form period vs loop delay, with our independent\n"
             "DDE integration overlaid (formula verified in Python)",
             fontsize=10)
ax.legend(fontsize=8, frameon=False, loc="upper left")

ax = axes[1]
cases = ["Belgacem\ncalib\n(tau=1.0)", "K3 measured\n(tau=2.0)",
         "K3 decomposed\n(tau=1.70)", "stress regime\n(fast p53)"]
Ts4 = [T1, T2, T3, T4]
colors = ["#888888", "#4C72B0", "#55A868", "#DD8452"]
ax.bar(range(4), Ts4, color=colors, width=0.55)
for i, v in enumerate(Ts4):
    ax.text(i, v + 0.12, f"{v:.2f} h", ha="center", fontsize=9)
ax.axhspan(5.2, 5.8, color="#C44E52", alpha=0.15)
ax.axhline(5.5, color="#C44E52", lw=1.0, ls="--")
ax.text(-0.45, 5.62, "measured 5.5 h", fontsize=8, color="#C44E52",
        ha="left")
ax.set_xticks(range(4))
ax.set_xticklabels(cases, fontsize=8)
ax.set_ylabel("closed-form period (h)")
ax.set_ylim(0, max(Ts4) * 1.25)
ax.set_title("b  Same closed form, four constant sets:\n"
             "which molecular inputs land on 5.5 h?", fontsize=10)

fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig4_closedform.png"), dpi=300,
            bbox_inches="tight")
print("wrote fig4_closedform.png  T =", [round(t, 3) for t in Ts4])
