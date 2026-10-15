# -*- coding: utf-8 -*-
"""paper_p53 fig2_census.png — rebuild: panel d title 94.5%, panel e redrawn
on the literature-matched dispersion series, inset cleanup.
Archives:
  代码51r3_p53_Fisher审计_结果.json (M=400 central h=2%/5%: evals, evecs, shares,
      capacity, gain dispersion)
  代码51r4B_Fisher_M4000.json (M=4000 central)
  代码51r4A_噪声标定与稳健性.json (lit_matched dispersion series)
  FHN Eq. with locked parameters s=20, tau_r=20, a=0.7, b=0.8, eps=0.08,
  I0=0.6, sigma=0.06 (reference pulse train).
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WS = Path(r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911")
RES = WS / "03_细胞线3/结果"
OUT = WS / "05_主线纲领与设计/论文_细胞通讯审计_2026-09-23/拆稿_LaTeX/paper_p53/fig2_census.png"

r3 = json.loads((RES / "代码51r3_p53_Fisher审计_结果.json").read_text(encoding="utf-8"))
r4b = json.loads((RES / "代码51r4B_Fisher_M4000.json").read_text(encoding="utf-8"))
r4a = json.loads((RES / "代码51r4A_噪声标定与稳健性.json").read_text(encoding="utf-8"))

fig, axes = plt.subplots(2, 3, figsize=(14.85, 8.19))
fig.subplots_adjust(hspace=0.55, wspace=0.32, left=0.06, right=0.99, top=0.90, bottom=0.10)

# ---- a: reference pulse train ------------------------------------------------
ax = axes[0, 0]
s, tau_r, a_, b_, eps, I0, sig = 20.0, 20.0, 0.7, 0.8, 0.08, 0.6, 0.06
dt, T = 0.02, 50.0
rng = np.random.default_rng(3)
n = int(T / dt)
v, w = -1.0, -0.6
traj = np.empty(n)
for i in range(n):
    t = i * dt
    I = I0 * np.exp(-t / tau_r)
    dv = s * (v - v**3 / 3 - w) + sig / np.sqrt(dt) * rng.standard_normal()
    dw = s * eps * (v + a_ - b_ * w - I)
    v += dv * dt; w += dw * dt
    traj[i] = v
tt = np.arange(n) * dt
npeaks = int(np.sum((traj[1:-1] > 1.0) & (traj[1:-1] > traj[:-2]) & (traj[1:-1] >= traj[2:])))
ax.plot(tt, traj, lw=0.9, color="#1f77b4")
ax.axhline(0, color="0.6", ls=":", lw=0.8)
ax.set_xlabel("time (au)"); ax.set_ylabel("v (FHN fast variable)")
ax.set_title(f"FHN pulse train, synchronous integration, N = {npeaks}", fontsize=10)
ax.text(-0.14, 1.02, "a", transform=ax.transAxes, fontsize=13, fontweight="bold")

# ---- b: Fisher spectra ---------------------------------------------------------
ax = axes[0, 1]
cfgs = [("M=400 central, h=2%", r3["fisher_h2"]["evals"], "#1f77b4"),
        ("M=400 central, h=5%", r3["fisher_h5"]["evals"], "#6baed6"),
        ("M=4000 central, h=2%", r4b["evals"], "#d95f02")]
x = np.arange(1, 8)
wd = 0.27
for k, (lab, ev, c) in enumerate(cfgs):
    ev = np.asarray(ev, float)
    ev[ev <= 0] = 1e-16
    ax.bar(x + (k - 1) * wd, np.log10(ev), width=wd, color=c, label=lab)
ax.set_xlabel("mode index"); ax.set_ylabel("log$_{10}$ eigenvalue")
ax.set_title("Fisher spectrum: $\\lambda_7/\\lambda_1 \\leq 4.7\\times 10^{-19}$", fontsize=10)
ax.legend(fontsize=7.5, loc="lower left", frameon=False)
ax.text(-0.14, 1.02, "b", transform=ax.transAxes, fontsize=13, fontweight="bold")

# ---- c: eigenvector heatmap -----------------------------------------------------
ax = axes[0, 2]
pnames = ["a", "b", "eps", "s", "tau_r", "I0", "sigma"]
E = np.abs(np.asarray(r3["fisher_h2"]["evecs"], float)).T  # 7 params x 7 modes
im = ax.imshow(E, aspect="auto", cmap="viridis", vmin=0, vmax=1)
ax.set_xticks(range(7)); ax.set_xticklabels([f"v{i+1}" for i in range(7)], fontsize=8)
ax.set_yticks(range(7)); ax.set_yticklabels(pnames, fontsize=8)
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
ax.set_title("eigenvectors (softest v7: sigma 0.996)", fontsize=10)
ax.text(-0.14, 1.02, "c", transform=ax.transAxes, fontsize=13, fontweight="bold")

# ---- d: channel shares -----------------------------------------------------------
ax = axes[1, 0]
obs = ["t1", "T", "A", "N", "w", "IPI_sd"]
sh = np.asarray(r3["fisher_h2"]["share"], float) * 100
colors = ["#1b9e77"] * 6; colors[2] = "#d95f02"
ax.bar(range(6), sh, color=colors)
for i, v_ in enumerate(sh):
    ax.text(i, v_ + 1.2, f"{v_:.1f}", ha="center", fontsize=8.5)
ax.set_xticks(range(6)); ax.set_xticklabels(obs, fontsize=9)
ax.set_ylabel("Fisher information share (%)")
timing = 100 - sh[2]
ax.set_title(f"channel shares: timing-type {timing:.1f}%, amplitude {sh[2]:.1f}%", fontsize=10)
ax.set_ylim(0, 72)
ax.text(-0.14, 1.02, "d", transform=ax.transAxes, fontsize=13, fontweight="bold")

# ---- e: dispersion ratio, literature-matched series --------------------------------
ax = axes[1, 1]
lm = r4a["lit_matched"]
xd = np.array([p["mu_d_over_T"] for p in lm])
yd = np.array([p["ratio_lit"] for p in lm])
ax.axhspan(0.03, 0.12, color="0.85", zorder=0)
ax.text(0.04, 0.045, "bare oscillator 0.03-0.12", fontsize=8, color="0.35")
ax.axhline(5.8, color="k", ls="--", lw=1.2)
ax.text(0.30, 7.3, "measured 5.8", fontsize=8.5)
ax.plot(xd, yd, "o-", color="#cc79a7", label="one upstream Exp delay\n(literature-matched)")
aq = np.sum(xd**2 * yd) / np.sum(xd**4)  # quadratic extrapolation a*x^2
ax.plot([1.0], [aq], "o", ms=10, mfc="none", mec="#cc79a7", mew=2,
        label=f"quadratic extrapolation: {aq:.1f} at $\\mu_d$ = T")
ax.axvline(0.73, color="#b8860b", ls=":", lw=1.2)
ax.text(0.75, 0.3, "literature calibration\n$\\mu_d/T$ = 0.73 (240/330 min)", fontsize=8, color="#b8860b")
ax.set_yscale("log")
ax.set_xlabel("$\\mu_d$ / T"); ax.set_ylabel("Var($t_1$) / Var(T)")
ax.set_title("stochastic sensing stage flips the dispersion ratio", fontsize=10)
ax.legend(fontsize=7.5, loc="upper left", frameon=True, framealpha=1.0, edgecolor="none")
ax.set_ylim(0.02, 60)
ax.text(-0.14, 1.02, "e", transform=ax.transAxes, fontsize=13, fontweight="bold")

# ---- f: counting law + gain-dispersion inset ---------------------------------------
ax = axes[1, 2]
cap = r3["capacity"]
ld = np.log(np.asarray(cap["doses"], float))
Nm = np.asarray(cap["N_means"], float); Nsd = np.asarray(cap["N_sds"], float)
ax.errorbar(ld, Nm, yerr=Nsd, fmt="o-", color="#1f77b4", capsize=3, label="N(I0), MC ensemble")
sl, ic = np.polyfit(ld, Nm, 1)
ax.plot(ld, sl * ld + ic, "--", color="0.35", label=f"log law, slope = {sl:.1f}")
ax.set_xlabel("ln $I_0$"); ax.set_ylabel("pulse count N")
ax.set_title("counting law; I(D;N) = 3.07 / 3.17 bits (96.8%)", fontsize=10)
ax.legend(fontsize=8, loc="upper left", frameon=False)
inset = ax.inset_axes([0.60, 0.10, 0.38, 0.42])
g = r3["gain"]
sk = [p["sigma_k"] for p in g]; cva = [p["cv_A"] for p in g]; cvt = [p["cv_T"] for p in g]
inset.plot(sk, cva, "s-", color="#d95f02", label="CV(A)")
inset.plot(sk, cvt, "o-", color="#1b9e77", label="CV(T)")
inset.set_xlabel("$\\sigma_k$", fontsize=8); inset.tick_params(labelsize=7)
inset.legend(fontsize=7, frameon=False, loc="upper left")
ax.text(-0.14, 1.02, "f", transform=ax.transAxes, fontsize=13, fontweight="bold")

fig.savefig(OUT, dpi=100)
print("written", OUT)
