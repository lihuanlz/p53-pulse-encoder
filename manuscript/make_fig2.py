# -*- coding: utf-8 -*-
"""paper_p53 fig2_census.png  -  drawn at print size (6.27 in, 7 pt class).

Archives:
  code51r3_p53_Fisher_audit_results.json (M=400 central h=2%/5%: evals, evecs,
      shares, capacity, gain dispersion)
  code51r4B_Fisher_M4000.json (M=4000 central)
  code51r4A_noise_calibration_robustness.json (lit_matched dispersion series)
  FHN Eq. with locked parameters s=20, tau_r=20, a=0.7, b=0.8, eps=0.08,
  I0=0.6, sigma=0.06 (reference pulse train).
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
RES = REPO / "results"
OUT = HERE / "fig2_census.png"

r3 = json.loads((RES / "code51r3_p53_Fisher_audit_results.json").read_text(encoding="utf-8"))
r4b = json.loads((RES / "code51r4B_Fisher_M4000.json").read_text(encoding="utf-8"))
r4a = json.loads((RES / "code51r4A_noise_calibration_robustness.json").read_text(encoding="utf-8"))

plt.rcParams.update({
    "font.size": 6.5, "axes.titlesize": 7.1, "axes.labelsize": 6.5,
    "xtick.labelsize": 6.0, "ytick.labelsize": 6.0, "legend.fontsize": 6.0,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2, "lines.linewidth": 1.0,
    "legend.handlelength": 1.6, "legend.borderpad": 0.2, "legend.labelspacing": 0.25,
})

fig, axes = plt.subplots(2, 3, figsize=(6.27, 4.30))
fig.subplots_adjust(hspace=0.52, wspace=0.45, left=0.085, right=0.975, top=0.915, bottom=0.115)


def tag(ax, letter):
    ax.text(-0.22, 1.06, letter, transform=ax.transAxes, fontsize=8.7,
            fontweight="bold", va="top", ha="left")


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
    v += dv * dt
    w += dw * dt
    traj[i] = v
tt = np.arange(n) * dt
npeaks = int(np.sum((traj[1:-1] > 1.0) & (traj[1:-1] > traj[:-2]) & (traj[1:-1] >= traj[2:])))
ax.plot(tt, traj, lw=0.8, color="#1f77b4")
ax.axhline(0, color="0.6", ls=":", lw=0.7)
ax.set_xlabel("time (au)")
ax.set_ylabel("$v$ (fast variable)")
ax.set_title(f"FHN pulse train, N = {npeaks}")
tag(ax, "a")

# ---- b: Fisher spectra -------------------------------------------------------
ax = axes[0, 1]
cfgs = [("M400, h=2%", r3["fisher_h2"]["evals"], "#1f77b4"),
        ("M400, h=5%", r3["fisher_h5"]["evals"], "#6baed6"),
        ("M4000, h=2%", r4b["evals"], "#d95f02")]
x = np.arange(1, 8)
wd = 0.27
for k, (lab, ev, c) in enumerate(cfgs):
    ev = np.asarray(ev, float)
    ev[ev <= 0] = 1e-16
    ax.bar(x + (k - 1) * wd, np.log10(ev), width=wd, color=c, label=lab)
ax.set_xlabel("mode index")
ax.set_ylabel(r"$\log_{10}$ eigenvalue")
ax.set_title(r"spectrum: $\lambda_7/\lambda_1$ = 4.7$\times10^{-19}$")
ax.legend(loc="lower left", frameon=False)
tag(ax, "b")

# ---- c: eigenvector heatmap -----------------------------------------------------
ax = axes[0, 2]
pnames = [r"$\alpha$", r"$\beta$", r"$\varepsilon$", r"$s$", r"$\tau_r$", r"$I_0$", r"$\sigma$"]
E = np.abs(np.asarray(r3["fisher_h2"]["evecs"], float)).T
im = ax.imshow(E, aspect="auto", cmap="viridis", vmin=0, vmax=1)
ax.set_xticks(range(7))
ax.set_xticklabels([f"v{i+1}" for i in range(7)], fontsize=5.5)
ax.set_yticks(range(7))
ax.set_yticklabels(pnames, fontsize=6.0)
cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
cb.ax.tick_params(labelsize=5)
ax.set_title(r"eigenvectors (softest $v_7$: $\sigma$)")
tag(ax, "c")

# ---- d: channel shares -----------------------------------------------------------
ax = axes[1, 0]
obs = ["$t_1$", "$T$", "$A$", "$N$", "$w$", "IPI$_{sd}$"]
sh = np.asarray(r3["fisher_h2"]["share"], float) * 100
colors = ["#1b9e77"] * 6
colors[2] = "#d95f02"
ax.bar(range(6), sh, color=colors, width=0.72)
for i, v_ in enumerate(sh):
    ax.text(i, v_ + 1.5, f"{v_:.1f}", ha="center", fontsize=6.0)
ax.set_xticks(range(6))
ax.set_xticklabels(obs, fontsize=6.5)
ax.set_ylabel("information share (%)")
ax.set_title("shares, representative calibration")
ax.set_ylim(0, 74)
ax.text(0.02, 0.94, "timing-type (green)\namplitude (orange)", transform=ax.transAxes,
        fontsize=5.7, color="0.3", va="top")
tag(ax, "d")

# ---- e: dispersion ratio, literature-matched series --------------------------------
ax = axes[1, 1]
lm = r4a["lit_matched"]
xd = np.array([p["mu_d_over_T"] for p in lm])
yd = np.array([p["ratio_lit"] for p in lm])
ax.axhspan(0.03, 0.12, color="0.85", zorder=0)
ax.text(0.05, 0.042, "bare 0.03-0.12", fontsize=5.7, color="0.35")
ax.axhline(5.8, color="k", ls="--", lw=0.8)
ax.text(0.52, 6.4, "measured 5.8", fontsize=6.0)
ax.plot(xd, yd, "o-", ms=2.6, color="#cc79a7", label="one Exp delay\n(lit.-matched)")
aq = np.sum(xd**2 * yd) / np.sum(xd**4)
ax.plot([1.0], [aq], "o", ms=4.5, mfc="none", mec="#cc79a7", mew=1.2,
        label=f"extrapolation {aq:.1f}")
ax.axvline(0.73, color="#b8860b", ls=":", lw=0.8)
ax.text(0.75, 0.16, "$\\mu_d/T$ = 0.73", fontsize=5.7, color="#b8860b")
ax.set_yscale("log")
ax.set_xlabel(r"$\mu_d / T$")
ax.set_ylabel(r"Var($t_1$) / Var($T$)")
ax.set_title("sensing stage flips the ratio")
ax.legend(loc="upper left", frameon=False, handletextpad=0.3)
ax.set_ylim(0.02, 70)
tag(ax, "e")

# ---- f: counting law + gain-dispersion inset ---------------------------------------
ax = axes[1, 2]
cap = r3["capacity"]
ld = np.log(np.asarray(cap["doses"], float))
Nm = np.asarray(cap["N_means"], float)
Nsd = np.asarray(cap["N_sds"], float)
ax.errorbar(ld, Nm, yerr=Nsd, fmt="o-", ms=2.6, color="#1f77b4",
            capsize=1.6, elinewidth=0.7, label="N($I_0$), MC")
sl, ic = np.polyfit(ld, Nm, 1)
ax.plot(ld, sl * ld + ic, "--", lw=0.8, color="0.35", label=f"slope {sl:.1f}")
ax.set_xlabel(r"$\ln I_0$")
ax.set_ylabel("pulse count $N$")
ax.set_title("counting law: 3.07/3.17 bits")
ax.set_ylim(0, 19)
ax.legend(loc="upper left", frameon=False, handletextpad=0.3)
inset = ax.inset_axes([0.55, 0.09, 0.42, 0.30])
g = r3["gain"]
sk = [p["sigma_k"] for p in g]
cva = [p["cv_A"] for p in g]
cvt = [p["cv_T"] for p in g]
inset.plot(sk, cva, "s-", ms=1.8, lw=0.7, color="#d95f02", label="CV(A)")
inset.plot(sk, cvt, "o-", ms=1.8, lw=0.7, color="#1b9e77", label="CV(T)")
inset.set_xlabel(r"$\sigma_k$", fontsize=5.5, labelpad=1)
inset.tick_params(labelsize=4.5, length=1.5, width=0.4)
inset.legend(fontsize=4.9, frameon=False, loc="upper left", handlelength=1.0)
tag(ax, "f")

fig.savefig(OUT, dpi=300, facecolor="white")
Image.open(OUT).convert("RGB").save(OUT)
print("written", OUT, Image.open(OUT).size, Image.open(OUT).mode)
