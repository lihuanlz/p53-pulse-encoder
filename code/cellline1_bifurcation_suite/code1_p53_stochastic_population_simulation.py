# -*- coding: utf-8 -*-
"""
code 1: stochastic simulation of an excitable p53 population
================================
Purpose: reproduce the four experimental signatures of Lahav 2004 / Batchelor 2008
  (1) pulse number increases with dose (digital encoding)
  (2) amplitude independent of dose (all-or-none)
  (3) inter-pulse interval (IPI) precise, ~5.5 h
  (4) first-pulse timing dispersed (phase noise-triggered)

Model: FitzHugh-Nagumo excitable reduction (Monke 2017 type)
  dv/dt = (v − v³/3 − w) × sf + noise        <- fast variable (read as p53)
  dw/dt = ε(v + a − b·w) × sf                <- slow recovery variable (read as Mdm2 etc.)
damage input: I(D,t) = (0.40+0.09D)*t/(t+0.5)*exp(-t/20)   <- rises then decays as repair proceeds

noise split in two (key point of the paper):
  - dynamical noise sigma tuned low -> preserves the beat (precise IPI)
  - extrinsic noise (per-cell gain/onset differ) -> creates phase dispersion (variable first-pulse timing)
Run: python3 code1_p53_stochastic_population_simulation.py   (approx. 1-2 min, outputs a four-panel PNG)
"""
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

rng = np.random.default_rng(7)          # fixed seed, reproducible results

# ---------- model parameters ----------
a, b, eps = 0.7, 0.8, 0.08              # FHN standard parameters (excitable regime)
sf   = 20.0                             # timescale acceleration factor -> period ~ 2 arbitrary units
sig  = 0.10                             # dynamical-system noise (kept low to preserve the beat)

T, dt = 24.0, 0.001                     # simulation window and step (arbitrary time units)
tt = np.arange(0, T, dt)
NC, doses = 200, [0.1, 0.3, 1.0, 3.0]   # 200 cells per dose, 4 doses

results = {}
for D in doses:
    # extrinsic (cell-to-cell) variation: input gain ~ log-normal, onset time ~ normal
    gain = rng.lognormal(0, 0.22, NC)
    ton  = np.clip(rng.normal(0.5, 0.18, NC), 0.1, None)
    v = np.full(NC, -1.0); w = np.full(NC, 0.0)   # starting from the resting state
    traj = np.empty((NC, len(tt)), dtype=np.float32)
    for i, t in enumerate(tt):
        traj[:, i] = v
        I = gain * (0.40 + 0.09*D) * (t/(t+ton)) * np.exp(-t/20.0)
        v += sf*(v - v**3/3 - w + I)*dt + sig*np.sqrt(dt)*rng.standard_normal(NC)
        w += sf*eps*(v + a - b*w)*dt
    # per-cell pulse detection and statistics
    stats = []
    for c in range(NC):
        pk, prop = find_peaks(traj[c], prominence=1.0, distance=int(1.0/dt), width=0.1/dt)
        tpk = tt[pk]
        stats.append(dict(n=len(pk),
                          amp=traj[c][pk].mean() if len(pk) else np.nan,
                          width=prop['widths'].mean()*dt if len(pk) else np.nan,
                          first=tpk[0] if len(pk) else np.nan,
                          ipi=np.diff(tpk) if len(pk) > 1 else np.array([])))
    results[D] = (traj, stats)

# ---------- time calibration: set median IPI = 5.5 h ----------
all_ipi = np.concatenate([s['ipi'] for _, ss in results.values() for s in ss if len(s['ipi'])])
hour = 5.5 / np.median(all_ipi)
print(f'median IPI = {np.median(all_ipi):.2f} arbitrary units -> 1 unit = {hour:.2f} h\n')
print(f'{"dose":>5} {"pulses":>11} {"ampl(CV)":>13} {"width(h)":>11} {"first(h)":>13} {"IPI(h)":>13}')
for D in doses:
    ss = results[D][1]
    ns   = np.array([s['n'] for s in ss]); amps = np.array([s['amp'] for s in ss])
    wid  = np.array([s['width'] for s in ss])*hour
    fst  = np.array([s['first'] for s in ss])*hour
    ipi  = np.concatenate([s['ipi'] for s in ss if len(s['ipi'])])*hour
    print(f'{D:5.1f} {ns.mean():4.1f}±{ns.std():<4.1f} {np.nanmean(amps):5.2f} '
          f'({np.nanstd(amps)/np.nanmean(amps)*100:3.1f}%) {np.nanmean(wid):5.1f}±{np.nanstd(wid):<4.1f} '
          f'{np.nanmean(fst):5.1f}±{np.nanstd(fst):<4.1f} {np.nanmean(ipi):5.2f}±{np.nanstd(ipi):<4.2f}')

# ---------- four-panel figure ----------
fig, axes = plt.subplots(2, 2, figsize=(13, 8))
t_h = tt*hour
ax = axes[0, 0]                                   # (a) example trajectories
for j, D in enumerate([0.3, 3.0]):
    for c in range(3):
        ax.plot(t_h, results[D][0][c]+j*3.5, lw=0.7, alpha=0.85,
                color=['#888', '#1f77b4'][j], label=f'D={D}' if c == 0 else None)
ax.set_xlabel('time (h)'); ax.set_ylabel('p53 level (offset display)'); ax.legend()
ax.set_title('(a) single-cell pulse trajectories: low vs high dose'); ax.set_xlim(0, 24)

ax = axes[0, 1]                                   # (b) count proportional to dose
bp = ax.boxplot([[s['n'] for s in results[D][1]] for D in doses],
                positions=doses, widths=0.15, patch_artist=True)
for p in bp['boxes']: p.set_facecolor('#aec7e8')
ax.set_xlabel('damage dose'); ax.set_ylabel('pulse count N'); ax.set_title('(b) digital encoding: N proportional to dose')

ax = axes[1, 0]                                   # (c) amplitude dose-independent
bp = ax.boxplot([[s['amp'] for s in results[D][1] if not np.isnan(s['amp'])] for D in doses],
                positions=doses, widths=0.15, patch_artist=True)
for p in bp['boxes']: p.set_facecolor('#ffbb78')
ax.set_xlabel('damage dose'); ax.set_ylabel('pulse amplitude'); ax.set_ylim(1.5, 2.1)
ax.set_title('(c) amplitude dose-independent (all-or-none)')

ax = axes[1, 1]                                   # (d) phase dispersion vs precise beat
ss = results[3.0][1]
fst = np.array([s['first'] for s in ss if not np.isnan(s['first'])])*hour
ipi = np.concatenate([s['ipi'] for s in ss if len(s['ipi'])])*hour
ax.hist(fst, bins=30, alpha=0.6, color='#d62728', label=f'first-pulse time (SD={fst.std():.1f} h)')
ax.hist(ipi, bins=30, alpha=0.6, color='#2ca02c', label=f'inter-pulse interval (SD={ipi.std():.2f} h)')
ax.set_xlabel('time (h)'); ax.set_ylabel('cells'); ax.legend()
ax.set_title('(d) D=3.0: phase dispersion vs precise beat')

plt.tight_layout()
plt.savefig('p53_pulse_simulation.png', dpi=200)
print('\nsaved p53_pulse_simulation.png')
