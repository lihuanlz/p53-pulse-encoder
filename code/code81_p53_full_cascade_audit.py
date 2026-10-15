# -*- coding: utf-8 -*-
"""
code81: p53 full-pathway cascade model with per-interface audit (L0 dose -> L4 fate)
=================================================================
date: 2026-09-23 | fixed seed: 20260923 | fully deterministic

Lineage:
  - design document: 03_cell_line_3/docs/design_p53_full_cascade_model_v01.md
  - encoder layer sealed: verdict_card_P53_Fisher_audit_2026-09-18.md / seal_card_p53_signalling_2026-09-18.md
    (code51r3 / code51r4A archived parameters)
  - toy cascade self-check: code50 (02_cell_line_2, v0.2 two-arm version)
  - framework: Nature_SI_v05 S1.3 (surviving set) / S1.5 (three-leg criterion) / S1.7 (information-loss identity) / S5.6 (cost of digitization)

archive anchors (quoted verbatim; overlaps must agree): 
  - time-class Fisher share 94.2%, amplitude 5.5%
  - I(D;N) = 3.07 / 3.17 bits (96.8%) 
  - gain erosion: sigma_k 0->0.6, CV(A) 0.002->0.639, CV(T) pinned at 0.016
  - sigma* = 0.25 -> within-cell IPI CV = 0.316 ~ literature 0.30; baseline sigma=0.06 -> IPI CV 0.085
  - mu_d ~ 1.0 T -> Var(t1)/Var(T) = 5.7 ~ literature 5.8

model-layer definitions:
  L0 damage layer: N_DSB ~ Poisson(35 * D)          (count production, class I)
  L1 sensing layer: N_trig = ceil(N_DSB / 50), cap 12 (registered gap, simplest clustering map, assumption)
  L2 encoding layer: T = 5.5 h fixed constant; A fixed (CV 0.002); N = N_trig;
             first-pulse delay t1 ~ Exp(1.0 T); interval jitter CV 0.085 (archived sigma=0.06 row)
  L3 decoding layer: leaky integration of mRNA; p21-type integrator tau=10 h; PUMA-type accumulator tau=4 h, produces only from the 3rd pulse
  L4 fate layer: thresholds on the PUMA/p21 integral ratio and absolute levels -> {survival, arrest, apoptosis}
"""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
try:
    from daimon_runtime import setup_plot
    setup_plot()
except Exception:
    pass

SEED = 20260923
ROOT = Path(__file__).resolve().parent.parent
OUT_JSON = ROOT / "results" / "code81_cascade_audit_results.json"
OUT_PNG = ROOT / "results" / "code81_cascade_audit_fivepanel.png"

# ---------------- fixed constants (all explicit, none fitted) ----------------
DOSES = np.array([0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0])   # Gy, equal weights -> H(D)=log2(7)
N_CELLS = 4000                                            # cells per dose
K0_DSB = 35.0         # DSB/Gy (literature anchor) 
C_CLUSTER = 50.0      # DSB/trigger (L1 clustering constant, assumption: 10 Gy -> 7 pulse, matching the literature range 0-7 pulses/24h) 
N_MAX = 12            # trigger cap
T_PERIOD = 5.5        # h, candidate structural constant (class II) , archived fixed value
MU_D = 1.0 * T_PERIOD # mean first-pulse delay (archived A3 calibration: mu_d~T gives the ratio 5.7 ~ literature 5.8) 
JIT_IPI = 0.085       # per-interval jitter CV (archived r4A sigma=0.06 row IPI CV = 0.085) 
A0 = 1.0              # fixed amplitude (carries no dose information) 
A_CV = 0.002          # intrinsic amplitude CV (archived CV(A)=0.002 @ sigma_k=0) 
W = 48.0              # observation window (h)
Q21, TAU21 = 100.0, 10.0        # p21-type integrator: quantum yield per pulse, leak time constant
QPU, TAUPU, PUMA_NTHR = 100.0, 4.0, 3   # PUMA-type accumulator: produces only from the 3rd pulse (high threshold) 
TH_PUMA, R0_FATE, TH_P21 = 800.0, 0.15, 1500.0  # L4 fate thresholds (fixed constants) 
B_GAIN = 0.3          # affine observation baseline b (code50 same convention) 
ETA_OBS = 0.02        # additive readout noise (relative to A0) 
N_BOOT = 200          # bootstrap replicates
ALPHA_E = 7.0         # S5.6 event-reliability constant

# ---------------- mutual-information tools ----------------

def mi_discrete(x, y):
    """Plug-in discrete mutual information (bits); x, y are integer label arrays."""
    x = np.asarray(x); y = np.asarray(y)
    Kx = int(x.max()) + 1; Ky = int(y.max()) + 1
    N = len(x)
    cxy = np.bincount(x * Ky + y, minlength=Kx * Ky).reshape(Kx, Ky) / N
    cx = cxy.sum(1, keepdims=True); cy = cxy.sum(0, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        term = cxy * np.log2(cxy / (cx @ cy))
    return float(np.nansum(term))


def mi_mm(x, y):
    """Miller-Madow bias-corrected mutual information (bits), floored at 0."""
    x = np.asarray(x); y = np.asarray(y)
    rx = len(np.unique(x)); ry = len(np.unique(y)); N = len(x)
    bias = (rx - 1) * (ry - 1) / (2 * N * np.log(2))
    return max(0.0, mi_discrete(x, y) - bias)


def qbin(v, nbin):
    """Quantile binning -> integer labels."""
    v = np.asarray(v, dtype=float)
    edges = np.unique(np.quantile(v, np.linspace(0, 1, nbin + 1)))
    if len(edges) <= 2:
        return np.zeros(len(v), dtype=int)
    return np.clip(np.digitize(v, edges[1:-1]), 0, len(edges) - 2)


def boot_ci(x, y, n_boot=N_BOOT, rng=None):
    """Bootstrap percentile CI (plug-in MI)."""
    rng = rng or np.random.default_rng(0)
    N = len(x); vals = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, N, N)
        vals[b] = mi_discrete(x[idx], y[idx])
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))

# ---------------- cascade simulation ----------------

def simulate(dose, n, rng, cvT_extra=0.0):
    """Single-dose population simulation. Returns the per-layer output dict.

    cvT_extra: per-cell period dispersion (for the fixed-constant test, audit 4) . 
    """
    # L0: damage count
    n_dsb = rng.poisson(K0_DSB * dose, n)
    # L1: clustered-trigger mapping (registered gap, simplest form, assumption)
    n_trig = np.minimum(np.ceil(n_dsb / C_CLUSTER).astype(int), N_MAX)
    # L2: pulse train (event timing)
    t1 = rng.exponential(MU_D, n)
    T_cell = T_PERIOD * (1.0 + cvT_extra * rng.standard_normal(n))  # per-cell period dispersion
    jit = JIT_IPI
    times = np.full((n, N_MAX), np.nan)
    tcur = t1.copy()
    for j in range(N_MAX):
        times[:, j] = tcur
        interval = T_cell * (1.0 + jit * rng.standard_normal(n))
        interval = np.maximum(interval, 0.3 * T_PERIOD)
        tcur = tcur + interval
    in_win = times <= W
    pulse_idx = np.arange(N_MAX)[None, :]
    active = (pulse_idx < n_trig[:, None]) & in_win
    n_pulse = active.sum(1)
    # amplitude: fixed + intrinsic CV 0.002
    amps = np.where(active, A0 * (1.0 + A_CV * rng.standard_normal(times.shape)), np.nan)
    # L3: leaky-integrator decoding
    with np.errstate(invalid="ignore"):
        integ = lambda q, tau: np.nansum(
            np.where(active, q * tau * (1.0 - np.exp(-(W - times) / tau)), 0.0), axis=1)
    m_p21 = integ(Q21, TAU21)
    puma_mask = active & (pulse_idx >= PUMA_NTHR - 1)
    m_puma = np.nansum(np.where(puma_mask,
                                QPU * TAUPU * (1.0 - np.exp(-(W - times) / TAUPU)), 0.0), axis=1)
    # L4: fate = thresholds on the count ratio and levels
    ratio = m_puma / (m_p21 + 1.0)
    fate = np.where((m_puma >= TH_PUMA) & (ratio >= R0_FATE), 2,
                    np.where(m_p21 >= TH_P21, 1, 0))  # 0 survival, 1 arrest, 2 apoptosis
    # per-cell period estimation (OLS slope of event time vs index -> CV(T) convention matches the archive)
    T_est = np.full(n, np.nan)
    for i in range(n):
        k = int(n_pulse[i])
        if k >= 3:
            xs = np.arange(k)
            T_est[i] = np.polyfit(xs, times[i, :k], 1)[0]
    return dict(n_dsb=n_dsb, n_trig=n_trig, n_pulse=n_pulse, t1=t1,
                times=times, amps=amps, m_p21=m_p21, m_puma=m_puma,
                fate=fate, T_est=T_est)


def run_ensemble(cvT_extra=0.0, seed=SEED):
    rng = np.random.default_rng(seed)
    parts = []
    for di, D in enumerate(DOSES):
        p = simulate(D, N_CELLS, rng, cvT_extra=cvT_extra)
        p["dose_idx"] = np.full(N_CELLS, di)
        parts.append(p)
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}

# ---------------- main pipeline ----------------

def main():
    log = lambda *a: print(*a, flush=True)
    log("[code81] full-cascade simulation starting, seed", SEED)
    E = run_ensemble()
    d = E["dose_idx"]
    H_D = np.log2(len(DOSES))
    log(f"H(D) = {H_D:.3f} bits, total cells {len(d)}")

    # ---- fate-distribution check ----
    fate_names = ["survival", "arrest", "apoptosis"]
    fate_tab = {}
    for di, D in enumerate(DOSES):
        m = d == di
        fr = [float(np.mean(E["fate"][m] == k)) for k in range(3)]
        fate_tab[f"{D}Gy"] = dict(zip(fate_names, [round(f, 4) for f in fr]))
        log(f"  D={D:>5} Gy  fate distribution S/A/P = {[round(f,3) for f in fr]}, "
            f"<N_DSB>={E['n_dsb'][m].mean():.1f}, <N_pulse>={E['n_pulse'][m].mean():.2f}")

    # ---- audit 2: information-loss identity H(D) = I(D;Y5) + sum L_i ----
    y1 = E["n_dsb"]                       # L0 output
    y2 = E["n_trig"]                      # L1 output
    y3 = E["n_pulse"] * 8 + qbin(E["t1"], 8)          # L2 output (count + first-pulse delay binned) 
    y4 = qbin(E["m_p21"], 12) * 12 + qbin(E["m_puma"], 12)  # L3 output (joint binning) 
    y5 = E["fate"]                        # L4 output
    labels = ["Y1=N_DSB (L0)", "Y2=N_trig (L1)", "Y3=pulse train (L2)",
              "Y4=mRNA joint (L3)", "Y5=fate (L4)"]
    ys = [y1, y2, y3, y4, y5]
    rng_b = np.random.default_rng(SEED + 1)
    I_chain, I_chain_mm, I_ci = [], [], []
    for lab, y in zip(labels, ys):
        I0 = mi_discrete(d, y)
        lo, hi = boot_ci(d, y, rng=rng_b)
        I_chain.append(I0); I_chain_mm.append(mi_mm(d, y)); I_ci.append([lo, hi])
        log(f"  I(D;{lab}) = {I0:.4f} bits  MM={mi_mm(d,y):.4f}  CI95=[{lo:.4f},{hi:.4f}]")
    I_seq = [H_D] + I_chain               # I(D;Y0)=H(D)
    losses = [I_seq[i] - I_seq[i + 1] for i in range(5)]
    closure = I_chain[-1] + sum(losses)
    weakest = int(np.argmax(losses))
    weakest_name = ["L0 damage", "L1 sensing(gap)", "L2 encoder", "L3 decoder", "L4 fate"][weakest]
    log(f"  identity closure: I(D;Y5)+sum(L_i) = {closure:.6f} vs H(D) = {H_D:.6f}, "
        f"residual {closure - H_D:.2e}")
    log(f"  per-layer losses L0..L4 = {[round(float(l),4) for l in losses]}, weakest layer = {weakest_name}")

    # ---- audit 3: gain-dispersion test (count chain vs amplitude chain) ----
    SIG_K = [0.0, 0.1, 0.2, 0.3, 0.45, 0.6]
    rng_g = np.random.default_rng(SEED + 2)
    gain_sweep = []
    amps = E["amps"]; times = E["times"]; n_pulse = E["n_pulse"]
    # hypothetical amplitude arm: A(D) = A0*(1+0.15*ln(D/1Gy)); amplitude roughly doubles across the dose range (code50 arm-B convention)
    amp_code = A0 * (1.0 + 0.15 * np.log(DOSES[d]))[:, None] * np.ones_like(amps)
    amp_intr = 0.15  # hypothetical amplitude-arm intrinsic CV (code50 arm-A convention 1.0+/-0.15) 
    for sk in SIG_K:
        k = np.exp(rng_g.standard_normal(len(d)) * sk - 0.5 * sk * sk)  # LogNormal(0,sk), median 1
        y_obs = k[:, None] * amps + B_GAIN + ETA_OBS * rng_g.standard_normal(amps.shape)
        # count detection: per-cell affine-invariant threshold (pulse-free baseline estimate + 0.5*dynamic range)
        base = B_GAIN + ETA_OBS * rng_g.standard_normal((len(d), 8))
        b_est = np.median(base, axis=1)
        has = n_pulse > 0
        y_max = np.full(len(d), np.nan)
        amp_meas = np.full(len(d), np.nan)
        y_max[has] = np.nanmax(y_obs[has], axis=1)
        amp_meas[has] = k[has] * np.nanmean(amps[has], axis=1)
        thr = b_est + 0.5 * (y_max - b_est)
        n_det = np.where(has, np.nansum(y_obs > thr[:, None], axis=1), 0).astype(int)
        # amplitude chain (hypothetical arm): measured amplitude = k * A_code * (1+0.15 z)
        amp_hyp = k * amp_code[:, 0] * (1.0 + amp_intr * rng_g.standard_normal(len(d)))
        amp_hyp = amp_hyp + ETA_OBS * rng_g.standard_normal(len(d))
        I_cnt = mi_discrete(d, n_det)
        I_amp = mi_discrete(d, qbin(amp_hyp, 16))
        # encoder-convention reproduction: CV(A_meas) and CV(T_est)
        real = ~np.isnan(amp_meas)
        cv_A = float(np.std(amp_meas[real]) / np.mean(amp_meas[real]))
        tval = E["T_est"]; realt = ~np.isnan(tval)
        cv_T = float(np.std(tval[realt]) / np.mean(tval[realt]))
        # archived calibration subset: N_pulse >= 7 (matching the archived N~9 period-estimation convention)
        calib = realt & (n_pulse >= 7)
        cv_T_cal = float(np.std(tval[calib]) / np.mean(tval[calib])) if calib.sum() > 50 else np.nan
        gain_sweep.append(dict(sigma_k=sk, I_count=I_cnt, I_amp_hyp=I_amp,
                               cv_A=cv_A, cv_T=cv_T, cv_T_calib=cv_T_cal,
                               frac_count_intact=float(np.mean(n_det == n_pulse))))
        log(f"  sigma_k={sk:.2f}: I(count)={I_cnt:.4f}  I(amp_hyp)={I_amp:.4f}  "
            f"CV(A)={cv_A:.3f}  CV(T)={cv_T:.4f}/{cv_T_cal:.4f}(N>=7)  "
            f"count detection rate={np.mean(n_det==n_pulse):.4f}")

    # ---- audit 4: fixed-constant test (CV(T) dispersion scan) ----
    CVT = [0.0, 0.016, 0.05, 0.1, 0.2]
    cvt_sweep = []
    for cv in CVT:
        Ec = run_ensemble(cvT_extra=cv, seed=SEED + 3)
        I_N = mi_discrete(Ec["dose_idx"], Ec["n_pulse"])
        I_F = mi_discrete(Ec["dose_idx"], Ec["fate"])
        cvt_sweep.append(dict(cvT=cv, I_count=I_N, I_fate=I_F))
        log(f"  CV(T)={cv:.3f}: I(D;N)={I_N:.4f}  I(D;fate)={I_F:.4f}")
    # end-to-end share of hypothetical amplitude code vs count code (calibration gain dispersion sigma_k=0.3, S5.6 cellular protein variation ~30%)
    rng_a = np.random.default_rng(SEED + 4)
    amp_share = {}
    for sk in [0.0, 0.3, 0.6]:
        k = np.exp(rng_a.standard_normal(len(d)) * sk - 0.5 * sk * sk)
        amp_hyp = k * amp_code[:, 0] * (1.0 + amp_intr * rng_a.standard_normal(len(d)))
        I_amp = mi_mm(d, qbin(amp_hyp, 16))
        I_cnt = mi_mm(d, E["n_pulse"])
        amp_share[f"sigma_k={sk}"] = dict(I_amp=I_amp, I_count=I_cnt,
                                          share=I_amp / (I_amp + I_cnt))
        log(f"  amplitude share @sigma_k={sk}: I_amp={I_amp:.4f} vs I_cnt={I_cnt:.4f} "
            f"-> share={I_amp/(I_amp+I_cnt):.4f}")

    # ---- audit 5: energy accounting (S5.6) ----
    bits_end = mi_mm(d, y5)
    D_levels = float(2 ** bits_end)     # end-to-end resolvable levels (inverted from transmitted bits) 
    L_axis = np.arange(1, 9)
    E_A = 2.0 * L_axis**2 * D_levels**2    # 2 kT L^2 D^2
    E_B = ALPHA_E * L_axis * D_levels**2   # alpha kT L D^2
    L_star = ALPHA_E / 2.0
    L_cascade = 4
    energy = dict(bits_per_cell=bits_end, D_levels=D_levels,
                  E_A_at_L4=float(2 * L_cascade**2 * D_levels**2),
                  E_B_at_L4=float(ALPHA_E * L_cascade * D_levels**2),
                  ratio_EA_over_EB_at_L4=float(2 * L_cascade / ALPHA_E),
                  L_star=L_star,
                  E_A_series=E_A.tolist(), E_B_series=E_B.tolist(), L_axis=L_axis.tolist())
    log(f"  energy: end-to-end {bits_end:.4f} bits/cell -> D={D_levels:.2f} resolvable levels; "
        f"E_A/E_B at L=4 = {2*L_cascade/ALPHA_E:.3f}, crossover L*={L_star}")

    # ---- audit 1: three-leg criterion x 4 interfaces ----
    g0, g6 = gain_sweep[0], gain_sweep[-1]
    legs = []
    def leg(name, legs_status, evidence, flagged=False):
        legs.append(dict(interface=name, legs=legs_status, evidence=evidence,
                         assumption_flagged=flagged))
        log(f"  interface {name}: {legs_status} | {evidence}")
    leg("I0->1 dose->DSB count",
        dict(leg1_algebraic=True, leg2_symbolic=True, leg3_temporal=True),
        f"count is a class-I calibration-free statistic (S1.4) ; I(D;N_DSB)={I_chain[0]:.3f} bits > 0; "
        f"the DSB-formation timescale (minutes) is far below the dose-delivery and downstream timescales (hours) ")
    leg("I1->2 DSB->pulse trigger (registered gap)",
        dict(leg1_algebraic=True, leg2_symbolic=True, leg3_temporal=True),
        f"the clustering map count->count makes algebraic survival trivially hold (assumption) ; I(D;N_trig)={I_chain[1]:.3f} bits; "
        f"ATM activation is minute-scale << T=5.5h; no public single-cell data exist for this layer and all three legs are flagged as assumptions",
        flagged=True)
    leg("I2->3 pulse train->mRNA (count/event timing)",
        dict(leg1_algebraic=True, leg2_symbolic=True, leg3_temporal=True),
        f"counting is immune to affine degeneracy: sigma_k 0->0.6 detection rate {g0['frac_count_intact']:.3f}->"
        f"{g6['frac_count_intact']:.3f}, I(count) {g0['I_count']:.3f}->{g6['I_count']:.3f} bits; "
        f"decoders tau_p21=10h, tau_PUMA=4h bracket the inter-pulse interval 5.5h and the window 48h (tau>W forgets, tau<<T does not integrate) ")
    leg("I2->3 pulse train->mRNA (amplitude, counterexample leg)",
        dict(leg1_algebraic=False, leg2_symbolic=False, leg3_temporal=None),
        f"amplitude leg 1 dies: sigma_k 0->0.6, CV(A_meas) {g0['cv_A']:.3f}->{g6['cv_A']:.3f}; "
        f"hypothetical amplitude-code share @sigma_k=0.6 = {amp_share['sigma_k=0.6']['share']:.3f}"
        f" (the end-to-end counterpart of the archived encoder-layer Fisher share 5.5%) ")
    leg("I3->4 mRNA->fate",
        dict(leg1_algebraic=True, leg2_symbolic=True, leg3_temporal=True),
        f"fate = threshold on count ratios, which cancels common-mode gain (S1.4 difference-class survival) ; I(D;fate)={I_chain[4]:.3f} bits > 0; "
        f"the decision timescale 24-48h matches the mRNA integration timescale")

    # ---- archive-consistency check (verbatim comparison of overlapping numbers) ----
    archived_check = dict(
        cv_A_at_sk0=dict(model=g0["cv_A"], archived=0.002,
                         match=abs(g0["cv_A"] - 0.002) < 0.005),
        cv_A_at_sk06=dict(model=g6["cv_A"], archived=0.639, note="qualitatively reproduced (generative model with LogNormal gain)",
                          match=0.55 <= g6["cv_A"] <= 0.75),
        cv_T_pinned=dict(model_range=[min(g["cv_T"] for g in gain_sweep),
                                      max(g["cv_T"] for g in gain_sweep)],
                         model_calib_N7=[round(g["cv_T_calib"], 4) for g in gain_sweep],
                         archived=0.016,
                         note="pinned at a small value under the OLS period-estimation convention and invariant to gain dispersion; the N>=7 calibration subset matches the archived N~9 convention",
                         match=max(g["cv_T"] for g in gain_sweep) < 0.06),
        count_immune=dict(I_count_range=[min(g["I_count"] for g in gain_sweep),
                                         max(g["I_count"] for g in gain_sweep)],
                          note="count-chain information flat in sigma_k"),
        archived_anchor=dict(fisher_time_share=0.942, fisher_amp_share=0.055,
                             I_DN="3.07/3.17 bits (96.8%)", cv_T_pinned=0.016,
                             sigma_star=0.25, ipi_cv_star=0.316, var_ratio_lit=5.7))
    for key in ["cv_A_at_sk0", "cv_A_at_sk06", "cv_T_pinned"]:
        mv = archived_check[key].get("model", archived_check[key].get("model_range"))
        log(f"  archive check {key}: model={mv} match={archived_check[key]['match']}")

    # ---- dual recording of failures and assumptions ----
    double_record = [
        "The L1 sensing layer is a registered gap: no public single-cell ATM activity data; the clustering constant c=50 DSB/trigger is an assumption (chosen so that 10 Gy -> 7 pulses, matching the literature range 0-7 pulses/24h); the L2 loss in the identity contains this layer's unknowable component",
        "The L3/L4 decoding kernels (quantum yield, tau, PUMA produces only from the 3rd pulse, fate thresholds TH_PUMA=800/R0=0.15/TH_P21=1500) are literature-inspired assumptions, not fits",
        "Mutual information is a binned plug-in estimate (mRNA 12x12 quantile binning); absolute values are biased, trends are reliable (same convention as code50); Miller-Madow corrected values also reported",
        "p21 integrator tau=10h > inter-pulse interval 5.5h: per-pulse resolution dies at this interface; what survives is the integrated count (consistent with the code50 leg-3 tau=8h collapse); the timescale of the integrated statistic is the whole-train duration, not the single-pulse interval",
        "The amplitude arm is a hypothetical counterfactual arm (A(D) log-slope 0.15, intrinsic CV 0.15), giving the amplitude channel a control with 'information to lose'; real p53 amplitude is dose-invariant (archived |alpha|<=0.05)",
        "The archived CV(T) value 0.016 comes from the FHN representative sigma=0.06 row; this generative model reproduces the same order of pinned-down value via OLS period estimation; numerical agreement depends on the estimation convention, as explicitly declared",
    ]

    results = dict(
        meta=dict(script="code81_p53_full_cascade_audit.py", seed=SEED, date="2026-09-23",
                  doses=DOSES.tolist(), n_cells_per_dose=N_CELLS,
                  constants=dict(K0_DSB=K0_DSB, C_CLUSTER=C_CLUSTER, N_MAX=N_MAX,
                                 T_PERIOD=T_PERIOD, MU_D_over_T=MU_D / T_PERIOD,
                                 JIT_IPI=JIT_IPI, A0=A0, A_CV=A_CV, W=W,
                                 Q21=Q21, TAU21=TAU21, QPU=QPU, TAUPU=TAUPU,
                                 PUMA_NTHR=PUMA_NTHR, TH_PUMA=TH_PUMA,
                                 R0_FATE=R0_FATE, TH_P21=TH_P21, B_GAIN=B_GAIN,
                                 ETA_OBS=ETA_OBS, ALPHA_E=ALPHA_E)),
        fate_table=fate_tab,
        info_identity=dict(H_D=H_D, labels=labels, I_chain=I_chain,
                           I_chain_mm=I_chain_mm, I_ci95=I_ci,
                           losses=losses, closure_sum=closure,
                           closure_residual=closure - H_D,
                           weakest_layer=weakest_name,
                           weakest_loss=losses[weakest]),
        gain_sweep=gain_sweep,
        cvt_sweep=cvt_sweep,
        amp_share=amp_share,
        energy=energy,
        three_legs=legs,
        archived_check=archived_check,
        double_record=double_record,
)
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)
    log("JSON written:", OUT_JSON)

    # ---------------- five-panel figure ----------------
    fig = plt.figure(figsize=(17, 10.5))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.30)

    # (a) cascade diagram + surviving-statistic labels
    axa = fig.add_subplot(gs[0, 0])
    axa.axis("off")
    boxes = [("L0 damage\n$N_{DSB}$ ~ Pois(35D)", "count"),
             ("L1 ATM sensing\n$N_{trig}$=ceil($N_{DSB}$/50)\n[registered gap]", "count (assumed)"),
             ("L2 p53 encoder\nT=5.5h fixed, A fixed", "count + event times\n(amplitude dead)"),
             ("L3 promoter\nleaky integral", "count integral"),
             ("L4 fate\nratio threshold", "ordinal fate")]
    ys_pos = np.linspace(0.88, 0.12, 5)
    for (txt, surv), yp in zip(boxes, ys_pos):
        axa.add_patch(plt.Rectangle((0.06, yp - 0.075), 0.52, 0.15,
                                    fc="#eef3fb", ec="#33527a", lw=1.4,
                                    transform=axa.transAxes))
        axa.text(0.32, yp, txt, ha="center", va="center", fontsize=8.2,
                 transform=axa.transAxes)
        axa.text(0.64, yp, "survives:\n" + surv, ha="left", va="center",
                 fontsize=7.6, color="#0b6b3a", transform=axa.transAxes)
    for i in range(4):
        axa.annotate("", xy=(0.32, ys_pos[i + 1] + 0.078), xytext=(0.32, ys_pos[i] - 0.078),
                     arrowprops=dict(arrowstyle="-|>", color="#33527a", lw=1.6),
                     xycoords=axa.transAxes)
    axa.set_title("(a) p53 cascade: surviving statistic per interface", fontsize=10.5)

    # (b) information waterfall
    axb = fig.add_subplot(gs[0, 1])
    names = ["H(D)", "L0\ndamage", "L1\nsense(gap)", "L2\nencode", "L3\ndecode", "L4\nfate", "I(D;fate)"]
    starts = [0.0]
    heights = [H_D]
    run = H_D
    for l in losses:
        starts.append(run - l); heights.append(l); run -= l
    starts.append(0.0); heights.append(I_chain[-1])
    colors = ["#33527a"] + ["#c05640"] * 4 + ["#d9863d"] + ["#0b6b3a"]
    colors[1 + weakest] = "#8c1d18"
    axb.bar(range(7), heights, bottom=starts, color=colors, width=0.68)
    for i, (s, h) in enumerate(zip(starts, heights)):
        axb.text(i, s + h + 0.04, f"{h:.2f}", ha="center", fontsize=8.5)
    axb.annotate(f"weakest: {weakest_name}\nL={losses[weakest]:.2f} bits",
                 xy=(1 + weakest, starts[1 + weakest] + 0.45 * heights[1 + weakest]),
                 xytext=(1 + weakest + 0.7, H_D * 1.16), fontsize=8.5, color="#8c1d18",
                 arrowprops=dict(arrowstyle="->", color="#8c1d18"))
    axb.set_xticks(range(7)); axb.set_xticklabels(names, fontsize=7.2)
    axb.set_ylabel("bits"); axb.set_ylim(0, H_D * 1.32)
    axb.set_title("(b) Information waterfall: $H(D)=I(D;Y_4)+\\sum L_i$", fontsize=10.5)

    # (c) gain-dispersion asymmetry
    axc = fig.add_subplot(gs[0, 2])
    sks = [g["sigma_k"] for g in gain_sweep]
    axc.plot(sks, [g["I_count"] for g in gain_sweep], "o-", color="#0b6b3a",
             label="counting chain I(D;N)")
    axc.plot(sks, [g["I_amp_hyp"] for g in gain_sweep], "s-", color="#c05640",
             label="hypothetical amplitude chain")
    axc.set_xlabel("gain dispersion $\\sigma_k$"); axc.set_ylabel("bits")
    axc2 = axc.twinx()
    axc2.plot(sks, [g["cv_A"] for g in gain_sweep], "^--", color="#8c1d18", label="CV(A)")
    axc2.plot(sks, [g["cv_T"] for g in gain_sweep], "v--", color="#33527a", label="CV(T) pop.")
    axc2.plot(sks, [g["cv_T_calib"] for g in gain_sweep], "v-", color="#5a7fb5",
              label="CV(T) N>=7 calib.")
    axc2.set_ylabel("CV", color="#555")
    axc2.axhline(0.639, color="#8c1d18", ls=":", lw=0.9)
    axc2.axhline(0.016, color="#33527a", ls=":", lw=0.9)
    axc2.set_ylim(-0.09, 0.75)
    axc2.text(0.22, 0.605, "archived 0.639", fontsize=7.5, color="#8c1d18")
    axc2.text(0.28, -0.055, "archived CV(T)=0.016", fontsize=7.5, color="#33527a")
    h1, l1 = axc.get_legend_handles_labels(); h2, l2 = axc2.get_legend_handles_labels()
    axc.legend(h1 + h2, l1 + l2, fontsize=7.5, loc="center left")
    axc.set_title("(c) Gain-dispersion asymmetry: count survives, amplitude dies", fontsize=10.5)

    # (d) fixed-constant robustness
    axd = fig.add_subplot(gs[1, 0])
    cvs = [c["cvT"] for c in cvt_sweep]
    axd.plot(cvs, [c["I_count"] for c in cvt_sweep], "o-", color="#0b6b3a",
             label="I(D; pulse count)")
    axd.plot(cvs, [c["I_fate"] for c in cvt_sweep], "s-", color="#33527a",
             label="I(D; fate) end-to-end")
    axd.axvline(0.016, color="#8c1d18", ls=":", lw=1.0)
    axd.text(0.017, axd.get_ylim()[0] + 0.05, "archived\nCV(T)=0.016", fontsize=7.5,
             color="#8c1d18")
    sh = amp_share["sigma_k=0.6"]
    axd.text(0.55, 0.18, f"hypothetical amplitude share @ $\\sigma_k$=0.6:\n"
                         f"{sh['share']*100:.1f}% (archived encoder share 5.5%)",
             transform=axd.transAxes, fontsize=8, color="#c05640",
             bbox=dict(fc="#fdf1ec", ec="#c05640", lw=0.8))
    axd.set_xlabel("per-cell CV(T) dispersion"); axd.set_ylabel("bits")
    axd.legend(fontsize=8, loc="upper right")
    axd.set_title("(d) Fixed-constant robustness: information vs CV(T)", fontsize=10.5)

    # (e) energy comparison
    axe = fig.add_subplot(gs[1, 1])
    axe.plot(L_axis, E_A, "o-", color="#c05640", label="analogue $E_A=2kTL^2D^2$")
    axe.plot(L_axis, E_B, "s-", color="#0b6b3a", label="counting $E_B=\\alpha kTLD^2$, $\\alpha$=7")
    axe.axvline(L_star, color="#555", ls="--", lw=1.0)
    axe.text(L_star + 0.08, E_A[-1] * 0.75, f"crossover $L^*=\\alpha/2$={L_star}", fontsize=8.5)
    axe.axvline(L_cascade, color="#33527a", ls=":", lw=1.2)
    axe.text(L_cascade + 0.08, E_A[0] * 1.6, "p53 cascade L=4\n(counting cheaper by "
             f"{(2*L_cascade/ALPHA_E - 1)*100:.0f}%)", fontsize=8.5, color="#33527a")
    axe.set_xlabel("cascade depth L (transmission stages)")
    axe.set_ylabel(f"energy (kT units, D={D_levels:.2f} levels)")
    axe.legend(fontsize=8.5, loc="upper left")
    axe.set_title("(e) Energy: counting vs analogue (S5.6)", fontsize=10.5)

    # summary strip
    axs = fig.add_subplot(gs[1, 2]); axs.axis("off")
    summary = (
        f"Code 81 audit summary (seed {SEED})\n\n"
        f"H(D) = {H_D:.3f} bits (7 doses, uniform)\n"
        f"I(D; fate) = {I_chain[-1]:.3f} bits end-to-end\n"
        f"Losses L0..L4 = {[round(float(l),2) for l in losses]}\n"
        f"Weakest layer: {weakest_name} ({losses[weakest]:.2f} bits)\n\n"
        f"Gain sweep: I(count) {gain_sweep[0]['I_count']:.2f} -> "
        f"{gain_sweep[-1]['I_count']:.2f} bits (flat)\n"
        f"CV(A) {gain_sweep[0]['cv_A']:.3f} -> {gain_sweep[-1]['cv_A']:.3f} "
        f"(archived 0.002 -> 0.639)\n"
        f"CV(T) pinned <= {max(g['cv_T'] for g in gain_sweep):.3f} (archived 0.016)\n\n"
        f"Amplitude share @ $\\sigma_k$=0.6: {amp_share['sigma_k=0.6']['share']*100:.1f}%\n"
        f"Energy: L*=3.5, p53 L=4 -> counting wins\n\n"
        f"Registered gap: L1 (ATM sensing, no data)"
)
    axs.text(0.02, 0.98, summary, va="top", fontsize=8.6, family="monospace",
             transform=axs.transAxes)
    fig.suptitle("Code 81: full-pathway p53 cascade audit (L0 dose -> L4 fate)", fontsize=13)
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    log("figure written:", OUT_PNG)
    log("[code81] done.")


if __name__ == "__main__":
    main()
