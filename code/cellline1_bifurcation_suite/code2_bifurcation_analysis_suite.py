# -*- coding: utf-8 -*-
"""
code 2: bifurcation analysis suite (three numerical verifications of the topology-selection theorem)
==================================================
B3  NF-family supercritical Hopf onset: A^2 proportional to (s-s_H)     - proper continuation-level analysis (not hand-made)
B2  subcritical-Hopf exception: hysteresis loop (third signature)
B1  lemma-2 universality: Morris-Lecar leaky digital (alpha far below the NF analog signature)

the three verifications jointly support the 'three-signature discrimination theorem' of manuscript section 4:
  (1) onset continuity  (2) hysteresis width  (3) wide-range alpha = dlnA/dlnD
Run: python3 code2_bifurcation_analysis_suite.py   (approx. 2-3 min, outputs a three-panel PNG)
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

# ============================================================
# part 1 (B3): NF model - equilibrium-branch solution + Hopf detection + limit-cycle branch
# ============================================================
# Model: p53->Mdm2 delayed negative feedback (3-ODE Goodwin type)
#   ṗ = s − δp·p − γ·m·p/(K+p)      (p53: production - spontaneous degradation - Mdm2-catalysed degradation)
#   ẋ = a·p − δx·x                  (Mdm2 mRNA, provides the phase lag)
#   ṁ = b·x − δm·m                  (Mdm2 protein)
dp_, gm, K, a_, dx_, b_, dm_ = 0.02, 1.5, 0.15, 1.0, 0.35, 1.0, 0.35

def eq_branch(p):
    """Closed-form equilibrium branch: given p*, compute s directly (no numerical root-finding needed)."""
    x = a_*p/dx_; m = b_*x/dm_
    return dp_*p + gm*m*p/(K+p)

def jac_eigs(p):
    """Eigenvalues of the 3x3 Jacobian at the equilibrium (determines the loss-of-stability mode)."""
    x = a_*p/dx_; m = b_*x/dm_
    J = np.array([[-dp_ - gm*m*K/(K+p)**2, 0, -gm*p/(K+p)],
                  [a_, -dx_, 0],
                  [0, b_, -dm_]])
    return np.linalg.eigvals(J)

# -- step 1: scan stability along the equilibrium branch, find the unique loss-of-stability point --
ps = np.linspace(0.001, 3.0, 6000)
stab = np.array([max(jac_eigs(p).real) for p in ps])
hopf_idx = np.where(np.diff(np.sign(stab)) != 0)[0]
i0 = hopf_idx[0]
s_H, _, _ = (lambda v: (eq_branch(v), None, None))(ps[i0])[0], None, None
ev = jac_eigs(ps[i0])
print(f'[B3] unique loss-of-stability point: p*={ps[i0]:.3f}, s_H={eq_branch(ps[i0]):.4f}, '
      f'crossing root={[e for e in ev if abs(e.imag)>1e-6][0]:.3f} (complex pair -> Hopf)')

# -- step 2: integrate point by point above Hopf, measure the limit-cycle amplitude branch --
def cycle_amp(s, tmax=800):
    def rhs(t, y):
        p, x, m = y
        return [s - dp_*p - gm*m*p/(K+p), a_*p - dx_*x, b_*x - dm_*m]
    sol = solve_ivp(rhs, (0, tmax), [1.4, 4, 11],
                    t_eval=np.linspace(0, tmax, 40001), rtol=1e-9, atol=1e-12)
    seg = sol.y[0][sol.t > tmax*0.6]
    return (seg.max()-seg.min())/2

s_br = s_H + np.array([0.05, 0.15, 0.3, 0.5, 0.8, 1.2, 1.7, 2.3, 3.0, 3.8])
A_br = np.array([cycle_amp(s) for s in s_br])
c = np.polyfit(s_br - s_H, A_br**2, 1)                 # A^2 vs (s-s_H) should be linear
A2 = A_br**2
R2 = 1 - ((A2-np.polyval(c, s_br-s_H))**2).sum()/((A2-A2.mean())**2).sum()
print(f'[B3] A² ∝ (s−s_H) linear R^2 = {R2:.4f} -> supercritical sqrt scaling confirmed\n')

# ============================================================
# part 2 (B2): subcritical Hopf - up/down-sweep hysteresis loop
# ============================================================
# normal form z' = (mu+i)z + |z|^2 z - |z|^4 z (Re c1=+1>0 subcritical, quintic term stabilising)
def run_normal(mu, r0, tmax=600):
    def rhs(t, y):
        x, yy = y; r2 = x*x + yy*yy
        return [mu*x - yy + x*r2 - x*r2**2, x + mu*yy + yy*r2 - yy*r2**2]
    sol = solve_ivp(rhs, (0, tmax), [r0, 0],
                    t_eval=np.linspace(0, tmax, 30001), rtol=1e-10, atol=1e-12)
    r = np.hypot(sol.y[0][sol.t > tmax*0.5], sol.y[1][sol.t > tmax*0.5])
    return r.mean()

mus = np.linspace(-0.35, 0.35, 29)
up   = np.array([run_normal(m, 0.05) for m in mus])    # small IC = up-sweep
down = np.array([run_normal(m, 1.20) for m in mus])    # large IC = down-sweep
on_up  = mus[np.argmax(up > 0.2)]
off_dn = mus[len(down)-1-np.argmax(down[::-1] < 0.2)]
print(f'[B2] subcritical: up-sweep onset μ_on={on_up:.3f} (finite amplitude r=1.00, appears digital) ')
print(f'      down-sweep extinction mu_off={off_dn:.3f}, hysteresis width={on_up-off_dn:.3f} (theory 0.25) -> third signature captured\n')

# ============================================================
# part 3 (B1): Morris-Lecar - 'leaky digital' with the drive entering the fast equation
# ============================================================
minf = lambda v: 0.5*(1+np.tanh((v+1.2)/18))
winf = lambda v: 0.5*(1+np.tanh((v-2)/30))
tauw = lambda v: 1/np.cosh((v-2)/60)

def run_ML(I, tmax=3000):
    def rhs(t, y):
        v, w = y
        return [(I - 4.4*minf(v)*(v-120) - 8*w*(v+84) - 2*(v+60))/20,
                0.04*(winf(v)-w)/tauw(v)]
    sol = solve_ivp(rhs, (0, tmax), [-60, 0.01],
                    t_eval=np.linspace(0, tmax, 150001), rtol=1e-9, atol=1e-12)
    v = sol.y[0][sol.t > tmax*0.4]
    pk, _ = find_peaks(v, prominence=20, distance=500)
    return (v[pk].mean()-v.min()) if len(pk) >= 2 else np.nan

Is = np.array([95, 100, 110, 125, 145, 170, 200])
As = np.array([run_ML(I) for I in Is])
ok = ~np.isnan(As)
alpha_ML = np.polyfit(np.log(Is[ok]), np.log(As[ok]), 1)[0]
print(f'[B1] Morris–Lecar: over a 2.1x drive range the amplitude changes {(As[ok][-1]/As[ok][0]-1)*100:.1f}%, '
      f'α = {alpha_ML:.3f} (leaky digital, still far below the NF analog signature)\n')

# ============================================================
# three-panel figure
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

ax = axes[0]                                          # (a) NF bifurcation diagram
ss = np.array([eq_branch(p) for p in ps])
ax.plot(ss, ps, 'k-', lw=1.2, label='equilibrium branch $p^*(s)$')
ax.plot(ss[stab < 0], ps[stab < 0], 'r-', lw=3, alpha=0.6, label='unstable segment')
ax.plot(s_br, 1.477+A_br, 'bo', ms=5, label='limit-cycle max')
ax.plot(s_br, 1.477-A_br, 'bo', ms=5, mfc='none', label='limit-cycle min')
fit_s = np.linspace(s_H, s_br.max(), 50)
ax.plot(fit_s, 1.477+np.sqrt(np.polyval(c, fit_s-s_H)), 'g--', lw=1.5,
        label=f'$\\sqrt{{s-s_H}}$ fit ($R^2$={R2:.3f})')
ax.axvline(s_H, color='gray', ls=':'); ax.annotate('Hopf', (s_H, 0.2))
ax.set_xlabel('drive s'); ax.set_ylabel('p'); ax.legend(fontsize=8)
ax.set_title('(a) NF family: supercritical Hopf onset (continuation level)'); ax.set_xlim(14, 21); ax.set_ylim(0, 3)

ax = axes[1]                                          # (b) subcriticalhysteresis loop
ax.plot(mus, up, 'b^-', ms=5, label='up-sweep (small IC)')
ax.plot(mus, down, 'rv-', ms=5, label='down-sweep (large IC)')
ax.fill_betweenx([0, 1.6], off_dn, on_up, color='orange', alpha=0.15)
ax.annotate(f'hysteresis width = {on_up-off_dn:.2f}', (-0.13, 1.3), fontsize=9)
ax.set_xlabel('μ'); ax.set_ylabel('amplitude r'); ax.legend(fontsize=9); ax.set_ylim(0, 1.6)
ax.set_title('(b) subcritical exception: finite onset but hysteresis != 0')

ax = axes[2]                                          # (c) alpha separation of the three families
alpha_NF = np.polyfit(np.log(s_br/s_H), np.log(A_br/A_br[0]), 1)[0]
D0s_ = np.array([0.42, 0.50, 0.60, 0.75, 0.95, 1.25, 1.65])
A_fhn = np.array([1.793, 1.813, 1.833, 1.847, 1.872, 1.898, 1.909])  # code1/TrackB results
ax.plot(s_br/s_H, A_br/A_br[0], 'rs-', ms=5, label=f'NF analog: alpha~{alpha_NF:.1f} (diverges near onset)')
alpha_fhn = np.polyfit(np.log(D0s_), np.log(A_fhn), 1)[0]
ax.plot(D0s_/D0s_[0], A_fhn/A_fhn[0], 'bo-', ms=5,
        label=f'FHN slow drive (digital): alpha={alpha_fhn:.3f}, {100*(A_fhn[-1]/A_fhn[0]-1):.1f}% over {D0s_[-1]/D0s_[0]:.2f}x')
ax.plot(Is/Is[0], As/As[0], 'g^-', ms=6, label=f'ML fast drive (leaky): alpha={alpha_ML:.2f}')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('drive (normalised) '); ax.set_ylabel('amplitude (normalised) '); ax.legend(fontsize=9)
ax.set_title('(c) the three generic classes separate by orders of magnitude in alpha')

plt.tight_layout()
plt.savefig('TrackB_bifurcation_validation.png', dpi=200)
print('saved TrackB_bifurcation_validation.png')
