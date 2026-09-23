#!/usr/bin/env python3
"""Phase-0 analysis of the AthenaK TOV baselines (progress/plan_athenak_bdnk_collapse.md §1):
max density and min lapse histories, the Cowling radial spectrum against the linear modes of
this star (F = 2.686, H1 = 4.550, H2 = 6.341, H3 = 8.108 kHz from radial_cowling_spectrum), and
the DG-FD hybrid of VALIDATION 7.13 on the same star for comparison when its series exists.
Usage: python p0_history.py [run tags...]   (default: the three 64^3 runs); writes
analysis/p0_history.png and prints the summary table."""
import os, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__)); runs = os.path.join(here, "..", "runs")
Msun_s = 4.925490947e-6
F_lin = [2.6861, 4.5495, 6.3414, 8.1082]     # Cowling (radial_cowling_spectrum, this star)
F_gr  = [1.4425, 3.9541, 5.9149]              # full GR (chandrasekhar_radial_omega2, same star)
tags = sys.argv[1:] or ["p0a_tov_cowling_K64_static", "p0a_tov_cowling_K64_kick", "p0b_tov_z4c_K64_static"]
def hst(tag):
    f = os.path.join(runs, tag, "tov.user.hst")
    if not os.path.exists(f): return None
    d = np.loadtxt(f); return d[:, 0], d[:, 2], d[:, 3]
def spectrum(t, q, fmin=1.0, fmax=10.0, n=4000):
    q = q - q.mean(); w = 0.5 - 0.5*np.cos(2*np.pi*np.arange(len(q))/(len(q)-1)); q = q*w
    f = np.linspace(fmin, fmax, n); om = 2*np.pi*f*1e3*Msun_s
    return f, np.array([abs(np.sum(q*np.exp(-1j*o*t)))**2 for o in om])
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
print(f"{'run':34s} {'t_end':>7s} {'rhomax/rhomax0-1':>17s} {'alpha_min':>10s} {'peaks [kHz]':>28s}")
for i, tag in enumerate(tags):
    d = hst(tag)
    if d is None: print(f"{tag:34s} (no history yet)"); continue
    t, rm, am = d; c = f"C{i}"
    ax[0].plot(t*Msun_s*1e3, rm/rm[0], c, lw=0.7, label=tag)
    ax[1].plot(t, am, c, lw=0.7, label=tag)
    f, P = spectrum(t, rm/rm[0]); P /= P.max()
    ax[2].semilogy(f, 10.0**(-2*i)*P, c, lw=0.8, label=tag + ("" if i == 0 else rf" ($\times10^{{-{2*i}}}$)"))
    pk = [f[j] for j in range(1, len(f)-1) if P[j] > P[j-1] and P[j] >= P[j+1] and P[j] > 1e-3]
    pk = sorted(pk, key=lambda x: -P[np.argmin(abs(f-x))])[:4]
    print(f"{tag:34s} {t[-1]:7.0f} {rm[-1]/rm[0]-1:17.3e} {am[-1]:10.4f} {str(np.round(sorted(pk),3)):>28s}")
hyb = os.path.join(here, "..", "..", "BDNKStar", "repro", "data", "dgcart3dfd_tov_series_K6_p5_unseeded.csv")
if os.path.exists(hyb):
    h = np.genfromtxt(hyb, delimiter=",", names=True)
    ax[0].plot(h["t_ms"], h["rhomax"]/h["rhomax"][0], "k--", lw=0.6, label="DG-FD hybrid K=6 (VALIDATION 7.13)")
for fl, name in zip(F_lin, ["F", r"H$_1$", r"H$_2$", r"H$_3$"]):
    ax[2].axvline(fl, color="k", ls=":", lw=0.8); ax[2].text(fl, 2.0, name, ha="center", fontsize=8)
for fl, name in zip(F_gr, ["F", r"H$_1$", r"H$_2$"]):
    ax[2].axvline(fl, color="C3", ls="-.", lw=0.8); ax[2].text(fl, 8e-9, name, ha="center", fontsize=8, color="C3")
ax[2].plot([], [], "k:", label="Cowling modes"); ax[2].plot([], [], "C3-.", label="full-GR modes")
ax[0].set_xlabel("time [ms]"); ax[0].set_ylabel(r"$\max\rho/\max\rho(0)$"); ax[0].set_title("AthenaK TOV: max density"); ax[0].legend(fontsize=7)
ax[1].set_xlabel(r"$t\ [M_\odot]$"); ax[1].set_ylabel(r"$\min\alpha$"); ax[1].set_title("minimum lapse (gauge settling with Z4c)"); ax[1].legend(fontsize=7)
ax[2].set_xlabel("f [kHz]"); ax[2].set_ylabel("power (arb.)"); ax[2].set_xlim(1, 10); ax[2].set_ylim(1e-9, 5)
ax[2].set_title("radial spectrum vs linear Cowling modes"); ax[2].legend(fontsize=7, loc="lower left")
plt.tight_layout(); out = os.path.join(here, "p0_history.png"); plt.savefig(out, dpi=150); print("figure:", out)
