#!/usr/bin/env python3
"""Figure for the DG/finite-difference hybrid on the radial Cowling star (DGStarFD): density
error, central density, and the radial spectrum converging to linear theory.
Reads repro/data/dgstarfd_hybrid_series_*.csv; writes paper/figs/dgstarfd_hybrid.png/.pdf."""
import numpy as np, os, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__)); data = os.path.join(here, "data")
out = os.path.join(here, "..", "paper", "figs"); os.makedirs(out, exist_ok=True)
Msun_s = 4.925490947e-6
F_lin = [2.6861, 4.5495, 6.3414, 8.1082]
def load(tag):
    f = os.path.join(data, f"dgstarfd_hybrid_series_{tag}.csv")
    return np.genfromtxt(f, delimiter=",", names=True) if os.path.exists(f) else None
def spectrum(t, q, fmin=1.0, fmax=10.0, n=4000):
    q = q - q.mean(); w = 0.5 - 0.5*np.cos(2*np.pi*np.arange(len(q))/(len(q)-1)); q = q*w
    f = np.linspace(fmin, fmax, n); om = 2*np.pi*f*1e3*Msun_s
    return f, np.array([abs(np.sum(q*np.exp(-1j*o*t)))**2 for o in om])
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
Ks = [(10, "C0"), (14, "C1"), (20, "C2"), (28, "C3")]
for K, c in Ks:
    d = load(f"uni_p5_K{K}_unseeded")
    if d is None: continue
    ax[0].semilogy(d["t"], d["errD"], c, lw=0.9, label=f"hybrid, K={K} per side")
d = load("uni_p5_K20_purefv")
if d is not None:
    ax[0].semilogy(d["t"], d["errD"], "k--", lw=0.9, label="pure finite volume, K=20")
ax[0].set_xlabel(r"$t\ [M_\odot]$"); ax[0].set_ylabel(r"err$[\tilde D]$")
ax[0].set_title("Density error"); ax[0].legend(fontsize=7); ax[0].set_ylim(1e-4, 1)
for K, c in Ks:
    d = load(f"uni_p5_K{K}_seeded")
    if d is None: continue
    ax[1].plot(d["t"], (d["rhoc"]/d["rhoc"][0]-1)*1e3, c, lw=0.7, label=f"K={K}")
ax[1].set_xlabel(r"$t\ [M_\odot]$"); ax[1].set_ylabel(r"$(\rho_c-\rho_{c,0})/\rho_{c,0}\ [10^{-3}]$")
ax[1].set_title(r"Central density, seeded $v=10^{-3}\sin(\pi r/R)$"); ax[1].legend(fontsize=7); ax[1].set_xlim(0, 4000)
for i, (K, c) in enumerate(Ks):
    d = load(f"uni_p5_K{K}_seeded")
    if d is None: continue
    f, P = spectrum(d["t"], d["rhoc"]/d["rhoc"][0])
    ax[2].semilogy(f, 10.0**(-2*i)*P/P.max(), c, lw=0.8, label=f"K={K}" + ("" if i == 0 else rf" ($\times10^{{-{2*i}}}$)"))
for fl, name in zip(F_lin, ["F", r"H$_1$", r"H$_2$", r"H$_3$"]):
    ax[2].axvline(fl, color="k", ls=":", lw=0.8); ax[2].text(fl, 2.0, name, ha="center", fontsize=8)
ax[2].set_xlabel("f [kHz]"); ax[2].set_ylabel("power (arbitrary units)")
ax[2].set_title("Radial spectrum vs linear theory"); ax[2].set_xlim(1, 10); ax[2].set_ylim(1e-11, 5); ax[2].legend(fontsize=7, loc="lower left")
plt.tight_layout()
for ext in ("png", "pdf"): plt.savefig(os.path.join(out, f"dgstarfd_hybrid.{ext}"), dpi=160)
print("figure written:", os.path.join(out, "dgstarfd_hybrid.png"))
