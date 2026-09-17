#!/usr/bin/env python3
"""Figure for the 3D cubed-sphere DG star (DGBall3D): density error, central density, radial and
ℓ=2 spectra. Reads repro/data/dgball3d_hkt_series_*.csv; writes paper/figs/dgball3d_hkt.png/.pdf."""
import numpy as np, matplotlib, os
matplotlib.use("Agg"); import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__)); data = os.path.join(here, "data"); out = os.path.join(here, "..", "paper", "figs")
Msun_s = 4.925490947e-6
def load(tag):
    f = os.path.join(data, f"dgball3d_hkt_series_{tag}.csv"); return np.genfromtxt(f, delimiter=",", names=True) if os.path.exists(f) else None
def spectrum(t, q, fmin=0.3, fmax=10.0, n=4000):
    q = q - q.mean(); w = 0.5 - 0.5*np.cos(2*np.pi*np.arange(len(q))/(len(q)-1)); q = q*w
    f = np.linspace(fmin, fmax, n); om = 2*np.pi*f*1e3*Msun_s
    return f, np.array([abs(np.sum(q*np.exp(-1j*o*t)))**2 for o in om])
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
for tag, lab, c in [("static_noWB", "static, no subtraction", "C1"), ("static_WB", "static, subtraction", "C0"),
                    ("l0_WB", "radial seed, subtraction", "C2"), ("l2_WB", "ℓ=2 seed, subtraction", "C3"), ("l2_noWB", "ℓ=2 seed, no subtraction", "C4")]:
    d = load(tag)
    if d is None: continue
    ax[0].semilogy(d["t"], np.maximum(d["errD"], 1e-17), c, lw=0.9, label=lab)
    ax[1].plot(d["t"], (d["rhoc"]/d["rhoc"][0]-1)*1e3, c, lw=0.8, label=lab)
ax[0].set_xlabel(r"$t\ [M_\odot]$"); ax[0].set_ylabel(r"err$[\tilde D]$"); ax[0].set_title("Density error"); ax[0].legend(fontsize=7)
ax[1].set_xlabel(r"$t\ [M_\odot]$"); ax[1].set_ylabel(r"$(\rho_c-\rho_{c,0})/\rho_{c,0}\ [10^{-3}]$"); ax[1].set_title("Central density"); ax[1].legend(fontsize=7)
d = load("l0_WB")
if d is not None:
    f, P = spectrum(d["t"], d["rhoc"]/d["rhoc"][0]); ax[2].semilogy(f, P/P.max(), "C2", lw=0.8, label=r"$\rho_c$, radial seed")
    for fl, name in [(2.6861, "F"), (4.5497, r"H$_1$"), (6.3414, r"H$_2$")]:
        ax[2].axvline(fl, color="k", ls=":", lw=0.8); ax[2].text(fl, 2.0, name, ha="center", fontsize=8)
for tag, lab, c, off in [("l2_WB", r"$q_{20}$, $\ell=2$ seed, nt=2 ($\times10^{-3}$)", "C3", 1e-3),
                         ("l2_WB_nt2pt5", r"$q_{20}$, nt=2, $p_t=5$ ($\times10^{-5}$)", "C4", 1e-5),
                         ("l2_WB_nt3s6", r"$q_{20}$, nt=3, shell filter s=6 ($\times10^{-7}$)", "C5", 1e-7)]:
    d = load(tag)
    if d is None: continue
    f, P = spectrum(d["t"], d["q20"]/abs(d["q00"][0])); ax[2].semilogy(f, off*P/P.max(), c, lw=0.8, label=lab)
ax[2].axvline(1.88291, color="C3", ls="--", lw=0.8); ax[2].text(1.88291, 2.0, "f", ha="center", fontsize=8, color="C3")
ax[2].set_xlabel("f [kHz]"); ax[2].set_ylabel("power (arbitrary units)"); ax[2].set_title("Spectra against linear theory"); ax[2].set_xlim(0.3, 10); ax[2].set_ylim(1e-12, 3); ax[2].legend(fontsize=7)
plt.tight_layout()
for ext in ("png", "pdf"): plt.savefig(os.path.join(out, f"dgball3d_hkt.{ext}"), dpi=160)
print("figure written")
