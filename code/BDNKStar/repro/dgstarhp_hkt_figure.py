#!/usr/bin/env python3
"""Figure for the 1D hp-adapted DG star (DGStarHP) — the analogues of Figs. 9, 12 and 13 of
Hébert, Kidder & Teukolsky (2018): density error, central density and radial spectrum.
Reads repro/data/dgstarhp_hkt_series_*.csv (from repro/dgstarhp_hkt.jl); writes
paper/figs/dgstarhp_hkt.png/.pdf. Simple titles and labels."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, os, glob
here = os.path.dirname(os.path.abspath(__file__)); data = os.path.join(here, "data")
out = os.path.join(here, "..", "paper", "figs"); os.makedirs(out, exist_ok=True)
Msun_s = 4.925490947e-6                      # seconds per M_sun
def load(tag):
    f = os.path.join(data, f"dgstarhp_hkt_series_{tag}.csv")
    return np.genfromtxt(f, delimiter=",", names=True) if os.path.exists(f) else None
def spectrum(t, q, fmax=15.0, n=6000):
    q = q - q.mean(); w = 0.5 - 0.5*np.cos(2*np.pi*np.arange(len(q))/(len(q)-1)); q = q*w
    f = np.linspace(0.3, fmax, n)                                    # kHz
    om = 2*np.pi*f*1e3*Msun_s                                        # rad per M_sun
    P = np.array([abs(np.sum(q*np.exp(-1j*o*t)))**2 for o in om])
    return f, P
F_lin = [2.6861, 4.5497, 6.3414, 8.1081, 9.8633, 11.6122]           # Cowling radial modes, kHz
cases = [("I1_minmod_noWB_static",  "I1, minmod, no subtraction (the paper's scheme)", "C1"),
         ("I1R_minmod_noWB_static", "I1R, minmod, no subtraction", "C4"),
         ("I2_wb_noWB_static",      "I2, wb limiter, no subtraction", "C2"),
         ("I2R_wb_noWB_static",     "I2R, wb limiter, no subtraction", "C5"),
         ("I2_wb_WB_static",        "I2, wb limiter, subtraction", "C0")]
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
for tag, lab, c in cases:
    d = load(tag)
    if d is None: continue
    ax[0].semilogy(d["t"], np.maximum(d["errD"], 1e-16), c, lw=1, label=lab)
    ax[1].plot(d["t"], d["rhoc_rel"]*1e3, c, lw=0.8, label=lab)
ax[0].set_xlabel(r"$t\ [M_\odot]$"); ax[0].set_ylabel(r"err$[\tilde D]$"); ax[0].set_title("Density error")
ax[0].legend(fontsize=7, loc="lower right"); ax[0].set_xlim(0, 10000)
ax[1].set_xlabel(r"$t\ [M_\odot]$"); ax[1].set_ylabel(r"$(\rho_c-\rho_{c,0})/\rho_{c,0}\ [10^{-3}]$"); ax[1].set_title("Central density")
ax[1].set_xlim(0, 10000); ax[1].legend(fontsize=7, loc="best")
for tag, lab, c, off in [("I1_minmod_noWB_static", "I1, minmod, unseeded (paper protocol)", "C1", 1.0),
                         ("I1R_minmod_noWB_seeded", "I1R, minmod, seeded", "C4", 1e-2),
                         ("I2_wb_WB_seeded", "I2, wb limiter, subtraction, seeded", "C0", 1e-4),
                         ("I2R_wb_WB_seeded", "I2R, wb limiter, subtraction, seeded", "C5", 1e-6)]:
    d = load(tag)
    if d is None: continue
    m = d["t"] <= 4000; f, P = spectrum(d["t"][m], d["rhoc_rel"][m])
    ax[2].semilogy(f, off*P/P.max(), c, lw=0.8, label=lab)
for k, fl in enumerate(F_lin):
    ax[2].axvline(fl, color="k", ls=":", lw=0.8)
ax[2].set_xlabel("f [kHz]"); ax[2].set_ylabel(r"$|\hat\rho_c|^2$ (arbitrary units)"); ax[2].set_title("Radial spectrum of the central density, t < 4000")
ax[2].set_xlim(0.3, 15); ax[2].set_ylim(1e-13, 3); ax[2].legend(fontsize=7, loc="upper right")
ax[2].text(F_lin[0], 2.0, "F", ha="center", fontsize=8); ax[2].text(F_lin[1], 2.0, r"H$_1$", ha="center", fontsize=8)
ax[2].text(F_lin[2], 2.0, r"H$_2$", ha="center", fontsize=8)
plt.tight_layout()
for ext in ("png", "pdf"):
    plt.savefig(os.path.join(out, f"dgstarhp_hkt.{ext}"), dpi=160)
print("figure written:", os.path.join(out, "dgstarhp_hkt.png"))
