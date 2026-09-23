#!/usr/bin/env python3
"""Evidence figure for DynGR1D well-balancing (VALIDATION.md 7.15 and 7.16).

(a) Raw static momentum residual against resolution. The plain operator's residual is spread
    through the bulk and falls as dr^2; hydrostatic reconstruction leaves the bulk at ROUND-OFF
    (flat, ~1e-13, eight decades lower) and what remains sits in the one cell where the star
    meets the constant-density artificial atmosphere. The dotted lines are the pre-2026-09-22
    operator, whose residual did not fall with dr at all.
(b) The unsubtracted static star. With the plain operator the centre starts moving within one
    M_sun, because the imbalance is everywhere; with hydrostatic reconstruction it sits at
    round-off for ~45 M_sun -- the time a signal needs to travel in from the surface -- and that
    delay GROWS with resolution as the unbalanced surface shell thins.
(c) The Phase-1 verdicts: the unstable star migrates when kicked outward and collapses when
    kicked inward, while the stable star merely oscillates.

Usage: python dyngr1d_wb_figure.py ; writes analysis/dyngr1d_wb.png"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
runs = os.path.join(here, "..", "runs")
fig, ax = plt.subplots(1, 3, figsize=(15.4, 4.5))

# ---- (a) residual vs resolution -------------------------------------------------------
r = np.genfromtxt(os.path.join(runs, "dyngr1d_wb_residual.csv"), delimiter=",", names=True,
                  dtype=None, encoding="utf-8")
for star, c in (("stable", "C0"), ("unstable", "C1")):
    ms = r["star"] == star
    p = ms & (r["scheme"] == "plain")
    h = ms & (r["scheme"] == "hydrostatic")
    ax[0].loglog(r["dr"][p], r["residual"][p], c + "o-", lw=1.2, ms=4,
                 label=f"{star}: plain operator")
    ax[0].loglog(r["dr"][h], r["residual"][h], c + "s--", lw=1.2, ms=4, alpha=0.75,
                 label=f"{star}: hydrostatic, whole star")
    ax[0].loglog(r["dr"][h], r["residual_bulk"][h], c + "^-", lw=1.8, ms=5,
                 label=f"{star}: hydrostatic, $r<0.95R$")
    ax[0].axhline(r["before_fix"][ms][0], color=c, ls=":", lw=1.2,
                  label=f"{star}: before the fix (any $\\Delta r$)")
dr = r["dr"][(r["star"] == "stable") & (r["scheme"] == "plain")]
y0 = r["residual"][(r["star"] == "stable") & (r["scheme"] == "plain")][0]
ax[0].loglog(dr, y0*(dr/dr[0])**2, "k--", lw=0.8, label=r"$\Delta r^2$")
ax[0].set_xlabel(r"$\Delta r\ [M_\odot]$")
ax[0].set_ylabel(r"$|\partial_t \tilde S|\,/\,|{\rm gravity}|$")
ax[0].set_title("Raw static momentum residual (nothing subtracted)")
ax[0].legend(fontsize=6.2, ncol=1, loc="center", framealpha=0.95,
             bbox_to_anchor=(0.5, 0.40), bbox_transform=ax[0].transAxes)
ax[0].grid(alpha=0.3, which="both")
ax[0].set_xticks([4e-3, 1e-2, 3e-2])
ax[0].set_xticklabels(["0.004", "0.01", "0.03"])
ax[0].set_xticks([], minor=True)

# ---- (b) the unsubtracted static star --------------------------------------------------
st = np.genfromtxt(os.path.join(runs, "dyngr1d_wb_static.csv"), delimiter=",", names=True,
                   dtype=None, encoding="utf-8")
NS = sorted(set(st["N"].astype(int)))
for scheme, c in (("plain", "C3"), ("hydrostatic", "C2")):
    for j, N in enumerate(NS):
        m = (st["scheme"] == scheme) & (st["N"] == N)
        if not m.any():
            continue
        d = np.abs(st["rhoc_over_rhoc0"][m] - 1)
        ax[1].semilogy(st["t"][m], np.maximum(d, 1e-16), c, lw=0.9, alpha=0.45 + 0.25*j,
                       label=f"{scheme}, N={N}")
        k = np.argmax(d > 1e-10) if (d > 1e-10).any() else None
        if scheme == "hydrostatic" and k:
            ax[1].plot(st["t"][m][k], d[k], c + "o", ms=5)
ax[1].axhline(0.867, color="k", ls=":", lw=1.2, label="before the fix (N=1600): 0.87")
ax[1].annotate("hydrostatic: centre still at round-off,\nno signal has reached it yet",
               xy=(22, 1.2e-16), xytext=(60, 4e-15), fontsize=7.5, color="C2",
               arrowprops=dict(arrowstyle="->", lw=0.9, color="C2"))
ax[1].annotate("plain: the centre moves\nwithin one $M_\\odot$", xy=(2.5, 2e-6),
               xytext=(18, 2e-9), fontsize=7.5, color="C3",
               arrowprops=dict(arrowstyle="->", lw=0.9, color="C3"))
ax[1].set_xlabel(r"$t\ [M_\odot]$")
ax[1].set_ylabel(r"$|\rho_c/\rho_c(0)-1|$")
ax[1].set_ylim(3e-17, 3.0)
ax[1].set_title("Static star, NO equilibrium subtraction")
ax[1].legend(fontsize=6.4, ncol=1, loc="center right", framealpha=0.95)

# ---- (c) the Phase-1 verdicts -----------------------------------------------------------
lab = {"p1a_collapse_eps0003": ("stable star, kick $-0.03$", "C0"),
       "p1b_collapse_font": ("unstable star, kick $-0.01$", "C3"),
       "p1b_migration_font": ("unstable star, kick $+0.05$", "C2")}
for tag, (l, c) in lab.items():
    f = os.path.join(runs, f"dyngr1d_p1_{tag}.csv")
    if not os.path.exists(f):
        continue
    d = np.genfromtxt(f, delimiter=",", names=True)
    ax[2].semilogy(d["tau_c"], d["rhoc"], c, lw=1.0, label=l)
ax[2].axhline(1.3e-3, color="k", ls="--", lw=0.8,
              label=r"stable-branch $\rho_c$ (migration target)")
ax[2].set_xlabel(r"central proper time $\tau=\int\alpha_c dt\ [M_\odot]$")
ax[2].set_ylabel(r"$\rho_c$")
ax[2].set_title("Phase-1 verdicts")
ax[2].legend(fontsize=7)

plt.tight_layout()
out = os.path.join(here, "dyngr1d_wb.png")
plt.savefig(out, dpi=150)
print("figure:", out)
