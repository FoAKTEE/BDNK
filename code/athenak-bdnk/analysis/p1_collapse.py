#!/usr/bin/env python3
"""Phase-1 analysis: ideal-fluid collapse / migration of the TOV star in AthenaK
(progress/plan_athenak_bdnk_collapse.md §1). Plots, for each run, the maximum density and the
minimum lapse against BOTH coordinate time and the central proper time tau = int(alpha_min dt)
— the gauge-invariant clock that lets the puncture-gauge AthenaK runs be compared with the
polar-sliced DynGR1D ones — plus the baryon mass. The DynGR1D curves are overlaid when their
CSVs exist, dashed. Its radial-mode sector is now validated to 0.05% against two independent
frequency-domain eigensolvers (VALIDATION.md 7.16), but the collapse and horizon numbers here
have no independent confirmation — obtaining one is what this comparison is for.
Usage: python p1_collapse.py [tags...]; writes analysis/p1_collapse.png and a summary table."""
import os, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
here = os.path.dirname(os.path.abspath(__file__)); runs = os.path.join(here, "..", "runs")
Msun_ms = 4.925490947e-3
TAGS = sys.argv[1:] or ["p1a_collapse_eps0003_K64", "p1b_collapse_font_K64", "p1b_migration_font_K64"]
# rho_c of each run's star, for normalisation and for the "did it migrate" test
RHOC = {"p1a_collapse_eps0003_K64": 2.4162e-3, "p1b_collapse_font_K64": 7.993e-3,
        "p1b_migration_font_K64": 7.993e-3}
DYN = {"p1a_collapse_eps0003_K64": "dyngr1d_p1_p1a_collapse_eps0003",
       "p1b_collapse_font_K64": "dyngr1d_p1_p1b_collapse_font",
       "p1b_migration_font_K64": "dyngr1d_p1_p1b_migration_font"}

# Resolution variants are named <base>_L1, _L2, ... for the centre-refinement sequence; they are
# the same star, so they share the base run's normalisation and DynGR1D reference curve.
def _base_tag(tag):
    for b in RHOC:
        if tag == b or tag.startswith(b + "_"):
            return b
    return None

def _rhoc_for(tag, fallback):
    b = _base_tag(tag)
    return RHOC[b] if b else fallback

def _dyn_for(tag):
    b = _base_tag(tag)
    return DYN.get(b, "") if b else ""

# A Z4c run that loses the star reports rho-max pinned at the density floor and alpha-min as
# DBL_MAX (the minimum of an empty set), and every Z4c norm and the diagnostic volume go to
# exactly zero. Integrating tau = int(alpha dt) through that gives inf and makes a run whose
# lapse collapsed read out as "migration/expansion" -- the opposite of what happened. Every
# series is therefore cut at the last physical sample and the breakdown time is reported.
DFLOOR = 1.0e-10

def _last_sane(rm, am):
    """Index one past the last sample that is still physical, or len() if the run is clean."""
    bad = ~np.isfinite(rm) | ~np.isfinite(am) | (rm <= 1.01*DFLOOR) | (am > 1.0) | (am <= 0.0)
    return int(np.argmax(bad)) if bad.any() else len(rm)

def load(tag):
    f = os.path.join(runs, tag, "tov.user.hst")
    if not os.path.exists(f): return None
    d = np.loadtxt(f)
    if d.ndim < 2 or len(d) < 3: return None
    t, rm, am = d[:, 0], d[:, 2], d[:, 3]
    k = _last_sane(rm, am)
    tbreak = t[k] if k < len(t) else None          # None => the run never broke down
    t, rm, am = t[:k], rm[:k], am[:k]
    if len(t) < 3: return None
    tau = np.concatenate(([0.0], np.cumsum(0.5*(am[1:] + am[:-1])*np.diff(t))))
    mb = None
    g = os.path.join(runs, tag, "tov.mhd.hst")
    if os.path.exists(g):
        e = np.loadtxt(g)
        if e.ndim == 2 and len(e) > 2:
            j = min(len(e), k)
            mb = (e[:j, 0], e[:j, 2])
    return t, rm, am, tau, mb, tbreak

fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
print(f"{'run':30s} {'t_end':>7s} {'tau_end':>8s} {'max rho/rho_c0':>15s} {'min alpha':>10s} {'Mb drift':>10s} {'breakdown':>10s}  verdict")
for i, tag in enumerate(TAGS):
    d = load(tag)
    if d is None: print(f"{tag:30s} (no history yet)"); continue
    t, rm, am, tau, mb, tbreak = d; c = f"C{i}"
    rho0 = _rhoc_for(tag, rm[0])
    ax[0].plot(tau, rm/rho0, c, lw=0.9, label=tag)
    ax[1].plot(tau, am, c, lw=0.9, label=tag)
    if mb is not None: ax[2].plot(mb[0], mb[1]/mb[1][0] - 1, c, lw=0.9, label=tag)
    dyn = os.path.join(runs, _dyn_for(tag) + ".csv")
    if os.path.exists(dyn):
        e = np.genfromtxt(dyn, delimiter=",", names=True)
        ax[0].plot(e["tau_c"], e["rhoc"]/rho0, c, lw=0.7, ls="--", alpha=0.6)
        ax[1].plot(e["tau_c"], e["alphac"], c, lw=0.7, ls="--", alpha=0.6)
    # The verdict is read off the WHOLE usable window, not off the last sample: a collapsing
    # lapse is the signature, and after it the matter leaves the grid and the last sample says
    # nothing. Lapse down by >2x with the density rising is collapse; density down by >2x with
    # the lapse recovered is migration.
    amin, rmax = am.min(), rm.max()
    if amin < 0.5*am[0] and rmax > 1.2*rho0:
        v = "COLLAPSE"
    elif rm[-1] < 0.5*rho0 and am[-1] > am[0]:
        v = "migration/expansion"
    else:
        v = "in progress"
    if tbreak is not None:
        v += f" (run broke down at t={tbreak:.1f})"
    mbd = f"{mb[1][-1]/mb[1][0]-1:+.1e}" if mb is not None else "  --"
    brk = f"{tbreak:10.1f}" if tbreak is not None else "      none"
    print(f"{tag:30s} {t[-1]:7.1f} {tau[-1]:8.2f} {rmax/rho0:15.4f} {amin:10.4f} {mbd:>10s} {brk}  {v}")
ax[0].plot([], [], "k--", alpha=0.6, label="DynGR1D (collapse not cross-checked, §7.16)")
ax[0].set_xlabel(r"central proper time $\tau=\int\alpha_{\min}dt\ [M_\odot]$")
ax[0].set_ylabel(r"$\max\rho/\rho_c(0)$"); ax[0].set_yscale("log")
ax[0].set_title("collapse / migration"); ax[0].legend(fontsize=7)
ax[1].set_xlabel(r"$\tau\ [M_\odot]$"); ax[1].set_ylabel(r"$\min\alpha$")
ax[1].set_title("lapse collapse (horizon formation)"); ax[1].legend(fontsize=7)
ax[2].set_xlabel(r"$t\ [M_\odot]$"); ax[2].set_ylabel(r"$M_b/M_b(0)-1$")
ax[2].set_title("baryon mass conservation"); ax[2].legend(fontsize=7)
plt.tight_layout(); out = os.path.join(here, "p1_collapse.png"); plt.savefig(out, dpi=150)
print("figure:", out)
