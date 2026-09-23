#!/usr/bin/env python3
"""Carter-Penrose (conformal) diagram of the spherical collapse we simulate, with what the runs
actually cover marked on it.

The conformal structure is the standard Oppenheimer-Snyder-type one for spherically symmetric
collapse to a Schwarzschild black hole. It is NOT computed from the simulation: that would need
a double-null / characteristic construction of the conformal factor, which neither engine
provides. What comes from the runs are the annotated numbers and the qualitative shape of the
slicing -- DynGR1D uses polar slicing + radial gauge, AthenaK the Z4c puncture gauge (1+log),
and BOTH are singularity-avoiding, so their slices pile up outside the horizon and the runs end
on a collapsing lapse rather than on the singularity. That is exactly why the horizon mass is
read off the areal-gauge 2m/r -> 1 crossing instead of from the geometry near r=0.

Radial null rays are at 45 degrees everywhere, which fixes the layout: i^- = (0,-1),
i^0 = (1,0), and i^+ = (0.3, 0.7) where H^+, the singularity and scri^+ all meet.

Numbers: runs/dyngr1d_p1_p1b_collapse_font.csv (Font et al. 2002 unstable TOV, rho_c=7.993e-3,
M=1.448, R=5.838, cubic inward kick -0.01; VALIDATION.md 7.15).
Usage: python collapse_penrose.py ; writes analysis/collapse_penrose.png/.pdf
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

here = os.path.dirname(os.path.abspath(__file__))

# ---- conformal layout -------------------------------------------------------------
YM, YS = -1.0, 0.70     # i^-  ;  conformal time of the singularity
X0, Y0 = 1.0, 0.0       # i^0
YEH = 0.40              # H^+ forms on the axis here
XP, YP = YS - YEH, YS   # i^+  (end of H^+ = right end of the singularity = end of scri^+)


def scri_p(x):
    """future null infinity: from i^0 (1,0) up-left to i^+ (0.3,0.7)"""
    return 1.0 - x


# ---- star surface: leaves i^-, static-like bulge, then infall onto the singularity --
u = np.linspace(0, 1, 500)
sx = 0.46 * np.sin(np.pi * u ** 0.72) ** 0.9 * (1 - 0.80 * u ** 3)
sy = YM + (YS - YM) * u
# a timelike worldline stays strictly inside the diamond (left of scri^- and scri^+)
sx = np.minimum(sx, 0.90 * np.minimum(sy - YM, 1.0 - sy))
sx[0], sy[0] = 0.0, YM
sx[-1], sy[-1] = 0.11, YS
icross = np.argmin(np.abs(sx - (sy - YEH)))      # the surface crosses H^+ here
xc, yc = sx[icross], sy[icross]

fig, ax = plt.subplots(figsize=(7.0, 9.0))

# ---- regions ----------------------------------------------------------------------
ax.fill(np.concatenate([[0], sx, [0]]), np.concatenate([[YM], sy, [YS]]),
        color="#cfe3f7", zorder=1)                                        # the star
ax.add_patch(mpatches.Polygon([[0, YEH], [XP, YP], [0, YS]], closed=True,
                              facecolor="#f6dede", edgecolor="none", zorder=0))  # trapped

# ---- constant-time slices: singularity-avoiding, every one ends at i^0 --------------
xs = np.linspace(0, 0.995, 300)
for ya in np.linspace(-0.90, 0.52, 9):
    yy = ya * (1 - xs) ** 0.9 + 0.14 * max(ya, 0.0) * np.sin(np.pi * xs ** 0.6)
    ax.plot(xs, yy, color="#8a8a00", lw=0.75, alpha=0.9, zorder=3)

# ---- radial null rays (45 degrees) --------------------------------------------------
for y_em, col in ((-0.30, "#1b7a33"), (0.54, "#b22222")):
    t = np.linspace(0, 1.4, 400)
    x, y = t, y_em + t
    keep = (y <= YS) & (y <= scri_p(x))
    ax.plot(x[keep], y[keep], color=col, lw=1.4, ls=(0, (4, 2)), zorder=6)

# ---- boundaries and horizons --------------------------------------------------------
ax.plot([0, 0], [YM, YS], "k-", lw=2.2, zorder=5)                 # r = 0 (regular centre)
ax.plot([0, X0], [YM, Y0], "k-", lw=2.2, zorder=5)                # scri^-
ax.plot([X0, XP], [Y0, YP], "k-", lw=2.2, zorder=5)               # scri^+
xz = np.linspace(0, XP, 200)                                      # singularity (spacelike)
ax.plot(xz, YS + 0.010 * np.sign(np.sin(xz * 110)), color="k", lw=2.4, zorder=7)
ax.plot([0, XP], [YEH, YP], color="k", lw=2.2, zorder=6)          # H^+ (null)
ah_y = np.linspace(0.50, yc, 80)                                  # apparent horizon,
ax.plot((ah_y - 0.50) / (yc - 0.50) * xc, ah_y,
        color="#b22222", lw=2.6, zorder=7)                        # spacelike inside the matter
ax.plot([xc, XP], [yc, YP], color="#b22222", lw=2.6,
        ls=(0, (1, 1.6)), zorder=7)                               # = H^+ = r=2M in vacuum
ax.plot(sx, sy, color="#1f5fa8", lw=2.6, zorder=6)                # star surface

# ---- labels: only the conformal infinities on the diagram, everything else in a legend
ax.text(-0.055, -0.30, "$r=0$  (regular centre)", rotation=90,
        ha="center", va="center", fontsize=9)
ax.text(XP / 2 - 0.04, YS + 0.07, "singularity   $r=0$", ha="center", fontsize=10.5)
ax.text(X0 + 0.03, Y0, "$i^0$", fontsize=12, va="center")
ax.text(XP + 0.05, YP + 0.02, "$i^+$", fontsize=12)
ax.text(-0.03, YM - 0.09, "$i^-$", fontsize=12, ha="center")
ax.text(0.60, -0.47, "$\\mathscr{I}^-$", fontsize=13, rotation=44)
ax.text(0.72, 0.335, "$\\mathscr{I}^+$", fontsize=13, rotation=-44)
ax.text(0.150, -0.40, "star", fontsize=11, color="#14487f")
ax.annotate("$\\alpha_c\\to10^{-3}$:\nthe runs stop\nhere", xy=(0.075, 0.395), xytext=(0.255, 0.30),
            fontsize=8.5, color="#6b6b00", ha="left",
            arrowprops=dict(arrowstyle="->", lw=0.9, color="#6b6b00"))

handles = [
    mlines.Line2D([], [], color="#1f5fa8", lw=2.4, label="star surface"),
    mpatches.Patch(facecolor="#cfe3f7", label="matter (the star)"),
    mlines.Line2D([], [], color="k", lw=2.0, label="event horizon $H^+$ (null)"),
    mlines.Line2D([], [], color="#b22222", lw=2.4,
                  label="apparent horizon ($2m/r\\to1$): $M_{\\rm AH}=1.273$"),
    mpatches.Patch(facecolor="#f6dede", label="trapped region"),
    mlines.Line2D([], [], color="#8a8a00", lw=0.9,
                  label="constant-time slices (polar slicing / 1+log)"),
    mlines.Line2D([], [], color="#1b7a33", lw=1.3, ls=(0, (4, 2)),
                  label="radial null ray: escapes to $\\mathscr{I}^+$"),
    mlines.Line2D([], [], color="#b22222", lw=1.3, ls=(0, (4, 2)),
                  label="radial null ray: trapped, ends on $r=0$"),
]
ax.legend(handles=handles, fontsize=8.2, loc="upper left",
          bbox_to_anchor=(0.0, 0.14), bbox_transform=ax.transAxes, framealpha=0.95)

note = ("DynGR1D, Font et al. 2002 unstable star\n"
        "$\\rho_c=7.993\\times10^{-3}$, $M=1.448$, $R=5.838$, kick $-1\\%$\n"
        "horizon at $t=53.7\\,M_\\odot$;  central proper time $\\tau_c=10.0$\n"
        "$M_{\\rm AH}=1.273$  (88% of $M$);   max $2m/r=0.96$\n\n"
        "Both gauges are singularity-avoiding: nothing\n"
        "above the red curve is covered by either code.")
ax.text(0.995, 0.14, note, fontsize=8.2, va="bottom", ha="right", transform=ax.transAxes,
        bbox=dict(boxstyle="round,pad=0.45", fc="#fffbe6", ec="#9a9a9a", lw=0.8))

ax.set_xlim(-0.17, 1.26)
ax.set_ylim(YM - 0.70, YS + 0.28)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title("Carter–Penrose diagram of the simulated collapse", fontsize=12.5, pad=10)
plt.tight_layout()
for ext in ("png", "pdf"):
    plt.savefig(os.path.join(here, f"collapse_penrose.{ext}"), dpi=170)
print("figure:", os.path.join(here, "collapse_penrose.png"))
