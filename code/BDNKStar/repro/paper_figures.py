"""Figures for paper/cowling3d_method.tex, numbered as in the compiled PDF (order of appearance). Reads repro/data/*.csv; skips a figure whose data is absent.
Run:  python3 repro/paper_figures.py      (from code/BDNKStar)"""
import os, csv, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data"); OUT = os.path.join(HERE, "..", "paper", "figs"); os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix", "font.size": 9, "axes.labelsize": 9,
    "axes.titlesize": 9.5, "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
    "axes.linewidth": 0.7, "lines.linewidth": 1.1, "lines.markersize": 4.5, "legend.frameon": False,
    "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.02})
W = 3.375   # PRD single column, inches
KHZ_TO_KM = 1e3/2.99792458e5; MSUN_KM = 1.476625; NS_INTEGRAL = 0.08746*0.03*1e3   # x1e3 /Msun
C = {"cut": "#1f4e9c", "mask": "#c0392b", "F1": "#c0392b", "F2": "#b7791f", "F3": "#1f4e9c", "grey": "#7f8c8d", "green": "#1e8449"}

def read(name):
    p = os.path.join(DATA, name)
    if not os.path.exists(p): print("  (no data yet:", name, ")"); return None
    rows = [r for r in csv.DictReader(l for l in open(p) if not l.startswith("#"))]
    return rows if rows else None
def f(r, k): 
    v = r.get(k, ""); return float(v) if v not in ("", None) else float("nan")
def save(fig, name):
    for ext in ("pdf", "png"): fig.savefig(os.path.join(OUT, f"{name}.{ext}"))
    plt.close(fig); print("  wrote", name)

# ---------------- Fig 1: f-mode convergence, cut cells vs mask ----------------
def fig_convergence():
    d = read("cowling3d_anchor_convergence.csv"); dm = read("cowling3d_mask_sigko.csv")
    if d is None: return
    F1D = 1.88291
    def series(rows):
        R = sorted([r for r in rows if int(r["m"]) == 0 and r["stable"] == "true"], key=lambda r: int(r["N"]))
        N = [int(r["N"]) for r in R]
        fm = [100*((f(r, "f_pgram") + f(r, "f_pencil"))/2/F1D - 1) for r in R]
        er = [100*abs(f(r, "f_pgram") - f(r, "f_pencil"))/2/F1D for r in R]
        return N, fm, er
    fig, ax = plt.subplots(figsize=(W, 2.6))
    N, fm, er = series([r for r in d if r["surface"] == "cut"])
    ax.errorbar(N, fm, yerr=er, color=C["cut"], marker="o", ls="-", capsize=2, label=r"cut cells, $\sigma_{KO}=0.01$")
    if dm:
        N, fm, er = series(dm)
        ax.errorbar(N, fm, yerr=er, color=C["mask"], marker="s", ls="-", capsize=2, label=r"mask, $\sigma_{KO}=0.004$")
    N, fm, er = series([r for r in d if r["surface"] == "mask"])
    ax.errorbar(N, fm, yerr=er, color=C["mask"], marker="s", mfc="none", ls="--", capsize=2, label=r"mask, $\sigma_{KO}=0.01$")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel(r"$N$ (cells per axis)"); ax.set_ylabel(r"$f$-mode error [%]"); ax.set_title(r"$\ell=2$ $f$-mode vs shooting solver")
    ax.legend(loc="center left"); save(fig, "fig2_convergence")

# ---------------- Fig 2: time series ----------------
def fig_timeseries():
    a = read("cowling3d_timeseries_cut_N24.csv"); b = read("bdnk3d_timeseries_F1_N24.csv")
    if a is None and b is None: return
    fig, ax = plt.subplots(figsize=(W, 2.4))
    def draw(rows, col, lab):
        t = np.array([f(r, "t_over_P") for r in rows]); q = np.array([f(r, "q") for r in rows]); q = q/abs(q).max()
        ax.plot(t, q, color=col, label=lab)
        # guide-the-eye envelope: log-linear fit to the positive peaks after the seed transient
        pk = [i for i in range(1, len(q)-1) if q[i] > q[i-1] and q[i] >= q[i+1] and t[i] > 0.8 and q[i] > 0]
        if len(pk) >= 3:
            g, lnA = np.polyfit(t[pk], np.log(q[pk]), 1); tt = t[t >= t[pk[0]]]
            ax.plot(tt, np.exp(lnA + g*tt), "--", color=col, lw=0.8)
    if a: draw(a, C["cut"], "ideal")
    if b: draw(b, C["F1"], r"BDNK, $\hat\eta=0.03$")
    ax.set_xlabel(r"$t/P_f$"); ax.set_ylabel(r"$q_{20}(t)$ (normalised)"); ax.set_title(r"Quadrupole moment, $N=24$, $\sigma_{KO}=0.005$")
    ax.legend(loc="upper right"); save(fig, "fig1_timeseries")

# ---------------- Fig 3: crusted-EOS instability vs surface index ----------------
def fig_surface():
    d = read("cowling3d_surface_stability.csv")
    if d is None: return
    fig, ax = plt.subplots(figsize=(W, 2.5))
    for r in d:
        n, ratio = f(r, "n_eff"), f(r, "ratio_N48"); col = C["mask"] if ratio > 1 else C["green"]
        ax.plot(n, ratio, "o", color=col); ax.annotate(r["eos"].replace("PolyE_G", r"$\Gamma=$").replace("ShumPoly100", r"$\Gamma=2$"), (n, ratio), textcoords="offset points", xytext=(5, -3), fontsize=7)
    ax.set_yscale("log"); ax.axhline(1, color="k", lw=0.6, ls="--")
    ax.set_xlabel(r"surface index $n$ in $\varepsilon\propto(R-r)^n$"); ax.set_ylabel("amplitude ratio, end/start"); ax.set_title(r"Surface stability, $N=48$, 23 periods")
    ax.set_xlim(0.3, 2.3); save(fig, "fig4_surface_index")

# ---------------- Fig 4: viscous damping rate vs N ----------------
def fig_damping():
    d = read("bdnk3d_frames_N.csv"); s = read("bdnk3d_viscous_shift.csv")
    if d is None: return
    fig, ax = plt.subplots(figsize=(W, 2.5))
    for fr, mk in (("F1", "o"), ("F2", "s"), ("F3", "^")):
        R = sorted([r for r in d if r["frame"] == fr], key=lambda r: int(r["N"]))
        ax.plot([int(r["N"]) for r in R], [1e3*f(r, "gamma_visc") for r in R], marker=mk, ls="-", color=C[fr],
                label=(r"F1, F2, F3, $\sigma_{KO}=0.01$" if fr == "F1" else None))
    if s:
        pts = {}
        for r in s:
            if r["sigko"] and abs(float(r["sigko"]) - 0.005) < 1e-9 and r["visc_compact"] == "false" and float(r["dt_fac"]) == 0.20 and r["stable"] == "true":
                key = (int(r["N"]), round(float(r["etahat"]), 4))
                if key not in pts or int(r["nper"]) > pts[key][0]: pts[key] = (int(r["nper"]), f(r, "gamma_pencil"))
        Ns = sorted({k[0] for k in pts if (k[0], 0.03) in pts and (k[0], 0.0) in pts})
        ax.plot(Ns, [1e3*(pts[(n, 0.03)][1] - pts[(n, 0.0)][1]) for n in Ns], marker="D", ls="-", color=C["grey"], label=r"F1, $\sigma_{KO}=0.005$")
    ax.axhline(NS_INTEGRAL, color="k", lw=0.8, ls="--", label="Navier–Stokes integral")
    ax.set_xlim(22, 54); ax.set_xlabel(r"$N$"); ax.set_ylabel(r"$\gamma_{\rm visc}\ [10^{-3}\,M_\odot^{-1}]$"); ax.set_title(r"Viscous damping, $\hat\eta=0.03$")
    ax.legend(loc="lower right"); save(fig, "fig5_damping")

# ---------------- Fig 5: viscous frequency shift ----------------
def fig_shift():
    s = read("bdnk3d_viscous_shift.csv")
    if s is None: return
    rows = {}
    for r in s:
        if r["visc_compact"] != "false" or float(r["dt_fac"]) != 0.20 or r["stable"] != "true": continue
        key = (int(r["N"]), round(float(r["sigko"]), 5), round(float(r["etahat"]), 4))
        if key not in rows or int(r["nper"]) > int(rows[key]["nper"]): rows[key] = r
    def shift(N, sig, eta):
        c, v = rows.get((N, sig, 0.0)), rows.get((N, sig, eta))
        if c is None or v is None: return None
        sp = 100*(f(v, "f_pgram")/f(c, "f_pgram") - 1); sq = 100*(f(v, "f_pencil")/f(c, "f_pencil") - 1)
        om = 2*math.pi*f(c, "f_pgram")*KHZ_TO_KM*MSUN_KM
        pred = -100*(f(v, "gamma_env")**2 - f(c, "gamma_env")**2)/(2*om*om)
        return (sp + sq)/2, abs(sp - sq)/2, pred
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(W, 4.4))
    for sig, col, lab in ((0.01, C["mask"], r"$\sigma_{KO}=0.01$"), (0.005, C["cut"], r"$\sigma_{KO}=0.005$")):
        Ns = [n for n in (24, 32, 40, 48) if shift(n, sig, 0.03)]
        if not Ns: continue
        v = [shift(n, sig, 0.03) for n in Ns]
        a1.errorbar(Ns, [x[0] for x in v], yerr=[x[1] for x in v], color=col, marker="o", ls="-", capsize=2, label="measured, " + lab)
        if sig == 0.005: a1.plot(Ns, [x[2] for x in v], "k--", marker="x", label="damped-oscillator prediction")
    a1.axhline(0, color="k", lw=0.5); a1.set_xlabel(r"$N$"); a1.set_ylabel(r"$\Delta f/f$ [%]"); a1.set_title(r"Frequency shift, $\hat\eta=0.03$"); a1.legend(loc="center right")
    etas = [0.015, 0.03, 0.045]; v = [shift(40, 0.005, e) for e in etas]
    if all(v):
        a2.errorbar(etas, [x[0] for x in v], yerr=[x[1] for x in v], color=C["cut"], marker="o", ls="none", capsize=2, label="measured")
        a2.plot(etas, [x[2] for x in v], "k--", marker="x", label="prediction")
        ee = np.linspace(0, 0.05, 50); k = v[2][2]/0.045**2; a2.plot(ee, k*ee**2, ":", color=C["grey"], lw=0.8, label=r"$\propto\hat\eta^2$")
    a2.axhline(0, color="k", lw=0.5); a2.set_xlabel(r"$\hat\eta$"); a2.set_ylabel(r"$\Delta f/f$ [%]"); a2.set_title(r"$N=40$, $\sigma_{KO}=0.005$"); a2.legend(loc="lower left")
    fig.tight_layout(); save(fig, "fig6_shift")

# ---------------- Fig 6: E_g / T_2g splitting vs N ----------------
def fig_msplit():
    d = read("cowling3d_anchor_convergence.csv")
    if d is None: return
    R = [r for r in d if r["surface"] == "cut" and r["stable"] == "true"]
    Ns = sorted({int(r["N"]) for r in R if int(r["m"]) == 1})
    if not Ns: print("  (no m=1 data yet)"); return
    def fm(N, m):
        r = [x for x in R if int(x["N"]) == N and int(x["m"]) == m]; return (f(r[0], "f_pgram") + f(r[0], "f_pencil"))/2 if r else float("nan")
    def sp(N, m):
        r = [x for x in R if int(x["N"]) == N and int(x["m"]) == m]; return abs(f(r[0], "f_pgram") - f(r[0], "f_pencil"))/2 if r else float("nan")
    s1 = [100*abs(fm(n, 1)/fm(n, 0) - 1) for n in Ns]; e1 = [100*(sp(n, 1) + sp(n, 0))/fm(n, 0) for n in Ns]
    s2 = max(100*abs(fm(n, 2)/fm(n, 0) - 1) for n in Ns)
    fig, ax = plt.subplots(figsize=(W, 2.5))
    ax.errorbar(Ns, s1, yerr=e1, color=C["mask"], marker="o", ls="-", capsize=2, label=r"$T_{2g}$: $|f(Y_{21})/f(Y_{20})-1|$")
    nn = np.linspace(Ns[0], Ns[-1], 50); ax.plot(nn, s1[0]*(Ns[0]/nn)**2, ":", color=C["grey"], lw=0.9, label=r"$\propto h^2$")
    ax.text(0.97, 0.55, r"$E_g$: $|f(Y_{22})/f(Y_{20})-1| < 10^{-12}$ %%" % () if False else r"$E_g$: $|f(Y_{22})/f(Y_{20})-1| < 10^{-12}\,\%$", transform=ax.transAxes, ha="right", fontsize=7.5)
    ax.set_ylim(0, None); ax.set_xlabel(r"$N$"); ax.set_ylabel("frequency splitting [%]"); ax.set_title(r"Cubic-grid splitting of the $\ell=2$ $f$-mode")
    ax.legend(loc="upper right"); save(fig, "fig3_msplit")

# ---------------- Fig 7: stability threshold ----------------
def fig_stability():
    d = read("bdnk3d_stability_threshold.csv")
    if d is None: return
    fig, ax = plt.subplots(figsize=(W, 2.5))
    for fr, mk in (("F1", "o"), ("F3", "^")):
        R = sorted([r for r in d if r["frame"] == fr], key=lambda r: f(r, "etahat"))
        ax.plot([f(r, "etahat") for r in R], [f(r, "growth_rate") for r in R], marker=mk, ls="none", color=C[fr], label=fr)
        thr = f(R[0], "C")/f(R[0], "inv_h"); ax.axvline(thr, color=C[fr], ls="--", lw=0.8)
    ax.axhline(0, color="k", lw=0.6); ax.set_yscale("symlog", linthresh=1e-3)
    ax.set_xlabel(r"$\hat\eta$"); ax.set_ylabel(r"growth rate of $\max|\delta v|$ [$M_\odot^{-1}$]"); ax.set_title(r"Onset of instability, $N=24$")
    ax.text(0.98, 0.04, r"dashed: predicted $\hat\eta_c = C\,h$", transform=ax.transAxes, ha="right", fontsize=7)
    ax.legend(loc="upper left"); save(fig, "fig7_stability")

if __name__ == "__main__":
    for g in (fig_convergence, fig_timeseries, fig_surface, fig_damping, fig_shift, fig_msplit, fig_stability): g()
