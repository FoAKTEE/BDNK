#!/usr/bin/env python3
"""Analysis of the unseeded B1-protocol run (repro/data/dgball3d_hkt_series_B1_unseeded_nt6_p3.csv):
the analogue of Figs. 15-16 of Hebert, Kidder & Teukolsky 2018 — err[D], rho_c(t), M_b(t) and the
Hanning-windowed spectrum of rho_c against the linear radial Cowling modes. Works on a partial file."""
import numpy as np, os, sys, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
here=os.path.dirname(os.path.abspath(__file__)); f=os.path.join(here,"data","dgball3d_hkt_series_B1_unseeded_nt6_p3.csv")
d=np.genfromtxt(f,delimiter=",",names=True)
# analyse only the pre-blow-up record: stop where err[D] first exceeds 0.1 (or at --tmax)
tmax=float(sys.argv[1]) if len(sys.argv)>1 else np.inf
ok=(d["errD"]<0.1)&np.isfinite(d["errD"])&(d["t"]<=tmax); last=np.argmin(ok) if (~ok).any() else len(ok)
d=d[:last] if last>0 else d
t=d["t"]; r=d["rhoc"]/d["rhoc"][0]-1; e=d["errD"]; Mb=d["Mb"]/d["Mb"][0]-1
Msun_s=4.925490947e-6; F_lin=[2.6861,4.5497,6.3414,8.1081,9.8633]
print(f"record: t_end={t[-1]:.0f}  errD end={e[-1]:.2e}  max after t>200 {e[t>200].max():.2e}  rhoc drift end={r[-1]:+.2e}  Mb drift={Mb[-1]:+.2e}  max|v_atm|={d['vatm'].max():.2f}")
for a in range(0,int(t[-1])+1,200):
    m=(t>=a)&(t<a+200)
    if m.sum()>10: print(f"  t in [{a},{a+200}): errD mean {e[m].mean():.2e}  rhoc mean {r[m].mean():+.2e} amp {np.abs(r[m]-r[m].mean()).max():.2e}  Mb {Mb[m][-1]:+.2e}")
def spectrum(t,q,fmin=0.5,fmax=12,n=4000):
    q=q-q.mean(); w=0.5-0.5*np.cos(2*np.pi*np.arange(len(q))/(len(q)-1)); q=q*w
    f=np.linspace(fmin,fmax,n); om=2*np.pi*f*1e3*Msun_s; return f,np.array([abs(np.sum(q*np.exp(-1j*o*t)))**2 for o in om])
m=t<=min(4000,t[-1]); fq,P=spectrum(t[m],r[m]); pk=[(fq[i],P[i]/P.max()) for i in range(1,len(P)-1) if P[i]>P[i-1] and P[i]>=P[i+1] and P[i]>1e-3*P.max()]
print("spectrum peaks (kHz, rel. power): "+", ".join(f"{a:.3f}({b:.2f})" for a,b in pk[:8])+"   linear F=2.686 H1=4.550 H2=6.341 H3=8.108 (paper: F,H1 resolved at B1)")
fig,ax=plt.subplots(1,3,figsize=(15,4.2))
ax[0].semilogy(t,e,lw=0.8); ax[0].set_xlabel(r"$t\ [M_\odot]$"); ax[0].set_ylabel(r"err$[\tilde D]$"); ax[0].set_title("Density error (paper B1: ~6e-4 at 4000)")
ax[1].plot(t,r*1e3,lw=0.6); ax[1].set_xlabel(r"$t\ [M_\odot]$"); ax[1].set_ylabel(r"$(\rho_c-\rho_{c,0})/\rho_{c,0}\ [10^{-3}]$"); ax[1].set_title("Central density (paper B1: settles at +0.25)")
ax[2].semilogy(fq,P/P.max(),lw=0.8); [ax[2].axvline(x,color="k",ls=":",lw=0.8) for x in F_lin]; ax[2].set_xlim(0.5,12); ax[2].set_ylim(1e-8,2)
ax[2].set_xlabel("f [kHz]"); ax[2].set_ylabel("power (arb.)"); ax[2].set_title(rf"Spectrum of $\rho_c$, $t<{t[-1]:.0f}$ (dotted: linear F, H$_1$, ...)")
plt.tight_layout(); out=os.path.join(here,"..","paper","figs","dgball3d_b1.png"); plt.savefig(out,dpi=150); print("figure:",out)
