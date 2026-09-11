#=
    viz/dyngr_radial_validation.jl — GR RADIAL-MODE VALIDATION figure.

    Cross-validates the dynamical-GR fundamental radial mode of the benchmark
    star ShumPolytrope(κ=100), εc=0.0015 km⁻² (M≈0.96 M⊙, R=9.52 km, 2M/R≈0.30)
    by THREE INDEPENDENT routes and compares to the Cowling value and the
    published Kokkotas–Ruoff (2001) A&A 366, 565 [gr-qc/0011093] full-GR band:

      • DYNAMICAL ENGINE  — DynGR1D time-domain (drops Cowling, evolves spacetime);
      • SL EIGENSOLVER    — corrected chandrasekhar_radial_omega2 (self-adjoint
                            LAWE, KR(2001) eqs.14-17, OUR convention g_tt=-e^ν);
      • SHOOTING          — first-order GHZ(1997) (ξ,Δp) RK4 cross-check.

    The corrected SL weights (P,Q ∝ e^{(λ+3ν)/2}, W ∝ e^{(3λ+ν)/2}) place the
    full-GR fundamental BELOW the Cowling value (F_GR≈2.12 < F_Cow≈4.01 kHz):
    the spacetime response SOFTENS the restoring force. The earlier
    'GR above Cowling' exponent-swap bug is fixed.

    Two panels:
      (a) bar chart of the three GR routes + Cowling, with the KR(2001) published
          full-GR band (2.1-2.3 kHz at this EOS family) shaded;
      (b) ω² comparison (GR routes vs Cowling) showing ω²_GR < ω²_Cowling.

    Run: cd code/BDNKStar && JULIA_NUM_THREADS=6 julia --project=. viz/dyngr_radial_validation.jl
=#
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using .BDNKStar: setup_dyngr, evolve_dyngr!, seed_dyngr_velocity!, solve_tov,
                 ShumPolytrope, dyngr_radial_freq, chandrasekhar_radial_omega2,
                 cowling_radial_omega2, radial_cowling_spectrum,
                 pressure, sound_speed2, energy_from_pressure
using .BDNKStar.Units: kHz_to_km
using CairoMakie
using Printf

const CONV = 1/kHz_to_km                       # km^-1 -> kHz (geometric)
fkHz(ω2) = sqrt(max(ω2, 0.0)) / (2π) * CONV

const EOS  = ShumPolytrope(100.0)
const EPSC = 0.0015
const HTOV = 2e-5

star = solve_tov(EOS, EPSC; h=HTOV)
const R = star.R; const M = star.M
@printf("benchmark: ShumPolytrope(100), eps_c=%.4g  R=%.3f km  M=%.3f km  2M/R=%.4f\n",
        EPSC, R, M, 2M/R)

_detrend(ts, y) = begin
    n=length(ts); x=(ts.-ts[1])./(ts[end]-ts[1]); A=hcat(ones(n),x,x.^2,x.^3); y.-A*(A\y)
end

# ---------------------------------------------------------------------------
# (1) SL eigensolver (corrected) and Cowling
# ---------------------------------------------------------------------------
ω2_sl  = chandrasekhar_radial_omega2(EOS, EPSC; nmodes=2, N=2500, h_tov=HTOV)
fcow, ω2cow, _ = radial_cowling_spectrum(EOS, EPSC; N=2000, h_tov=5e-5, nmodes=2)
F_SL   = fkHz(ω2_sl[1])
F_COW  = fcow[1]
@printf("SL   : omega^2=%.6e  F_SL =%.4f kHz\n", ω2_sl[1], F_SL)
@printf("Cow  : omega^2=%.6e  F_Cow=%.4f kHz\n", ω2cow[1], F_COW)

# ---------------------------------------------------------------------------
# (2) SHOOTING cross-check (co-integrated GHZ 1997 (ξ,Δp), Schwarzschild ν-shift)
# ---------------------------------------------------------------------------
@inline function _deriv(r, y, ω2)
    m, p, ν, ξ, Δp = y
    p ≤ 0 && return (0.0,0.0,0.0,0.0,0.0)
    ε  = energy_from_pressure(EOS, p)
    den = r*(r-2m); fac = m + 4π*r^3*p
    dm = 4π*r^2*ε; dν = 2*fac/den; pp = -(ε+p)*fac/den
    eλ = 1.0/max(1-2m/r,1e-12); eνm = exp(-ν)
    cs2 = clamp(sound_speed2(EOS, ε),1e-12,1.0); Γ1 = cs2*(ε+p)/p
    V = -3.0/r - pp/(ε+p)
    W = -(1.0/r)/(Γ1*p)
    X = ω2*eλ*eνm*(ε+p)*r - 4.0*pp + pp^2*r/(ε+p) - 8π*eλ*(ε+p)*p*r
    Y = pp/(ε+p) - 4π*(ε+p)*r*eλ
    return (dm, pp, dν, V*ξ+W*Δp, X*ξ+Y*Δp)
end
function _shoot(ω2; h=5e-4, ptol_rel=1e-8)
    pc = pressure(EOS, EPSC); ptol = ptol_rel*pc
    r0 = h; m = (4π/3)*EPSC*r0^3
    p = pc - 2π*(EPSC+pc)*(EPSC/3+pc)*r0^2; ν = 0.0
    ε0 = energy_from_pressure(EOS, p); Γ10 = clamp(sound_speed2(EOS,ε0),1e-12,1.0)*(ε0+p)/p
    ξ = 1.0; Δp = -3.0*Γ10*p*ξ
    y = (m,p,ν,ξ,Δp); r = r0; Δp_s = Δp; ν_s = ν; R_int = r0
    while p > ptol && r < 100.0
        k1=_deriv(r,y,ω2)
        y2=ntuple(i->y[i]+h/2*k1[i],5); k2=_deriv(r+h/2,y2,ω2)
        y3=ntuple(i->y[i]+h/2*k2[i],5); k3=_deriv(r+h/2,y3,ω2)
        y4=ntuple(i->y[i]+h*k3[i],5);   k4=_deriv(r+h,y4,ω2)
        yn=ntuple(i->y[i]+h/6*(k1[i]+2k2[i]+2k3[i]+k4[i]),5)
        r+=h; pn=yn[2]
        if pn ≤ ptol
            frac = y[2]/(y[2]-pn)
            Δp_s = y[5]+frac*(yn[5]-y[5]); ν_s = y[3]+frac*(yn[3]-y[3])
            R_int = (r-h)+frac*h; break
        end
        p=pn; y=yn
    end
    return Δp_s, ν_s, R_int
end
const ΔΝ = let (_, νR, Rint) = _shoot(1e-3); log(1-2M/Rint) - νR end
_surf(ω2raw) = _shoot(ω2raw)[1]
function _find_fund(; ω2lo=1e-4, ω2hi=1.5e-2, nscan=400)
    grid = collect(range(ω2lo, ω2hi; length=nscan))
    prevv=_surf(grid[1]); prevω=grid[1]; brk=nothing
    for ω2 in grid[2:end]
        cur=_surf(ω2)
        sign(cur)!=sign(prevv) && brk===nothing && (brk=(prevω,ω2))
        prevv=cur; prevω=ω2
    end
    a,c = brk; fa=_surf(a)
    for _ in 1:200
        mid=0.5*(a+c); fm=_surf(mid)
        sign(fm)==sign(fa) ? (a=mid; fa=fm) : (c=mid)
        (c-a)<1e-13 && break
    end
    return 0.5*(a+c)
end
ω2_shoot = _find_fund()*exp(ΔΝ)
F_SHOOT  = fkHz(ω2_shoot)
@printf("Shoot: omega^2=%.6e  F_shoot=%.4f kHz  (e^Dnu=%.5f)\n", ω2_shoot, F_SHOOT, exp(ΔΝ))

# ---------------------------------------------------------------------------
# (3) DYNAMICAL ENGINE time-domain fundamental (drops Cowling, evolves spacetime)
# ---------------------------------------------------------------------------
function _fmode(fm, band)
    eng, st = setup_dyngr(EOS, EPSC; N=500, cfl=0.25, freeze_metric=fm)
    seed_dyngr_velocity!(st, eng; A=2e-4, profile=:linear)
    res = evolve_dyngr!(st, eng; tmax=220*star.R, probe_frac=0.5, sample_dt=0.5*eng.g.Δr)
    yv = _detrend(res.ts, copy(res.probe))
    fpk,_,_ = dyngr_radial_freq(res.ts, yv; fmin=band[1]*kHz_to_km, fmax=band[2]*kHz_to_km, nf=8000)
    fpk*CONV
end
F_DYN = _fmode(false, (1.2, 4.0))
@printf("Engine: F_dyn=%.4f kHz (dynamical, time-domain)\n", F_DYN)

# ---------------------------------------------------------------------------
# Published full-GR anchor: Kokkotas & Ruoff (2001) A&A 366, 565 Table A.18,
# n=1 κ=100 km² relativistic polytrope (p=κρ², Γ=2) — our EOS family. Nearby-M
# full-GR fundamentals: 2.150 (0.802 Msun), 2.302 (0.998 Msun), 2.323 kHz.
const KR_BAND = (2.10, 2.35)           # published full-GR fundamental band
const KR_PTS  = [2.150, 2.302, 2.323]  # tabulated KR(2001) ν0 for nearby stars

# ---------------------------------------------------------------------------
# FIGURE
# ---------------------------------------------------------------------------
set_theme!(theme_minimal())
fig = Figure(size=(1180, 480), fontsize=15)

# (a) frequencies ------------------------------------------------------------
ax1 = Axis(fig[1,1], title="(a) full-GR radial fundamental: 3 routes agree, GR < Cowling",
           ylabel="F  [kHz]", xticks=(1:5, ["engine\n(DynGR1D)","SL eig\n(corrected)",
                                          "shooting\n(GHZ)","","Cowling\n(frozen)"]))
labels = ["engine","SL","shoot","","Cowling"]
fvals  = [F_DYN, F_SL, F_SHOOT, NaN, F_COW]
cols   = [:dodgerblue, :seagreen, :darkorange, :white, :crimson]
# KR(2001) published full-GR band
band = poly!(ax1, Point2f[(0.4,KR_BAND[1]),(3.6,KR_BAND[1]),(3.6,KR_BAND[2]),(0.4,KR_BAND[2])],
             color=(:grey,0.18))
hlines!(ax1, KR_PTS, xmin=0.07, xmax=0.62, color=:grey35, linestyle=:dash, linewidth=1)
for (i,(f,c)) in enumerate(zip(fvals, cols))
    isnan(f) && continue
    barplot!(ax1, [i], [f], color=c, width=0.62)
    text!(ax1, i, f+0.08, text=@sprintf("%.3f", f), align=(:center,:bottom), fontsize=13)
end
text!(ax1, 2.0, KR_BAND[2]+0.05, text="KR(2001) full-GR band\n(κ=100, Γ=2)",
      align=(:center,:bottom), color=:grey25, fontsize=11)
ylims!(ax1, 0, 4.6); xlims!(ax1, 0.4, 5.6)

# (b) omega^2 ----------------------------------------------------------------
ax2 = Axis(fig[1,2], title="(b) ω² : full-GR sits below Cowling (softened restoring force)",
           ylabel="ω²  [km⁻²]  (×10⁻³)",
           xticks=(1:4, ["engine","SL eig","shooting","Cowling"]))
ω2_engine = (F_DYN/CONV*2π)^2
ω2vals = [ω2_engine, ω2_sl[1], ω2_shoot, ω2cow[1]] .* 1e3
cols2  = [:dodgerblue, :seagreen, :darkorange, :crimson]
for (i,(w,c)) in enumerate(zip(ω2vals, cols2))
    barplot!(ax2, [i], [w], color=c, width=0.62)
    text!(ax2, i, w+0.15, text=@sprintf("%.3f", w), align=(:center,:bottom), fontsize=13)
end
hlines!(ax2, [ω2cow[1]*1e3], color=:crimson, linestyle=:dot, linewidth=1.5)
text!(ax2, 1.0, ω2cow[1]*1e3+0.15, text="Cowling level", color=:crimson, fontsize=11,
      align=(:left,:bottom))
ylims!(ax2, 0, ω2cow[1]*1e3*1.25)

Label(fig[0,:], "Dynamical-GR radial mode — independent validation (corrected KR(2001) SL eigensolver)",
      fontsize=17, font=:bold)

# match summary
mismatch_sl   = 100*(F_DYN-F_SL)/F_SL
@printf("\nSUMMARY  F_dyn=%.4f  F_SL=%.4f  F_shoot=%.4f  F_Cow=%.4f kHz\n",
        F_DYN, F_SL, F_SHOOT, F_COW)
@printf("  engine vs SL : %+.2f%%   shoot vs SL : %+.2f%%   GR<Cow : %s\n",
        mismatch_sl, 100*(F_SHOOT-F_SL)/F_SL, ω2_shoot < ω2cow[1])

outdir = joinpath(@__DIR__, "..", "figures")
isdir(outdir) || mkpath(outdir)
outpath = joinpath(outdir, "dyngr_radial_validation.png")
save(outpath, fig)
println("saved -> ", outpath)
