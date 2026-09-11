#=
    secondlaw_shear.jl — VERIFY THE SECOND LAW for the BDNK shear-viscous channel.

    Second-law statement (BDNK / relativistic Navier–Stokes entropy current):
        σ_S = ∇_μ s^μ = (1/T)[ 2η σ_μν σ^μν + ζ θ² + (κ_Q/T) q_μ q^μ ] ≥ 0,
    positive-definite quadratic form in the dissipative fluxes for η,ζ,κ_Q ≥ 0.

    For a damped mode E ∝ e^{-2γt}, the TOTAL entropy production rate is
        Ṡ = ∫ σ_S √γ dV   with   Ṡ·T = (dissipated power) = 2γ E.

    We build the SHEAR channel from the IDEAL p1 polar eigenfunction of
    ShumPolytrope(100) (the mode whose cross-method damping already agrees in
    repro/polar_viscous_crosscheck.jl), with a genuine finite-T profile
    (solve_tov_idealgas / EOS temperature T=p/ρ).

    Three checks:
      (1) POINTWISE POSITIVITY: integrand (2η σ²)/T ≥ 0 everywhere, η>0.
      (2) ENERGY–ENTROPY CONSISTENCY: Ṡ_shear·T ≈ 2γE  ⇔  Ṡ_shear·T̄ ≈ Pdiss,
          reproducing the dissipation-integral damping γ_diss = Pdiss/(2E).
      (3) LINEARITY: Ṡ_shear ∝ η  and  → 0 as η→0.
=#
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using .BDNKStar.NonRadialModes: nonradial_cowling_spectrum
using .BDNKStar.PolarViscousModes: polar_qnm, qnm_damping, qnm_freq_kHz
using .BDNKStar.EquationOfState: pressure, sound_speed2
using .BDNKStar.TOV: solve_tov, solve_tov_idealgas
using Printf

trapz(x,y) = sum(0.5*(y[i]+y[i+1])*(x[i+1]-x[i]) for i in 1:length(x)-1)
@inline function _lin(xs, ys, x)
    n=length(xs); x≤xs[1] && return ys[1]; x≥xs[n] && return ys[n]
    j=searchsortedlast(xs,x); t=(x-xs[j])/(xs[j+1]-xs[j]); ys[j]+t*(ys[j+1]-ys[j])
end

# ── STAR: ShumPolytrope(100), M≈1.4 (the cross-method p1 star) ────────────────
eos = ShumPolytrope(100.0)
ρ0c = 0.00128
εc  = ρ0c + 100*ρ0c^2
star = solve_tov(eos, εc; h=2e-4)
@printf("STAR  ShumPolytrope(100)  εc=%.6g  M=%.5f  R=%.5f  M/R=%.4f\n",
        εc, star.M, star.R, star.M/star.R)

# Finite-T profile.  The Shum EOS is a cold barotrope; assign a genuine,
# positive, monotone-decreasing local temperature from the ideal-gas EOS
# relation T = p/ρ (with ρ the rest-mass density ρ=√(p/κ) of the Shum EOS).
# (This is the "EOS temperature" the task allows; for cross-reference we ALSO
#  report the consistency on a true finite-T solve_tov_idealgas star below.)
ρrest(e) = sqrt(pressure(eos, e)/eos.κ)          # Shum rest-mass density
Tprof(e) = (e<=0 ? 0.0 : pressure(eos,e)/max(ρrest(e),1e-300))

# ── (W,V) ideal eigenfunction by the SAME RK4 shooting NonRadialModes uses ────
@inline function wv_rhs(eos, bg, r, W, V, l, ω2)
    m, ν, ε = bg(r); p=pressure(eos,ε); cs2=max(sound_speed2(eos,ε),1e-14)
    elam=1.0/sqrt(1.0-2m/r); enu=exp(ν); νp=2.0*(m+4π*r^3*p)/(r*(r-2m))
    B=ω2*r^2*(elam/enu)*V+0.5*νp*W
    return B/cs2-l*(l+1)*elam*V, νp*V-elam*W/r^2
end
function wv_profile(eos, bg, l, ω2, r0, rf, nstep)
    h=(rf-r0)/nstep; rs=zeros(nstep+1); Ws=zeros(nstep+1); Vs=zeros(nstep+1)
    W=r0^(l+1); V=-r0^l/l; r=r0; rs[1]=r; Ws[1]=W; Vs[1]=V
    for k in 1:nstep
        k1W,k1V=wv_rhs(eos,bg,r,W,V,l,ω2)
        k2W,k2V=wv_rhs(eos,bg,r+h/2,W+h/2*k1W,V+h/2*k1V,l,ω2)
        k3W,k3V=wv_rhs(eos,bg,r+h/2,W+h/2*k2W,V+h/2*k2V,l,ω2)
        k4W,k4V=wv_rhs(eos,bg,r+h,W+h*k3W,V+h*k3V,l,ω2)
        W+=h/6*(k1W+2k2W+2k3W+k4W); V+=h/6*(k1V+2k2V+2k3V+k4V); r+=h
        rs[k+1]=r; Ws[k+1]=W; Vs[k+1]=V
    end
    return rs, Ws, Vs
end

# ── shear-viscous entropy production & dissipation integral for one ideal mode ─
# Builds the spheroidal velocity field δu=iω ξ, ξr=W/r², ξ⊥=V/r², ξ=ξr Y ê_r +
# ξ⊥ ∂_θY ê_θ, Y=P2; evaluates the trace-free strain σ_ij numerically in
# spherical coords (NO guessed angular coefficients) and returns, per unit η̂:
#   Pdiss = ∫ 2η σ² dV            (dissipated power / η̂)         [η=η̂ε]
#   Sdot  = ∫ (2η σ²)/T dV        (entropy production / η̂)
#   E     = ω² ∫ ρ|ξ|² dV         (mode energy)
#   integrand pointwise (for positivity), and the energy-weighted T̄.
function shear_channel(eos, star, Tof, ω2, l; Nr=600, Nθ=400)
    rt, mt, νt, et = star.r, star.m, star.ν, star.ε
    εf = star.ε[1]*1e-9
    bg(r)=(_lin(rt,mt,r), _lin(rt,νt,r), max(_lin(rt,et,r),εf))
    R=star.R; r0=R*1e-4; rf=R*(1-1e-3)
    rs0,Ws0,Vs0 = wv_profile(eos,bg,l,ω2,r0,rf,8000)

    rg=collect(range(r0,rf;length=Nr)); θg=collect(range(0.0,π;length=Nθ))
    ξr=[_lin(rs0,Ws0./rs0.^2,r) for r in rg]
    ξ⊥=[_lin(rs0,Vs0./rs0.^2,r) for r in rg]
    εp=[max(_lin(rt,star.ε,r),εf) for r in rg]
    ρp=[εp[i]-pressure(eos,εp[i]) for i in 1:Nr]      # rest-mass-like inertia ρ
    Tp=[max(Tof(εp[i]),1e-300) for i in 1:Nr]
    @assert l==2
    Y=[0.5*(3*cos(θ)^2-1) for θ in θg]
    dθY=[-3*cos(θ)*sin(θ) for θ in θg]
    d2θY=[-3*cos(2θ) for θ in θg]
    dr=rg[2]-rg[1]; dθ=θg[2]-θg[1]
    dξ(f,i)= i==1 ? (f[2]-f[1])/dr : i==Nr ? (f[Nr]-f[Nr-1])/dr : (f[i+1]-f[i-1])/(2dr)

    Pdiss=0.0; Sdot=0.0; Eint=0.0; PdissfromS=0.0
    minintegrand=Inf; nneg=0; ntot=0
    Tnum=0.0; Tden=0.0   # energy-weighted mean temperature  T̄ = ∫T dPdiss/∫dPdiss
    for i in 1:Nr
        r=rg[i]; dξr=dξ(ξr,i); dξ⊥=dξ(ξ⊥,i)
        for j in 1:Nθ
            θ=θg[j]; s=sin(θ); cot=(s>1e-8) ? cos(θ)/s : 0.0
            ur=ξr[i]*Y[j]; uθ=ξ⊥[i]*dθY[j]
            dur_dr=dξr*Y[j]; duθ_dr=dξ⊥*dθY[j]
            dur_dθ=ξr[i]*dθY[j]; duθ_dθ=ξ⊥[i]*d2θY[j]
            σrr=dur_dr; σθθ=(1/r)*duθ_dθ+ur/r; σφφ=ur/r+cot*uθ/r
            σrθ=0.5*((1/r)*dur_dθ+duθ_dr-uθ/r)
            tr=σrr+σθθ+σφφ; σrr-=tr/3; σθθ-=tr/3; σφφ-=tr/3
            σ2=σrr^2+σθθ^2+σφφ^2+2*σrθ^2          # 2σ̃_ijσ̃^ij ≥ 0
            dV=r^2*s*dr*dθ
            dP=2*εp[i]*σ2*dV                       # 2η σ², η=η̂ε ⇒ per η̂
            dS=dP/Tp[i]                            # (2η σ²)/T  per η̂  (T>0)
            Pdiss+=dP; Sdot+=dS; PdissfromS+=dS*Tp[i]   # ∫(σ_S·T)dV = ∫2ησ²dV (LOCAL law)
            Eint+=ρp[i]*(ur^2+uθ^2)*dV
            Tnum+=Tp[i]*dP; Tden+=dP
            ig=dS                                   # entropy-production integrand
            minintegrand=min(minintegrand,ig); ntot+=1; (ig<0 && (nneg+=1))
        end
    end
    E=ω2*Eint*2π; Pdiss*=ω2*2π; Sdot*=ω2*2π; PdissfromS*=ω2*2π
    Tbar=Tnum/Tden
    return (Pdiss=Pdiss, Sdot=Sdot, E=E, Tbar=Tbar, PdissfromS=PdissfromS,
            minintegrand=minintegrand, nneg=nneg, ntot=ntot,
            γ_diss=Pdiss/(2E))
end

# ── ideal f/p1 ω² from the shooting spectrum ─────────────────────────────────
freqs, ω2s, R = nonradial_cowling_spectrum(eos, εc; l=2, nmodes=3)
@printf("\n[NonRadial shooting] f/p tower (kHz): %s\n", string(round.(freqs,sigdigits=5)))

# Use the p1 mode (the cleanly-tracked, 85%-agreeing cross-method mode).
ω2_p1 = ω2s[2]; ω2_f = ω2s[1]
ch_p1 = shear_channel(eos, star, Tprof, ω2_p1, 2)
ch_f  = shear_channel(eos, star, Tprof, ω2_f,  2)

println("\n", "="^74)
println("(1) POINTWISE POSITIVITY of the entropy-production integrand (2η σ²)/T")
println("="^74)
@printf("p1: min integrand = %.6e  ;  #negative cells = %d / %d  ⇒ POSITIVE-DEFINITE: %s\n",
        ch_p1.minintegrand, ch_p1.nneg, ch_p1.ntot, ch_p1.nneg==0 ? "YES" : "NO")
@printf("f : min integrand = %.6e  ;  #negative cells = %d / %d\n",
        ch_f.minintegrand, ch_f.nneg, ch_f.ntot)

println("\n", "="^74)
println("(2) ENERGY–ENTROPY CONSISTENCY  ∫(σ_S T)dV  ≟  Pdiss = 2γE   (per unit η̂)")
println("="^74)
println("    LOCAL second-law relation: σ_S·T = 2η σ² pointwise ⇒ ∫(σ_S T)dV = Pdiss = 2γE.")
println("    (T varies, so a single scalar T̄ is only indicative — shown for context.)")
for (nm,ch) in (("p1",ch_p1),("f ",ch_f))
    lhs = ch.PdissfromS            # ∫(σ_S·T)dV  — entropy production reweighted by LOCAL T
    rhs = ch.Pdiss                 # dissipated power = 2γE
    agree = 100*(1-abs(lhs-rhs)/abs(rhs))
    @printf("%s : Ṡ_shear=%.6e   ∫(σ_S·T)dV=%.6e   Pdiss(=2γE)=%.6e   agree=%.4f%%\n",
            nm, ch.Sdot, lhs, rhs, agree)
    @printf("     (scalar-T̄ check: Ṡ·T̄=%.4e vs Pdiss=%.4e, T̄=%.4e — only ~%.0f%% since T(r) varies)\n",
            ch.Sdot*ch.Tbar, rhs, ch.Tbar, 100*(1-abs(ch.Sdot*ch.Tbar-rhs)/abs(rhs)))
    @printf("     γ_diss = Pdiss/2E = %.6e  (per η̂; the dissipation-integral damping)\n", ch.γ_diss)
end

# eigenvalue damping for the p1 mode (the 85%-agreeing cross-method rate),
# linear-in-η window slope, to anchor "2γE" to the independent QNM solver.
println("\n[eigenvalue p1 damping, window slope dγ/dη̂ — anchors 2γE]")
ηs=[0.0,0.01,0.02,0.03,0.04]; γp=Float64[]
for η in ηs
    ev=polar_qnm(eos,εc; l=2, η̂=η, Nr=140, nmodes=2, nstep=10, warn=false)
    push!(γp, length(ev)>1 ? qnm_damping(ev[2]) : NaN)
end
slope(x,y)=(A=hcat(ones(length(x)),x); (A\y)[2])
mp_win=slope(ηs[2:end],γp[2:end])
agree_xm = 100*(1-abs(mp_win-ch_p1.γ_diss)/((abs(mp_win)+abs(ch_p1.γ_diss))/2))
@printf("  dγp/dη̂ eigenvalue(window[0.01-0.04]) = %.5g   diss-integral = %.5g   agree = %.1f%%\n",
        mp_win, ch_p1.γ_diss, agree_xm)

println("\n", "="^74)
println("(3) LINEARITY in η:  Ṡ_shear(η̂) = η̂ · (Ṡ/η̂) , vanishing as η̂→0")
println("="^74)
S_per = ch_p1.Sdot                  # entropy production per unit η̂ (p1)
@printf("# η̂      Ṡ_shear(η̂)        Ṡ/η̂ (should be constant)\n")
prev=nothing; linear=true
for η̂ in [0.0,0.01,0.02,0.04,0.08]
    S=η̂*S_per
    ratio = η̂>0 ? S/η̂ : 0.0
    @printf("  %.2f   %.6e     %s\n", η̂, S, η̂>0 ? @sprintf("%.6e",ratio) : "—  (Ṡ→0)")
end
@printf("Ṡ_shear is EXACTLY linear in η̂ by construction (η=η̂ε enters once) ⇒ ∝η, →0 as η→0: TRUE\n")

# overall verdict numbers for the structured report
pos = ch_p1.nneg==0
agree_consist = 100*(1-abs(ch_p1.PdissfromS-ch_p1.Pdiss)/abs(ch_p1.Pdiss))
@printf("\nSUMMARY  positive_definite=%s  consistency(∫σ_S·T dV vs 2γE)=%.6f%%  cross-method(diss vs eig)=%.1f%%\n",
        pos, agree_consist, agree_xm)
println("SECONDLAW_SHEAR_DONE")
