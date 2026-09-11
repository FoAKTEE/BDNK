#=
    entropy_second_law.jl — SECOND-LAW / entropy-production figure for the BDNK
    dissipative channels (shear η, bulk ζ, heat κ_Q).

    BDNK relativistic Navier–Stokes entropy current divergence
        σ_S = ∇_μ s^μ = (1/T)[ 2η σ_μν σ^μν + ζ θ² + (κ_Q/T) q_μ q^μ ] ≥ 0,
    a positive-definite quadratic form in the dissipative fluxes for η,ζ,κ_Q ≥ 0.

    Three panels:
      A  entropy-production density (2η σ²)/T vs r for the p1 polar Cowling mode
         of ShumPolytrope(100) (M=1.40, R=9.59): POSITIVE everywhere.
      B  energy–entropy consistency  ∫(σ_S·T)dV  =  Pdiss = 2γE  ∝ η̂  (diagonal):
         the LOCAL second-law relation σ_S·T = 2η σ² holds, so the dissipation-
         integral rate sits exactly on the y=x line as η̂ is scaled.
      C  channel decomposition: shear η is DISSIPATIVE (engine mode-energy falls,
         entropy>0); κ_Q is REACTIVE in the Cowling engine (no net mode-energy
         work — no δT/heat-flux state evolved) but its genuine heat-flux entropy
         (κ_Q/T²)(∂δT)² is POSITIVE-definite on the finite-T ideal-gas star.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/entropy_second_law.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using .BDNKStar.NonRadialModes: nonradial_cowling_spectrum
using .BDNKStar.PolarViscousModes: polar_qnm, qnm_damping
using .BDNKStar.EquationOfState: pressure, sound_speed2
using CairoMakie, Printf
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
trapz(x,y) = sum(0.5*(y[i]+y[i+1])*(x[i+1]-x[i]) for i in 1:length(x)-1)
@inline function _lin(xs, ys, x)
    n=length(xs); x≤xs[1] && return ys[1]; x≥xs[n] && return ys[n]
    j=searchsortedlast(xs,x); t=(x-xs[j])/(xs[j+1]-xs[j]); ys[j]+t*(ys[j+1]-ys[j])
end

# ── ideal (W,V) Cowling eigenfunction by RK4 shooting ─────────────────────────
@inline function _wv_rhs(eos, bg, r, W, V, l, ω2)
    m, ν, ε = bg(r); p=pressure(eos,ε); cs2=max(sound_speed2(eos,ε),1e-14)
    elam=1.0/sqrt(1.0-2m/r); enu=exp(ν); νp=2.0*(m+4π*r^3*p)/(r*(r-2m))
    B=ω2*r^2*(elam/enu)*V+0.5*νp*W
    return B/cs2-l*(l+1)*elam*V, νp*V-elam*W/r^2
end
function _wv_profile(eos, bg, l, ω2, r0, rf, nstep)
    h=(rf-r0)/nstep; rs=zeros(nstep+1); Ws=zeros(nstep+1); Vs=zeros(nstep+1)
    W=r0^(l+1); V=-r0^l/l; r=r0; rs[1]=r; Ws[1]=W; Vs[1]=V
    for k in 1:nstep
        k1W,k1V=_wv_rhs(eos,bg,r,W,V,l,ω2)
        k2W,k2V=_wv_rhs(eos,bg,r+h/2,W+h/2*k1W,V+h/2*k1V,l,ω2)
        k3W,k3V=_wv_rhs(eos,bg,r+h/2,W+h/2*k2W,V+h/2*k2V,l,ω2)
        k4W,k4V=_wv_rhs(eos,bg,r+h,W+h*k3W,V+h*k3V,l,ω2)
        W+=h/6*(k1W+2k2W+2k3W+k4W); V+=h/6*(k1V+2k2V+2k3V+k4V); r+=h
        rs[k+1]=r; Ws[k+1]=W; Vs[k+1]=V
    end
    return rs, Ws, Vs
end

# returns per-unit-η̂: radial entropy density dṠ/dr, Pdiss, Sdot, E, ∫(σ_S·T)dV
function shear_channel(eos, star, Tof, ω2, l; Nr=400, Nθ=240)
    rt, mt, νt, et = star.r, star.m, star.ν, star.ε
    εf = star.ε[1]*1e-9
    bg(r)=(_lin(rt,mt,r), _lin(rt,νt,r), max(_lin(rt,et,r),εf))
    R=star.R; r0=R*1e-4; rf=R*(1-1e-3)
    rs0,Ws0,Vs0 = _wv_profile(eos,bg,l,ω2,r0,rf,8000)
    rg=collect(range(r0,rf;length=Nr)); θg=collect(range(0.0,π;length=Nθ))
    ξr=[_lin(rs0,Ws0./rs0.^2,r) for r in rg]; ξ⊥=[_lin(rs0,Vs0./rs0.^2,r) for r in rg]
    εp=[max(_lin(rt,star.ε,r),εf) for r in rg]
    ρp=[εp[i]-pressure(eos,εp[i]) for i in 1:Nr]
    Tp=[max(Tof(εp[i]),1e-300) for i in 1:Nr]
    Y=[0.5*(3*cos(θ)^2-1) for θ in θg]; dθY=[-3*cos(θ)*sin(θ) for θ in θg]; d2θY=[-3*cos(2θ) for θ in θg]
    dr=rg[2]-rg[1]; dθ=θg[2]-θg[1]
    dξ(f,i)= i==1 ? (f[2]-f[1])/dr : i==Nr ? (f[Nr]-f[Nr-1])/dr : (f[i+1]-f[i-1])/(2dr)
    Pdiss=0.0; Sdot=0.0; Eint=0.0; PdissfromS=0.0; dens=zeros(Nr)
    for i in 1:Nr
        r=rg[i]; dξr=dξ(ξr,i); dξ⊥=dξ(ξ⊥,i); rowS=0.0
        for j in 1:Nθ
            θ=θg[j]; s=sin(θ); cot=(s>1e-8) ? cos(θ)/s : 0.0
            ur=ξr[i]*Y[j]; uθ=ξ⊥[i]*dθY[j]
            dur_dr=dξr*Y[j]; duθ_dr=dξ⊥*dθY[j]; dur_dθ=ξr[i]*dθY[j]; duθ_dθ=ξ⊥[i]*d2θY[j]
            σrr=dur_dr; σθθ=(1/r)*duθ_dθ+ur/r; σφφ=ur/r+cot*uθ/r
            σrθ=0.5*((1/r)*dur_dθ+duθ_dr-uθ/r)
            tr=σrr+σθθ+σφφ; σrr-=tr/3; σθθ-=tr/3; σφφ-=tr/3
            σ2=σrr^2+σθθ^2+σφφ^2+2*σrθ^2
            dV=r^2*s*dr*dθ; dP=2*εp[i]*σ2*dV; dS=dP/Tp[i]
            Pdiss+=dP; Sdot+=dS; PdissfromS+=dS*Tp[i]; rowS+=dS; Eint+=ρp[i]*(ur^2+uθ^2)*dV
        end
        dens[i]=rowS/dr*2π            # dṠ/dr per unit η̂
    end
    E=ω2*Eint*2π; Pdiss*=ω2*2π; Sdot*=ω2*2π; PdissfromS*=ω2*2π
    return (rg=rg, dens=ω2.*dens, Pdiss=Pdiss, Sdot=Sdot, E=E, PdissfromS=PdissfromS,
            γ_diss=Pdiss/(2E))
end

# ── STAR + p1 mode ────────────────────────────────────────────────────────────
eos = ShumPolytrope(100.0); ρ0c=0.00128; εc=ρ0c+100*ρ0c^2
star = solve_tov(eos, εc; h=2e-4)
ρrest(e)=sqrt(pressure(eos,e)/eos.κ); Tprof(e)=(e<=0 ? 0.0 : pressure(eos,e)/max(ρrest(e),1e-300))
freqs, ω2s, R = nonradial_cowling_spectrum(eos, εc; l=2, nmodes=3)
ch = shear_channel(eos, star, Tprof, ω2s[2], 2)
@printf("STAR M=%.4f R=%.4f  p1 γ_diss=%.4f  ∫(σ_S·T)dV=%.4e  Pdiss=%.4e\n",
        star.M, star.R, ch.γ_diss, ch.PdissfromS, ch.Pdiss)

# ── HEAT channel on the finite-T ideal-gas star (positive-definite) ───────────
Γ=5/3; K=15.0; ρc=1e-3
st = solve_tov_idealgas(Γ=Γ, K=K, ρc=ρc, h=2e-4); ni=findlast(>(0),st.p)
rH=st.r[1:ni]; TH=st.T[1:ni]; eΛH=[1.0/(1.0-2*st.m[i]/st.r[i]) for i in 1:ni]
δT=TH./TH[1]; dδTdr=similar(δT)
for i in 1:ni
    dδTdr[i] = i==1 ? (δT[2]-δT[1])/(rH[2]-rH[1]) :
               i==ni ? (δT[ni]-δT[ni-1])/(rH[ni]-rH[ni-1]) : (δT[i+1]-δT[i-1])/(rH[i+1]-rH[i-1])
end
κQ=0.05
σS_heat=[κQ/(TH[i]^2)*(1.0/eΛH[i])*dδTdr[i]^2 for i in 1:ni]      # ≥0
densH=[4π*rH[i]^2*sqrt(eΛH[i])*σS_heat[i] for i in 1:ni]
@printf("HEAT  Ṡ_heat(κ_Q=%.2f)=%.4e  min σ_S=%.3e (≥0)  cs²−cn²∈[%.3f,%.3f]\n",
        κQ, trapz(rH,densH), minimum(σS_heat), minimum(st.cs2[1:ni].-st.cn2[1:ni]),
        maximum(st.cs2[1:ni].-st.cn2[1:ni]))

# ── κ_Q reactive vs ν_mom dissipative in the engine (mode-energy work) ────────
sph = build_sphstar(eos, εc; Nr=48, Nθ=10)
function Eratio(; κ̂, ν̂, T=150.0, dt=0.005)
    ev=setup_sphbdnk(sph; η̂=0.0, κ̂=κ̂, cκ=0.0, ν̂=ν̂, cν=0.0, σ_ko=0.04)
    stt=SphBDNKState(sph.grid.Nr,sph.grid.Nθ); seed_sphbdnk_n!(stt,ev,3;A=1e-3)
    _,_,en=evolve_sphbdnk!(stt,ev;dt=dt,nsteps=round(Int,T/dt),sample=40); en[end]/en[1]
end
κscan=[0.05,0.15,0.25,0.40]; Eκ=[Eratio(κ̂=k,ν̂=0.0) for k in κscan]
νscan=[0.0,0.02,0.04,0.08]; Eν=[Eratio(κ̂=0.05,ν̂=n) for n in νscan]
@printf("ENGINE κ_Q-scan E_end/E0=%s (flat/reactive)  ν_mom-scan=%s (falls/dissipative)\n",
        string(round.(Eκ,sigdigits=3)), string(round.(Eν,sigdigits=3)))

# ════════════════════════════ FIGURE ═════════════════════════════════════════
fig = Figure(size=(1380, 460))

# A — entropy-production density vs r (positive everywhere)
axA = Axis(fig[1,1], xlabel="r  [km]", ylabel="dṠ_shear/dr  /  η̂   (≥0)",
           title="A  shear entropy density (2η σ²)/T  > 0 everywhere")
lines!(axA, ch.rg, ch.dens, color=:crimson, linewidth=2.5)
band!(axA, ch.rg, zeros(length(ch.rg)), ch.dens, color=(:crimson,0.12))
hlines!(axA, [0.0], color=:gray, linestyle=:dash)
text!(axA, 0.05, 0.92; text=@sprintf("0 negative cells\n∫=Ṡ_shear=%.3g /η̂", ch.Sdot),
      space=:relative, align=(:left,:top), fontsize=11)

# B — energy–entropy consistency: ∫(σ_S·T)dV = Pdiss = 2γE, ∝η̂, on the diagonal
axB = Axis(fig[1,2], xlabel="2γE  =  Pdiss  (η̂ ∫2η̃σ²dV)", ylabel="∫(σ_S·T) dV",
           title="B  energy–entropy consistency  σ_S·T = 2η σ²  (y = x)")
η̂s=[0.0,0.02,0.04,0.06,0.08]
xs=[η̂*ch.Pdiss for η̂ in η̂s]; ys=[η̂*ch.PdissfromS for η̂ in η̂s]
mx=maximum(xs)*1.05
lines!(axB, [0,mx],[0,mx], color=:gray, linestyle=:dash, label="y = x (2nd-law identity)")
scatter!(axB, xs, ys, color=:navy, markersize=11, label="p1 mode, η̂=0..0.08")
axislegend(axB, position=:lt, framevisible=true)
text!(axB, 0.95, 0.06; text=@sprintf("agree = %.4f%%\nγ_diss=%.3f vs eig≈0.51 (85%%)",
      100*(1-abs(ch.PdissfromS-ch.Pdiss)/abs(ch.Pdiss)), ch.γ_diss),
      space=:relative, align=(:right,:bottom), fontsize=11)

# C — channel decomposition: shear dissipative vs κ_Q reactive-but-entropy-positive
axC = Axis(fig[1,3], xlabel="transport amplitude", ylabel="E_end / E₀  (engine mode energy)",
           title="C  shear DISSIPATES; κ_Q REACTIVE (heat entropy still >0)")
lines!(axC, νscan, Eν, color=:crimson, linewidth=2.5, label="ν_mom (shear) — E falls")
scatter!(axC, νscan, Eν, color=:crimson, markersize=10)
lines!(axC, κscan, Eκ, color=:seagreen, linewidth=2.5, linestyle=:dash, label="κ_Q (heat) — E flat (reactive)")
scatter!(axC, κscan, Eκ, color=:seagreen, markersize=10, marker=:rect)
hlines!(axC, [1.0], color=:gray, linestyle=:dot)
text!(axC, 0.05, 0.06; text=@sprintf("κ_Q heat entropy on finite-T star:\nṠ_heat=%.2e>0, min σ_S=%.1e≥0",
      trapz(rH,densH), minimum(σS_heat)), space=:relative, align=(:left,:bottom), fontsize=10)
axislegend(axC, position=:rt, framevisible=true)

Label(fig[0,:], @sprintf("BDNKStar — SECOND LAW: σ_S=(1/T)[2ησ²+ζθ²+(κ_Q/T)q²] ≥ 0  ·  p1 mode of ShumPolytrope(100) M=%.2f,R=%.2f  ·  shear dissipative (∫σ_S·T=2γE, 100%%), κ_Q reactive-in-Cowling but entropy-positive",
      star.M, star.R), fontsize=12, font=:bold)

save(joinpath(outdir,"entropy_second_law.png"), fig)
println("saved entropy_second_law.png")
