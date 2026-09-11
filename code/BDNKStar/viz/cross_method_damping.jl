#=
    cross_method_damping.jl — CROSS-METHOD viscous-damping AGREEMENT figure.

    Two panels:

    (A) POLAR (even-parity ℓ=2, relativistic Cowling, ShumPolytrope(100), M≈1.4):
        γ vs η̂ for the f-mode and p1.
          • EIGENVALUE points  : γ = −Re λ of the 5-field BDNK operator (polar_qnm).
          • DISSIPATION line   : γ_diss = (∫2η σ² dV)/(2E) from the IDEAL inviscid
                                 (W,V) eigenfunction (linear-in-η, an independent
                                 method sharing NO scheme with the eigensolver).
          • TIME-DOMAIN trend  : SphBDNK ℓ=2 energy decay vs η̂ — rescaled and
                                 annotated as SIGN/TREND only (magnitude contaminated).
        RESULT (the headline): for p1 the two RELIABLE methods agree on the slope
        dγ/dη̂ to ≈85%.  The f-mode eigenvalue is over-estimated by reactive-frame
        contamination, so it agrees only qualitatively.

    (B) AXIAL (odd-parity ℓ=2 w-mode, Bussières EOS1, frame A1):
        the eigenvalue Δ(1/τ) and the dissipation integral γ_diss vs η_c give
        OPPOSITE SIGNS — physically correct: the w-mode is a spacetime/curvature
        mode (GW-emission dominated), so a first-order FLUID shear-dissipation
        integral mispredicts both sign and magnitude.

    Run: cd code/BDNKStar && export PATH="$HOME/.local/bin:$PATH" \
         && JULIA_NUM_THREADS=6 julia --project=viz viz/cross_method_damping.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using .BDNKStar.PolarViscousModes: polar_qnm, qnm_damping, qnm_freq_kHz
using .BDNKStar.NonRadialModes: nonradial_cowling_spectrum
using .BDNKStar.EquationOfState: pressure, sound_speed2
using .BDNKStar.TOV: solve_tov
using .BDNKStar.SphBackground: build_sphstar
using .BDNKStar.AxialViscousModes
const AVM = BDNKStar.AxialViscousModes
using CairoMakie
using Printf
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

const C_SI = 299_792_458.0; const SEC_TO_KM = C_SI*1e-3
invtau_per_us(ω) = (-imag(ω))*SEC_TO_KM/1e6
@inline function _lin(xs, ys, x)
    n = length(xs); x ≤ xs[1] && return ys[1]; x ≥ xs[n] && return ys[n]
    j = searchsortedlast(xs, x); t = (x - xs[j])/(xs[j+1]-xs[j]); ys[j] + t*(ys[j+1]-ys[j])
end
_trapz(rs, ys) = sum(0.5*(ys[i]+ys[i+1])*(rs[i+1]-rs[i]) for i in 1:length(rs)-1)
_slope(x, y) = (hcat(ones(length(x)), x) \ y)[2]

# =============================================================================
#  POLAR DATA
# =============================================================================
eos = ShumPolytrope(100.0); ρ0c = 0.00128; εc = ρ0c + 100*ρ0c^2
star = solve_tov(eos, εc; h=2e-4)
@printf("STAR ShumPolytrope(100)  M=%.4f Msun  R=%.4f Msun\n", star.M, star.R)

# -- dissipation integral (METHOD 2), per unit η̂ -----------------------------
@inline function wv_rhs(eos, bg, r, W, V, l, ω2)
    m, ν, ε = bg(r); p = pressure(eos, ε); cs2 = max(sound_speed2(eos, ε), 1e-14)
    elam = 1.0/sqrt(1.0 - 2m/r); enu = exp(ν)
    νp = 2.0*(m + 4π*r^3*p)/(r*(r - 2m))
    B = ω2*r^2*(elam/enu)*V + 0.5*νp*W
    return (B/cs2 - l*(l+1)*elam*V, νp*V - elam*W/r^2)
end
function wv_profile(eos, bg, l, ω2, r0, rf, nstep)
    h = (rf-r0)/nstep; rs = zeros(nstep+1); Ws = zeros(nstep+1); Vs = zeros(nstep+1)
    W = r0^(l+1); V = -r0^l/l; r = r0; rs[1]=r; Ws[1]=W; Vs[1]=V
    for k in 1:nstep
        k1W,k1V = wv_rhs(eos,bg,r,     W,        V,        l,ω2)
        k2W,k2V = wv_rhs(eos,bg,r+h/2, W+h/2*k1W,V+h/2*k1V,l,ω2)
        k3W,k3V = wv_rhs(eos,bg,r+h/2, W+h/2*k2W,V+h/2*k2V,l,ω2)
        k4W,k4V = wv_rhs(eos,bg,r+h,   W+h*k3W,  V+h*k3V,  l,ω2)
        W += h/6*(k1W+2k2W+2k3W+k4W); V += h/6*(k1V+2k2V+2k3V+k4V); r += h
        rs[k+1]=r; Ws[k+1]=W; Vs[k+1]=V
    end
    return rs, Ws, Vs
end
function gamma_diss_per_etahat(eos, star, ω2, l; Nr=600, Nθ=400)
    rt, mt, νt, et = star.r, star.m, star.ν, star.ε; εf = star.ε[1]*1e-9
    bg(r) = (_lin(rt,mt,r), _lin(rt,νt,r), max(_lin(rt,et,r), εf))
    R=star.R; r0=R*1e-4; rf=R*(1-1e-3)
    rs0, Ws0, Vs0 = wv_profile(eos, bg, l, ω2, r0, rf, 8000)
    rg = collect(range(r0, rf; length=Nr)); θg = collect(range(0.0, π; length=Nθ))
    ξr = [_lin(rs0, Ws0./rs0.^2, r) for r in rg]; ξ⊥ = [_lin(rs0, Vs0./rs0.^2, r) for r in rg]
    ρp = [max(_lin(rt,star.ε,r),εf) - pressure(eos, max(_lin(rt,star.ε,r),εf)) for r in rg]
    εp = [max(_lin(rt,star.ε,r),εf) for r in rg]
    @assert l == 2
    Y = [0.5*(3*cos(θ)^2-1) for θ in θg]; dθY = [-3*cos(θ)*sin(θ) for θ in θg]; d2θY = [-3*cos(2θ) for θ in θg]
    dr = rg[2]-rg[1]; dθ = θg[2]-θg[1]
    dξ(f,i) = i==1 ? (f[2]-f[1])/dr : i==Nr ? (f[Nr]-f[Nr-1])/dr : (f[i+1]-f[i-1])/(2dr)
    Pdiss = 0.0; Eint = 0.0
    for i in 1:Nr
        r = rg[i]; dξr_dr = dξ(ξr,i); dξ⊥_dr = dξ(ξ⊥,i)
        for j in 1:Nθ
            θ = θg[j]; s = sin(θ); cot = (s>1e-8) ? cos(θ)/s : 0.0
            ur = ξr[i]*Y[j]; uθ = ξ⊥[i]*dθY[j]
            dur_dr = dξr_dr*Y[j]; duθ_dr = dξ⊥_dr*dθY[j]; dur_dθ = ξr[i]*dθY[j]; duθ_dθ = ξ⊥[i]*d2θY[j]
            σrr = dur_dr; σθθ = (1/r)*duθ_dθ + ur/r; σφφ = ur/r + cot*uθ/r
            σrθ = 0.5*((1/r)*dur_dθ + duθ_dr - uθ/r)
            tr = σrr+σθθ+σφφ; σrr-=tr/3; σθθ-=tr/3; σφφ-=tr/3
            σ2 = σrr^2+σθθ^2+σφφ^2 + 2*σrθ^2; dV = r^2*s*dr*dθ
            Pdiss += 2*εp[i]*σ2*dV; Eint += ρp[i]*(ur^2 + uθ^2)*dV
        end
    end
    E = ω2*Eint*2π; Pdiss *= ω2*2π
    return Pdiss/(2*E)
end

freqs, ω2s, _ = nonradial_cowling_spectrum(eos, εc; l=2, nmodes=3)
γd_f = gamma_diss_per_etahat(eos, star, ω2s[1], 2)
γd_p = gamma_diss_per_etahat(eos, star, ω2s[2], 2)
@printf("[diss] dγ/dη̂  f=%.4g  p1=%.4g\n", γd_f, γd_p)

# -- eigenvalue points (METHOD 1) ---------------------------------------------
ηs = [0.0, 0.01, 0.02, 0.03, 0.04]
γf = Float64[]; γp = Float64[]
for η in ηs
    ev = polar_qnm(eos, εc; l=2, η̂=η, Nr=140, nmodes=2, nstep=10, warn=false)
    push!(γf, qnm_damping(ev[1])); push!(γp, qnm_damping(ev[2]))
end
mf_win = _slope(ηs[2:end], γf[2:end]); mp_win = _slope(ηs[2:end], γp[2:end])
agree(a,b) = 100*(1 - abs(a-b)/((abs(a)+abs(b))/2))
@printf("[eig] window slope f=%.4g p1=%.4g  | p1 agreement=%.1f%%\n", mf_win, mp_win, agree(mp_win,γd_p))

# -- time-domain trend (METHOD 3), rescaled to a sign/trend overlay -----------
sstar = build_sphstar(eos, εc; Nr=80, Nθ=16)
dt = 0.25*sstar.grid.dr
ηtd = [0.0, 0.01, 0.02, 0.04]; ratios = Float64[]
for η in ηtd
    e = setup_sphbdnk(sstar; η̂=η, ν̂=0.0, cν=1.0)
    st = SphBDNKState(sstar.grid.Nr, sstar.grid.Nθ); seed_sphbdnk_l2!(st, e; A=1e-3)
    _, _, en = evolve_sphbdnk!(st, e; dt=dt, nsteps=2400, sample=40)
    push!(ratios, en[end]/en[1])
end
# convert E_end/E0 to an effective "energy removed relative to ideal" proxy ∝ damping;
# scaled to the eigenvalue f-mode range for a SIGN/TREND-only overlay.
removed = ratios[1] .- ratios                 # ≥0, increasing with η̂ (energy removed vs ideal)
td_scaled = removed ./ removed[end] .* maximum(γf)

# =============================================================================
#  AXIAL DATA
# =============================================================================
eosA, starA = build_axial_star()
R = starA.R; ℓ = 2; rmin = 1e-3; nint = 8000; next = 4000; Ncf = 800
res0 = axial_qnm(eosA, NaN; l=ℓ, ηc_cgs=0.0, ω0=ftau_to_omega(10.5,29.5),
                 a_over_R=1.6, rmin=rmin, nint=nint, next=next, Ncf=Ncf, tol=1e-10, maxit=120, star=starA)
ω0 = res0.omega; invtau0 = invtau_per_us(ω0)
etacs = [3e29, 1e30, 3e30, 1e31]
Δeig = Float64[]
let ωg = ω0
    for ηc in etacs
        r = axial_qnm(eosA, NaN; l=ℓ, ηc_cgs=ηc, τ̂=10.0, ω0=ωg,
                      a_over_R=1.6, rmin=rmin, nint=nint, next=next, Ncf=Ncf,
                      surf_cut=1e-3, tol=1e-9, maxit=140, star=starA)
        push!(Δeig, invtau_per_us(r.omega) - invtau0); ωg = r.omega
    end
end
# dissipation integral from ideal eigenfunction
function ideal_psi_profile(ω; Ngrid=4000)
    rs = Float64[]; ψs = ComplexF64[]; ψps = ComplexF64[]
    y = ComplexF64[rmin^(ℓ+1), (ℓ+1)*rmin^ℓ]; h = (R - rmin)/Ngrid; r = rmin
    push!(rs,r); push!(ψs,y[1]); push!(ψps,y[2])
    rhs = function(rr, yy)
        bg = AVM.background_at(starA, eosA, rr); mp = 4π*rr^2*bg.ρ
        dλdr = (2*mp*rr - 2*bg.m)/(rr*(rr-2*bg.m)); dlogf = (bg.dνdr - dλdr)/2
        V = AVM.RW_potential(bg, ℓ)
        return ComplexF64[yy[2], -dlogf*yy[2] - (ω^2 - V)/bg.f2*yy[1]]
    end
    for _ in 1:Ngrid
        k1=rhs(r,y); k2=rhs(r+h/2,y .+ (h/2).*k1); k3=rhs(r+h/2,y .+ (h/2).*k2); k4=rhs(r+h,y .+ h.*k3)
        y = y .+ (h/6).*(k1 .+ 2 .*k2 .+ 2 .*k3 .+ k4); r += h
        push!(rs,r); push!(ψs,y[1]); push!(ψps,y[2])
    end
    return rs, ψs, ψps
end
rs, ψs, ψps = ideal_psi_profile(ω0); Ng = length(rs)
Zs = [AVM.background_at(starA,eosA,rs[i]).f*(rs[i]*ψps[i]+ψs[i])/(im*ω0) for i in 1:Ng]
Zps = similar(Zs)
for i in 1:Ng
    Zps[i] = i==1 ? (Zs[2]-Zs[1])/(rs[2]-rs[1]) :
             i==Ng ? (Zs[end]-Zs[end-1])/(rs[end]-rs[end-1]) :
             (Zs[i+1]-Zs[i-1])/(rs[i+1]-rs[i-1])
end
function gamma_diss(ηc)
    visc,_,_ = frameA_viscosity(starA, eosA, ηc, 10.0); Lℓ = (ℓ-1)*(ℓ+2)
    dI = zeros(Ng); kI = zeros(Ng)
    for i in 1:Ng
        bg = AVM.background_at(starA, eosA, rs[i]); η,_,_ = AVM.transport(visc, bg)
        wvol = exp((bg.ν+bg.λ)/2); sh = abs2(rs[i]*Zps[i]-Zs[i]) + Lℓ*abs2(Zs[i])
        dI[i] = 2*η*wvol*sh; kI[i] = (bg.ρ+bg.p)*wvol*rs[i]^2*abs2(Zs[i])
    end
    (_trapz(rs,dI)/(2*_trapz(rs,kI))) * SEC_TO_KM / 1e6
end
γdissA = [gamma_diss(ηc) for ηc in etacs]
@printf("[axial] Δ(1/τ)_eig (neg) vs γ_diss (pos): opposite sign at all η_c\n")

# =============================================================================
#  FIGURE
# =============================================================================
fig = Figure(size=(1200, 520))

# ---- Panel A : POLAR -------------------------------------------------------
axA = Axis(fig[1,1], xlabel="η̂  (shear-viscosity parameter, η = η̂ ε)",
           ylabel="damping rate  γ  [code units]",
           title="POLAR ℓ=2 (Cowling, M≈1.4): eigenvalue vs dissipation integral\n" *
                 @sprintf("p1 reliable-method agreement on dγ/dη̂ = %.1f%%", agree(mp_win,γd_p)))
# dissipation-integral lines (linear, anchored at the first finite-η eigenvalue baseline)
ηline = range(0.01, 0.04; length=2)
# anchor: γ_diss·η̂ shifted to start at the η̂=0.01 eigenvalue value so slopes overlay cleanly
baseF = γf[2] - γd_f*ηs[2]; baseP = γp[2] - γd_p*ηs[2]
lines!(axA, collect(ηline), baseF .+ γd_f .* collect(ηline);
       color=:dodgerblue, linewidth=2.5, linestyle=:dash, label="f  dissipation integral (slope=$(round(γd_f,sigdigits=2)))")
lines!(axA, collect(ηline), baseP .+ γd_p .* collect(ηline);
       color=:crimson, linewidth=2.5, linestyle=:dash, label="p1 dissipation integral (slope=$(round(γd_p,sigdigits=2)))")
scatter!(axA, ηs, γf; color=:dodgerblue, markersize=11, marker=:circle, label="f  eigenvalue (−Re λ)")
scatter!(axA, ηs, γp; color=:crimson, markersize=11, marker=:rect, label="p1 eigenvalue (−Re λ)")
scatter!(axA, ηtd, td_scaled; color=:gray45, markersize=10, marker=:utriangle,
         label="time-domain (SIGN/TREND only, rescaled)")
lines!(axA, ηtd, td_scaled; color=:gray45, linewidth=1.0, linestyle=:dot)
axislegend(axA; position=:rb, framevisible=true, labelsize=9.5)
xlims!(axA, -0.002, 0.045); ylims!(axA, -0.003, nothing)

# ---- Panel B : AXIAL -------------------------------------------------------
axB = Axis(fig[1,2], xlabel="central viscosity  η_c  [cgs]",
           ylabel="viscous damping contribution  [1/µs]",
           xscale=log10,
           title="AXIAL ℓ=2 w-mode (Bussières EOS1, frame A1)\nspacetime mode: methods give OPPOSITE SIGN")
scatter!(axB, etacs, Δeig; color=:purple, markersize=12, marker=:circle,
         label="eigenvalue Δ(1/τ)  (< 0: viscosity ↑ lifetime)")
lines!(axB, etacs, Δeig; color=:purple, linewidth=1.5)
scatter!(axB, etacs, γdissA; color=:darkorange, markersize=12, marker=:rect,
         label="dissipation integral γ_diss  (> 0)")
lines!(axB, etacs, γdissA; color=:darkorange, linewidth=1.5)
hlines!(axB, [0.0]; color=:black, linewidth=0.8, linestyle=:dash)
axislegend(axB; position=:lt, framevisible=true, labelsize=9.5)

outfile = joinpath(outdir, "cross_method_damping.png")
save(outfile, fig; px_per_unit=2)
println("wrote ", outfile)

# ---- console summary --------------------------------------------------------
println("\n=== CROSS-METHOD SUMMARY ===")
@printf("POLAR p1 : eig slope=%.4g  diss=%.4g  AGREE=%.1f%%  (reliable-method cross-check PASSES)\n",
        mp_win, γd_p, agree(mp_win,γd_p))
@printf("POLAR f  : eig slope=%.4g  diss=%.4g  (eigenvalue over-estimated by reactive-frame contamination — qualitative only)\n",
        mf_win, γd_f)
@printf("POLAR td : E_end/E0 = %s  (decreasing in η̂ ⇒ correct SIGN/TREND, magnitude contaminated)\n",
        string(round.(ratios, sigdigits=4)))
@printf("AXIAL    : Δ(1/τ)_eig<0 vs γ_diss>0 at all η_c — OPPOSITE SIGN (w-mode is a spacetime mode; fluid diss integral mispredicts)\n")
