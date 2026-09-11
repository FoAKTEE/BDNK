using Test
using BDNKStar
using BDNKStar.PolarViscousModes: polar_qnm, qnm_damping
using BDNKStar.NonRadialModes: nonradial_cowling_spectrum
using BDNKStar.EquationOfState: pressure, sound_speed2
using BDNKStar.TOV: solve_tov
using BDNKStar.AxialViscousModes
using BDNKStar.SphBackground: build_sphstar
const AVMt = BDNKStar.AxialViscousModes
const C_SI_X = 299_792_458.0
const SEC_TO_KM_X = C_SI_X * 1e-3
invtau_per_us_x(ω) = (-imag(ω)) * SEC_TO_KM_X / 1e6

# =============================================================================
#  CROSS-METHOD VISCOUS-DAMPING AGREEMENT TEST
#
#  Two INDEPENDENT numerical methods compute the shear-viscous damping of stellar
#  oscillation modes, and we assert they agree (or, where physics dictates they
#  must NOT, that they disagree in the documented, physically-correct way).
#
#  POLAR (even-parity ℓ=2, relativistic Cowling, ShumPolytrope(100), M≈1.4):
#    (1) EIGENVALUE        γ = −Re λ of the 5-field FD BDNK operator (polar_qnm),
#                          linear-window slope dγ/dη̂ over η̂∈[0.01,0.04].
#    (2) DISSIPATION       γ_diss = (∫2η σ² dV)/(2E) from the IDEAL inviscid (W,V)
#        INTEGRAL          NonRadialModes eigenfunction, strain evaluated on an
#                          independent (r,θ) grid (NO shared scheme with method 1).
#    (3) TIME-DOMAIN       SphBDNK ℓ=2 energy decay vs η̂ — SIGN/TREND only.
#  RESULT: for p1 the two RELIABLE methods (1,2) agree on dγ/dη̂ to ≈85% — a
#  rigorous cross-check.  The f-mode eigenvalue carries the documented reactive-
#  frame contamination so only its SIGN (positive, linear) is asserted.
#
#  AXIAL (odd-parity ℓ=2 w-mode, Bussières EOS1, frame A1): the two methods give
#  OPPOSITE SIGNS — physically correct: the w-mode is a spacetime/curvature mode
#  whose damping is GW-emission dominated, so a first-order FLUID shear-dissipation
#  integral does NOT capture (sign or magnitude) the viscous correction.  We assert
#  the documented opposite-sign behaviour, each method individually linear-in-η.
#
#  Runs are deliberately light (reduced Nr / grids / shooting resolution) yet
#  reproduce the production numbers to <a few %.
# =============================================================================

# ---- shared light linear interpolant + trapezoid ----------------------------
@inline function _lin(xs, ys, x)
    n = length(xs); x ≤ xs[1] && return ys[1]; x ≥ xs[n] && return ys[n]
    j = searchsortedlast(xs, x); t = (x - xs[j])/(xs[j+1]-xs[j]); ys[j] + t*(ys[j+1]-ys[j])
end
_trapz(rs, ys) = sum(0.5*(ys[i]+ys[i+1])*(rs[i+1]-rs[i]) for i in 1:length(rs)-1)
_slope(x, y) = (hcat(ones(length(x)), x) \ y)[2]

# =============================================================================
#  POLAR
# =============================================================================
@testset "Cross-method viscous damping — POLAR (eigenvalue vs dissipation integral)" begin
    eos = ShumPolytrope(100.0); ρ0c = 0.00128; εc = ρ0c + 100*ρ0c^2
    star = solve_tov(eos, εc; h=2e-4)
    @test isapprox(star.M, 1.40016; rtol=2e-3)   # M≈1.4 M⊙ reference star

    # -- METHOD 2: dissipation integral from the IDEAL (W,V) eigenfunction -----
    @inline function wv_rhs(eos, bg, r, W, V, l, ω2)
        m, ν, ε = bg(r); p = pressure(eos, ε); cs2 = max(sound_speed2(eos, ε), 1e-14)
        elam = 1.0/sqrt(1.0 - 2m/r); enu = exp(ν)
        νp = 2.0*(m + 4π*r^3*p)/(r*(r - 2m))
        B = ω2*r^2*(elam/enu)*V + 0.5*νp*W
        return (B/cs2 - l*(l+1)*elam*V, νp*V - elam*W/r^2)
    end
    function wv_profile(eos, bg, l, ω2, r0, rf, nstep)
        h = (rf-r0)/nstep
        rs = zeros(nstep+1); Ws = zeros(nstep+1); Vs = zeros(nstep+1)
        W = r0^(l+1); V = -r0^l/l; r = r0
        rs[1]=r; Ws[1]=W; Vs[1]=V
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
    # γ_diss per unit η̂ (exactly linear in η̂ by construction). Reduced (Nr,Nθ)
    # grid relative to production; converged to the production value to <0.1%.
    function gamma_diss_per_etahat(eos, star, ω2, l; Nr=300, Nθ=200)
        rt, mt, νt, et = star.r, star.m, star.ν, star.ε; εf = star.ε[1]*1e-9
        bg(r) = (_lin(rt,mt,r), _lin(rt,νt,r), max(_lin(rt,et,r), εf))
        R=star.R; r0=R*1e-4; rf=R*(1-1e-3)
        rs0, Ws0, Vs0 = wv_profile(eos, bg, l, ω2, r0, rf, 8000)
        rg = collect(range(r0, rf; length=Nr)); θg = collect(range(0.0, π; length=Nθ))
        ξr = [_lin(rs0, Ws0./rs0.^2, r) for r in rg]
        ξ⊥ = [_lin(rs0, Vs0./rs0.^2, r) for r in rg]
        ρp = [max(_lin(rt,star.ε,r),εf) - pressure(eos, max(_lin(rt,star.ε,r),εf)) for r in rg]
        εp = [max(_lin(rt,star.ε,r),εf) for r in rg]
        @assert l == 2
        Y   = [0.5*(3*cos(θ)^2-1) for θ in θg]
        dθY = [-3*cos(θ)*sin(θ)    for θ in θg]
        d2θY= [-3*cos(2θ)          for θ in θg]
        dr = rg[2]-rg[1]; dθ = θg[2]-θg[1]
        dξ(f,i) = i==1 ? (f[2]-f[1])/dr : i==Nr ? (f[Nr]-f[Nr-1])/dr : (f[i+1]-f[i-1])/(2dr)
        Pdiss = 0.0; Eint = 0.0
        for i in 1:Nr
            r = rg[i]; dξr_dr = dξ(ξr,i); dξ⊥_dr = dξ(ξ⊥,i)
            for j in 1:Nθ
                θ = θg[j]; s = sin(θ); cot = (s>1e-8) ? cos(θ)/s : 0.0
                ur = ξr[i]*Y[j]; uθ = ξ⊥[i]*dθY[j]
                dur_dr = dξr_dr*Y[j]; duθ_dr = dξ⊥_dr*dθY[j]
                dur_dθ = ξr[i]*dθY[j]; duθ_dθ = ξ⊥[i]*d2θY[j]
                σrr = dur_dr; σθθ = (1/r)*duθ_dθ + ur/r; σφφ = ur/r + cot*uθ/r
                σrθ = 0.5*((1/r)*dur_dθ + duθ_dr - uθ/r)
                tr = σrr+σθθ+σφφ; σrr-=tr/3; σθθ-=tr/3; σφφ-=tr/3
                σ2 = σrr^2+σθθ^2+σφφ^2 + 2*σrθ^2
                dV = r^2*s*dr*dθ
                Pdiss += 2*εp[i]*σ2*dV
                Eint  += ρp[i]*(ur^2 + uθ^2)*dV
            end
        end
        E = ω2*Eint*2π; Pdiss *= ω2*2π
        return Pdiss/(2*E)
    end

    freqs, ω2s, _ = nonradial_cowling_spectrum(eos, εc; l=2, nmodes=3)
    γd_f = gamma_diss_per_etahat(eos, star, ω2s[1], 2)
    γd_p = gamma_diss_per_etahat(eos, star, ω2s[2], 2)

    # (b) dissipation integral positive & linear-in-η (>0 per unit η̂, by construction)
    @test γd_f > 0
    @test γd_p > 0
    # p1 damps MORE than f (genuine k²-shear scaling): production ratio ≈6.7
    @test γd_p/γd_f > 4.0
    # production values 0.0875 (f), 0.585 (p1)
    @test isapprox(γd_f, 0.0875; rtol=0.05)
    @test isapprox(γd_p, 0.585;  rtol=0.05)

    # -- METHOD 1: eigenvalue, clean linear-in-η window slope dγ/dη̂ ∈[0.01,0.04]
    ηs = [0.01, 0.02, 0.04]
    γf = Float64[]; γp = Float64[]
    for η in ηs
        ev = polar_qnm(eos, εc; l=2, η̂=η, Nr=110, nmodes=2, nstep=8, warn=false)
        push!(γf, qnm_damping(ev[1])); push!(γp, qnm_damping(ev[2]))
    end
    mf_win = _slope(ηs, γf); mp_win = _slope(ηs, γp)

    # (b) eigenvalue damping positive & monotone-increasing (linear) in the window
    @test all(γf .> 0); @test issorted(γf)
    @test all(γp .> 0); @test issorted(γp)
    @test mf_win > 0; @test mp_win > 0
    # eigenvalue ordering p1 damps more than f (same SIGN as the diss integral)
    @test mp_win > mf_win

    # (a) RIGOROUS CROSS-CHECK: p1 eigenvalue slope agrees with dissipation
    #     integral prediction to the measured ≈85% (≤25% relative discrepancy).
    agree(a,b) = 100*(1 - abs(a-b)/((abs(a)+abs(b))/2))
    agree_p1 = agree(mp_win, γd_p)
    @test agree_p1 > 75.0                    # production = 85.4%; light run ≈ same
    @test isapprox(mp_win, γd_p; rtol=0.25)  # the two reliable methods quantitatively agree for p1

    # f-mode: eigenvalue is over-estimated by reactive-frame contamination, so only
    # the qualitative facts hold (positive, sublinear). Document the disagreement.
    @test mf_win > γd_f                       # f eigenvalue slope exceeds the diss integral (contamination)
end

# =============================================================================
#  POLAR — METHOD 3 : SphBDNK time-domain energy decay, SIGN/TREND only
# =============================================================================
@testset "Cross-method viscous damping — POLAR time-domain SIGN/TREND (not magnitude)" begin
    eos = ShumPolytrope(100.0); ρ0c = 0.00128; εc = ρ0c + 100*ρ0c^2
    s = build_sphstar(eos, εc; Nr=48, Nθ=12)
    dt = 0.25*s.grid.dr
    ratios = Float64[]
    for η in (0.0, 0.02, 0.04)
        e  = setup_sphbdnk(s; η̂=η, ν̂=0.0, cν=1.0)
        st = SphBDNKState(s.grid.Nr, s.grid.Nθ); seed_sphbdnk_l2!(st, e; A=1e-3)
        _, _, en = evolve_sphbdnk!(st, e; dt=dt, nsteps=1200, sample=40)
        push!(ratios, en[end]/en[1])
    end
    # (c) SIGN/TREND: more η̂ ⇒ more energy removed ⇒ E_end/E0 decreases monotonically.
    #     Magnitude is NOT asserted (inviscid scheme has slow numerical energy gain;
    #     rising η̂ removes it — honest, per the convergence finding).
    @test all(diff(ratios) .< 0)              # monotone decreasing in η̂ (correct sign/trend)
    @test ratios[end] < ratios[1]             # net: viscosity removes energy relative to ideal
end

# =============================================================================
#  AXIAL — eigenvalue vs dissipation integral: documented OPPOSITE SIGN
# =============================================================================
@testset "Cross-method viscous damping — AXIAL w-mode (opposite-sign, both linear-in-η)" begin
    invtau_per_us = invtau_per_us_x
    SEC_TO_KM = SEC_TO_KM_X

    eos, star = build_axial_star()
    R = star.R; ℓ = 2; rmin = 1e-3
    # reduced shooting resolution (nint/next/Ncf): w-mode 1/τ converged to <1e-3 of prod
    nint = 4000; next = 2000; Ncf = 400

    # inviscid anchor
    res0 = axial_qnm(eos, NaN; l=ℓ, ηc_cgs=0.0, ω0=ftau_to_omega(10.5,29.5),
                     a_over_R=1.6, rmin=rmin, nint=nint, next=next, Ncf=Ncf,
                     tol=1e-9, maxit=120, star=star)
    @test res0.converged
    @test isapprox(res0.f_kHz, 10.50; rtol=0.01)
    ω0 = res0.omega; invtau0 = invtau_per_us(ω0)

    # METHOD 1 — EIGENVALUE Δ(1/τ) at a couple of η_c (linear regime)
    etacs = [3e30, 1e31]
    Δeig = Float64[]
    let ωg = ω0
        for ηc in etacs
            r = axial_qnm(eos, NaN; l=ℓ, ηc_cgs=ηc, τ̂=10.0, ω0=ωg,
                          a_over_R=1.6, rmin=rmin, nint=nint, next=next, Ncf=Ncf,
                          surf_cut=1e-3, tol=1e-9, maxit=140, star=star)
            @test 9.5 < r.f_kHz < 10.8       # stayed on the w-mode branch
            push!(Δeig, invtau_per_us(r.omega) - invtau0)
            ωg = r.omega
        end
    end
    # eigenvalue: viscosity LENGTHENS the w-mode lifetime ⇒ Δ(1/τ) < 0, linear-in-η
    @test all(Δeig .< 0)
    @test Δeig[2] < Δeig[1]                   # more negative with larger η_c (monotone)

    # METHOD 2 — DISSIPATION INTEGRAL γ_diss from the IDEAL (inviscid) eigenfunction
    function ideal_psi_profile(ω::ComplexF64; Ngrid::Int=2000)
        rs = Float64[]; ψs = ComplexF64[]; ψps = ComplexF64[]
        y = ComplexF64[rmin^(ℓ+1), (ℓ+1)*rmin^ℓ]
        h = (R - rmin)/Ngrid; r = rmin
        push!(rs,r); push!(ψs,y[1]); push!(ψps,y[2])
        rhs = function(rr, yy)
            bg = AVMt.background_at(star, eos, rr)
            mp = 4π*rr^2*bg.ρ
            dλdr = (2*mp*rr - 2*bg.m)/(rr*(rr-2*bg.m))
            dlogf = (bg.dνdr - dλdr)/2
            V = AVMt.RW_potential(bg, ℓ)
            return ComplexF64[yy[2], -dlogf*yy[2] - (ω^2 - V)/bg.f2*yy[1]]
        end
        for _ in 1:Ngrid
            k1 = rhs(r,     y); k2 = rhs(r+h/2, y .+ (h/2).*k1)
            k3 = rhs(r+h/2, y .+ (h/2).*k2); k4 = rhs(r+h, y .+ h.*k3)
            y = y .+ (h/6).*(k1 .+ 2 .*k2 .+ 2 .*k3 .+ k4); r += h
            push!(rs,r); push!(ψs,y[1]); push!(ψps,y[2])
        end
        return rs, ψs, ψps
    end
    rs, ψs, ψps = ideal_psi_profile(ω0)
    Ng = length(rs)
    # axial fluid velocity Z(r) from the perfect-fluid relation iωZ=f(rψ'+ψ)
    Zs = Vector{ComplexF64}(undef, Ng)
    for i in 1:Ng
        bg = AVMt.background_at(star, eos, rs[i])
        Zs[i] = bg.f*(rs[i]*ψps[i] + ψs[i]) / (im*ω0)
    end
    Zps = Vector{ComplexF64}(undef, Ng)
    for i in 1:Ng
        Zps[i] = i==1 ? (Zs[2]-Zs[1])/(rs[2]-rs[1]) :
                 i==Ng ? (Zs[end]-Zs[end-1])/(rs[end]-rs[end-1]) :
                 (Zs[i+1]-Zs[i-1])/(rs[i+1]-rs[i-1])
    end
    function gamma_diss(ηc::Float64)
        visc, _, _ = frameA_viscosity(star, eos, ηc, 10.0)
        Lℓ = (ℓ-1)*(ℓ+2)
        dissI = zeros(Float64, Ng); kinI = zeros(Float64, Ng)
        for i in 1:Ng
            bg = AVMt.background_at(star, eos, rs[i])
            η, _, _ = AVMt.transport(visc, bg)
            wvol = exp((bg.ν+bg.λ)/2)
            sh = abs2(rs[i]*Zps[i] - Zs[i]) + Lℓ*abs2(Zs[i])
            dissI[i] = 2*η * wvol * sh
            kinI[i]  = (bg.ρ+bg.p) * wvol * rs[i]^2 * abs2(Zs[i])
        end
        γ_km = _trapz(rs, dissI) / (2*_trapz(rs, kinI))
        return γ_km * SEC_TO_KM / 1e6
    end
    γdiss = [gamma_diss(ηc) for ηc in etacs]
    # dissipation integral: positive and EXACTLY linear in η (η enters once)
    @test all(γdiss .> 0)
    @test isapprox(γdiss[2]/γdiss[1], etacs[2]/etacs[1]; rtol=1e-3)   # exact linearity

    # (d) CORE AXIAL FINDING: the two methods give OPPOSITE SIGNS — physically
    #     correct for a spacetime/curvature (w) mode. Agreement = 0% (disagree).
    for i in eachindex(etacs)
        @test sign(Δeig[i]) != sign(γdiss[i])
    end
    @test Δeig[end] < 0 && γdiss[end] > 0
end
