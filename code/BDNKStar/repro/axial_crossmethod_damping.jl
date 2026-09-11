#=
    repro/axial_crossmethod_damping.jl

    CROSS-METHOD viscous damping of the AXIAL (odd-parity) sector, computed two
    INDEPENDENT ways in the validated AxialViscousModes framework, and compared
    in the linear-in-η regime.

    METHOD 1 — EIGENVALUE (QNM root find).
        AxialViscousModes.axial_qnm:  the viscosity-shifted ℓ=2 fundamental
        w-mode 1/τ at fixed physical central viscosity η_c (cgs), several values.
        The viscous DAMPING contribution is the η-induced shift
            Δ(1/τ)|_eig(η_c) = (1/τ)(η_c) − (1/τ)(0).
        (Frame A1, τ̂=10 — the frame validated <0.05% vs Bussières Table II.)

    METHOD 2 — DISSIPATION INTEGRAL (first-order perturbation theory).
        Independent of the complex root find.  Reconstruct the AXIAL FLUID
        VELOCITY eigenfunction from the *inviscid* (ideal) w-mode ψ via the
        perfect-fluid axial relation [Bussières line 322]
            iω Z = f (r ψ' + ψ) ,           f = √(e^{ν−λ}),
        Z(r) being the amplitude of the toroidal velocity perturbation.  The
        shear-viscous power dissipated by this velocity field, over the mode
        (kinetic) energy, gives the standard viscous damping rate
            γ_diss = ( ∫ 2 η σ_{ab}σ^{ab} √(-g) d³x ) / ( 2 E_kin ),
        E_kin = ∫ (ρ+p) |δv|² √(-g) d³x   [Cutler–Lindblom-type formula].
        For the axial toroidal field the angular integrals reduce the volume
        integrals to radial integrals; the (l-dependent) angular factor cancels
        between numerator and denominator up to the shear-tensor structure
            S_diss ∝ ∫ η e^{(ν+λ)/2} [ (r Z' − Z)² + (ℓ−1)(ℓ+2) Z² ] dr
            S_kin  ∝ ∫ (ρ+p) e^{(λ−ν)/2} r² |Z|² dr           (kinetic norm)
        with the common angular normalisation N_ang = ℓ(ℓ+1)/(2ℓ+1) cancelling.
        γ_diss is LINEAR in η by construction (η appears once, in S_diss), so it
        is the linear-in-η viscous damping coefficient; we report γ_diss(η_c)
        evaluated with the SAME physical η(r) profile (frame A) used by method 1.

    Both methods share the SAME star, SAME η(r) profile (frame A, η_c calibrated
    from cgs), SAME ℓ=2 — the only difference is eigen-root-find vs energy
    integral, so the comparison is genuinely cross-method.

    η-MODE NOTE.  The axial η-mode is a non-hydrodynamic second-sound resonance
    (frame-dependent per the frame-independence study); its "damping" is not a
    perturbative viscous correction to a hydro mode, so the dissipation-integral
    (perturbation-theory) comparison applies to the w-mode.  We still REPORT the
    η-mode eigenvalue 1/τ for completeness and flag it as frame-dependent.

    RUN:
      cd code/BDNKStar && export PATH="$HOME/.local/bin:$PATH" \
        && JULIA_NUM_THREADS=6 julia --project=. repro/axial_crossmethod_damping.jl
=#

using BDNKStar
using BDNKStar.AxialViscousModes
const AVM = BDNKStar.AxialViscousModes
using BDNKStar.EquationOfState
using BDNKStar.TOV
using Printf

# ----------------------------------------------------------------------
# unit constants (same as the module)
# ----------------------------------------------------------------------
const C_SI      = 299_792_458.0
const SEC_TO_KM = C_SI * 1e-3
const KHZ_TO_KM = 1e3 / SEC_TO_KM
# 1/τ in 1/µs:  τ[µs] = (-1/Imω)/SEC_TO_KM·1e6  =>  1/τ[1/µs] = -Imω·SEC_TO_KM/1e6
invtau_per_us(ω) = (-imag(ω)) * SEC_TO_KM / 1e6

# ======================================================================
# STAR  (Bussières EOS1 reference; same builder as the validated solver)
# ======================================================================
eos, star = build_axial_star()        # κ=100, n=1, ρc=3e15 g/cc
R = star.R; M = star.M
ℓ = 2
a = 1.6R; rmin = 1e-3; nint = 8000; next = 4000; Ncf = 800; surf_cut = 1e-3

@printf("TOV star: M=%.5f km (%.4f M⊙ approx), R=%.5f km, M/R=%.4f\n",
        M, M/1.4766, R, M/R)

# ======================================================================
# METHOD 1 — EIGENVALUE QNM, viscous w-mode
# ======================================================================
println("\n", "="^74)
println("METHOD 1 — EIGENVALUE (axial_qnm) : w-mode 1/τ vs η_c  [frame A1, τ̂=10]")
println("="^74)

# inviscid anchor
res0 = axial_qnm(eos, NaN; l=ℓ, ηc_cgs=0.0, ω0=ftau_to_omega(10.5,29.5),
                 a_over_R=1.6, rmin=rmin, nint=nint, next=next, Ncf=Ncf,
                 tol=1e-10, maxit=120, star=star)
ω0 = res0.omega
invtau0 = invtau_per_us(ω0)        # 1/τ [1/µs] inviscid
@printf("INVISCID w-mode: f=%.4f kHz, τ=%.4f µs, 1/τ=%.6e 1/µs  |Δ|=%.1e\n",
        res0.f_kHz, res0.tau_us, invtau0, res0.residual)

# central-viscosity values (cgs).  Stay in the LINEAR regime (small η_c) plus a
# couple larger to show where linearity breaks.
etacs = [3e29, 1e30, 3e30, 1e31]

# w-mode branch defined by |Re ω| within a window of the inviscid w-mode (avoid
# spurious off-branch roots).  Continuation guess = nearest previous w-mode root.
eig_rows = NamedTuple[]
let ωg = ω0
    for ηc in etacs
        r = axial_qnm(eos, NaN; l=ℓ, ηc_cgs=ηc, τ̂=10.0, ω0=ωg,
                      a_over_R=1.6, rmin=rmin, nint=nint, next=next, Ncf=Ncf,
                      surf_cut=surf_cut, tol=1e-9, maxit=140, star=star)
        f_in_window = 9.5 < r.f_kHz < 10.8
        invtau = invtau_per_us(r.omega)
        Δinvtau = invtau - invtau0
        push!(eig_rows, (ηc=ηc, f=r.f_kHz, tau=r.tau_us, invtau=invtau,
                         Δinvtau=Δinvtau, res=r.residual, ok=r.converged && f_in_window))
        @printf("η_c=%.2e : f=%.4f kHz τ=%.4f µs  1/τ=%.6e  Δ(1/τ)=%.6e 1/µs  |Δ|=%.1e %s\n",
                ηc, r.f_kHz, r.tau_us, invtau, Δinvtau, r.residual,
                (r.converged && f_in_window) ? "" : "[!off-branch]")
        if r.converged && f_in_window
            ωg = r.omega   # continue from this w-mode root
        end
    end
end

# ======================================================================
# METHOD 2 — DISSIPATION INTEGRAL (first-order perturbation theory)
# ======================================================================
println("\n", "="^74)
println("METHOD 2 — DISSIPATION INTEGRAL γ_diss from the IDEAL eigenfunction")
println("="^74)

# --- reconstruct the IDEAL (inviscid) eigenfunction ψ(r), ψ'(r) on a grid ---
# Integrate the decoupled inviscid RW/w-mode equation with the regular seed,
# at the converged inviscid eigenfrequency ω0, recording ψ,ψ' along the way.
function ideal_psi_profile(ω::ComplexF64; Ngrid::Int=4000)
    rs  = Float64[]
    ψs  = ComplexF64[]
    ψps = ComplexF64[]
    rℓ  = rmin^(ℓ+1); drℓ = (ℓ+1)*rmin^ℓ
    y = ComplexF64[rℓ, drℓ]      # (ψ, ψ')
    h = (R - rmin)/Ngrid
    r = rmin
    push!(rs, r); push!(ψs, y[1]); push!(ψps, y[2])
    # inviscid ψ:  ψ'' = -dlogf ψ' - (ω²-V)/f² ψ   (decoupled RW, interior)
    rhs = function(rr, yy)
        bg = AVM.background_at(star, eos, rr)
        mp = 4π*rr^2*bg.ρ
        dλdr = (2*mp*rr - 2*bg.m)/(rr*(rr-2*bg.m))
        dlogf = (bg.dνdr - dλdr)/2
        V = AVM.RW_potential(bg, ℓ)
        ψ, ψp = yy[1], yy[2]
        return ComplexF64[ψp, -dlogf*ψp - (ω^2 - V)/bg.f2*ψ]
    end
    for _ in 1:Ngrid
        k1 = rhs(r,      y)
        k2 = rhs(r+h/2,  y .+ (h/2).*k1)
        k3 = rhs(r+h/2,  y .+ (h/2).*k2)
        k4 = rhs(r+h,    y .+ h.*k3)
        y = y .+ (h/6).*(k1 .+ 2 .*k2 .+ 2 .*k3 .+ k4)
        r += h
        push!(rs, r); push!(ψs, y[1]); push!(ψps, y[2])
    end
    return rs, ψs, ψps
end

rs, ψs, ψps = ideal_psi_profile(ω0)
Ngrid = length(rs)

# --- axial fluid VELOCITY amplitude Z(r) from the perfect-fluid relation ---
#   iω Z = f (r ψ' + ψ)   [Bussières line 322]  =>  Z = f(rψ'+ψ)/(iω)
Zs  = Vector{ComplexF64}(undef, Ngrid)
for i in 1:Ngrid
    bg = AVM.background_at(star, eos, rs[i])
    Zs[i] = bg.f*(rs[i]*ψps[i] + ψs[i]) / (im*ω0)
end
# radial derivative Z'(r) by central finite differences
Zps = Vector{ComplexF64}(undef, Ngrid)
for i in 1:Ngrid
    if i == 1
        Zps[i] = (Zs[2]-Zs[1])/(rs[2]-rs[1])
    elseif i == Ngrid
        Zps[i] = (Zs[end]-Zs[end-1])/(rs[end]-rs[end-1])
    else
        Zps[i] = (Zs[i+1]-Zs[i-1])/(rs[i+1]-rs[i-1])
    end
end

# --- trapezoid integrator over rs ---
function trapz(rs, ys)
    s = 0.0
    for i in 1:length(rs)-1
        s += 0.5*(ys[i]+ys[i+1])*(rs[i+1]-rs[i])
    end
    return s
end

# Build the dissipation- and kinetic-norm integrands using the SAME η(r) profile
# as method 1 (frame A, η_c calibrated from cgs).  γ_diss is linear in η, so we
# compute the per-unit-η shape once and scale by each η_c's calibration.
function gamma_diss(ηc::Float64)
    visc, _, _ = frameA_viscosity(star, eos, ηc, 10.0)
    Lℓ = (ℓ-1)*(ℓ+2)        # angular shear factor for the toroidal field
    dissI = zeros(Float64, Ngrid)
    kinI  = zeros(Float64, Ngrid)
    for i in 1:Ngrid
        bg = AVM.background_at(star, eos, rs[i])
        η, _, _ = AVM.transport(visc, bg)
        eλ = exp(bg.λ); eν = exp(bg.ν)
        # proper-volume / metric weights for axial shear & kinetic norm
        #   √(-g) = e^{(ν+λ)/2} r² ;  shear ~ |rZ'-Z|²+Lℓ|Z|²  ;  kinetic ~ |Z|²
        # (η enters the dissipation integrand linearly)
        sh = abs2(rs[i]*Zps[i] - Zs[i]) + Lℓ*abs2(Zs[i])
        wvol = exp((bg.ν+bg.λ)/2)
        dissI[i] = 2*η * wvol * sh
        kinI[i]  = (bg.ρ+bg.p) * wvol * rs[i]^2 * abs2(Zs[i])
    end
    Sdiss = trapz(rs, dissI)
    Skin  = trapz(rs, kinI)
    # γ_diss = (dissipated power)/(2·E_kin).  Units: [1/km].
    γ_km = Sdiss / (2*Skin)
    # convert 1/km -> 1/µs :  γ[1/µs] = γ[1/km]·SEC_TO_KM/1e6
    γ_per_us = γ_km * SEC_TO_KM / 1e6
    return γ_per_us, Sdiss, Skin
end

println("\n# η_c[cgs]   γ_diss[1/µs]   (linear-in-η dissipation-integral damping)")
diss_rows = NamedTuple[]
for ηc in etacs
    γ, Sd, Sk = gamma_diss(ηc)
    push!(diss_rows, (ηc=ηc, γ=γ))
    @printf("η_c=%.2e : γ_diss=%.6e 1/µs   (Sdiss=%.3e, Skin=%.3e)\n",
            ηc, γ, Sd, Sk)
end

# ======================================================================
# COMPARISON  (linear-in-η regime)
# ======================================================================
println("\n", "="^74)
println("CROSS-METHOD COMPARISON : Δ(1/τ)_eigen  vs  γ_diss   (both 1/µs)")
println("="^74)
function slope(xs, ys)
    n=length(xs); sx=sum(xs); sy=sum(ys); sxx=sum(xs.^2); sxy=sum(xs.*ys)
    (n*sxy - sx*sy)/(n*sxx - sx^2)
end

function do_comparison()
    @printf("%-12s %-16s %-16s %-14s %-10s\n",
            "η_c[cgs]", "Δ(1/τ)_eig", "γ_diss", "|ratio| diss/|eig|", "agree")
    best_agree = 0.0
    for i in 1:length(etacs)
        de = eig_rows[i].Δinvtau
        gd = diss_rows[i].γ
        ratio = gd/abs(de)
        agree = 100*(1 - abs(gd-abs(de))/abs(de))   # magnitude agreement
        if i <= 2
            best_agree = max(best_agree, agree)
        end
        signnote = sign(de)==sign(gd) ? "same-sign" : "OPP-sign"
        @printf("%-12.2e % -16.6e % -16.6e %-14.4f  %.2f%%  [%s]\n",
                etacs[i], de, gd, ratio, agree, signnote)
    end

    # linear-regime slope check vs η_c (all branch points are linear-ish)
    nlin = min(3, length(etacs))
    xs = collect(etacs[1:nlin])
    se = slope(xs, [eig_rows[i].Δinvtau for i in 1:nlin])   # signed
    sd = slope(xs, [diss_rows[i].γ for i in 1:nlin])         # positive
    @printf("\nlinear-in-η slope (first %d η_c):  d/dη_c[Δ(1/τ)_eig]=%.4e (signed) ; d/dη_c[γ_diss]=%.4e\n",
            nlin, se, sd)
    @printf("|slope| magnitude ratio diss/eig = %.4f  -> |slope| agreement = %.2f%%\n",
            abs(sd)/abs(se), 100*(1-abs(abs(sd)-abs(se))/abs(se)))
    @printf("sign(slope) eig=%+d  diss=%+d  -> %s\n",
            Int(sign(se)), Int(sign(sd)), sign(se)==sign(sd) ? "SAME sign" : "OPPOSITE sign")
    @printf("\nbest MAGNITUDE agreement (2 smallest η_c): %.2f%%\n", best_agree)
    return best_agree
end
best_agree = do_comparison()

# ======================================================================
# η-MODE (non-hydro / second-sound; frame-dependent) — eigenvalue only
# ======================================================================
println("\n", "="^74)
println("η-MODE (non-hydro second-sound; FRAME-DEPENDENT) — eigenvalue 1/τ only")
println("="^74)
println("Reported for completeness; NOT subject to the perturbative dissipation")
println("integral (it is not a viscous correction to a hydro mode).")
# locate at a large η_c via a guess from the avoided-crossing data (frame A here).
# Use frame A1 large η_c; the η-mode sits at lower Re ω, larger |Im ω|.
ηc_eta = 5e31
# guess near Re ω ~ 0.13 /km (from avoided_crossing E-branch); Im negative (decaying)
ωη_guess = complex(0.13, -0.02)
visc_eta, _, _ = frameA_viscosity(star, eos, ηc_eta, 10.0)
gη = ω -> AVM.matching_residual(star, eos, ω, ℓ, visc_eta;
                                a=a, rmin=rmin, nint=nint, next=next, Ncf=Ncf,
                                surf_cut=surf_cut)
ωη, resη, okη = AVM.find_qnm(gη, ωη_guess; tol=1e-8, maxit=160)
if okη && imag(ωη) < 0
    fη, τη = omega_to_ftau(ωη)
    @printf("η-mode @η_c=%.0e (frame A1): ω=%.5f%+.5fi  f=%.3f kHz  τ=%.4f µs  1/τ=%.4e 1/µs  |Δ|=%.1e\n",
            ηc_eta, real(ωη), imag(ωη), fη, τη, invtau_per_us(ωη), resη)
else
    @printf("η-mode: not cleanly isolated from guess (ω=%.4f%+.4fi |Δ|=%.1e ok=%s) — frame-dependent, hard to bracket\n",
            real(ωη), imag(ωη), resη, okη)
end

println("\nAXIAL_CROSSMETHOD_DONE")
