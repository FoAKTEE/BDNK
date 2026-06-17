# ======================================================================
# repro/axial_qnm_verify.jl
#
# R4 / stage 3 — ADVERSARIAL VERIFICATION of repro/axial_qnm.jl.
#
# Independent cross-checks of the Bussieres 2604.13208 EOS1 ell=2 w-mode:
#   target inviscid (commented tex l.543): (10.50 kHz, 29.54 us)
#   A1 eta_c=3e29 (Table II, l.495-499) : (10.4884, 29.5870)
#   A1 eta_c=1e31 (Table II, l.530-534) : (10.0898, 30.8857)
#
# Checks performed here:
#   (V1) Unit conversion omega(km^-1) <-> (f,tau) FROM FIRST PRINCIPLES, with an
#        independent constant chain.  convention omega = 2pi f - i/tau [tex l.560].
#   (V2) INVISCID w-mode by an INDEPENDENT exterior method:
#        instead of the Leaver continued fraction, integrate the vacuum
#        Regge-Wheeler eq from R outward to a LARGE radius and impose the
#        purely-outgoing condition psi ~ e^{i omega r_*} directly (asymptotic
#        log-derivative).  If this reproduces the same omega as axial_qnm.jl's
#        Leaver match, then interior shooting + Leaver + unit conversion are all
#        cross-validated, independently of each other.
#   (V3) Leaver up_logderiv vs direct vacuum-RW outgoing log-derivative at the
#        SAME radius a (pure exterior consistency, no interior).
#   (V4) C1-SIGN sensitivity: compute the viscous A1 modes with BOTH signs of the
#        Z->psi coupling (+C1 as written in tex eq.17 l.280, and -C1 as the file
#        uses) and report which reproduces Table II, and the magnitude of the
#        difference.  This isolates whether the published match hinges on the
#        sign flip.
#   (V5) Convergence: vary nint/next/Ncf and matching radius a; confirm the modes
#        are grid-independent (genuine eigenvalues, not numerical artefacts).
#
# REUSE: we include the stage-2 solver file directly (which itself includes the
# stage-1 wave-eq file and the package).  We re-use its integrate_interior,
# integrate_vacuum, up_logderiv, find_qnm, build_star, frameA_viscosity,
# omega_to_ftau, ftau_to_omega, etc.  We DO NOT edit it.
# ======================================================================

include("/data/haiyangw/claude/BDNK/code/BDNKStar/repro/axial_qnm.jl")

using LinearAlgebra
using Printf

# ----------------------------------------------------------------------
# (V1) Independent unit-conversion chain.
# omega is geometric (km^-1).  In geometric units length=time, so an angular
# frequency in km^-1 is converted to s^-1 by multiplying by c[km/s].
#   omega[s^-1]  = omega[km^-1] * c[km/s]              (c = 2.99792458e5 km/s)
#   f[Hz]        = omega_real[s^-1] / (2 pi)
#   f[kHz]       = f[Hz] / 1e3
#   tau[s]       = -1 / Im(omega[s^-1])                (omega=2pi f - i/tau)
#   tau[us]      = tau[s] * 1e6
# This is an independent derivation (no precomputed KHZ_TO_KM/SEC_TO_KM constants).
# ----------------------------------------------------------------------
const C_KM_PER_S = 299_792.458   # speed of light in km/s

function omega_to_ftau_indep(ω::ComplexF64)
    ω_re_s = real(ω) * C_KM_PER_S          # s^-1
    ω_im_s = imag(ω) * C_KM_PER_S          # s^-1
    f_kHz  = ω_re_s / (2π) / 1e3
    τ_us   = (-1.0 / ω_im_s) * 1e6
    return (f_kHz, τ_us)
end

function ftau_to_omega_indep(f_kHz::Real, τ_us::Real)
    ω_re_s = 2π * (f_kHz * 1e3)             # s^-1
    ω_im_s = -1.0 / (τ_us * 1e-6)           # s^-1
    return complex(ω_re_s, ω_im_s) / C_KM_PER_S   # km^-1
end

# ----------------------------------------------------------------------
# (V2/V3) INDEPENDENT exterior outgoing log-derivative by direct integration.
# Vacuum Regge-Wheeler eq (tex eq.33/l.395-403):
#   psi'' = -(f'/f) psi' - (omega^2 - V)/f^2 psi,  f=1-2M/r,
#   V = f[ l(l+1)/r^2 - 6M/r^3].
# Outgoing BC: at large r, psi ~ e^{i omega r_*}, r_* = r + 2M ln(r/2M - 1),
#   so psi'/psi -> i omega (1 + 2M/(r-2M))^... -> i omega dr_*/dr = i omega / f.
# We start FAR OUT (r_far) with the asymptotic outgoing seed (a few-term
# large-r expansion) and integrate INWARD to the matching radius a, returning
# psi'(a)/psi(a).  This is fully independent of Leaver's continued fraction.
# ----------------------------------------------------------------------
function up_logderiv_direct(M::Float64, ω::ComplexF64, ℓ::Int, a::Float64;
                            r_far::Float64=4000.0, nsteps::Int=400_000)
    # asymptotic outgoing solution: psi = e^{i omega r_*} * (1 + a1/(omega r) + ...).
    # Use the leading asymptotic log-derivative at r_far:
    #   d/dr ln psi ~ i omega / f  + correction.  We seed with the exact
    #   tortoise-phase log-derivative plus the leading 1/r centrifugal term.
    f_far  = 1 - 2M/r_far
    # leading-order outgoing: psi'/psi = i*omega/f  (since r_* ' = 1/f).
    # Build psi, psi' at r_far from this; the inward RK4 then resolves the
    # subleading structure self-consistently (the QNM eigenvalue is set by the
    # INTERIOR match, the exterior just needs to be the outgoing solution).
    Lfar = im*ω/f_far + ℓ*(ℓ+1)/(2*im*ω*r_far^2)   # +1st centrifugal correction
    ψ  = ComplexF64(1.0)
    ψp = Lfar*ψ
    y  = ComplexF64[ψ, ψp]
    h  = (a - r_far) / nsteps   # negative (integrate inward)
    r  = r_far
    for _ in 1:nsteps
        k1 = _vac_rhs(M, r,     ω, ℓ, y)
        k2 = _vac_rhs(M, r+h/2, ω, ℓ, y .+ (h/2).*k1)
        k3 = _vac_rhs(M, r+h/2, ω, ℓ, y .+ (h/2).*k2)
        k4 = _vac_rhs(M, r+h,   ω, ℓ, y .+ h.*k3)
        y = y .+ (h/6).*(k1 .+ 2 .*k2 .+ 2 .*k3 .+ k4)
        r += h
    end
    return y[2]/y[1]
end

# ----------------------------------------------------------------------
# (V4) C1-sign-explicit interior log-derivative.
# Re-implements in_logderiv but with an EXPLICIT C1 sign multiplier so we can
# compare +C1 (as literally written in tex eq.17, l.280: "+ C1 Z") vs -C1.
# Everything else identical to the stage-2 in_logderiv (interior RK4 + surface
# BC eq.24 + vacuum continuation).
# ----------------------------------------------------------------------
@inline function _interior_rhs_sign(star, eos, r, ω, ℓ, visc, y, c1sign::Float64)
    P, Q, _ = axial_linear_system(star, eos, r, ω, ℓ, visc)
    ψ, Z, ψp, Zp = y[1], y[2], y[3], y[4]
    Q12 = c1sign * Q[1, 2]
    dψ  = ψp
    dZ  = Zp
    dψp = Q[1,1]*ψ + Q12*Z + P[1,1]*ψp + P[1,2]*Zp
    dZp = Q[2,1]*ψ + Q[2,2]*Z + P[2,1]*ψp + P[2,2]*Zp
    return ComplexF64[dψ, dZ, dψp, dZp]
end

function integrate_interior_sign(star, eos, ω, ℓ, visc, seed::Symbol;
                                 rmin, rsurf, nsteps, c1sign::Float64)
    rℓ  = rmin^(ℓ+1); drℓ = (ℓ+1)*rmin^ℓ
    y = seed === :psi ? ComplexF64[rℓ,0,drℓ,0] : ComplexF64[0,rℓ,0,drℓ]
    h = (rsurf - rmin)/nsteps; r = rmin
    for _ in 1:nsteps
        k1 = _interior_rhs_sign(star,eos,r,     ω,ℓ,visc,y,           c1sign)
        k2 = _interior_rhs_sign(star,eos,r+h/2, ω,ℓ,visc,y.+(h/2).*k1,c1sign)
        k3 = _interior_rhs_sign(star,eos,r+h/2, ω,ℓ,visc,y.+(h/2).*k2,c1sign)
        k4 = _interior_rhs_sign(star,eos,r+h,   ω,ℓ,visc,y.+h.*k3,    c1sign)
        y = y .+ (h/6).*(k1 .+ 2 .*k2 .+ 2 .*k3 .+ k4); r += h
    end
    return y
end

function in_logderiv_sign(star, eos, ω, ℓ, visc; a, c1sign::Float64,
                          rmin=1e-3, nint=8000, next=4000, surf_cut=1e-3)
    R = star.R; M = star.M
    r_s = (1 - surf_cut)*R
    y1 = integrate_interior_sign(star,eos,ω,ℓ,visc,:psi; rmin=rmin,rsurf=r_s,nsteps=nint,c1sign=c1sign)
    y2 = integrate_interior_sign(star,eos,ω,ℓ,visc,:Z;   rmin=rmin,rsurf=r_s,nsteps=nint,c1sign=c1sign)
    B1,B2,B3,B4 = surface_condition(star, ω, visc)
    num = B1*y1[2]+B2*y1[4]+B3*y1[1]+B4*y1[3]
    den = B1*y2[2]+B2*y2[4]+B3*y2[1]+B4*y2[3]
    K = -num/den
    ψR  = y1[1]+K*y2[1]; ψpR = y1[3]+K*y2[3]
    ψa,ψpa = integrate_vacuum(M, ω, ℓ, ψR, ψpR; r0=r_s, r1=a, nsteps=next)
    return ψpa/ψa
end

function matching_residual_sign(star, eos, ω, ℓ, visc; a, c1sign, Ncf=800,
                                rmin=1e-3, nint=8000, next=4000, surf_cut=1e-3)
    Lin = in_logderiv_sign(star,eos,ω,ℓ,visc; a=a, c1sign=c1sign,
                           rmin=rmin, nint=nint, next=next, surf_cut=surf_cut)
    Lup = up_logderiv(star.M, ω, ℓ, a; Ncf=Ncf)
    return Lin - Lup
end

# ======================================================================
function verify()
    println("="^72)
    println("ADVERSARIAL VERIFY — axial_qnm.jl  (Bussieres 2604.13208, EOS1, l=2)")
    println("="^72)

    # ----- (V1) unit conversion -----
    println("\n[V1] UNIT CONVERSION cross-check (independent constant chain)")
    for (lbl,f,τ) in [("inviscid",10.50,29.54),("A1 3e29",10.4884,29.5870),
                      ("A1 1e31",10.0898,30.8857)]
        ωf = ftau_to_omega(f,τ)            # file's converter
        ωi = ftau_to_omega_indep(f,τ)      # independent
        fi,τi = omega_to_ftau_indep(ωf)    # round-trip via independent
        ff,τf = omega_to_ftau(ωi)          # round-trip via file
        @printf("  %-9s  file_omega=%.8f%+.8fi  indep_omega=%.8f%+.8fi  |dω|=%.2e\n",
                lbl, real(ωf),imag(ωf), real(ωi),imag(ωi), abs(ωf-ωi))
        @printf("            roundtrip (f,τ): indep(%.4f,%.4f) file(%.4f,%.4f)\n",
                fi,τi, ff,τf)
    end

    eos, star = build_star()
    ℓ = 2; R = star.R; M = star.M
    a = 1.6*R
    @printf("\nTOV: M=%.5f Msun  R=%.5f km  M/R=%.4f   matching a=%.4f km\n",
            mass_solar(star), R, M/R, a)

    # ----- (V2) INVISCID w-mode via INDEPENDENT exterior (direct outgoing) -----
    println("\n[V2] INVISCID w-mode — independent exterior (direct outgoing RW, NO Leaver)")
    visc0 = inviscid()
    g_indep = ω -> begin
        Lin = in_logderiv(star, eos, ω, ℓ, visc0; a=a, rmin=1e-3, nint=8000, next=4000)
        Lup = up_logderiv_direct(M, ω, ℓ, a)
        Lin - Lup
    end
    ω_guess = ftau_to_omega(10.5, 29.5)
    ω_indep, res_i, ok_i = find_qnm(g_indep, ω_guess; tol=1e-10, maxit=100)
    f_i,τ_i = omega_to_ftau(ω_indep)
    @printf("  INDEP-exterior inviscid:  ω=%.8f%+.8fi  |Δ|=%.2e\n",
            real(ω_indep),imag(ω_indep),res_i)
    @printf("    -> (f,τ)=(%.4f kHz, %.4f μs)   target (10.50, 29.54)\n", f_i,τ_i)

    # Leaver version for comparison
    g_leaver = ω -> matching_residual(star, eos, ω, ℓ, visc0; a=a, rmin=1e-3,
                                      nint=8000, next=4000, Ncf=800)
    ω_lv, res_l, _ = find_qnm(g_leaver, ω_guess; tol=1e-10, maxit=100)
    f_l,τ_l = omega_to_ftau(ω_lv)
    @printf("  LEAVER-exterior inviscid: ω=%.8f%+.8fi  -> (f,τ)=(%.4f, %.4f)\n",
            real(ω_lv),imag(ω_lv), f_l,τ_l)
    @printf("  |ω_indep - ω_leaver| = %.3e  (cross-validates Leaver + interior + units)\n",
            abs(ω_indep - ω_lv))

    # ----- (V3) pure exterior log-derivative consistency at a -----
    println("\n[V3] EXTERIOR-only: Leaver vs direct outgoing log-deriv at a (no interior)")
    for ωt in [ω_lv, complex(0.21,-0.11), complex(0.25,-0.05)]
        Llv = up_logderiv(M, ωt, ℓ, a; Ncf=800)
        Ldr = up_logderiv_direct(M, ωt, ℓ, a)
        @printf("  ω=%.4f%+.4fi  Leaver L=%.6f%+.6fi  direct L=%.6f%+.6fi  |Δ|=%.2e\n",
                real(ωt),imag(ωt), real(Llv),imag(Llv), real(Ldr),imag(Ldr), abs(Llv-Ldr))
    end

    # ----- (V4) C1-SIGN sensitivity on the viscous A1 modes -----
    println("\n[V4] C1-SIGN test on viscous A1 modes (+C1 = literal tex eq.17 ; -C1 = file)")
    for (lbl, ηc, ftgt, τtgt) in [("3e29",3e29,10.4884,29.5870),
                                  ("1e31",1e31,10.0898,30.8857)]
        visc,η̂,ηg = frameA_viscosity(star, eos, ηc, 10.0)
        for sgn in (+1.0, -1.0)
            gv = ω -> matching_residual_sign(star, eos, ω, ℓ, visc; a=a, c1sign=sgn)
            ωv, resv, okv = find_qnm(gv, ω_lv; tol=1e-9, maxit=150)
            fv,τv = omega_to_ftau(ωv)
            tag = sgn>0 ? "+C1(tex)" : "-C1(file)"
            @printf("  η_c=%-5s %-9s -> (f,τ)=(%.4f, %.4f)  target(%.4f, %.4f)  Δf=%+.3f%% Δτ=%+.3f%%  |Δ|=%.1e\n",
                    lbl, tag, fv, τv, ftgt, τtgt,
                    100*(fv-ftgt)/ftgt, 100*(τv-τtgt)/τtgt, resv)
        end
    end

    # ----- (V5) Convergence (grid + matching-radius independence) -----
    println("\n[V5] CONVERGENCE — vary grid & matching radius a (inviscid + A1 3e29)")
    visc329,_,_ = frameA_viscosity(star, eos, 3e29, 10.0)
    for (nint,next,Ncf) in [(4000,2000,400),(8000,4000,800),(16000,8000,1200)]
        g0 = ω -> matching_residual(star,eos,ω,ℓ,visc0; a=a,rmin=1e-3,nint=nint,next=next,Ncf=Ncf)
        ω0,_,_ = find_qnm(g0, ω_guess; tol=1e-10, maxit=100)
        f0,τ0 = omega_to_ftau(ω0)
        gV = ω -> matching_residual(star,eos,ω,ℓ,visc329; a=a,rmin=1e-3,nint=nint,next=next,Ncf=Ncf,surf_cut=1e-3)
        ωV,_,_ = find_qnm(gV, ω0; tol=1e-9, maxit=120)
        fV,τV = omega_to_ftau(ωV)
        @printf("  grid(%5d,%5d,%4d): inviscid(%.4f,%.4f)  A1_3e29(%.4f,%.4f)\n",
                nint,next,Ncf, f0,τ0, fV,τV)
    end
    for afac in (1.3, 1.6, 1.9)
        aa = afac*R
        g0 = ω -> matching_residual(star,eos,ω,ℓ,visc0; a=aa,rmin=1e-3,nint=8000,next=4000,Ncf=800)
        ω0,_,_ = find_qnm(g0, ω_guess; tol=1e-10, maxit=100)
        f0,τ0 = omega_to_ftau(ω0)
        @printf("  a=%.2fR=%.3f km (4M=%.2f,2R=%.2f): inviscid (f,τ)=(%.4f, %.4f)\n",
                afac, aa, 4M, 2R, f0,τ0)
    end

    println("\n" * "="^72)
    println("VERIFY COMPLETE")
    println("="^72)
end

verify()
