# ======================================================================
# repro/axial_qnm_verify2.jl
#
# Focused fix of the V2/V3 independent-exterior check.  Goal: an INDEPENDENT
# (non-Leaver) determination of the outgoing log-derivative psi'/psi at the
# matching radius a, to cross-validate up_logderiv and hence the inviscid mode.
#
# Method: build the outgoing solution by a high-order ASYMPTOTIC series at large
# r and integrate INWARD to a.  We then (a) compare to Leaver directly, and
# (b) verify it IS outgoing by checking the asymptotic phase.  We also test the
# Leaver solution by integrating it OUTWARD and confirming outgoing behaviour.
# ======================================================================

include("/data/haiyangw/claude/BDNK/code/BDNKStar/repro/axial_qnm.jl")
using Printf

# Asymptotic outgoing series for vacuum Regge-Wheeler.
# psi = e^{i omega r_*} * F(r),  r_* = r + 2M ln(r/2M - 1),
#   F(r) = 1 + sum_{k>=1} a_k / (omega r)^k ... but easier: build psi, psi'
# directly from the tortoise form to high order.  We use the known result:
#   psi ~ e^{i omega r_*} [ 1 + (l(l+1) - ... )/(2 i omega r) + ... ].
# To avoid hand-coding many terms we instead seed at very large r where the
# potential V ~ l(l+1)/r^2 -> 0 and f -> 1, so the EXACT outgoing solution is
# psi = e^{i omega r} (Coulomb-like).  We seed with the leading e^{i omega r_*}
# log-derivative and push r_far HUGE so corrections are negligible, integrating
# inward with very fine steps.  KEY FIX vs verify.jl: use r_* phase log-deriv
#   psi'/psi = i*omega * (dr_*/dr) = i*omega / f   (this part was right)
# but seed at r_far large enough AND check we are on the OUTGOING branch by
# comparing |psi| growth direction.  The previous failure was r_far/steps too
# small for Im(omega) large => exponential mis-seeding.  Here we test r_far.
# ----------------------------------------------------------------------
function out_logderiv_tortoise(M, ω, ℓ, a; r_far, nsteps)
    f_far = 1 - 2M/r_far
    # leading + first two subleading corrections to the outgoing log-derivative.
    # For RW, psi = e^{i w r_*} sum_k b_k/r^k. Leading log-deriv:
    #   L = i w / f + (correction).  First correction from V=f*l(l+1)/r^2:
    #   d/dr ln psi = i w/f  -  V/(2 i w f) + O(1/r^3-ish)   (WKB-ish).
    V_far = f_far*(ℓ*(ℓ+1)/r_far^2 - 6M/r_far^3)
    L = im*ω/f_far - V_far/(2*im*ω*f_far)
    y = ComplexF64[1.0, L]
    h = (a - r_far)/nsteps; r = r_far
    for _ in 1:nsteps
        k1 = _vac_rhs(M, r,     ω, ℓ, y)
        k2 = _vac_rhs(M, r+h/2, ω, ℓ, y .+ (h/2).*k1)
        k3 = _vac_rhs(M, r+h/2, ω, ℓ, y .+ (h/2).*k2)
        k4 = _vac_rhs(M, r+h,   ω, ℓ, y .+ h.*k3)
        y = y .+ (h/6).*(k1 .+ 2 .*k2 .+ 2 .*k3 .+ k4); r += h
    end
    return y[2]/y[1]
end

# Outward test: take Leaver L at a, integrate OUTWARD to r_far, check the
# solution matches the outgoing phase e^{i w r_*} (NOT the ingoing e^{-i w r_*}).
function leaver_is_outgoing(M, ω, ℓ, a; Ncf=800, r_far=200.0, nsteps=200_000)
    L0 = up_logderiv(M, ω, ℓ, a; Ncf=Ncf)
    y = ComplexF64[1.0, L0]
    h = (r_far - a)/nsteps; r = a
    for _ in 1:nsteps
        k1 = _vac_rhs(M, r,     ω, ℓ, y)
        k2 = _vac_rhs(M, r+h/2, ω, ℓ, y .+ (h/2).*k1)
        k3 = _vac_rhs(M, r+h/2, ω, ℓ, y .+ (h/2).*k2)
        k4 = _vac_rhs(M, r+h,   ω, ℓ, y .+ h.*k3)
        y = y .+ (h/6).*(k1 .+ 2 .*k2 .+ 2 .*k3 .+ k4); r += h
    end
    # at r_far compute log-deriv and compare to outgoing (+) and ingoing (-)
    f = 1 - 2M/r_far
    L_num = y[2]/y[1]
    L_out = im*ω/f
    L_in  = -im*ω/f
    return (L_num, L_out, L_in)
end

function main()
    eos, star = build_star()
    ℓ=2; M=star.M; R=star.R; a=1.6*R
    @printf("M=%.5f km  R=%.5f km  a=%.5f km\n\n", M,R,a)

    # test outgoing direction of Leaver at the known inviscid mode
    ω_mode = complex(0.22009212, -0.11294334)
    println("[A] Is the Leaver up-solution OUTGOING? (integrate it outward to r_far)")
    for rf in (60.0, 120.0, 240.0)
        Ln,Lo,Li = leaver_is_outgoing(M, ω_mode, ℓ, a; r_far=rf, nsteps=300_000)
        @printf("  r_far=%6.1f:  L_num=%.5f%+.5fi  outgoing iw/f=%.5f%+.5fi (|Δout|=%.2e)  ingoing=%.5f%+.5fi (|Δin|=%.2e)\n",
                rf, real(Ln),imag(Ln), real(Lo),imag(Lo), abs(Ln-Lo),
                real(Li),imag(Li), abs(Ln-Li))
    end

    println("\n[B] Independent tortoise-seed outgoing L at a vs Leaver L at a")
    for rf in (300.0, 1000.0, 3000.0)
        Ldir = out_logderiv_tortoise(M, ω_mode, ℓ, a; r_far=rf, nsteps=600_000)
        Llv  = up_logderiv(M, ω_mode, ℓ, a; Ncf=800)
        @printf("  r_far=%7.1f:  direct L=%.6f%+.6fi   Leaver L=%.6f%+.6fi   |Δ|=%.3e\n",
                rf, real(Ldir),imag(Ldir), real(Llv),imag(Llv), abs(Ldir-Llv))
    end

    # Now the full INDEPENDENT inviscid root-find using the tortoise outgoing
    # exterior (NOT Leaver), interior from in_logderiv (shooting).
    println("\n[C] INVISCID w-mode with INDEPENDENT (tortoise) exterior, root-find:")
    visc0 = inviscid()
    g = ω -> in_logderiv(star,eos,ω,ℓ,visc0; a=a, rmin=1e-3, nint=8000, next=4000) -
             out_logderiv_tortoise(M, ω, ℓ, a; r_far=1000.0, nsteps=400_000)
    ω_guess = ftau_to_omega(10.5, 29.5)
    ωm, res, ok = find_qnm(g, ω_guess; tol=1e-9, maxit=80)
    f,τ = omega_to_ftau(ωm)
    @printf("  ω=%.8f%+.8fi  |Δ|=%.2e  -> (f,τ)=(%.4f kHz, %.4f μs)  target (10.50, 29.54)\n",
            real(ωm),imag(ωm),res, f,τ)
    @printf("  vs Leaver mode ω=0.22009212-0.11294334i : |Δω|=%.3e\n",
            abs(ωm - ω_mode))
end

main()
