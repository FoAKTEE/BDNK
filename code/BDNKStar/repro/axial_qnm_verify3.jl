# ======================================================================
# repro/axial_qnm_verify3.jl
#
# FULLY INDEPENDENT inviscid w-mode determination that NEVER uses Leaver.
#
# Method (classic direct-integration QNM, e.g. Chandrasekhar-Detweiler):
#   1. interior: integrate the regular psi-seed from r_min to R (decoupled
#      inviscid wave eq), giving (psi(R), psi'(R)).
#   2. exterior: continue with the vacuum Regge-Wheeler eq OUTWARD from R to a
#      large radius r_far (>> 1/|omega|, outside the potential).
#   3. outgoing condition: at r_far the solution must be purely outgoing,
#      psi ~ e^{i omega r_*}.  We decompose the numerical (psi,psi') at r_far
#      into outgoing/ingoing pieces and require the INGOING amplitude to vanish.
#      Equivalently, the QNM residual is
#         g(omega) = (psi'/psi)|_{r_far} - (i omega/f)|_{r_far}  -> 0 .
#   This is independent of the Leaver continued fraction entirely.
#
# Because the outgoing wave grows for Im(omega)<0, integrating OUTWARD keeps the
# growing (outgoing) solution dominant and well-conditioned (opposite to the
# inward integration in verify.jl which failed).  This is the textbook reason
# direct outward integration works for QNMs while inward seeding does not.
# ======================================================================

include("/data/haiyangw/claude/BDNK/code/BDNKStar/repro/axial_qnm.jl")
using Printf

# outward-integrated outgoing residual at r_far
function outgoing_residual(star, eos, ω, ℓ; r_far, rmin=1e-3, nint=8000,
                           next_in=4000, next_out=200_000)
    R = star.R; M = star.M
    # interior (inviscid): regular psi seed to R
    y = integrate_interior(star, eos, ω, ℓ, inviscid(), :psi; rmin=rmin, rsurf=R, nsteps=nint)
    ψR, ψpR = y[1], y[3]
    # exterior outward R -> r_far (vacuum RW)
    ψf, ψpf = integrate_vacuum(M, ω, ℓ, ψR, ψpR; r0=R, r1=r_far, nsteps=next_out)
    f = 1 - 2M/r_far
    # purely-outgoing log-derivative at r_far (leading + 1st WKB correction)
    V = f*(ℓ*(ℓ+1)/r_far^2 - 6M/r_far^3)
    L_out = im*ω/f - V/(2*im*ω*f)
    return ψpf/ψf - L_out
end

function main()
    eos, star = build_star()
    ℓ=2; M=star.M; R=star.R
    @printf("M=%.5f km (=%.5f Msun)  R=%.5f km\n\n", M, mass_solar(star), R)

    println("[D] INVISCID w-mode by DIRECT OUTWARD integration + outgoing BC (NO Leaver)")
    ω_guess = ftau_to_omega(10.5, 29.5)
    for r_far in (60.0, 100.0, 150.0)
        g = ω -> outgoing_residual(star, eos, ω, ℓ; r_far=r_far, nint=8000, next_out=300_000)
        ωm, res, ok = find_qnm(g, ω_guess; tol=1e-8, maxit=80)
        f,τ = omega_to_ftau(ωm)
        @printf("  r_far=%6.1f km:  ω=%.8f%+.8fi  |Δ|=%.2e -> (f,τ)=(%.4f kHz, %.4f μs)\n",
                r_far, real(ωm),imag(ωm), res, f, τ)
    end

    # Leaver reference (from axial_qnm.jl machinery)
    a = 1.6*R
    gL = ω -> matching_residual(star, eos, ω, ℓ, inviscid(); a=a, rmin=1e-3, nint=8000, next=4000, Ncf=800)
    ωL, resL, _ = find_qnm(gL, ω_guess; tol=1e-10, maxit=100)
    fL,τL = omega_to_ftau(ωL)
    @printf("\n  LEAVER reference:  ω=%.8f%+.8fi -> (f,τ)=(%.4f, %.4f)   target (10.50, 29.54)\n",
            real(ωL),imag(ωL), fL,τL)
    println("\n  => If the two methods agree, the inviscid w-mode (and hence interior")
    println("     shooting + unit conversion) is confirmed by a Leaver-independent route.")
end

main()
