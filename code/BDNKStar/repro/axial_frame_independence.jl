# ======================================================================
# repro/axial_frame_independence.jl
#
# BDNK FRAME-INDEPENDENCE TEST via the validated AXIAL viscous QNM solver
# (Bussières et al. 2604.13208; the repo's best-validated viscous solver,
#  <0.05% vs Table II).
#
# PROTOCOL
#   Hold the PHYSICAL transport fixed:
#     * central shear viscosity  η_c = 1e31 cgs  (calibrated into η(0)_geom),
#     * the star fixed (Bussières EOS1: κ=100, n=1, ρc=3e15 g/cc; M/R≈0.21),
#     * EOS / ℓ=2 fixed.
#   SWEEP the FRAME relaxation parameter  τ̂ ∈ {5,10,20,40,80}  (stays causal /
#   well-posed) and the frame-A/B parametrization selector.
#
#   Track:
#     (a) the inviscid w-mode (continued from η_c=0; spacetime mode),
#     (b) the VISCOUS w-mode  (continued from the inviscid w-mode at fixed η_c),
#     (c) the VISCOSITY-DRIVEN η-mode (second-sound resonance; located by
#         complex scan at fixed η_c and tracked across τ̂).
#
#   A PHYSICAL mode must be ~frame-INVARIANT.  A genuine FRAME / non-hydro mode
#   moves with τ̂ (gauge artifact, expected).
#
#   NOTE ON c_η:  the second-sound speed c_η² = η/(τ(p+ρ)) with τ = τ̂ L0 (...).
#   At FIXED physical η, increasing τ̂ DECREASES c_η, so the η-mode frequency
#   (∝ c_η) is EXPECTED to drift with τ̂ — τ̂ sets the relaxation TIME, a genuine
#   transport scale in this parametrization, not a pure gauge.  The frame-
#   invariance claim is sharpest for the spacetime w-mode (which should barely
#   move) and for the η_c→0 viscous shift of the w-mode.
#
# Run: cd code/BDNKStar && JULIA_NUM_THREADS=6 julia --project=. \
#         repro/axial_frame_independence.jl
# ======================================================================
include(joinpath(@__DIR__, "axial_ultracompact.jl"))   # main-guarded; pulls in
                                                       # all stage-1/2 machinery
                                                       # + frameA/B + η-mode tools
using Printf

# ----------------------------------------------------------------------
# fixed physics
# ----------------------------------------------------------------------
const ETA_C_CGS = 1e31            # PHYSICAL central shear viscosity (held fixed)
const TAUHATS   = [5.0, 10.0, 20.0, 40.0, 80.0]   # frame relaxation sweep
const ELL       = 2

eos, star = eta_mode_setup()
R = star.R; M = star.M
a = 1.6*R; rmin = 1e-3; nint = 7000; next = 3500; Ncf = 800; surf_cut = 1e-3
ρc = star.ε[1]; pc = star.p[1]

gof(v) = (ω -> matching_residual(star, eos, ω, ELL, v;
                                 a=a, rmin=rmin, nint=nint, next=next, Ncf=Ncf,
                                 surf_cut=surf_cut))

const RULE = "="^74
const SUB  = "-"^74
println(RULE)
@printf("AXIAL BDNK FRAME-INDEPENDENCE  (Bussières EOS1, ℓ=%d)\n", ELL)
println(RULE)
@printf("star: M=%.5f M⊙  R=%.5f km  M/R=%.4f\n", mass_solar(star), R, M/R)
@printf("FIXED physical η_c = %.2e cgs ;  sweep τ̂ ∈ %s ;  frames A & B\n\n",
        ETA_C_CGS, string(TAUHATS))

# ----------------------------------------------------------------------
# inviscid w-mode (frame-independent by construction: η=0 → no viscosity)
# ----------------------------------------------------------------------
ω_inv, _, ok_inv = find_qnm(gof(inviscid()), ftau_to_omega(10.5, 29.5);
                            tol=1e-10, maxit=120)
f_inv, τ_inv = omega_to_ftau(ω_inv)
@printf("[inviscid w-mode]  f=%.5f kHz  τ=%.5f μs  (1/τ=%.6f 1/μs)  ok=%s\n\n",
        f_inv, τ_inv, 1/τ_inv, ok_inv)

# storage rows: (frame, kind, τ̂, f_kHz, invτ_per_us, τ_us, ok)
rows = NamedTuple[]
push!(rows, (frame="-", kind="inviscid-w", tauhat=NaN, f=f_inv,
             invtau=1/τ_inv, tau=τ_inv, ok=ok_inv))

# second-sound centre speed for the η-mode scaling diagnostic
ceta_c_A(η̂, τ̂) = sqrt(η̂*(ρc+pc)*sound_speed2(eos, ρc) / (τ̂*star.R*(pc+ρc)))  # param A
ceta_c_B(η̂, τ̂) = sqrt(η̂*ρc / (τ̂*(pc+ρc)))                                   # param B (η̂ p L0 form)

# ======================================================================
# (b) VISCOUS w-mode: continue from the inviscid w-mode, fixed η_c, sweep τ̂
# ======================================================================
for (frame, mkvisc) in (("A", frameA_viscosity), ("B", frameB_viscosity))
    println(SUB)
    @printf("[VISCOUS w-mode]  frame %s,  η_c=%.1e cgs (FIXED)\n", frame, ETA_C_CGS)
    println(SUB)
    for τ̂ in TAUHATS
        v, η̂, ηg = mkvisc(star, eos, ETA_C_CGS, τ̂)
        ωr, res, ok = find_qnm(gof(v), ω_inv; tol=1e-9, maxit=140)
        f, τ = omega_to_ftau(ωr)
        @printf("  τ̂=%5.1f : f=%.5f kHz  τ=%.5f μs  1/τ=%.6f 1/μs  |Δ|=%.1e ok=%s\n",
                τ̂, f, τ, 1/τ, res, ok)
        push!(rows, (frame=frame, kind="viscous-w", tauhat=τ̂, f=f,
                     invtau=1/τ, tau=τ, ok=ok))
    end
    @printf("\n")
end

# ======================================================================
# (c) η-mode: locate at fixed η_c (large enough to split into long-lived strip)
#     for each τ̂, then track f & 1/τ vs τ̂.  Frame B (Table II nears an η-mode).
# ======================================================================
"""
locate_eta(v) -> (ω, ok) ; scan the long-lived strip and pick the longest-lived
genuinely-viscous root (smallest |Im ω|), the η-mode signature [main.tex 588].
"""
function locate_eta(v; seed=nothing)
    g = gof(v)
    if seed !== nothing
        ωr, r, ok = _local_root(g, seed; w=0.02, n=13, tol=1e-8, maxit=130)
        if ok && real(ωr) > 0.01 && -0.15 < imag(ωr) < -1e-6 && real(ωr) < 0.6
            return ωr, true
        end
    end
    seeds = scan_complex_for_root(g; re_lo=0.02, re_hi=0.55,
                                  im_lo=-0.10, im_hi=-0.0008, nre=90, nim=46)
    cands = ComplexF64[]
    for (s, _) in seeds
        ωr, r, ok = find_qnm(g, s; tol=1e-8, maxit=120)
        if ok && real(ωr) > 0.01 && -0.12 < imag(ωr) < -1e-6 && real(ωr) < 0.6 &&
           all(abs(ωr - z) > 3e-3 for z in cands)
            push!(cands, ωr)
        end
    end
    isempty(cands) && return (ComplexF64(0), false)
    sort!(cands, by=z->abs(imag(z)))   # longest-lived ⇒ most η-mode-like
    return (cands[1], true)
end

for frame in ("B", "A")
    mkvisc = frame == "B" ? frameB_viscosity : frameA_viscosity
    cetac  = frame == "B" ? ceta_c_B : ceta_c_A
    println(SUB)
    @printf("[η-mode]  frame %s,  η_c=%.1e cgs (FIXED), sweep τ̂\n", frame, ETA_C_CGS)
    println(SUB)
    seed = nothing
    for τ̂ in TAUHATS
        v, η̂, _ = mkvisc(star, eos, ETA_C_CGS, τ̂)
        ωη, ok = locate_eta(v; seed=seed)
        if !ok
            @printf("  τ̂=%5.1f : [no η-mode isolated in long-lived strip]\n", τ̂)
            push!(rows, (frame=frame, kind="eta-mode", tauhat=τ̂, f=NaN,
                         invtau=NaN, tau=NaN, ok=false))
            continue
        end
        seed = ωη
        f, τ = omega_to_ftau(ωη)
        ce = cetac(η̂, τ̂)
        @printf("  τ̂=%5.1f : f=%.5f kHz  τ=%.5f μs  1/τ=%.6f 1/μs  c_η(0)=%.4f /km  f·√τ̂=%.4f\n",
                τ̂, f, τ, 1/τ, ce, f*sqrt(τ̂))
        push!(rows, (frame=frame, kind="eta-mode", tauhat=τ̂, f=f,
                     invtau=1/τ, tau=τ, ok=true))
    end
    @printf("\n")
end

# ======================================================================
# INVARIANCE SUMMARY
# ======================================================================
spread(xs) = isempty(xs) ? NaN : (maximum(xs) - minimum(xs))
relspread(xs) = isempty(xs) ? NaN : 100*(maximum(xs) - minimum(xs))/abs(sum(xs)/length(xs))

println(RULE)
@printf("INVARIANCE SUMMARY  (variation of physical observable across τ̂ sweep)\n")
println(RULE)

function summarise(frame, kind)
    sel = [r for r in rows if r.frame == frame && r.kind == kind && r.ok && isfinite(r.f)]
    isempty(sel) && (return @printf("  %-12s frame %s : [no data]\n", kind, frame))
    fs   = [r.f for r in sel]
    invs = [r.invtau for r in sel]
    @printf("  %-12s frame %s : f∈[%.4f,%.4f] kHz  Δf=%.4f kHz (%.3f%%) | " *
            "1/τ∈[%.5f,%.5f] Δ(1/τ)=%.5f (%.3f%%)\n",
            kind, frame, minimum(fs), maximum(fs), spread(fs), relspread(fs),
            minimum(invs), maximum(invs), spread(invs), relspread(invs))
end

summarise("A", "viscous-w")
summarise("B", "viscous-w")
summarise("A", "eta-mode")
summarise("B", "eta-mode")
@printf("\ninviscid w-mode reference: f=%.5f kHz  1/τ=%.6f 1/μs\n", f_inv, 1/τ_inv)

# machine-readable dump
@printf("\n# CSV frame,kind,tauhat,f_kHz,invtau_per_us,tau_us,ok\n")
for r in rows
    @printf("CSV,%s,%s,%.4f,%.6f,%.6f,%.6f,%s\n",
            r.frame, r.kind, r.tauhat, r.f, r.invtau, r.tau, r.ok)
end
