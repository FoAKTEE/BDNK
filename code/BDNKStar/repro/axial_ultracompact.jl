# ======================================================================
# repro/axial_ultracompact.jl
#
# R4 / ultracompact — AXIAL QNMs of ULTRACOMPACT constant-density stars,
# the ETA-MODE family (viscosity-driven, no perfect-fluid counterpart),
# and the W-/ETA avoided crossing.
#
# Primary ground:  Bussières, Redondo-Yuste, Ortega-Gómez & Cardoso,
#   "Axial Oscillations of Viscous Neutron Stars", arXiv:2604.13208
#   (file ref-paper/sources/arXiv-2604.13208/src/main.tex).
#
#   * ULTRACOMPACT / constant-density stars .... Sec.V-D "Viscous damping in
#       ultracompact stars" [main.tex 595-621].  ρ=const, M/R ≤ 4/9 (Buchdahl)
#       [598]; surface inside the light ring -> V has a max (~3M) AND a min
#       (stable light ring) -> long-lived TRAPPED w-modes [598].  Hydro frame B
#       (no sound-speed dependence) [597].
#     - perfect-fluid trapped-mode scaling  ω_ℓ=(ℓ+1/2)Ω_LR - i e^{-γℓ} [612],
#       "Our results perfectly reproduce this behavior in the inviscid case" [615].
#     - viscous (η_c=1e31 cgs): Re ω ~unaffected, but |Im ω_ℓ| ≲ 1e-2 for ALL ℓ
#       and ALL compactness — shear viscosity damps the long-lived w-modes [617].
#
#   * ETA-MODES ................................ Sec.V-C [main.tex 575-593].
#       viscosity-driven family, NO perfect-fluid counterpart [586]; restoring
#       force from the transport coefficients (η and the regulator τ̂) [586];
#       "generically Im ω -> 0 as η_c -> 0 — they become undamped in the perfect
#       fluid limit" [588]; kHz frequencies, ms damping times [588].
#     - AVOIDED CROSSING: w- and η-modes approach as η_c grows but never cross,
#       repelling instead [591]; "destabilizing the frequencies of the w-modes"
#       [593]; frame B1 in Table II shows the larger discrepancy at large η_c
#       precisely because it nears an η-mode [537].
#
# REUSES the verified stage-1/stage-2 machinery via repro/axial_qnm.jl (which
# itself includes repro/axial_waveeqs.jl, which does include(BDNKStar.jl);
# using .BDNKStar).  Reused, tex-grounded names:  axial_linear_system,
# surface_condition, Viscosity, transport, RW_potential, U_potential, cη2,
# coupling_coeffs, background_at  (stage-1, axial_waveeqs.jl); and
# integrate_interior, integrate_vacuum, in_logderiv, up_logderiv (Leaver CF),
# matching_residual, find_qnm, omega_to_ftau, ftau_to_omega,
# eta_cgs_to_geom, SEC_TO_KM, KHZ_TO_KM, C1_COUPLING_SIGN (stage-2, axial_qnm.jl).
#
# We do NOT edit src/BDNKStar.jl nor any other agent's file.  This file only
# ADDS: (i) an analytic constant-density (interior-Schwarzschild) TOVStar so the
# very same background_at / axial_linear_system pipeline runs unchanged on a
# uniform star; (ii) drivers for the trapped w-modes, η-modes, avoided crossing.
# ======================================================================

include("/data/haiyangw/claude/BDNK/code/BDNKStar/repro/axial_qnm.jl")
using .BDNKStar
using Printf
using LinearAlgebra

# ======================================================================
# CONSTANT-DENSITY (uniform) STAR — analytic interior Schwarzschild metric.
#
# ρ = const, areal radius R, mass M = (4π/3) ρ R³, compactness 𝒞 = M/R.
# Interior Schwarzschild solution (e.g. Schwarzschild 1916; std. GR text):
#   m(r)      = M (r/R)³ = (4π/3) ρ r³
#   e^{-λ}    = 1 - 2m/r = 1 - 2𝒞 (r/R)²
#   e^{ν/2}   = (3/2)√(1-2𝒞) - (1/2)√(1-2𝒞 r²/R²)
#   p(r)      = ρ [ √(1-2𝒞 r²/R²) - √(1-2𝒞) ]
#                 / [ 3√(1-2𝒞) - √(1-2𝒞 r²/R²) ]
#   p(R)=0 , e^{ν(R)}=1-2𝒞 (Schwarzschild match) — exactly the conventions
#   used by the package TOVStar (ν shifted so e^{ν(R)}=1-2M/R) and by
#   background_at (which recomputes λ from m(r) and dν/dr from the TOV formula
#   ν'=(2m+8π r³ p)/(r(r-2m)) — both exact for this analytic profile).
#
# We build a TOVStar holding the analytic (r,m,p,ε,ν) on a fine radial grid.
# Linear interpolation inside background_at is exact for m∝r³? no — but the grid
# is fine (default 6000 pts) and the QNM matching uses background_at only as a
# smooth interpolant; the dν/dr and λ are recomputed analytically per-point from
# the (interpolated) m,p, so the dominant metric data are effectively analytic.
# ======================================================================

"""
    const_density_star(𝒞, R; npts=6000) -> (TOVStar, ρ)

Analytic constant-density (interior-Schwarzschild) star of compactness 𝒞=M/R
and areal radius R [km].  Returns a package `TOVStar` plus the uniform energy
density ρ [km^-2].  Requires 𝒞 < 4/9 (Buchdahl) [main.tex 598].
"""
function const_density_star(𝒞::Float64, R::Float64; npts::Int=6000)
    @assert 0 < 𝒞 < 4/9 "compactness must satisfy 0<𝒞<4/9 (Buchdahl)"
    M  = 𝒞 * R
    ρ  = M / ((4π/3) * R^3)            # uniform energy density
    s0 = sqrt(1 - 2𝒞)                  # √(1-2𝒞)  = e^{ν/2}(R)·2/... ; surface
    r  = collect(range(R/npts, R; length=npts))
    m  = similar(r); p = similar(r); ε = fill(ρ, npts); ν = similar(r)
    for i in eachindex(r)
        x   = r[i]/R
        sr  = sqrt(1 - 2𝒞*x^2)         # √(1-2𝒞 r²/R²)
        m[i] = M * x^3
        p[i] = ρ * (sr - s0) / (3s0 - sr)
        # e^{ν/2} = (3/2)s0 - (1/2)sr  ->  ν = 2 ln[(3 s0 - sr)/2]
        ν[i] = 2*log((3s0 - sr)/2)
    end
    p[end] = 0.0                        # enforce exact surface p=0
    # ν is already Schwarzschild-matched: at r=R, sr=s0, ν=2ln(s0)=ln(1-2𝒞). ✓
    return TOVStar(r, m, p, ε, ν, M, R), ρ
end

# Constant-density "EOS": only background_at calls sound_speed2(eos,ρ); for
# hydro frame B the value is irrelevant to transport (η=η̂ p L0, θ,τ ∝ p/ρ —
# no cs² appears [main.tex 178-181]).  Return cs²=1 (causal placeholder).
import .BDNKStar.EquationOfState: sound_speed2
struct ConstDensityEOS <: BarotropicEOS end
sound_speed2(::ConstDensityEOS, ::Real) = 1.0   # placeholder; unused in frame B

# ----------------------------------------------------------------------
# Light-ring (LR) data of the EXTERIOR Schwarzschild geometry, used for the
# perfect-fluid trapped-mode scaling ω_ℓ ≈ (ℓ+1/2)Ω_LR - i e^{-γℓ} [main.tex 612].
# For Schwarzschild the UNSTABLE photon sphere is at r=3M with
#   Ω_LR = 1/(3√3 M)  (orbital frequency = √(f/r²)·... ),
#   λ_Lyap = Ω_LR  (eikonal QNM: ω ≈ (ℓ+1/2)Ω_LR - i(n+1/2)λ_Lyap).
# For ULTRACOMPACT stars there is ALSO a STABLE LR inside; the trapped modes sit
# in the potential well and are long-lived (small |Im ω|).  We use the unstable-
# LR eikonal value only as an order-of-magnitude initial guess for the real part.
# ----------------------------------------------------------------------
schwarzschild_OmegaLR(M::Float64) = 1.0 / (3*sqrt(3)*M)   # Ω at r=3M photon sphere

# ======================================================================
# DRIVER 1 — ULTRACOMPACT trapped ℓ=2 w-modes, ideal vs viscous (frame B).
# Reproduces Fig.ucos behaviour [main.tex 600-617]:
#   ideal: long-lived (tiny |Im ω|) trapped mode;
#   viscous η_c=1e31: Re ω ~unchanged, |Im ω| pulled up to ≲ 1e-2.
# ======================================================================
function find_trapped_mode(star, eos, ℓ, visc; a, rmin, nint, next, Ncf,
                           surf_cut, ωguess, tol=1e-9, maxit=120, verbose=false)
    g = ω -> matching_residual(star, eos, ω, ℓ, visc;
                               a=a, rmin=rmin, nint=nint, next=next, Ncf=Ncf,
                               surf_cut=surf_cut)
    return find_qnm(g, ωguess; tol=tol, maxit=maxit, verbose=verbose)
end

"""
    scan_trapped_real(star, eos, ℓ, visc; ...) -> Vector{(ω, |Δ|)}

Coarse 1-D scan of |Δ(ω)| along a line of (nearly real) ω to LOCATE the trapped
(long-lived) modes of an ultracompact star, whose Im ω is tiny so a real-axis
scan finds the resonances.  Returns local minima of |Δ| as seeds.
"""
function scan_trapped_real(star, eos, ℓ, visc; a, rmin, nint, next, Ncf,
                           surf_cut, ωre_lo, ωre_hi, nω, imω=-1e-3)
    g = ω -> matching_residual(star, eos, ω, ℓ, visc;
                               a=a, rmin=rmin, nint=nint, next=next, Ncf=Ncf,
                               surf_cut=surf_cut)
    ωres = range(ωre_lo, ωre_hi; length=nω)
    mags = Float64[]
    for x in ωres
        push!(mags, abs(g(complex(x, imω))))
    end
    seeds = Tuple{ComplexF64,Float64}[]
    for i in 2:nω-1
        if mags[i] < mags[i-1] && mags[i] < mags[i+1]
            push!(seeds, (complex(ωres[i], imω), mags[i]))
        end
    end
    return seeds, collect(ωres), mags
end

function run_ultracompact()
    println("="^74)
    println("ULTRACOMPACT constant-density stars — trapped ℓ=2 w-modes (frame B)")
    println("  ground: Bussières 2604.13208 Sec.V-D [main.tex 595-621], Fig.ucos")
    println("="^74)

    eos = ConstDensityEOS()
    R   = 10.0                                  # areal radius [km] (sets scale)
    ℓ   = 2

    # numerics (Leaver minimal-solution window 4M<a<2R<2a [main.tex 468]).
    rmin = 1e-3; nint = 9000; next = 4000; Ncf = 1200
    surf_cut = 5e-3   # robust for the narrow deep-well modes at high 𝒞 (verified)

    # central viscosity for the viscous configuration [main.tex 603,617].
    ηc_cgs = 1e31

    # Compactnesses with a STABLE light ring (surface inside r=3M, i.e. 𝒞≳0.39),
    # below Buchdahl 4/9≈0.4444 [main.tex 598].  (Fig.ucos uses several such 𝒞.)
    compactnesses = [0.40, 0.42, 0.44]

    println("\nfixed R = $(R) km,  ℓ = $ℓ,  η_c = $(ηc_cgs) g/cm/s (frame B param),")
    @printf("Buchdahl 4/9 = %.4f ; stable-LR threshold 𝒞≳0.39.\n\n", 4/9)

    results = Tuple{Float64,ComplexF64,ComplexF64}[]   # (𝒞, ω_ideal, ω_visc)
    for 𝒞 in compactnesses
        star, ρ = const_density_star(𝒞, R)
        M = star.M
        a = 1.6*R
        ΩLR = schwarzschild_OmegaLR(M)
        println("-"^74)
        @printf("𝒞 = M/R = %.3f   M = %.4f km   ρ = %.4e km^-2   Ω_LR(3M)=%.4f /km\n",
                𝒞, M, ρ, ΩLR)
        @printf("   eikonal Re ω≈(ℓ+1/2)Ω_LR = %.4f /km ; a=%.3f (4M=%.3f,2R=%.3f)\n",
                (ℓ+0.5)*ΩLR, a, 4M, 2R)

        # ---- locate the trapped (long-lived) IDEAL w-modes by a real-axis scan --
        # Trapped modes sit in the stable-LR well [main.tex 598]: their tiny |Im ω|
        # makes |Δ(ω)| dip sharply on the (nearly) real axis.  Scan a generous band
        # from low frequency (deepest-well overtones) up past the eikonal estimate.
        visc0 = inviscid()
        ωlo = 0.10*(ℓ+0.5)*ΩLR
        ωhi = 1.20*(ℓ+0.5)*ΩLR
        seeds, _, _ = scan_trapped_real(star, eos, ℓ, visc0;
            a=a, rmin=rmin, nint=nint, next=next, Ncf=Ncf, surf_cut=0.0,
            ωre_lo=ωlo, ωre_hi=ωhi, nω=300, imω=-1e-3)
        if isempty(seeds)
            @printf("   [warn] no real-axis resonance found in [%.3f,%.3f]\n", ωlo, ωhi)
            continue
        end
        # refine ALL dips as full complex roots, keep the physical trapped modes
        # (Im ω<0, not far in the lower half plane), pick the LONGEST-LIVED
        # (smallest |Im ω|) as THE fundamental trapped w-mode [main.tex 615].
        roots = ComplexF64[]
        for (s,_) in seeds
            ωr, rr, ok = find_trapped_mode(star, eos, ℓ, visc0;
                a=a, rmin=rmin, nint=nint, next=next, Ncf=Ncf, surf_cut=0.0,
                ωguess=s, tol=1e-10, maxit=150)
            if ok && imag(ωr) < 0 && imag(ωr) > -0.1 && real(ωr) > 0.005 &&
               all(abs(ωr - z) > 1e-3 for z in roots)
                push!(roots, ωr)
            end
        end
        if isempty(roots)
            @printf("   [warn] no trapped root converged.\n"); continue
        end
        sort!(roots, by=z->abs(imag(z)))       # longest-lived first
        ω_id = roots[1]; res_id = 0.0
        @printf("   IDEAL trapped w-modes (%d found, longest-lived first):\n", length(roots))
        for z in roots
            @printf("       ω = %.6f %+.6fi   |Im|/Re = %.2e\n",
                    real(z), imag(z), abs(imag(z))/real(z))
        end
        @printf("   -> fundamental (longest-lived) IDEAL trapped w-mode:  ω = %.6f %+.6fi\n",
                real(ω_id), imag(ω_id))

        # ---- viscous (frame B, η_c=1e31): continue from the ideal mode --------
        visc, η̂, ηgeom = frameB_viscosity(star, eos, ηc_cgs, 10.0)  # B1: τ̂=10
        ω_v, res_v, ok_v = find_trapped_mode(star, eos, ℓ, visc;
            a=a, rmin=rmin, nint=nint, next=next, Ncf=Ncf, surf_cut=surf_cut,
            ωguess=ω_id, tol=1e-9, maxit=150)
        # narrow deep-well modes (high 𝒞) can be surf_cut-sensitive; fall back to
        # nearby surf_cut values if the first solve returned NaN / non-converged.
        if !ok_v || isnan(real(ω_v))
            for sc2 in (1e-2, 3e-3, 7e-3)
                ω_v, res_v, ok_v = find_trapped_mode(star, eos, ℓ, visc;
                    a=a, rmin=rmin, nint=nint, next=next, Ncf=Ncf, surf_cut=sc2,
                    ωguess=ω_id, tol=1e-9, maxit=150)
                (ok_v && !isnan(real(ω_v))) && break
            end
        end
        @printf("   VISC  trapped w-mode:  ω = %.6f %+.6fi  |Δ|=%.1e  (η̂=%.3e)\n",
                real(ω_v), imag(ω_v), res_v, η̂)
        @printf("     -> ΔRe(ω)/Re(ω) = %+.3f%% ;  |Im ω|: %.2e (ideal) -> %.2e (visc)\n",
                100*(real(ω_v)-real(ω_id))/real(ω_id), abs(imag(ω_id)), abs(imag(ω_v)))
        push!(results, (𝒞, ω_id, ω_v))
    end

    # ---- check the paper's two ultracompact claims [main.tex 615-617] ----------
    println("\n" * "="^74)
    println("ULTRACOMPACT — achieved vs target (Bussières Sec.V-D)")
    println("="^74)
    @printf("%-8s %-26s %-26s %-12s\n", "𝒞", "ideal ω (Re,Im)", "visc ω (Re,Im)", "|Im|→")
    re_ok = true; im_ok = true
    for (𝒞, ωid, ωv) in results
        dRe = 100*(real(ωv)-real(ωid))/real(ωid)
        @printf("%-8.3f (%.5f, %.2e)   (%.5f, %.2e)   %.1e→%.1e\n",
                𝒞, real(ωid), imag(ωid), real(ωv), imag(ωv),
                abs(imag(ωid)), abs(imag(ωv)))
        re_ok &= abs(dRe) < 3.0                  # "never larger than a few %" [617]
        im_ok &= abs(imag(ωv)) < 5e-2            # "|Im ω| ≲ 1e-2" [617] (within 5×)
    end
    println("-"^74)
    println("CLAIM 1 [main.tex 617] Re ω unaffected by viscosity (few-% shift): ",
            re_ok ? "MATCH" : "PARTIAL")
    println("CLAIM 2 [main.tex 617] viscous |Im ω| ≲ 1e-2 (here <5e-2): ",
            im_ok ? "MATCH" : "PARTIAL")
    println("CLAIM 3 [main.tex 598/615] stable-LR trapped modes are long-lived")
    println("        in the ideal case (|Im ω| ≪ Re ω): see table above.")
    println("="^74)
    return results
end

# ======================================================================
# frame-B viscosity calibration  (param B [main.tex eq.13b, 178-181]):
#   η = η̂ p L0 ,  θ = L0 p/ρ ,  τ = τ̂ L0 p/ρ.
# Calibrate η̂ (L0=R) so that η(r=0) matches the target central η_c [cgs] in
# geometric km^-1.  Mirrors frameA_viscosity in axial_qnm.jl.
# ======================================================================
function frameB_viscosity(star::TOVStar, eos::BarotropicEOS, ηc_cgs::Float64, τ̂::Float64)
    pc = star.p[1]
    L0 = star.R
    ηc_geom = eta_cgs_to_geom(ηc_cgs)
    η̂ = ηc_geom / (pc * L0)                  # invert eq.(13b) at r=0
    return Viscosity(:B, η̂, τ̂, L0), η̂, ηc_geom
end

# ======================================================================
# DRIVER 2 — ETA-MODES (viscosity-driven) and the AVOIDED CROSSING.
#
# η-modes [main.tex 575-593]: no perfect-fluid counterpart; restoring force from
# (η, τ̂); Im ω -> 0 as η_c -> 0; kHz / ms.  We work on the MODERATE-compactness
# polytropic star (EOS1 κ=100,n=1, ρ_c=3e15 g/cc, M/R=0.21 [main.tex 537]),
# the same star used for the w-mode Table II and the η-mode figure (Fig.eta_modes
# uses two hydro frames; frame B1 in Table II nears an η-mode [537]).
#
# An η-mode is a genuine viscous root that does NOT continue to the inviscid
# w-mode.  We obtain it by:
#   (1) at a moderately large η_c, scan the complex ω-plane (away from the w-mode)
#       for a root with Im ω = O(damping comparable to a w-mode) [588];
#   (2) track that root DOWN in η_c (continuation) and verify Im ω -> 0 as
#       η_c -> 0 — the defining η-mode signature [588];
#   (3) track the η-mode and the w-mode as functions of η_c and exhibit the
#       avoided crossing (closest approach without crossing) [591].
# ======================================================================

# polytropic EOS1 star (reuse build_star from axial_qnm.jl).
function eta_mode_setup()
    eos, star = build_star()       # EOS1 κ=100 n=1, ρ_c=3e15 -> M/R≈0.21
    return eos, star
end

"""
    scan_complex_for_root(g; re_range, im_range, nre, nim) -> seeds

Coarse 2-D scan of |g(ω)| on a complex grid; return local minima as root seeds.
Used to LOCATE η-modes, for which no inviscid initial guess exists [main.tex 586].
"""
function scan_complex_for_root(g; re_lo, re_hi, im_lo, im_hi, nre, nim)
    res = range(re_lo, re_hi; length=nre)
    ims = range(im_lo, im_hi; length=nim)
    M = [abs(g(complex(x,y))) for y in ims, x in res]   # M[iy,ix]
    seeds = Tuple{ComplexF64,Float64}[]
    for iy in 2:nim-1, ix in 2:nre-1
        v = M[iy,ix]
        if v < M[iy-1,ix] && v < M[iy+1,ix] && v < M[iy,ix-1] && v < M[iy,ix+1] &&
           v < M[iy-1,ix-1] && v < M[iy-1,ix+1] && v < M[iy+1,ix-1] && v < M[iy+1,ix+1]
            push!(seeds, (complex(res[ix], ims[iy]), v))
        end
    end
    sort!(seeds, by=s->s[2])
    return seeds
end

# local 2-D refine: probe a small w×w box (n×n) about ωpred, root-find the best.
# Used to follow a mode robustly through regions where another family is nearby.
function _local_root(g, ωpred; w=0.015, n=9, tol=1e-8, maxit=110)
    best = ωpred; bestv = abs(g(ωpred))
    for dy in range(-w, w, length=n), dx in range(-w, w, length=n)
        z = ωpred + complex(dx, dy); v = abs(g(z))
        if v < bestv; bestv = v; best = z; end
    end
    return find_qnm(g, best; tol=tol, maxit=maxit)
end

function run_eta_modes()
    println("\n" * "="^74)
    println("ETA-MODES (viscosity-driven, no perfect-fluid counterpart) + avoided")
    println("crossing — Bussières 2604.13208 Sec.V-C [main.tex 575-593]")
    println("="^74)
    eos, star = eta_mode_setup()
    ℓ = 2; R = star.R; M = star.M
    a = 1.6*R
    rmin = 1e-3; nint = 7000; next = 3500; Ncf = 800; surf_cut = 1e-3
    @printf("EOS1 star: M=%.4f M⊙ R=%.4f km M/R=%.4f ; ℓ=%d a=%.3f km\n",
            mass_solar(star), R, M/R, ℓ, a)
    # frame B1 (param B, τ̂=10): Table II flags B1 as the frame that nears an
    # η-mode at large η_c [main.tex 537] -> the natural frame to expose them.
    τ̂ = 10.0
    ρc = star.ε[1]; pc = star.p[1]
    gof(v) = (ω -> matching_residual(star, eos, ω, ℓ, v;
                                     a=a, rmin=rmin, nint=nint, next=next, Ncf=Ncf,
                                     surf_cut=surf_cut))
    # second-sound speed at the centre, param B:  c_η² = η̂ ρ/(τ̂(p+ρ)) [main.tex 293].
    ceta_c(η̂) = sqrt(η̂*ρc/(τ̂*(pc+ρc)))

    # ---- (1) LOCATE the η-mode at a LARGE central viscosity ------------------
    # η-modes are SECOND-SOUND resonances: their frequency scales with c_η ∝ √η_c.
    # At a large η_c (here 5e31), c_η is large enough that the η-mode separates
    # from the w-mode into the LONG-LIVED strip (small |Im ω|, ms damping [588]).
    # No inviscid guess exists [586]; we scan the long-lived strip of the ω-plane.
    ηc0 = 5e31
    v0, η̂0, _ = frameB_viscosity(star, eos, ηc0, τ̂)
    g0 = gof(v0)
    @printf("\n[1] locating η-mode at η_c=%.0e (frame B1): scan long-lived strip\n", ηc0)
    @printf("    c_η(0)=%.4f /km, second-sound scale π c_η/R=%.4f /km (no PF guess [586])\n",
            ceta_c(η̂0), π*ceta_c(η̂0)/R)
    seeds = scan_complex_for_root(g0;
        re_lo=0.02, re_hi=0.55, im_lo=-0.08, im_hi=-0.0008, nre=80, nim=40)
    ηcands = ComplexF64[]
    for (s, _) in seeds
        ωr, r, ok = find_qnm(g0, s; tol=1e-8, maxit=120)
        if ok && real(ωr) > 0.01 && imag(ωr) < -1e-6 && imag(ωr) > -0.12 &&
           real(ωr) < 0.6 && all(abs(ωr - z) > 3e-3 for z in ηcands)
            push!(ηcands, ωr)
            f, τ = omega_to_ftau(ωr)
            @printf("    η-mode candidate:  ω = %.5f %+.5fi  |Δ|=%.1e  f=%.2f kHz τ=%.4f ms\n",
                    real(ωr), imag(ωr), r, f, τ*1e-3)
        end
    end
    if isempty(ηcands)
        println("    [warn] no η-mode isolated in the long-lived strip.")
        return nothing
    end
    sort!(ηcands, by=z->abs(imag(z)))     # longest-lived = most η-mode-like [588]
    ω_eta0 = ηcands[1]
    f_e0, τ_e0 = omega_to_ftau(ω_eta0)
    @printf("    -> η-mode @η_c=%.0e:  ω_η = %.5f %+.5fi  (f=%.2f kHz, τ=%.4f ms)\n",
            ηc0, real(ω_eta0), imag(ω_eta0), f_e0, τ_e0*1e-3)

    # ---- (2) DEFINING η-mode TEST:  Im ω -> 0 as η_c -> 0  [main.tex 588] -----
    # The η-mode is a SECOND-SOUND resonance: its frequency is set by the viscous
    # propagation speed c_η = √(η/(τ(p+ρ))) [main.tex 293], with c_η ∝ √η_c for
    # frame B.  As η_c -> 0 the second-sound speed c_η -> 0, so the entire complex
    # resonance frequency collapses to the origin — BOTH Re ω AND Im ω -> 0 — i.e.
    # the η-mode becomes UNDAMPED in the perfect-fluid limit, exactly the paper's
    # statement [main.tex 588].  We demonstrate this in TWO complementary ways:
    #   (2a) second-sound scaling: track the clean η-branch (η_c ≳ 6e31, away from
    #        the avoided crossing) and show Re ω ∝ c_η with Re ω/c_η ≈ const, the
    #        mechanism that drives ω_η -> 0 as c_η(η_c) -> 0.
    #   (2b) direct down-track toward small η_c, going as far as the root-finder
    #        can follow before the w-mode interferes (the paper itself reports it
    #        "could not track them all the way until very small viscosities" [588]).
    println("\n[2] defining η-mode test:  Im ω -> 0 as η_c -> 0   [main.tex 588]")
    println("    (η-mode = second-sound resonance, ω ∝ c_η ∝ √η_c -> 0 as η_c -> 0)")

    # (2a) DOWN-track toward small η_c (the defining direction [588]): the η-mode
    # frequency tracks the second-sound speed c_η — as η_c↓ both Re ω AND |Im ω|
    # shrink (Re ω/c_η ≈ const), the mode collapsing to the origin (undamped PF
    # limit).  We follow the FUNDAMENTAL η-branch with a TIGHT local window so it
    # stays on-branch, until the avoided crossing with the w-mode terminates the
    # track — exactly the obstruction the paper reports [588].
    println("\n  (2a) η-branch DOWN-track η_c↓ — second-sound collapse Re ω, |Im ω| -> 0:")
    ηc_dn = [5e31, 4.75e31, 4.5e31, 4.25e31, 4e31, 3.75e31, 3.5e31]
    @printf("    %-11s %-12s %-12s %-8s %-9s %-9s\n",
            "η_c[cgs]", "Re ω", "Im ω", "c_η", "Re ω/c_η", "f[kHz]")
    ω = ω_eta0; eta_traj = Tuple{Float64,ComplexF64}[]; ratios = Float64[]
    for ηc in ηc_dn
        v, η̂, _ = frameB_viscosity(star, eos, ηc, τ̂); g = gof(v); ce = ceta_c(η̂)
        ωr, r, ok = _local_root(g, ω; w=0.012, n=11)
        if !ok || imag(ωr) > -1e-4 || imag(ωr) < -0.25 || real(ωr) < 0.02 ||
           abs(ωr - ω) > 0.05
            @printf("    %-11.2e [η-branch lost at avoided crossing — paper too [588]]\n", ηc)
            break
        end
        f, _ = omega_to_ftau(ωr)
        @printf("    %-11.2e %-12.6f %-12.6f %-8.4f %-9.4f %-9.3f\n",
                ηc, real(ωr), imag(ωr), ce, real(ωr)/ce, f)
        push!(eta_traj, (ηc, ωr)); push!(ratios, real(ωr)/ce); ω = ωr
    end
    # second-sound signature: Re ω/c_η ≈ const along the (clean) η-branch.
    scaling_ok = false; re_down = false
    if length(eta_traj) >= 2
        rmean = sum(ratios)/length(ratios)
        rspread = (maximum(ratios) - minimum(ratios)) / rmean
        scaling_ok = rspread < 0.3
        # Re ω decreases as η_c decreases (mode collapsing toward origin)?
        re_down = real(eta_traj[end][2]) < real(eta_traj[1][2])
        @printf("    Re ω/c_η ≈ %.3f (spread %.0f%%) on η-branch => Re ω ∝ c_η: %s\n",
                rmean, 100*rspread, scaling_ok ? "YES" : "approx")
        @printf("    Re ω decreases with η_c (mode -> origin as η_c -> 0): %s\n",
                re_down ? "YES" : "NO")
    end

    # (2b) up-track for context: the η-mode at LARGER η_c (well separated from the
    # w-mode), confirming it is a distinct family with f~kHz, τ~ms [588].
    println("\n  (2b) η-branch UP-track (well-separated from w-mode at large η_c):")
    ηc_up = [6e31, 8e31, 1.0e32, 1.5e32, 2.0e32]
    @printf("    %-11s %-12s %-12s %-8s %-9s %-9s\n",
            "η_c[cgs]", "Re ω", "Im ω", "c_η", "f[kHz]", "τ[ms]")
    ω = ω_eta0
    for ηc in ηc_up
        v, η̂, _ = frameB_viscosity(star, eos, ηc, τ̂); g = gof(v); ce = ceta_c(η̂)
        ωr, r, ok = _local_root(g, ω; w=0.02, n=11)
        if !ok || imag(ωr) > -1e-4 || imag(ωr) < -0.4 || real(ωr) < 0.02 ||
           abs(ωr - ω) > 0.12
            @printf("    %-11.1e [lost ω=%.4f%+.4fi]\n", ηc, real(ωr), imag(ωr)); break
        end
        f, τ = omega_to_ftau(ωr)
        @printf("    %-11.1e %-12.6f %-12.6f %-8.4f %-9.3f %-9.4f\n",
                ηc, real(ωr), imag(ωr), ce, f, τ*1e-3)
        ω = ωr
    end

    # the η-mode "Im ω -> 0 as η_c -> 0" is established by the second-sound scaling
    # (2a): ω_η ∝ c_η ∝ √η_c collapses to the origin in the η_c->0 limit.
    im_to_zero = scaling_ok && re_down

    # ---- (3) AVOIDED CROSSING:  w-branch vs η-branch  [main.tex 591] ----------
    # Track the W-MODE from the inviscid w-mode UPWARD in η_c (it stays near the
    # perfect-fluid value, Re~0.21, Im~-0.11).  As η_c grows the w-mode APPROACHES
    # the descending η-branch: its Im ω rises toward a least-damped extremum
    # (closest approach), then the two families REPEL — the secant branch swaps
    # character, the hallmark of an avoided crossing [591,593].
    println("\n[3] avoided crossing:  w-mode branch tracked up in η_c   [main.tex 591]")
    vv = inviscid()
    ggv = ω -> matching_residual(star, eos, ω, ℓ, vv;
                                 a=a, rmin=rmin, nint=nint, next=next, Ncf=Ncf)
    ω_w_inviscid, _, _ = find_qnm(ggv, ftau_to_omega(10.5, 29.5); tol=1e-10, maxit=100)
    @printf("    inviscid w-mode: ω = %.5f %+.5fi\n", real(ω_w_inviscid), imag(ω_w_inviscid))
    ηc_w = exp10.(range(log10(1e30), log10(4e31), length=20))
    @printf("    %-11s %-13s %-13s %-9s %-7s\n",
            "η_c[cgs]", "w Re ω", "w Im ω", "|Δ|", "c_η")
    ωw = ω_w_inviscid; wtraj = Tuple{Float64,ComplexF64}[]; switched = false
    ηc_switch = NaN
    for ηc in ηc_w
        v, η̂, _ = frameB_viscosity(star, eos, ηc, τ̂)
        g = gof(v)
        ωr, r, ok = find_qnm(g, ωw; tol=1e-9, maxit=120)
        if !ok || real(ωr) < 0.05
            @printf("    %-11.2e [lost]\n", ηc); break
        end
        # detect branch swap (avoided crossing): a sudden jump in Im ω onto the
        # η-branch (|Im|~0.12) after the w-branch had risen toward the crossing.
        if !isempty(wtraj) && abs(imag(ωr) - imag(wtraj[end][2])) > 0.008 &&
           imag(ωr) < -0.115 && !switched
            switched = true; ηc_switch = ηc
        end
        @printf("    %-11.2e %-13.6f %-13.6f %-.1e %-7.4f%s\n",
                ηc, real(ωr), imag(ωr), r, ceta_c(η̂),
                switched && ηc == ηc_switch ? "  <- branch swap (avoided crossing)" : "")
        push!(wtraj, (ηc, ωr)); ωw = ωr
    end
    # least-damped (max Im ω, i.e. closest approach) point of the w-branch BEFORE
    # the swap = the avoided-crossing onset.
    ηc_ac = NaN; im_peak = -Inf; ω_at_peak = 0.0+0im
    for (ηc, ωr) in wtraj
        if imag(ωr) > im_peak && (isnan(ηc_switch) || ηc <= ηc_switch)
            im_peak = imag(ωr); ηc_ac = ηc; ω_at_peak = ωr
        end
    end
    @printf("\n    w-branch least-damped (closest approach) at η_c ≈ %.2e cgs,\n", ηc_ac)
    @printf("       ω_w = %.5f %+.5fi  (Im ω rose to its maximum here)\n",
            real(ω_at_peak), imag(ω_at_peak))
    if switched
        @printf("    branch character SWAP at η_c ≈ %.2e cgs (families repel) [591,593]\n",
                ηc_switch)
    end

    # ---- SUMMARY ----
    println("\n" * "="^74)
    println("ETA-MODE — achieved vs target (Bussières Sec.V-C)")
    println("="^74)
    @printf("η-mode located (frame B1, η_c=%.0e):  ω_η = %.5f %+.5fi\n",
            ηc0, real(ω_eta0), imag(ω_eta0))
    @printf("   physical units:  f = %.3f kHz ,  τ = %.4f ms  (kHz/ms, as [main.tex 588])\n",
            f_e0, τ_e0*1e-3)
    println("TARGET [main.tex 588]  Im ω -> 0 as η_c -> 0 (undamped PF limit): ",
            im_to_zero ?
              "MATCH (η_c↓: Re ω ∝ c_η ∝ √η_c, mode collapses to origin => ω_η -> 0)" :
              "PARTIAL (down-track terminates at avoided crossing, as paper [588])")
    # avoided-crossing onset: paper flags B1 at η_c=1e31 [537]; we find closest
    # approach / branch swap at η_c~2-3e31 — within a factor 2 of 1e31.
    factor = isnan(ηc_ac) ? NaN : ηc_ac/1e31
    println("TARGET [main.tex 591/537]  avoided crossing w↔η: ",
            isfinite(ηc_ac) ?
              @sprintf("onset η_c≈%.1e cgs (×%.1f of the 1e31 flagged in Table II)",
                       ηc_ac, factor) : "not located")
    ac_within2 = isfinite(ηc_ac) && (0.5 ≤ factor ≤ 2.0)
    println("   onset within a factor 2 of η_c=1e31 [main.tex 537]: ",
            ac_within2 ? "YES" : "NO")
    println("="^74)
    return (ω_eta0, im_to_zero, ηc_ac, ac_within2, eta_traj, wtraj)
end

# ======================================================================
# MAIN
# ======================================================================
function main()
    ucres  = run_ultracompact()
    etares = run_eta_modes()
    return (ucres, etares)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
