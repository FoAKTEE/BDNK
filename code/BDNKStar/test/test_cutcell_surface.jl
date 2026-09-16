using Test
using BDNKStar
using LinearAlgebra: svd, eigvals, pinv
using BDNKStar.NonRadialModes: nonradial_cowling_spectrum

# CUT-CELL SURFACE TREATMENT (`cut=true`) in BOTH Cartesian engines.
#
# WHY THIS TEST EXISTS. Everything the cut-cell path buys — the ideal f-mode at a
# few tenths of a percent, the converged viscous damping, and a FULL-STAR run with
# no excision — sits behind a default-off flag, so the rest of the suite never
# touches it. Without this file those results have no regression protection at all
# and a refactor could silently kill them while the suite stayed green.
#
# WHAT IT PINS.
#  1. GEOMETRY. κ (volume fraction) and the six face apertures both reduce to 1D
#     quadratures of one disc–rectangle area, so accuracy is set by quadrature
#     order, not by an interface reconstruction. Two independent checks: Σκ dx³ vs
#     4πR³/3, and the CLOSURE condition — for a closed cut cell the coordinate-face
#     area vectors must be balanced by the boundary face, so their vector-sum
#     magnitude is the embedded interface area and must sum to 4πR².
#  2. THE FIX IS REAL. The masked scheme must be several percent off while cut
#     cells are sub-percent, on the same star at the same resolution.
#  3. VISCOUS DAMPING against the independent dissipation integral (0.08746 per
#     unit η̂ for the ℓ=2 f-mode, test_cross_method), same η₀=η̂·ε₀ normalisation.
#  4. FULL STAR, NO EXCISION for the full-frame BDNK engine — the thing the
#     CowlingBDNK3D docstring long recorded as unachieved after six failed
#     regularizers.
#
# Damping is read with a MATRIX PENCIL, not `damping_rate`: the latter returns NaN
# below 3 peaks, and the f-mode here is weakly damped over a short record. The
# pencil also separates the f-mode from p-mode contamination, which biases a
# periodogram peak. Resolutions are deliberately modest to keep the suite fast;
# the converged numbers live in the module docstrings.
@testset "Cut-cell surface treatment (Cartesian engines)" begin
    eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2
    M2c = BDNKStar.Units.Msun_to_km * BDNKStar.Units.kHz_to_km
    feig, _, _ = nonradial_cowling_spectrum(eos, εc; l=2, nmodes=3, N=1200, nscan=300)
    CE = BDNKStar.CowlingEvolve3D
    CB = BDNKStar.CowlingBDNK3D

    # ---- 1. geometry ------------------------------------------------------
    s32 = build_star3d(eos, εc; N=32, Lfac=1.20)
    ecut = setup_evo3d(s32; σ_ko=0.004, cut=true)
    dx = s32.grid.dx
    @test all(0 .<= ecut.κ .<= 1)
    @test isapprox(sum(ecut.κ)*dx^3, 4π*s32.R^3/3; rtol=1e-3)   # volume
    area = 0.0
    for I in CartesianIndices(ecut.κ)
        0 < ecut.κ[I] < 1 || continue
        sx = (ecut.axp[I]-ecut.axm[I]); sy = (ecut.ayp[I]-ecut.aym[I]); sz = (ecut.azp[I]-ecut.azm[I])
        area += sqrt(sx^2+sy^2+sz^2)*dx^2
    end
    @test isapprox(area, 4π*s32.R^2; rtol=5e-3)                 # closure ⇒ interface area

    # ---- helpers ----------------------------------------------------------
    fpeak(ts,q) = begin
        wide = range(0.3, 9.0; length=6000)
        iw = argmax(CE.periodogram(ts, q, wide.*M2c))
        fb = range(max(0.3, wide[iw]-0.3), wide[iw]+0.3; length=4000)
        fb[argmax(CE.periodogram(ts, q, fb.*M2c))]
    end
    function pencil_f(q, dts)                     # damped-exponential fit
        M = length(q); L = div(M,3)
        H = [q[i+j] for i in 1:(M-L), j in 0:(L-1)]
        U,S,V = svd(H); p = min(8, count(S .> S[1]*1e-11)); Vp = V[:,1:p]
        z = eigvals(pinv(Vp[1:end-1,:])*Vp[2:end,:])
        out = Tuple{Float64,Float64}[]
        for zi in z
            abs(zi) < 1e-12 && continue
            sm = log(complex(zi))/dts; f = abs(imag(sm))/(2π*M2c)
            0.05 < f && push!(out, (f, -real(sm)))
        end
        out
    end
    function run_evo(N, ηh, cut; T=500.0)
        s = build_star3d(eos, εc; N=N, Lfac=1.20)
        # explicit on BOTH branches: cut=true is now the package default, so the masked
        # (rigid-wall) comparison run must ask for cut=false by name
        e = setup_evo3d(s; σ_ko=0.004, η̂=ηh, cut=cut)
        dt = 0.30*s.grid.dx
        st = EvolState(N); seed_l2!(st, e; A=1e-3)
        evolve3d!(st, e; dt=dt, nsteps=round(Int,T/dt), sample=6)
    end

    # ---- 2. cut cells fix the ideal f-mode --------------------------------
    ts_c, q_c, _ = run_evo(32, 0.0, true)
    ts_m, q_m, _ = run_evo(32, 0.0, false)
    @test all(isfinite, q_c) && all(isfinite, q_m)
    err_c = 100*(fpeak(ts_c,q_c) - feig[1])/feig[1]
    err_m = 100*(fpeak(ts_m,q_m) - feig[1])/feig[1]
    @test abs(err_c) < 2.0          # measured −0.58% at N=32 (−0.12% at N=44)
    @test abs(err_m) > 4.0          # measured −6.94%: the stair-stepped wall
    @test abs(err_c) < abs(err_m)/3 # the fix is large, not marginal

    # ---- 3. viscous damping vs the dissipation integral -------------------
    ηs = [0.0, 0.04, 0.08]; γs = Float64[]
    for ηh in ηs
        ts, q, _ = run_evo(32, ηh, true; T=700.0)
        @test all(isfinite, q)
        md = [m for m in pencil_f(q, ts[2]-ts[1]) if 1.0 < m[1] < 3.0]
        @test !isempty(md)
        push!(γs, md[argmin([abs(m[1]-1.88) for m in md])][2])
    end
    @test issorted(γs)                                   # viscosity damps, monotonically
    slope = (length(ηs)*sum(ηs.*γs) - sum(ηs)*sum(γs)) /
            (length(ηs)*sum(abs2,ηs) - sum(ηs)^2)
    @test 0.6 < slope/0.08746 < 1.5    # measured ~105% at N=32; window covers the
                                       # κ-clipping systematic and the Newtonian-vs-
                                       # relativistic reference-form ambiguity (~5–8%)

    # ---- 4. full-frame BDNK: FULL STAR, NO EXCISION -----------------------
    s24 = build_star3d(eos, εc; N=24, Lfac=1.20)
    eb = setup_bdnk3d(s24; η̂=0.05, ζ̂=0.0, τ̂=2.0, σ_ko=0.004,
                      cut=true, excise_frac=0.0)
    @test eb.cut && eb.wfloor == 0.0                     # nothing excised
    @test count(>(0), eb.κ) > count(s24.interior)        # cut cells reach past the mask
    stb = BDNKState(24); CB.seed_bdnk_l2!(stb, eb; A=1e-3)
    tsb, qb, _ = CB.evolve_bdnk3d!(stb, eb; dt=0.20*s24.grid.dx,
                                   nsteps=round(Int, 500.0/(0.20*s24.grid.dx)), sample=6)
    @test all(isfinite, qb)
    @test maximum(abs, qb) < 10*maximum(abs, qb[1:20])   # bounded, no surface blow-up
    @test abs(100*(fpeak(tsb,qb) - feig[1])/feig[1]) < 5.0
end
