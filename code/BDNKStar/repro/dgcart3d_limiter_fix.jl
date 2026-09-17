# ======================================================================
# repro/dgcart3d_limiter_fix.jl — THE EQUILIBRIUM-PRESERVING LIMITER OF THE 3+1D DG ENGINE
#
# RESULT (data: repro/data/dgcart3d_limiter_drift.csv, dgcart3d_fmode_scan.csv, 2026-09-17;
# write-up: VALIDATION.md §7.8, paper Sec. VI):
#   The DGCart3D background used to erode (ρ_c −10% in 4.6 f-mode periods, −31 to −43% by
#   3000 M⊙; absolute f −16%). Attribution: limiter OFF → unseeded star static to 1e-12;
#   limiter ON → the same unseeded star thrown into a ±3% radial oscillation in 50 M⊙. The
#   mean-based Zhang–Shu limiter flattens the equilibrium's own intra-element profile in the
#   surface elements whenever it engages. Fix: limit (and dissipate) on the DEVIATION from the
#   stored equilibrium — reference R = U_eq + (mean deviation shared ∝ background), carrying
#   the kinetic energy of its own momentum (src/dg/DGCart3D.jl, _limit_wb!, _rus5wb).
#
#   Fixed engine, K=6, p=2, seed v^r = A(r/R)Y20:
#     unseeded, limiter ON     : |ρ_c/ρ_c0−1| ≤ 2e-13, |v| ≤ 4e-11 over 250 M⊙   (was −10%)
#     A=1e-2, t=600 / 1000     : ρ_c −5.8e-4 / −7.3e-4 (decelerating), star elements never
#                                on the fallback path, max|v| 0.003 after the transient
#     A=1e-3, t=600            : ρ_c −1.3e-5  (∝ A²: thermalised mode energy)
#     f-mode, K=6/8/10, T=600–1300: 1.69–1.85 kHz (−2 to −10%, estimator scatter of the same
#                                size, non-monotone in K), envelope DECAYING (pre-fix: 1.57–1.61
#                                kHz, growing). The residual is the staircase/atmosphere surface
#                                of a nodal DG star, not the limiter — see VALIDATION.md §7.8.
#
# Runtime: part 1 ≈ 1 min; part 2 ≈ 8 min single-threaded; part 3 (K-scan) ≈ 1 h on 5 threads.
# Set BDNK_DG_FULL=1 for parts 2–3.
# ======================================================================
using BDNKStar, Printf
using BDNKStar: Msun_to_km, kHz_to_km
eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2
full = get(ENV, "BDNK_DG_FULL", "0") == "1"
datadir = joinpath(@__DIR__, "data")

# ---- part 1: fixed point --------------------------------------------------------------
eng, st = setup_dgcart3d(eos, εc; Kx=6, Ky=6, Kz=6, p=2, cfl=0.2)
D0 = copy(st.D); BDNKStar.DGCart3D._limit!(st, eng)
println("limiter leaves U_eq bitwise unchanged: ", st.D == D0)
ρc0 = dgcart3d_central_density(st, eng)
ts, q2, ρc = evolve_dgcart3d!(st, eng; tmax=250.0, sample_dt=50.0)
mn = dgcart3d_prim_minmax(st)
@printf("unseeded, limiter ON, t=250: max|ρc/ρc0−1| = %.1e   max|v| = %.1e\n", maximum(abs.(ρc./ρc0 .- 1)), max(abs(mn[5]),abs(mn[6])))

full || exit()
# ---- part 2: seeded drift + limiter census ---------------------------------------------
open(joinpath(datadir, "dgcart3d_limiter_drift.csv"), "w") do io
    println(io, "A,K,t,rhoc_drift,sumD_drift,vmax,elems_free,elems_scaled,elems_fallback,star_fallback")
    for (A, T) in ((1e-2, 1000.0), (1e-3, 600.0))
        eng, st = setup_dgcart3d(eos, εc; Kx=6, Ky=6, Kz=6, p=2, cfl=0.2)
        seed_dgcart3d_l2!(st, eng; A=A)
        ρc0 = dgcart3d_central_density(st, eng); ΣD0 = sum(st.D); t = 0.0
        while t < T - 1e-9
            ts, q2, ρc = evolve_dgcart3d!(st, eng; tmax=100.0, sample_dt=2.0); t += 100
            mn = dgcart3d_prim_minmax(st); c = dgcart3d_limiter_census(st, eng)
            @printf(io, "%.0e,6,%.0f,%.3e,%.2e,%.3f,%d,%d,%d,%d\n", A, t, ρc[end]/ρc0-1, sum(st.D)/ΣD0-1,
                    max(abs(mn[5]),abs(mn[6])), c.free, c.scaled, c.fallback, c.star_fallback); flush(io)
            @printf("A=%.0e t=%4.0f  ρc/ρc0−1=%+.3e  census free/scaled/fallback=%d/%d/%d  star_fallback=%d\n",
                    A, t, ρc[end]/ρc0-1, c.free, c.scaled, c.fallback, c.star_fallback)
        end
    end
end

# ---- part 3: f-mode resolution scan --------------------------------------------------
ν_ref = 1.88291*kHz_to_km*Msun_to_km
open(joinpath(datadir, "dgcart3d_fmode_scan.csv"), "w") do io
    println(io, "K,A,T,f_pgram,f_pencil,agree_pct,df_kHz,stable,ratio,rhoc_drift,gamma")
    for (K, A, T) in ((6, 1e-2, 600.0), (6, 1e-3, 600.0), (8, 1e-2, 600.0), (10, 1e-2, 600.0), (8, 1e-2, 1300.0))
        eng, st = setup_dgcart3d(eos, εc; Kx=K, Ky=K, Kz=K, p=2, cfl=0.2)
        seed_dgcart3d_l2!(st, eng; A=A)
        ρc0 = dgcart3d_central_density(st, eng)
        ts, q2, ρc = evolve_dgcart3d!(st, eng; tmax=T, sample_dt=2.0)
        a = analyze_qnm(ts, q2; ν_ref=ν_ref)
        @printf(io, "%d,%.1e,%.0f,%.6f,%.6f,%.3f,%.4f,%s,%.3g,%.4e,%.4e\n", K, A, T, a.f_kHz, a.f_pencil_kHz,
                a.agree_pct, a.df_kHz, a.stable, a.envelope_ratio, ρc[end]/ρc0-1, a.γ); flush(io)
        @printf("K=%2d A=%.0e T=%4.0f: f=%.5f/%.5f kHz (%+.2f%% vs 1D)  df/f=%.1f%%  ρc drift=%+.2e  stable=%s\n",
                K, A, T, a.f_kHz, a.f_pencil_kHz, 100*(a.f_kHz/1.88291-1), 100*a.df_kHz/a.f_kHz, ρc[end]/ρc0-1, a.stable)
    end
end
