# ======================================================================
# repro/dgball3d_hkt.jl — THE 3+1D CUBED-SPHERE DG STAR (DGBall3D), the 3D method of
# Hébert–Kidder–Teukolsky (2018) Sec. VI.B on the anchor star, coarse grid.
#
# Grid: ball_grid(R; nt=2, p_t=3, p_int=3, p_surf=2, p_ext=3): rounded cube + 18 cubed-sphere
# shells (2 centre, 7 interior, 5 quadratic surface shells with the surface on a shell boundary,
# 5 exterior to 3R), 464 elements, 27 776 nodes, Δt ≈ 0.10 M⊙, ≈0.03 s/step on 8 threads.
# Runs (all with the :wb limiter on the surface shells and the paper's momentum filter):
#   static_WB     unseeded, static-residual subtraction         → fixed point to 1e-14
#   static_noWB   unseeded, the paper's scheme (no subtraction) → settles, compare their B1 err[D̃]
#   l0_WB         v = 1e-3 sin(πr/R) n̂, subtraction             → radial spectrum vs F, H1, H2
#   l2_WB         v = 1e-3 sin(πr/R) Y20 n̂, subtraction         → f-mode vs 1.88291 kHz
#   l2_noWB       same without subtraction
#   l2_WB_nt3     nt=3 (999 elements)                           → resolution check
# Output: repro/data/dgball3d_hkt_summary.csv, repro/data/dgball3d_hkt_series_<tag>.csv.
# Runtime ≈ 2.5 h on 8 threads (set BDNK_BALL_FULL=1 for the nt=3 run).
# ======================================================================
using BDNKStar, Printf
using BDNKStar: Msun_to_km, kHz_to_km
eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2
datadir = joinpath(@__DIR__, "data")
ν_ref = 1.88291*kHz_to_km*Msun_to_km
full = get(ENV, "BDNK_BALL_FULL", "0") == "1"
runs = [ # tag, nt, seed, wellbalanced, T
    ("static_WB",   2, :none,   true,  200.0),
    ("static_noWB", 2, :none,   false, 2000.0),
    ("l0_WB",       2, :radial, true,  2000.0),
    ("l2_WB",       2, :l2,     true,  2000.0),
    ("l2_noWB",     2, :l2,     false, 2000.0),
]
full && push!(runs, ("l2_WB_nt3", 3, :l2, true, 2000.0))
open(joinpath(datadir, "dgball3d_hkt_summary.csv"), "w") do io
    println(io, "tag,nt,seed,wellbalanced,T,elements,nodes,dt,seconds,errD_end,rhoc_drift_end,Mb_drift_end,vatm_max,resid_interior,resid_surface,f_pgram_kHz,f_pencil_kHz,df_over_f,gamma,stable,peak1_kHz,peak2_kHz,peak3_kHz")
    for (tag, nt, seed, wb, T) in runs
        eng, st = setup_dgball3d(eos, εc; nt=nt, p_t=3, p_int=3, p_surf=2, p_ext=3, wellbalanced=wb, limiter=:wb)
        res = dgball3d_static_residual(st, eng)
        seed == :radial && seed_dgball3d_radial!(st, eng; A=1e-3)
        seed == :l2 && seed_dgball3d_l2!(st, eng; A=1e-3)
        rec = evolve_dgball3d!(st, eng; tmax=T, sample_dt=1.0)
        open(joinpath(datadir, "dgball3d_hkt_series_$(tag).csv"), "w") do ios
            println(ios, "t,rhoc,errD,Mb,q00,q20,q22,vatm")
            for i in eachindex(rec.t); @printf(ios, "%.3f,%.10e,%.6e,%.10e,%.10e,%.10e,%.10e,%.4e\n", rec.t[i], rec.ρc[i], rec.errD[i], rec.Mb[i], rec.q00[i], rec.q20[i], rec.q22[i], rec.vatm[i]); end
        end
        fp = NaN; fpen = NaN; dff = NaN; γ = NaN; stable = false; pk = Float64[]
        if seed == :l2
            a = analyze_qnm(rec.t, rec.q20; ν_ref=ν_ref); fp = a.f_kHz; fpen = a.f_pencil_kHz; dff = a.df_kHz/a.f_kHz; γ = a.γ; stable = a.stable
        elseif seed == :radial
            pk, _, _ = dgstarhp_spectrum(rec.t, rec.ρc; npeaks=3)
        end
        pk3 = vcat(pk, fill(NaN, 3-length(pk)))
        @printf(io, "%s,%d,%s,%s,%.0f,%d,%d,%.5f,%.0f,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.5f,%.5f,%.4f,%.3e,%s,%.4f,%.4f,%.4f\n", tag, nt, seed, wb, T, eng.K, eng.Ntot, rec.dt, rec.seconds,
                rec.errD[end], rec.ρc[end]/rec.ρc[1]-1, rec.Mb[end]/rec.Mb[1]-1, maximum(rec.vatm), res[:interior], res[:surface], fp, fpen, dff, γ, stable, pk3...); flush(io)
        @printf("%-12s K=%d N=%d dt=%.4f [%.0f s]: errD=%.2e ρc drift=%+.2e Mb drift=%+.2e vatm_max=%.2e  f=%.4f/%.4f kHz  peaks=%s\n", tag, eng.K, eng.Ntot, rec.dt, rec.seconds,
                rec.errD[end], rec.ρc[end]/rec.ρc[1]-1, rec.Mb[end]/rec.Mb[1]-1, maximum(rec.vatm), fp, fpen, join([@sprintf("%.3f",x) for x in pk], ",")); flush(stdout)
    end
end
println("BALL DONE")
