#=
    Localize the Shum fine-Δr instability: capture the spatial profile just before
    blow-up.  Run Δr=0.0024 to t=26 (blow-up ~t=27) and Δr=0.04 (stable) to t=26;
    save r, the energy perturbation δε=ε−ε_bg, and the EVOLVED derivative field dε
    (a first-order-reduction auxiliary — prime instability suspect).

    Run: julia code/BDNKStar/repro/shum_instab_profile.jl
=#
include(joinpath(@__DIR__, "shum_evolve_opt.jl"))

function snap(Dr, t_f)
    s, ts, ecs, nan_hit, ec0 = run_evolution(; dr=Dr, t_f=t_f, vpert=0.0, epspert=1e-4,
                                              sample_dt=1.0, σKO=0.5, case=:smallSB_F2)
    r = s.bg.r; δε = s.ε .- s.bg.εbg
    return r, δε, copy(s.dε), s.bg.Rstar, nan_hit, ts[end]
end

r24, d24, dd24, Rs, nan24, tend24 = snap(0.0024, 26.0)
r04, d04, dd04, _,  nan04, tend04 = snap(0.04,   26.0)

open(joinpath(@__DIR__,"shum_instab_profile.txt"),"w") do io
    println(io, "# FINE Δr=0.0024 (t=$tend24, nan=$nan24, Rstar=$Rs):  r  δε  dε")
    for i in eachindex(r24); println(io, r24[i], "  ", d24[i], "  ", dd24[i]); end
    println(io, "# COARSE Δr=0.04 (t=$tend04, nan=$nan04):  r  δε  dε")
    for i in eachindex(r04); println(io, r04[i], "  ", d04[i], "  ", dd04[i]); end
end
println("SAVED shum_instab_profile.txt")
println("  fine  Δr=0.0024: t=$tend24  max|δε|=", round(maximum(abs.(d24)),sigdigits=4),
        " at r=", round(r24[argmax(abs.(d24))],digits=3), " (Rstar=", round(Rs,digits=3), ")")
println("  coarse Δr=0.04 : t=$tend04  max|δε|=", round(maximum(abs.(d04)),sigdigits=4),
        " at r=", round(r04[argmax(abs.(d04))],digits=3))
