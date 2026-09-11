using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using Printf
using Statistics

eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2

# CFL: keep dt/dr fixed across resolutions. Reference (Nr=64): dr=R/64, dt=0.005.
# Build one star to get R for the CFL scale.
sref = build_sphstar(eos, εc; Nr=64, Nθ=12)
R = sref.R
cfl = 0.005 / (R/64)   # dt = cfl * dr  -> matched Courant number
@printf("R = %.4f ; CFL number dt/dr = %.4f\n", R, cfl)

# --- Γ extraction -----------------------------------------------------------
# en(t) ≈ E0 exp(-2 Γ_amp t). Γ_E = -d ln E / dt ; Γ_amp = Γ_E/2.
# Fit ln E over a clean late-time window: drop the first `skip` fraction
# (seed-dephasing transient) and the last small tail.
function gamma_amp(ts, en; skip=0.40, tail=0.0)
    n = length(en)
    @assert all(isfinite, en)
    @assert all(>(0), en)
    i0 = max(2, round(Int, skip*n))
    i1 = max(i0+3, round(Int, (1-tail)*n))
    t = ts[i0:i1]; y = log.(en[i0:i1])
    tb = mean(t); yb = mean(y)
    slope = sum((t .- tb).*(y .- yb)) / sum((t .- tb).^2)   # = -Γ_E
    ΓE = -slope
    # R^2 of fit
    yhat = yb .+ slope.*(t .- tb)
    ss_res = sum((y .- yhat).^2); ss_tot = sum((y .- yb).^2)
    r2 = ss_tot > 0 ? 1 - ss_res/ss_tot : 1.0
    (Γamp = ΓE/2, ΓE = ΓE, r2 = r2)
end

# --- one run ----------------------------------------------------------------
# returns (Γamp, Eend/E0, max/E0, r2, finite)
function run_case(; Nr, Nθ, η̂, σ_ko, dt, T=400.0, n_seed=3, kw...)
    s  = build_sphstar(eos, εc; Nr=Nr, Nθ=Nθ)
    e  = setup_sphbdnk(s; η̂=η̂, σ_ko=σ_ko, kw...)
    st = SphBDNKState(s.grid.Nr, s.grid.Nθ); seed_sphbdnk_n!(st, e, n_seed; A=1e-3)
    nsteps = round(Int, T/dt)
    sample = max(1, nsteps ÷ 300)
    ts, _, en = evolve_sphbdnk!(st, e; dt=dt, nsteps=nsteps, sample=sample)
    fin = all(isfinite, en) && all(>(0), en)
    if !fin
        return (Γamp=NaN, ratio=NaN, mx=NaN, r2=NaN, fin=false)
    end
    g = gamma_amp(ts, en)
    (Γamp=g.Γamp, ratio=en[end]/en[1], mx=maximum(en)/en[1], r2=g.r2, fin=true)
end

dt_for(Nr) = cfl * (R/Nr)

println("\n#########################################################")
println("# TEST 1 — η̂=0 BASELINE (numerical contribution)")
println("# Γ₀(σ_ko, resolution). Γamp<0 means GROWTH.")
println("#########################################################")
sko_list = [0.005, 0.01, 0.02, 0.04]
res_list = [(64,12),(96,18),(128,24)]
println(@sprintf("%-12s %-10s %-12s %-12s %-10s %-8s", "Nr,Nθ", "σ_ko", "Γamp", "Eend/E0", "max/E0", "R2"))
T1 = Dict{Tuple{Int,Int,Float64},Float64}()
for (Nr,Nθ) in res_list
    dt = dt_for(Nr)
    for σ in sko_list
        r = run_case(Nr=Nr, Nθ=Nθ, η̂=0.0, σ_ko=σ, dt=dt, T=400.0)
        T1[(Nr,Nθ,σ)] = r.Γamp
        @printf("%-12s %-10.3f % .4e  % .4e  % .3e  %.3f\n",
                "($Nr,$Nθ)", σ, r.Γamp, r.ratio, r.mx, r.r2)
    end
end

println("\n#########################################################")
println("# TEST 2 — σ_ko INDEPENDENCE at η̂=0.05 (decisive)")
println("#########################################################")
println(@sprintf("%-10s %-12s %-12s %-10s", "σ_ko", "Γamp", "Eend/E0", "R2"))
g2 = Float64[]
for σ in sko_list
    r = run_case(Nr=64, Nθ=12, η̂=0.05, σ_ko=σ, dt=dt_for(64), T=400.0)
    push!(g2, r.Γamp)
    @printf("%-10.3f % .4e  % .4e  %.3f\n", σ, r.Γamp, r.ratio, r.r2)
end
# slope dΓ/dσ_ko
let x=sko_list, y=g2
    xb=mean(x); yb=mean(y); sl=sum((x.-xb).*(y.-yb))/sum((x.-xb).^2)
    @printf("dΓ/dσ_ko (η̂=0.05) = % .4e ; spread max-min = %.4e\n", sl, maximum(y)-minimum(y))
end
# compare with η̂=0 σ-sensitivity at same resolution
let x=sko_list, y=[T1[(64,12,σ)] for σ in sko_list]
    xb=mean(x); yb=mean(y); sl=sum((x.-xb).*(y.-yb))/sum((x.-xb).^2)
    @printf("dΓ/dσ_ko (η̂=0  )  = % .4e ; spread max-min = %.4e\n", sl, maximum(y)-minimum(y))
end

println("\n#########################################################")
println("# TEST 3 — RESOLUTION CONVERGENCE at η̂=0.05, σ_ko=0.02")
println("#########################################################")
println(@sprintf("%-12s %-10s %-12s %-12s %-10s", "Nr,Nθ", "dt", "Γamp", "Eend/E0", "R2"))
res3 = [(48,9),(64,12),(96,18),(128,24)]
g3 = Float64[]
for (Nr,Nθ) in res3
    dt = dt_for(Nr)
    r = run_case(Nr=Nr, Nθ=Nθ, η̂=0.05, σ_ko=0.02, dt=dt, T=400.0)
    push!(g3, r.Γamp)
    @printf("%-12s %-10.5f % .4e  % .4e  %.3f\n", "($Nr,$Nθ)", dt, r.Γamp, r.ratio, r.r2)
end
@printf("Γ spread over resolution (η̂=0.05): max-min = %.4e\n", maximum(g3)-minimum(g3))
@printf("Γ₀ spread over resolution (η̂=0, σ=0.02): %s\n",
        string(round.([T1[(Nr,Nθ,0.02)] for (Nr,Nθ) in res_list], sigdigits=3)))

println("\n#########################################################")
println("# TEST 4 — η̂ SCALING + ANALYTIC MATCH (σ_ko=0.02, Nr=64)")
println("#########################################################")
ηs = [0.0, 0.02, 0.04, 0.06, 0.08, 0.10]
println(@sprintf("%-8s %-12s %-12s %-10s", "η̂", "Γamp", "Eend/E0", "R2"))
g4 = Float64[]
for η̂ in ηs
    r = run_case(Nr=64, Nθ=12, η̂=η̂, σ_ko=0.02, dt=dt_for(64), T=400.0)
    push!(g4, r.Γamp)
    @printf("%-8.2f % .4e  % .4e  %.3f\n", η̂, r.Γamp, r.ratio, r.r2)
end
# (Γ - Γ₀) linear fit through η̂>0
let x=ηs, y=g4
    Γ0=y[1]; dy=y .- Γ0
    xb=mean(x); db=mean(dy); sl=sum((x.-xb).*(dy.-db))/sum((x.-xb).^2)
    # R^2 of linearity
    dhat = db .+ sl.*(x.-xb); ss_res=sum((dy.-dhat).^2); ss_tot=sum((dy.-db).^2)
    r2lin = ss_tot>0 ? 1-ss_res/ss_tot : 1.0
    @printf("slope d(Γ-Γ₀)/dη̂ = % .4e ; linearity R2 = %.4f ; Γ₀=% .4e\n", sl, r2lin, Γ0)
    # analytic: Γ_amp_phys = ν_mom k²/2 = cν η̂ k²/2  (cν=1) -> slope_pred = k²/2
    # dominant k of n=3 seed: radial sin(3πr/R) -> k = 3π/R ; (angular ℓ=2 weaker)
    k_r = 3*π/R
    @printf("k_r (n=3 radial) = %.4f ; analytic slope k_r²/2 = %.4e\n", k_r, k_r^2/2)
    @printf("measured/analytic = %.3f\n", sl/(k_r^2/2))
end

println("\n#########################################################")
println("# TEST 5 — dt INDEPENDENCE (Nr=64, η̂=0.05, σ_ko=0.02)")
println("#########################################################")
println(@sprintf("%-10s %-12s %-12s %-10s", "dt", "Γamp", "Eend/E0", "R2"))
g5 = Float64[]
for dt in [0.01, 0.005, 0.0025]
    r = run_case(Nr=64, Nθ=12, η̂=0.05, σ_ko=0.02, dt=dt, T=400.0)
    push!(g5, r.Γamp)
    @printf("%-10.5f % .4e  % .4e  %.3f\n", dt, r.Γamp, r.ratio, r.r2)
end
@printf("Γ spread over dt = %.4e\n", maximum(g5)-minimum(g5))

println("\n#########################################################")
println("# TEST 6 — κ_Q REACTIVE CHECK (ν̂=cν=0, scan κ̂; Nr=64, σ_ko=0.02)")
println("# η̂=0 so νm=0 always; vary κ̂ (conduction). Should NOT dissipate.")
println("#########################################################")
println(@sprintf("%-10s %-12s %-12s %-10s", "κ̂", "Γamp", "Eend/E0", "R2"))
for κ in [0.0, 0.03, 0.1, 0.3]
    r = run_case(Nr=64, Nθ=12, η̂=0.0, σ_ko=0.02, dt=dt_for(64), T=400.0,
                 ν̂=0.0, cν=0.0, κ̂=κ, cκ=0.0)
    @printf("%-10.3f % .4e  % .4e  %.3f\n", κ, r.Γamp, r.ratio, r.r2)
end
println("# Now contrast: ν_mom channel at same magnitude (ν̂ scan, κ̂=0):")
println(@sprintf("%-10s %-12s %-12s %-10s", "ν̂", "Γamp", "Eend/E0", "R2"))
for ν in [0.0, 0.03, 0.1, 0.3]
    r = run_case(Nr=64, Nθ=12, η̂=0.0, σ_ko=0.02, dt=dt_for(64), T=400.0,
                 ν̂=ν, cν=0.0, κ̂=0.0, cκ=0.0)
    @printf("%-10.3f % .4e  % .4e  %.3f\n", ν, r.Γamp, r.ratio, r.r2)
end

println("\nDONE")
