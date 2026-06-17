#=
    PMP heat_stationary reproduction (2209.09265): the closed-form initial
    acceleration ε̈(x)=(κ T')'/τ_ε from the time-symmetric heat-flow data
    T=A e^{-x²/w²}+δ, P=const (eq:heat_ID_EOM). σ̂=0 ⇒ κ=0 ⇒ ε̈≡0 (stationary,
    top panel); σ̂=1/3 ⇒ κ≠0 ⇒ localized ε̈ (bottom panel). No evolution needed —
    this is the quantitative content of fig:heat_stationary. Reuses the verified
    engine repro/pmp_viscous_core.jl for the EOS / transport coefficients.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/pmp_heat.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "repro", "pmp_viscous_core.jl"))   # engine (self-test on load)
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

# closed-form helpers (copied from repro/pmp_heat_stationary.jl, eq:heat_ID_EOM)
kappa_heat(g, ε, n, σ) = (P=pressure_eos(g,ε,n); ρ=ε+P; T=P/n; σ*ρ^2/(n^2*T))   # κ=σρ²/(n²T)
function heat_prims(g, x; A=0.4, w=10.0, δ=1.0, P0=1.0)
    T = @. A*exp(-x^2/w^2) + δ
    ε = @. P0*(g.m/T + 1/(g.Γ-1)); n = @. P0/T
    return ε, n, T
end
function epsddot(fr, x, ε, n, T)
    N=length(x); dx=x[2]-x[1]; τε = fr.L*fr.Vhat*fr.τhat
    κ = similar(x)
    for i in 1:N
        σ = transport_coeffs(fr, ε[i], n[i]).σ
        κ[i] = kappa_heat(fr.g, ε[i], n[i], σ)
    end
    Tp = similar(x); for i in 2:N-1; Tp[i]=(T[i+1]-T[i-1])/(2dx); end; Tp[1]=Tp[2]; Tp[N]=Tp[N-1]
    fl = κ.*Tp; dd = zeros(N)
    for i in 2:N-1; dd[i]=(fl[i+1]-fl[i-1])/(2dx)/τε; end
    return dd
end

x = collect(range(-50.0, 50.0; length=2049))
fig = Figure(size=(720, 540))
for (j, σh) in enumerate((0.0, 1/3))
    fr = pmp_frame(; Γ=4/3, m=0.1, Vhat=2/15, σhat=σh, τhat=1.5)
    ε, n, T = heat_prims(fr.g, x)
    dd = epsddot(fr, x, ε, n, T)
    ax = Axis(fig[j,1], xlabel = j==2 ? "x" : "", ylabel="|ε̈|",
              title="σ̂ = $(σh==0 ? "0  (κ=0 ⇒ ε̈≡0, stationary)" : "1/3  (localized ε̈)")")
    lines!(ax, x, abs.(dd), color=:black, linewidth=1.8)
    @info "σ̂=$σh  max|ε̈|" maximum(abs.(dd))
end
Label(fig[0,:], "BDNKStar — reproduce PMP heat_stationary: closed-form ε̈=(κT')'/τ_ε  (σ̂=0 stationary; σ̂=1/3 localized)",
      fontsize=12.5, font=:bold)
save(joinpath(outdir, "pmp_heat_reproduction.png"), fig)
println("saved pmp_heat_reproduction.png")
