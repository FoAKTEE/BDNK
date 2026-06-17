#=
    PMP telegraphers_plot reproduction (2209.09265): the BDNK heat-flow linear
    mode obeys the telegrapher equation  T̈ = c_h² T'' − τ_ε⁻¹ Ṫ  (hybrid form,
    repro/pmp_telegrapher.jl). Evolve the heat-flow Gaussian T(x,0)=0.1 e^{-x²/25}+1,
    Ṫ=0, for the three frames (σ̂,τ̂)=(0.15,1.5),(1.5,15),(7.5,75) ⇒ same thermal
    speed c_h=0.20739 but τ_ε=(2/15)τ̂=0.2,2,10. Panels at t=16,39,312 show the
    heat→wave transition: small τ_ε stays diffusive/peaked, large τ_ε splits into
    two outgoing pulses at x=±c_h t (d'Alembert).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/pmp_telegrapher.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

const CH2 = 0.04301075                       # c_h² (thermal speed² = κ(Γ-1)/(n τ_ε), frame-invariant)
const τε  = (0.2, 2.0, 10.0)                 # (2/15)·τ̂ for τ̂=1.5,15,75  (σ̂=0.15,1.5,7.5)
labels    = ("σ̂=0.15", "σ̂=1.5", "σ̂=7.5")
greys     = (0.7, 0.45, 0.0)

x = collect(range(-100.0, 100.0; length=1601)); dx = x[2]-x[1]
T0 = @. 0.1*exp(-x^2/25) + 1.0
ch = sqrt(CH2)
dt = 0.5*dx/ch                               # CFL = 0.5 (stable, leapfrog)

# Standard STABLE explicit telegrapher scheme (central in time + space):
#   (T^{n+1}-2T^n+T^{n-1})/dt² = c_h²·lap(T^n) − (1/τ)(T^{n+1}-T^{n-1})/(2dt)
#   ⇒ T^{n+1} = [2a T^n + (b-a)T^{n-1} + c_h²·lap] / (a+b),  a=1/dt², b=1/(2τ dt)
function evolve(τ, tend)
    N = length(x)
    Tprev = copy(T0)
    a = 1/dt^2; b = 1/(2τ*dt)
    # first step (T_t=0): Tnow = T0 + ½ dt² c_h² lap(T0)
    Tnow = copy(T0)
    @inbounds for i in 2:N-1
        Tnow[i] = T0[i] + 0.5*dt^2*CH2*(T0[i-1]-2T0[i]+T0[i+1])/dx^2
    end
    nsteps = ceil(Int, tend/dt)
    Tnext = similar(Tnow)
    for _ in 1:nsteps
        @inbounds for i in 2:N-1
            lap = CH2*(Tnow[i-1]-2Tnow[i]+Tnow[i+1])/dx^2
            Tnext[i] = (2a*Tnow[i] + (b-a)*Tprev[i] + lap)/(a+b)
        end
        Tnext[1]=1.0; Tnext[N]=1.0
        Tprev, Tnow, Tnext = Tnow, Tnext, Tprev
    end
    return Tnow
end

times = (16.0, 39.0, 312.0)
fig = Figure(size=(1300, 430))
for (j, t) in enumerate(times)
    ax = Axis(fig[1,j], title="t = $(Int(t))", xlabel="x", ylabel = j==1 ? "T" : "")
    for (k, τ) in enumerate(τε)
        lines!(ax, x, evolve(τ, t), linewidth=2.0,
               color=RGBf(greys[k],greys[k],greys[k]),
               label = j==1 ? labels[k] : nothing)
    end
    xlims!(ax, -95, 95); ylims!(ax, 0.999, 1.103)
    j==1 && axislegend(ax, position=:rt, framevisible=false, labelsize=12)
end
Label(fig[0,:], "BDNKStar — reproduce PMP telegraphers_plot: heat→wave (d'Alembert split) vs relaxation τ_ε",
      fontsize=14, font=:bold)
outfile = joinpath(outdir, "pmp_telegrapher_reproduction.png")
save(outfile, fig); println("saved ", outfile)
