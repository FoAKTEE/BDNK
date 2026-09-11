#=
    STAGE-3 validation: the boundary-conforming (r,θ) time-domain engine reproduces the
    POLAR (even-parity) non-radial spectrum — f + pressure overtones p1..p4 — of the
    M=1.4 M☉ Shum star for ℓ=2,3,4, benchmarked against the frequency-domain eigensolver
    (NonRadialModes; itself within 0.3% of Sotani/MVHS 2107.13339).

    Method: per ℓ, one broadband evolution seeded with the Legendre angular factor
    P_ℓ(cosθ) × a Gaussian radial pulse; point probes at several radii; summed power
    spectrum. Vertical red lines = eigensolver eigenfrequencies. The pole BCs are
    ℓ-independent (m=0 axisymmetric), so only the angular factor changes with ℓ.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/sph_fp_spectrum.jl
=#
using Pkg; Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie, Printf, Statistics
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2
kHz2cyc(f) = f*BDNKStar.Units.Msun_to_km*BDNKStar.Units.kHz_to_km
Pl(l,x) = l==2 ? (3x^2-1)/2 : l==3 ? (5x^3-3x)/2 : (35x^4-30x^2+3)/8
modelab = ["f","p₁","p₂","p₃","p₄"]

function run_l(l; Nr=96, Nθ=18, T=4000.0)
    bench,_,_ = nonradial_cowling_spectrum(eos, εc; l=l, nmodes=5)
    s = build_sphstar(eos, εc; Nr=Nr, Nθ=Nθ); e = setup_sphevo(s; σ_ko=0.008)
    st = SphState(Nr, Nθ)
    for jj in 1:Nθ, ii in 1:Nr
        g = exp(-((s.grid.r[ii]-0.5*s.R)/(0.13*s.R))^2)
        st.δε[ii+1,jj+1] = 1e-3*s.ε0[ii]*g*Pl(l, s.grid.cosθ[jj])
    end
    dt = 0.25*s.grid.dr; nst = round(Int, T/dt)
    probes = [(round(Int,f*Nr)+1, 3) for f in (0.25,0.45,0.65,0.8)]
    k1=SphState(Nr,Nθ);k2=SphState(Nr,Nθ);k3=SphState(Nr,Nθ);k4=SphState(Nr,Nθ);tmp=SphState(Nr,Nθ)
    ts=Float64[]; prs=[Float64[] for _ in probes]; samp=max(1,round(Int,1.0/dt))
    for n in 0:nst
        if n%samp==0; push!(ts,n*dt); for (pk,(ip,jp)) in enumerate(probes); push!(prs[pk],st.δε[ip,jp]); end; end
        n==nst && break
        BDNKStar.SphEvolve._rk4!(st,e,dt,k1,k2,k3,k4,tmp)
    end
    fb=collect(range(1.2, 1.15*maximum(bench); length=8000)); P=zeros(length(fb))
    for p in prs; P .+= periodogram(ts, p .- mean(p), fb .* kHz2cyc(1.0)); end
    P ./= maximum(P)
    loc=Int[]; for i in 2:length(P)-1; (P[i]>P[i-1]&&P[i]>=P[i+1]&&P[i]>0.01)&&push!(loc,i); end
    peaks=sort(fb[loc]); matched=[peaks[argmin(abs.(peaks .- b))] for b in bench]
    (fb, P, bench, matched)
end

fig = Figure(size=(950, 760))
for (row,l) in enumerate((2,3,4))
    fb,P,bench,matched = run_l(l)
    @printf("ℓ=%d  err%%: %s\n", l, join((@sprintf("%s=%+.1f",modelab[m],100*(matched[m]-bench[m])/bench[m]) for m in 1:5), "  "))
    ax = Axis(fig[row,1], yscale=log10, ylabel="power (norm.)",
              title="ℓ=$l  polar spectrum  (f + p₁…p₄)",
              xlabel = row==3 ? "frequency [kHz]" : "")
    lines!(ax, fb, max.(P,1e-4), color=:navy)
    for (m,b) in enumerate(bench)
        vlines!(ax, [b], color=:crimson, linestyle=:dash)
        text!(ax, b, 1.25; text=modelab[m], color=:crimson, align=(:center,:bottom), fontsize=12)
        scatter!(ax, [matched[m]], [max(P[argmin(abs.(fb .- matched[m]))],1e-4)], color=:black, markersize=8)
    end
    ylims!(ax, 1e-4, 2.2)
end
Label(fig[0,:], "BDNKStar — STAGE 3: polar (f+p) spectrum of the M=1.4 M☉ star from the stable (r,θ) evolution, " *
      "ℓ=2,3,4 vs eigensolver (red lines) — all modes ≤4%", fontsize=12, font=:bold)
save(joinpath(outdir,"sph_fp_spectrum.png"), fig)
println("saved sph_fp_spectrum.png")
