using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using Printf, Statistics

eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2
sref = build_sphstar(eos, εc; Nr=64, Nθ=12); R = sref.R
cfl = 0.005/(R/64); dt_for(Nr) = cfl*(R/Nr)

function gamma_amp(ts, en; skip=0.40, tail=0.0)
    n=length(en); i0=max(2,round(Int,skip*n)); i1=max(i0+3,round(Int,(1-tail)*n))
    t=ts[i0:i1]; y=log.(en[i0:i1]); tb=mean(t); yb=mean(y)
    sl=sum((t.-tb).*(y.-yb))/sum((t.-tb).^2); yhat=yb.+sl.*(t.-tb)
    r2 = sum((y.-yb).^2)>0 ? 1-sum((y.-yhat).^2)/sum((y.-yb).^2) : 1.0
    (Γamp=-sl/2, r2=r2)
end
function run_case(; Nr,Nθ,η̂,σ_ko,dt,T=400.0,n_seed=3,kw...)
    s=build_sphstar(eos,εc;Nr=Nr,Nθ=Nθ); e=setup_sphbdnk(s;η̂=η̂,σ_ko=σ_ko,kw...)
    st=SphBDNKState(s.grid.Nr,s.grid.Nθ); seed_sphbdnk_n!(st,e,n_seed;A=1e-3)
    ns=round(Int,T/dt); sm=max(1,ns÷300)
    ts,_,en=evolve_sphbdnk!(st,e;dt=dt,nsteps=ns,sample=sm)
    fin=all(isfinite,en)&&all(>(0),en)
    fin || return (Γamp=NaN,ratio=NaN,r2=NaN,fin=false)
    g=gamma_amp(ts,en); (Γamp=g.Γamp,ratio=en[end]/en[1],r2=g.r2,fin=true)
end

# A) Convergence with KO turned essentially OFF (σ=0.001 min for stability) at η̂=0.05.
# If Γ is physical (ν_mom k²/2) it should be resolution-converged & equal to the σ=0.02 value.
println("# A — resolution convergence at η̂=0.05, MINIMAL KO (σ=0.002)")
println(@sprintf("%-12s %-10s %-12s %-12s %-8s","Nr,Nθ","dt","Γamp","Eend/E0","R2"))
for (Nr,Nθ) in [(48,9),(64,12),(96,18),(128,24)]
    r=run_case(Nr=Nr,Nθ=Nθ,η̂=0.05,σ_ko=0.002,dt=dt_for(Nr),T=400.0)
    @printf("%-12s %-10.5f % .4e  % .4e  %.3f\n","($Nr,$Nθ)",dt_for(Nr),r.Γamp,r.ratio,r.r2)
end

# B) Same resolution convergence but keep dt FIXED (decouple dt from Nr) at η̂=0.05,σ=0.02.
# Isolates pure spatial-resolution effect on Γ.
println("\n# B — resolution convergence at η̂=0.05, σ=0.02, FIXED dt=0.0025")
println(@sprintf("%-12s %-12s %-12s %-8s","Nr,Nθ","Γamp","Eend/E0","R2"))
for (Nr,Nθ) in [(48,9),(64,12),(96,18),(128,24)]
    r=run_case(Nr=Nr,Nθ=Nθ,η̂=0.05,σ_ko=0.02,dt=0.0025,T=400.0)
    @printf("%-12s % .4e  % .4e  %.3f\n","($Nr,$Nθ)",r.Γamp,r.ratio,r.r2)
end

# C) k-dependence: at fixed grid (96,18) & η̂=0.05,σ=0.02 vary radial overtone n.
# Physical viscous: Γ_phys-Γ0 ∝ k² ∝ n². Check the n-scaling.
println("\n# C — k-scaling: seed overtone n (Nr=96,Nθ=18,η̂=0.05,σ=0.02), and η̂=0 baseline per n")
println(@sprintf("%-6s %-12s %-12s %-12s","n","Γ(η̂=.05)","Γ0(η̂=0)","Γ-Γ0"))
for n in [2,3,4,5]
    rv=run_case(Nr=96,Nθ=18,η̂=0.05,σ_ko=0.02,dt=dt_for(96),T=400.0,n_seed=n)
    r0=run_case(Nr=96,Nθ=18,η̂=0.0,σ_ko=0.02,dt=dt_for(96),T=400.0,n_seed=n)
    @printf("%-6d % .4e  % .4e  % .4e   (k²∝n²: n²=%d)\n",n,rv.Γamp,r0.Γamp,rv.Γamp-r0.Γamp,n^2)
end

# D) Pure-viscous isolation: η̂=0.05 but FORCE κ̂=0,cκ=0,τ knobs default; subtract η̂=0 SAME-σ baseline.
# Report Γ_phys ≡ Γ(η̂)-Γ0 across resolution with matched-CFL — does the SUBTRACTED rate converge?
println("\n# D — SUBTRACTED Γ_phys = Γ(η̂=0.05)-Γ0(η̂=0) per resolution (matched CFL, σ=0.02)")
println(@sprintf("%-12s %-12s %-12s %-12s","Nr,Nθ","Γ(.05)","Γ0","Γ_phys"))
for (Nr,Nθ) in [(48,9),(64,12),(96,18),(128,24)]
    dt=dt_for(Nr)
    rv=run_case(Nr=Nr,Nθ=Nθ,η̂=0.05,σ_ko=0.02,dt=dt,T=400.0)
    r0=run_case(Nr=Nr,Nθ=Nθ,η̂=0.0,σ_ko=0.02,dt=dt,T=400.0)
    @printf("%-12s % .4e  % .4e  % .4e\n","($Nr,$Nθ)",rv.Γamp,r0.Γamp,rv.Γamp-r0.Γamp)
end

println("\nDONE2")
