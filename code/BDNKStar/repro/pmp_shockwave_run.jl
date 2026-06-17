#=
    PMP 2209.09265 shockwave_plot (fig:shockwave_profile) — steady planar BDNK
    shock structure ε(x), v(x), n(x).  Left state {ε_L,v_L,n_L}={1,0.8,0.1};
    Γ=4/3, m=0.1, V̂=2/15, σ̂=0, τ̂=1.5; RH right state {4.439,0.4143,0.2929}.
    Evolve from a smooth initial shock to the relaxed steady state.

    Run: julia code/BDNKStar/repro/pmp_shockwave_run.jl
=#
include(joinpath(@__DIR__, "pmp_viscous_core.jl"))

_sherf(x) = (t=1/(1+0.3275911*abs(x)); y=1-(((((1.061405429t-1.453152027)t)+1.421413741)t-0.284496736)t+0.254829592)t*exp(-x^2); x≥0 ? y : -y)
_uv(v) = v/sqrt(1-v^2)

function init_shock(fr; N=2049, xmin=-20.0, xmax=20.0, cfl=0.1, εL,vL,nL,εR,vR,nR,w)
    x=collect(range(xmin,xmax;length=N)); dx=x[2]-x[1]; z=zeros(N)
    ε=similar(x); n=similar(x); u=similar(x)
    for i in 1:N
        ξ=x[i]
        ε[i]=(εR-εL)/2*(_sherf(ξ/w)+1)+εL
        v=(vL-vR)/2*(1-_sherf(ξ/w))+vR
        n[i]=(nL-nR)/2*(1-_sherf(ξ/w))+nR
        u[i]=_uv(v)
    end
    s=VState(fr,x,dx,cfl*dx,ε,n,u,copy(z),copy(z),copy(z),copy(z),copy(z),copy(z),copy(z),copy(z),false,false)
    for i in 1:N
        Ttt,Ttx,_,Jt,_=ideal_stress(fr.g,ε[i],n[i],u[i]); s.Ttt[i]=Ttt; s.Ttx[i]=Ttx; s.Jt[i]=Jt
    end
    _update_aux!(s); return s
end

const IDsw = (εL=1.0, vL=0.8, nL=0.1, εR=4.439, vR=0.4143, nR=0.2929, w=2.0)
fr = pmp_frame(; Γ=4/3, m=0.1, Vhat=2/15, σhat=0.0, τhat=1.5)
s  = init_shock(fr; N=513, xmin=-15.0, xmax=15.0, IDsw...)

# save initial profile
xs = copy(s.x)
ε0 = copy(s.ε); v0 = [s.u[i]/sqrt(1+s.u[i]^2) for i in eachindex(s.x)]; n0 = copy(s.n)

T=12.0; nsteps=ceil(Int,T/s.dt); chunk=max(1,nsteps÷40); doneR=Ref(0)
while doneR[]<nsteps
    nb=min(chunk,nsteps-doneR[]); evolve!(s,nb); doneR[]+=nb
    any(!isfinite,s.ε) && break
end
done = doneR[]
εf = copy(s.ε); vf = [s.u[i]/sqrt(1+s.u[i]^2) for i in eachindex(s.x)]; nf = copy(s.n)

open(joinpath(@__DIR__,"pmp_shockwave.txt"),"w") do io
    println(io, "# x  eps0 v0 n0  epsF vF nF   (t=0 initial vs t=$(round(done*s.dt,digits=1)) relaxed; nan=$(any(!isfinite,s.ε)))")
    for i in eachindex(xs)
        println(io, xs[i],"  ",ε0[i]," ",v0[i]," ",n0[i],"  ",εf[i]," ",vf[i]," ",nf[i])
    end
end
println("SAVED pmp_shockwave.txt | t_f=", round(done*s.dt,digits=1),
        " εR(final center+)=", round(εf[findmin(abs.(xs.-3))[2]],sigdigits=5),
        " vL=", round(vf[1],sigdigits=4), " nan=", any(!isfinite,s.ε))
