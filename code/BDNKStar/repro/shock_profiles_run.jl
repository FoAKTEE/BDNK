#=
    Detached runner: evolve the viscous-BDNK shock for τ̂=3 (t=27, unstable) and
    τ̂=1.5 (t=100, stable) with the verified engine (pmp_viscous_core.jl), and
    save the v(x) + c₊(x) profiles to repro/shock_prof_t<τ>.txt for the figure.
    Run detached (nohup) so it survives the foreground limit.
=#
include(joinpath(@__DIR__, "pmp_viscous_core.jl"))   # engine + self-test

const ID1 = (εL=1.0, vL=0.9, nL=1.0, εR=11.5174, vR=0.354727, nR=5.44212, w=10.0)
_sherf(x) = (t=1/(1+0.3275911*abs(x)); y=1-(((((1.061405429t-1.453152027)t)+1.421413741)t-0.284496736)t+0.254829592)t*exp(-x^2); x≥0 ? y : -y)
_uv(v) = v/sqrt(1-v^2)

function init_shock(fr; N=1025, xmin=-100.0, xmax=100.0, cfl=0.1, εL,vL,nL,εR,vR,nR,w)
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

function run_save(τh, T, fname)
    fr=pmp_frame(; Γ=4/3, m=0.1, Vhat=4/3, σhat=0.0, τhat=τh)
    s=init_shock(fr; ID1...)
    nsteps=ceil(Int,T/s.dt); chunk=max(1,nsteps÷40); done=0
    while done<nsteps
        nb=min(chunk,nsteps-done); evolve!(s,nb); done+=nb
        any(!isfinite,s.ε) && break
    end
    open(fname,"w") do io
        println(io,"# x  v  cplus  (τ̂=$τh, t=$(round(done*s.dt,digits=1)))")
        for i in eachindex(s.x)
            v=s.u[i]/sqrt(1+s.u[i]^2); cp=sqrt(max(cpm2_closed(fr,s.ε[i],s.n[i])[1],0.0))
            println(io, s.x[i], "  ", v, "  ", cp)
        end
    end
    println("SAVED $fname  (t=$(round(done*s.dt,digits=1)), nans=$(any(!isfinite,s.ε)))")
end

run_save(3.0, 27.0,  joinpath(@__DIR__,"shock_prof_t3.txt"))
run_save(1.5, 80.0, joinpath(@__DIR__,"shock_prof_t1p5.txt"))
println("SHOCK_PROFILES_DONE")
