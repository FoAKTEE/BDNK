#=
    PMP 2209.09265 conv_plot — self-convergence factor Q_N(t) of the WENO5 BDNK
    shock PDE.  Evolve the {1,0.8,0.1}_L shock (Γ=4/3, m=0.1, V̂=2/15, τ̂=1.5) at
    aligned resolutions N=513,1025,2049 (Δx halving) and compute
       Q_N(t) = log2( ||u_513 − u_1025↓|| / ||u_1025 − u_2049↓|| )
    on the common coarse grid.  Expect Q≈4 (5th-order WENO + finite-volume → ~4 on
    this smooth-shock evolution), dropping near shock-steepening events.

    Run (detached): julia code/BDNKStar/repro/pmp_conv_run.jl
=#
include(joinpath(@__DIR__, "pmp_viscous_core.jl"))
_sherf(x) = (t=1/(1+0.3275911*abs(x)); y=1-(((((1.061405429t-1.453152027)t)+1.421413741)t-0.284496736)t+0.254829592)t*exp(-x^2); x≥0 ? y : -y)
_uv(v) = v/sqrt(1-v^2)

function init_shock(fr, N; xmin=-15.0, xmax=15.0, cfl=0.1, εL,vL,nL,εR,vR,nR,w)
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

# same physical dt across resolutions (use the finest dt so step counts are aligned)
Ns = [513, 1025, 2049]
states = [init_shock(fr, N; IDsw...) for N in Ns]
dtf = states[3].dt                    # finest dt
for s in states; s.dt = dtf; end       # lock all to the finest dt (CFL-safe: coarser had larger dt)
T = 12.0; nout = 24
tout = collect(range(T/nout, T; length=nout))
l1(a,b) = sum(abs.(a .- b))/length(a)

open(joinpath(@__DIR__,"pmp_conv.txt"),"w") do io
    println(io, "# t  Q_N  d_lo  d_hi   (Q=log2(||513-1025||/||1025-2049||) of ε)")
    doneR = Ref(0)
    for (k,t) in enumerate(tout)
        target = round(Int, t/dtf)
        for s in states; evolve!(s, target-doneR[]); end
        doneR[] = target
        e513 = states[1].ε
        e1025 = states[2].ε[1:2:end]      # restrict 1025->513
        e2049 = states[3].ε[1:4:end]      # restrict 2049->513
        dlo = l1(e513, e1025); dhi = l1(e1025, e2049)
        Q = (dhi>0 && dlo>0) ? log2(dlo/dhi) : NaN
        println(io, round(t,digits=3), "  ", round(Q,digits=4), "  ", dlo, "  ", dhi)
        flush(io)
    end
end
println("PMP_CONV_DONE  (finest dt=", round(dtf,sigdigits=4), ", T=", T, ")")
