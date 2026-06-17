#=
    PMP convergence — diagnostic: is the engine's low Q from the SCHEME order or
    from the under-resolved steepening shock?  Test on SMOOTH data (a wide, small
    ε bump, no shock) so no near-discontinuity forms.  Self-convergence factor Q
    on smooth data reveals the true spatial order.

    Run: julia code/BDNKStar/repro/pmp_conv_smooth.jl
=#
include(joinpath(@__DIR__, "pmp_viscous_core.jl"))

function init_smooth(fr, N; xmin=-50.0, xmax=50.0, cfl=0.1, A=0.2, w=400.0, n0=0.1)
    x=collect(range(xmin,xmax;length=N)); dx=x[2]-x[1]; z=zeros(N)
    ε=similar(x); n=fill(n0,N); u=zeros(N)
    for i in 1:N
        ε[i]=1.0 + A*exp(-x[i]^2/w)          # wide smooth bump, v=0
    end
    s=VState(fr,x,dx,cfl*dx,ε,n,u,copy(z),copy(z),copy(z),copy(z),copy(z),copy(z),copy(z),copy(z),false,false)
    for i in 1:N
        Ttt,Ttx,_,Jt,_=ideal_stress(fr.g,ε[i],n[i],u[i]); s.Ttt[i]=Ttt; s.Ttx[i]=Ttx; s.Jt[i]=Jt
    end
    _update_aux!(s); return s
end

fr = pmp_frame(; Γ=4/3, m=0.1, Vhat=2/15, σhat=0.0, τhat=1.5)
Ns = [513, 1025, 2049]
states = [init_smooth(fr, N) for N in Ns]
dtf = states[3].dt
for s in states; s.dt = dtf; end
l1(a,b)=sum(abs.(a .- b))/length(a)

doneR=Ref(0)
open(joinpath(@__DIR__,"pmp_conv_smooth.txt"),"w") do io
    println(io, "# t  Q   (smooth-bump self-convergence, locked dt; isolates spatial order)")
    for t in [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
        target=round(Int,t/dtf)
        for s in states; evolve!(s, target-doneR[]); end
        doneR[]=target
        e1=states[1].ε; e2=states[2].ε[1:2:end]; e3=states[3].ε[1:4:end]
        dlo=l1(e1,e2); dhi=l1(e2,e3)
        Q=(dhi>0 && dlo>0) ? log2(dlo/dhi) : NaN
        println(io, round(t,digits=2), "  ", round(Q,digits=4))
        println("  t=", t, "  Q=", round(Q,digits=3))
    end
end
println("SMOOTH_CONV_DONE")
