#=
    Pandya 2201.12317 Tab_cons — discrete conservation of ∫T^tt for the conformal-
    BDNK Gaussian clump.  The flux-form (finite-VOLUME) update of the conserved
    variable T^tt telescopes, so ∫T^tt is conserved to MACHINE PRECISION; the same
    T^tt RECONSTRUCTED from the (non-conservatively evolved) primitives ξ,u — the
    finite-DIFFERENCE quantity — drifts at the truncation level.  Both jump once
    the split pulses reach the outflow boundary.

      FV curve:  |∫ Ttt_conserved(t) − ∫ Ttt(0)|      (the evolved conserved var)
      FD curve:  |∫ T_tt(ξ,u,∂ξ,∂u)(t) − ∫(0)|        (rebuilt from primitives)

    Run: julia code/BDNKStar/repro/conf_tab_cons.jl
=#
include("/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl")
using .BDNKStar
using .BDNKStar.ConformalBDNK: ConformalFrame, T_tt
using Printf

eta0 = 10.0^0.25/(3π)
fr   = ConformalFrame(eta0, (25/7)*eta0, (25/4)*eta0)
s    = init_gaussian(fr; N=257, xmin=-200.0, xmax=200.0, A=1.0, x0=0.0, w=25.0, c=0.1, cfl=0.1)

NG = 3
function integrals(s)                       # interior FV and FD integrals of T^tt
    Ifv = 0.0; Ifd = 0.0
    for i in NG+1:length(s.x)-NG
        Ifv += s.Ttt[i]
        Ifd += T_tt(fr, s.ξ[i], s.u[i], s.ξt[i], s.ξx[i], s.ut[i], s.ux[i])
    end
    return Ifv*s.dx, Ifd*s.dx
end

Ifv0, Ifd0 = integrals(s)
ts=Float64[]; efv=Float64[]; efd=Float64[]
chunk=16; nout=128       # to t≈320 (captures the outflow-boundary jump near t~200)
for k in 0:nout
    if k>0; evolve!(s, chunk); end
    Ifv, Ifd = integrals(s)
    push!(ts, k*chunk*s.dt)
    push!(efv, abs(Ifv-Ifv0)); push!(efd, abs(Ifd-Ifd0))
end
open(joinpath(@__DIR__,"conf_tab_cons.txt"),"w") do io
    println(io, "# t  err_FV  err_FD   (|∫Ttt(t)-∫Ttt(0)|;  FV=conserved var, FD=primitive-reconstructed)")
    for i in eachindex(ts); println(io, ts[i], "  ", efv[i], "  ", efd[i]); end
end
@printf("Ifv0=%.6e  Ifd0=%.6e\n", Ifv0, Ifd0)
println("FV err: min=", minimum(filter(>(0),efv)), " plateau~", efv[nout÷4], " final=", efv[end])
println("FD err: plateau~", efd[nout÷4], " final=", efd[end])
println("SAVED conf_tab_cons.txt")
