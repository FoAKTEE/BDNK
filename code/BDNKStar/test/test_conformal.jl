using Test
using BDNKStar
using BDNKStar.ConformalBDNK: T_tt, T_tx, T_xx

@testset "Conformal BDNK: Rankine–Hugoniot steady shock (PMP/Pandya)" begin
    εR, vR = rankine_hugoniot(1.0, 0.8)        # εL=1, vL=0.8
    @test isapprox(εR, 4.4074; atol=1e-4)      # arXiv:2201.12317 benchmark
    @test isapprox(vR, 0.41667; atol=1e-5)     # = 1/(3·0.8)
end

@testset "Conformal BDNK: perfect-fluid limit (gradients = 0)" begin
    fr = pmp_luminal_frame(10.0)
    for ξ in (-0.5, 0.0, 0.7), u in (-0.6, 0.0, 0.9)
        e = exp(ξ); W = sqrt(1+u^2); p = e/3
        # PF conformal: T^{tt}=(e+p)W²-p, T^{tx}=(e+p)W²v, v=u/W
        @test isapprox(T_tt(fr, ξ, u, 0,0,0,0), (e+p)*W^2 - p; rtol=1e-12)
        @test isapprox(T_tx(fr, ξ, u, 0,0,0,0), (e+p)*W^2*(u/W); rtol=1e-12)
    end
end

@testset "Conformal BDNK: primitive recovery round-trip (CrossCheck vs solver.c)" begin
    fr = pmp_luminal_frame(10.0)
    maxerr = 0.0
    # pick state + frozen spatial gradients + TRUE time derivatives, go forward
    # through the BDNK stress, then invert and check we recover (ξ̇, u̇).
    for ξ in (-0.3, 0.0, 0.4), u in (-0.4, 0.0, 0.7)
        for (ξx, ux, ξt_true, ut_true) in
                ((0.05, -0.03, 0.02, -0.01), (-0.1, 0.07, -0.05, 0.04))
            ut = sqrt(1+u^2)
            T00 = T_tt(fr, ξ, u, ξx, ux, ξt_true, ut_true)
            T01 = T_tx(fr, ξ, u, ξx, ux, ξt_true, ut_true)
            ξt, ut_dot = recover_time_derivs(fr, ξ, u, ξx, ux, T00, T01)
            maxerr = max(maxerr, abs(ξt - ξt_true), abs(ut_dot - ut_true))
        end
    end
    @info "conformal BDNK recovery round-trip max error" maxerr
    @test maxerr ≤ 1e-9
end
