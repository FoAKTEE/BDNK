using Test
using BDNKStar
const MSUN_KM = BDNKStar.Units.Msun_to_km

# STAGE 3 eigensolver: relativistic Cowling NON-radial (polar, l>=2) f/p modes.
# Equations: McDermott-Van Horn-Scholl / Sotani W,V form (arXiv:2412.18569 Eqs 11-14).
# Benchmark: relativistic Cowling l=2 spectrum of the n=1, K=100, rho_c=0.00128
# polytrope (M=1.4, R~14.15 km) — Font/Stergioulas/Kokkotas, as tabulated in
# arXiv:2107.13339 Table VI.
@testset "Non-radial Cowling f/p modes (Sotani/MVHS 2412.18569; benchmark 2107.13339)" begin
    eos = ShumPolytrope(100.0)
    εc  = 0.00128 + 100 * 0.00128^2                 # n=1 K=100 polytrope ⇒ M=1.4
    bench = [1.8825, 4.1060, 6.0298, 7.8670, 9.6663]  # f, p1, p2, p3, p4 [kHz]

    f, ω2, R = nonradial_cowling_spectrum(eos, εc; l=2, nmodes=5, N=2500, nscan=450)
    @info "l=2 Cowling spectrum [kHz]" computed=f benchmark=bench R_km=R*MSUN_KM
    @test length(f) ≥ 5
    for i in 1:5
        @test isapprox(f[i], bench[i]; rtol=3e-3)   # every mode within 0.3% of the paper
    end
    @test all(ω2 .> 0)                               # real ω² ⇒ dynamically stable (Cowling)
    @test issorted(f)                                # f < p1 < p2 < ...
    @test isapprox(R * MSUN_KM, 14.15; atol=0.1)     # the M=1.4 benchmark star

    # f-mode frequency rises with angular degree l (l=3 above l=2)
    f3, _, _ = nonradial_cowling_spectrum(eos, εc; l=3, nmodes=1, N=2500,
                                          nscan=300, ω2lo=1.5e-3, ω2hi=8e-3)
    @test f3[1] > f[1]

    # shooting convergence: a coarse grid reproduces the production f-mode to <0.2%
    fc, _, _ = nonradial_cowling_spectrum(eos, εc; l=2, nmodes=1, N=1500,
                                          nscan=300, ω2lo=1.5e-3, ω2hi=8e-3)
    @test isapprox(fc[1], f[1]; rtol=2e-3)
end
