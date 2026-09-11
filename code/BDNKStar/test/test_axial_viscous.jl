using Test
using BDNKStar

# STAGE 1B/3 — AXIAL (odd-parity) VISCOUS quasi-normal modes, full GR (no Cowling):
# interior shooting + exterior vacuum (Leaver continued fraction) + Wronskian/log-derivative
# matching + complex-ω root finding, following Bussières–Redondo-Yuste–Ortega-Gómez–Cardoso
# (arXiv:2604.13208). Unlike the polar sector, the axial fluid sector is otherwise trivial,
# so viscosity DRIVES new dynamics (the η-mode family) and couples to the spacetime w-modes.
# Benchmark: Bussières Table II (EOS1 polytrope κ=100,n=1, ρc=3e15 g/cc → M≈1.27 M☉).
@testset "Axial viscous QNM (full GR): Bussières Table II w-mode + viscous shift" begin
    G = 6.6743015e-11; c = 299_792_458.0
    εc = (3e15 * 1e3) * G / c^2 * 1e6                  # ρc=3e15 g/cc → km^-2 (EOS1 reference)
    eos = PolytropeEnergy(100.0, 1.0)

    # inviscid fundamental w-mode (target: f=10.50 kHz, τ=29.54 μs)
    inv = axial_qnm(eos, εc; l=2, ηc_cgs=0.0)
    @test inv.converged
    @test isapprox(inv.f_kHz, 10.50; rtol=0.01)
    @test isapprox(inv.tau_us, 29.54; rtol=0.01)

    # viscous (BDNK frame A, η_c = 1e31 cgs): Bussières Table II → f≈10.090, τ≈30.886
    vis = axial_qnm(eos, εc; l=2, ηc_cgs=1e31, τ̂=10.0, ω0=inv.omega)
    @test vis.converged
    @test isapprox(vis.f_kHz, 10.0898; rtol=0.01)
    @test isapprox(vis.tau_us, 30.8857; rtol=0.01)
    @test vis.f_kHz < inv.f_kHz                        # shear viscosity lowers the w-mode frequency
end
