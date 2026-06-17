using Test
using BDNKStar

@testset "Kovtun dispersion (1907.08191): stability + analytic limits" begin
    # c_v(φ=0) = relativistic velocity addition (v0±c0)/(1±v0 c0)
    cp, cm = kovtun_cv(0.5, 0.0; c0=0.5)
    @test isapprox(cp, 0.8; atol=1e-12)
    @test isapprox(cm, 0.0; atol=1e-12)
    # subluminal for v0²<1 (paper: c_v(φ)²<1)
    for v0 in (-0.9, 0.3, 0.9), φ in range(0, π; length=7)
        a, b = kovtun_cv(v0, φ; c0=0.5)
        @test abs(a) ≤ 1 + 1e-9 && abs(b) ≤ 1 + 1e-9
    end
    # gapped shear mode at k=0: Im ω = √(1-v0²)/(η v0²-θ) (units w0=η=1)
    v0 = 0.9; θ = 2.0
    w1, w2 = kovtun_shear_modes(0.0, 0.0; v0=v0, θη=θ)
    gap = sqrt(1-v0^2)/(v0^2 - θ)        # = -0.366
    @test isapprox(min(imag(w1), imag(w2)), gap; atol=1e-9)
    # STABILITY (the paper's central claim for θ/η=2): Im ω ≤ 0 for all k, φ, v0
    maxim = -Inf
    for v0 in (0.0, 0.5, 0.9), φ in range(0, π/2; length=9), k in range(0, 3; length=80)
        w1, w2 = kovtun_shear_modes(k, φ; v0=v0, θη=2.0)
        maxim = max(maxim, imag(w1), imag(w2))
    end
    @test maxim ≤ 1e-10
end
