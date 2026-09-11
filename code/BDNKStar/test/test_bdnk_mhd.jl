using Test
using BDNKStar

# ===========================================================================
# BDNK viscoresistive relativistic MHD — FOUNDATION (Stage 1a)
# Lier, Armas, Porth 2026, arXiv:2606.22691.  Validates the constitutive
# tensors, the causality / front-velocity analysis in every sector, and the
# boosted-telegrapher analytic benchmark, against the paper's analytic formulas
# and its REPORTED numbers (KH/OT/Harris front velocities + the τ_X-essential
# Orszag–Tang result + the telegrapher 2nd-order convergence). Fast, no GPU.
# ===========================================================================

@testset "BDNK-MHD: ideal one-form constitutive tensors (kinematic identities)" begin
    s = MHDState([0.1, 0.2, 0.0], [0.3, 0.0, 0.4])      # E = B×v (ideal MHD)
    bμ, b2 = magnetic_fourvector(s)
    Γ = mhd_lorentz(s.v)
    uμ = Γ .* [1.0, s.v...]
    # u·u = -1, b·u = 0 (one-form-MHD constraints)
    udotu = -uμ[1]^2 + sum(uμ[2:4].^2)
    bdotu = -bμ[1]*uμ[1] + sum(bμ[2:4].*uμ[2:4])
    @test isapprox(udotu, -1.0; atol=1e-12)
    @test isapprox(bdotu, 0.0; atol=1e-12)
    @test b2 ≥ 0                                          # spacelike b
    # T^{μν} symmetric, J^{μν} antisymmetric
    T = ideal_Tmunu(1.0, s); J = ideal_Jmunu(s)
    @test maximum(abs.(T - transpose(T))) < 1e-12
    @test maximum(abs.(J + transpose(J))) < 1e-12
    # trace T^μ_μ = -ε+3p+ ... for ε=3p the conformal trace of the FLUID part is 0;
    # the magnetic part contributes -b² (traceless EM has T^μ_μ=0 in 4D ⇒ total trace 0)
    tr = -T[1,1] + T[2,2] + T[3,3] + T[4,4]
    @test isapprox(tr, 0.0; atol=1e-12)                  # conformal + EM both traceless
end

@testset "BDNK-MHD: coefficient maps + 2nd-law positivity (Eqs.10,11,17,18)" begin
    c = bdnk_coeffs_from_Dmaps(Du=1e-2, Dε=2e-3, rb=1e-2, τu=2e-1, τX=2e-1,
                               τb=8e-2, ε=30.0, b2=0.5)
    w = (4/3)*30.0
    @test isapprox(BDNKStar.BDNKMHD.sigma(c),     w*c.Dε;        rtol=1e-12)  # σ = w D_ε
    @test isapprox(BDNKStar.BDNKMHD.shear_eta(c), w*c.Du;        rtol=1e-12)  # η = w D_u
    @test isapprox(BDNKStar.BDNKMHD.bulk_zeta(c), (2/3)*w*c.Du;  rtol=1e-12)  # ζ = 2wD_u/3
    # ζ − ⅔η = 0
    @test isapprox(BDNKStar.BDNKMHD.bulk_zeta(c) - (2/3)*BDNKStar.BDNKMHD.shear_eta(c),
                   0.0; atol=1e-12)
    @test isapprox(BDNKStar.BDNKMHD.tau_eps(c), 2*c.τu; rtol=1e-12)           # τ_ε := 2τ_u
    # resistivity Eq.17: (w+b²)/w · r⊥ = r∥ = r_b
    @test isapprox((c.w+c.b2)/c.w * BDNKStar.BDNKMHD.r_perp(c),
                   BDNKStar.BDNKMHD.r_par(c); rtol=1e-12)
    @test isapprox(BDNKStar.BDNKMHD.r_par(c), c.rb; rtol=1e-12)
    # 2nd-law positivity: holds for non-negative D-maps, fails if a coefficient is negative
    @test second_law_ok(c)
    cbad = bdnk_coeffs_from_Dmaps(Du=-1e-2, Dε=2e-3, rb=1e-2, τu=2e-1, τX=2e-1,
                                  τb=8e-2, ε=30.0, b2=0.5)
    @test !second_law_ok(cbad)                                                # η = wD_u < 0
end

@testset "BDNK-MHD: SOUND subluminal window (Eq.15/16)" begin
    Du, Dε, τu = 1e-2, 2e-3, 2e-1                          # τ_ε = 2τ_u
    lo, hi = sound_subluminal_window(Du, Dε, τu)
    @test lo ≈ 2*Du                                        # lower edge = 2D_u
    # at τ_X = hi the larger sound branch is exactly luminal W² = 1 (Eq.16 boundary)
    Wp = sound_front_W2(Du, Dε, τu, hi; branch=:plus)
    Wm = sound_front_W2(Du, Dε, τu, hi; branch=:minus)
    @test isapprox(real(max(real(Wp), real(Wm))), 1.0; atol=1e-9)
    # inside the window: real and subluminal (W² ≤ 1)
    for τX in range(lo, hi; length=11)
        W2 = sound_front_W2(Du, Dε, τu, τX; branch=:plus)
        @test abs(imag(W2)) < 1e-12                        # real
        @test real(W2) ≤ 1.0 + 1e-9                        # subluminal
        @test sound_is_subluminal(Du, Dε, τu, τX)
    end
    # below the window (τ_X < 2D_u): reality lost (smaller branch W² < 0)
    @test !sound_is_subluminal(Du, Dε, τu, 0.0)
    @test real(sound_front_W2(Du, Dε, τu, 0.0; branch=:minus)) < 0
    # above the window: superluminal
    @test !sound_is_subluminal(Du, Dε, τu, hi + 0.05)
    @test real(sound_front_W2(Du, Dε, τu, hi + 0.05; branch=:plus)) > 1.0
end

@testset "BDNK-MHD: magnetic TELEGRAPHER causality + analytic + 2nd-order solver (Eq.9/21)" begin
    # causality predicate τ_b > r'_b (Eq.9)
    @test telegrapher_causal(0.5, 0.4)
    @test !telegrapher_causal(0.3, 0.4)
    # analytic mode satisfies the PDE residual ≈ 0 (finite-difference check)
    τb, rpb, k = 0.5, 0.4, 3.0
    Θ2 = (rpb/τb)*k^2 - 1/(4τb^2)
    @test Θ2 > 0                                            # propagating regime
    f(t,x) = telegrapher_analytic(t, x, k, τb, rpb)
    t0, x0, ht, hx = 0.7, 1.1, 1e-4, 1e-4
    btt = (f(t0+ht,x0) - 2f(t0,x0) + f(t0-ht,x0))/ht^2
    bt  = (f(t0+ht,x0) - f(t0-ht,x0))/(2ht)
    bxx = (f(t0,x0+hx) - 2f(t0,x0) + f(t0,x0-hx))/hx^2
    @test isapprox(τb*btt + bt - rpb*bxx, 0.0; atol=1e-5)   # PDE satisfied
    # 1D solver CONVERGES to the analytic solution at 2nd order (Fig.2)
    Ns = [64, 128, 256, 512]
    errs, orders = telegrapher_convergence(Ns, 1.0, k, τb, rpb)
    @test all(>(1.85), orders)                              # observed order ≈ 2
    @test errs[end] < errs[1]/50                            # error falls ≳ N⁻²
end

@testset "BDNK-MHD: ALFVÉN + MAGNETOSONIC branches (App.B) — grading & finiteness" begin
    c = bdnk_coeffs_from_Dmaps(Du=1e-2, Dε=2e-3, rb=1e-2, τu=2e-1, τX=2e-1,
                               τb=8e-2, ε=30.0, b2=1.0)
    # the magnetosonic quartic has NO surviving odd powers of √x (the √x-grading
    # of M_ms is exact ⇒ a genuine quartic in x) — at several angles
    for θ in range(0, π/2; length=9)
        coeffs, oddmax = magnetosonic_quartic(c, θ)
        @test oddmax < 1e-9                                 # grading witness
        @test all(isfinite, coeffs)
        @test length(magnetosonic_x(c, θ)) == 4             # four roots
    end
    # Alfvén roots real & finite on this (stable, τX large) background
    for θ in range(0, π/2; length=9)
        xp, xm = alfven_x(c, θ)
        @test isfinite(xp) && isfinite(xm)
    end
end

@testset "BDNK-MHD: KH initial v_max ≈ 0.910 / 0.950 (Table I, Eq.23)" begin
    b = kh_initial_b()                                      # ≈ 0.795
    @test 0.7 < b < 0.9
    khpar = [(1e-4, 5e-5, 1e-4, 1e-3, 5e-4, 5e-4),         # KH-a  → 0.9098
             (1e-4, 5e-5, 1e-3, 1e-3, 5e-4, 5e-3),         # KH-b  → 0.9098
             (1e-3, 5e-4, 1e-4, 8e-3, 5e-3, 5e-4),         # KH-c  → 0.9503
             (1e-3, 5e-4, 1e-3, 8e-3, 5e-3, 5e-3)]         # KH-d  → 0.9503
    targets = (0.9098, 0.9098, 0.9503, 0.9503)
    for (i, (Du,Dε,rb,τu,τX,τb)) in enumerate(khpar)
        c = bdnk_coeffs_from_Dmaps(Du=Du, Dε=Dε, rb=rb, τu=τu, τX=τX, τb=τb,
                                   ε=1.0, b2=b^2)
        v = front_velocity_max(c; nθ=721).vmax
        @test isapprox(v, targets[i]; atol=0.01)            # within 1% of the paper
        @test v < 1.0                                       # subluminal (causal)
    end
end

@testset "BDNK-MHD: Orszag–Tang — τ_X is ESSENTIAL (Eq.27): τ_X=0 COMPLEX, τ_X=0.2 REAL" begin
    b = ot_max_b()                                          # ≈ 1.404
    # τ_X = 0  ⇒  anti-diffusive / UNSTABLE: a COMPLEX front velocity.
    c0 = bdnk_coeffs_from_Dmaps(Du=1e-2, Dε=2e-3, rb=1e-2, τu=2e-1, τX=0.0,
                                τb=8e-2, ε=30.0, b2=b^2)
    r0 = front_velocity_max(c0; nθ=721)
    @test r0.max_imW > 0.05                                 # NON-zero imaginary part
    @test isapprox(r0.max_imW, 0.18834; atol=5e-4)          # paper's Im = 0.18834
    # τ_X = 0.2 ⇒ REAL, bounded, stable.
    c2 = bdnk_coeffs_from_Dmaps(Du=1e-2, Dε=2e-3, rb=1e-2, τu=2e-1, τX=0.2,
                                τb=8e-2, ε=30.0, b2=b^2)
    r2 = front_velocity_max(c2; nθ=721)
    @test r2.max_imW < 1e-6                                 # purely real
    @test isapprox(r2.vmax, 0.86717; atol=0.012)            # paper's 0.86717 (within ~1.3%)
    @test r2.vmax < 1.0                                     # subluminal
    # the decisive contrast: τ_X=0 unstable (Im>0), τ_X=0.2 stable (Im≈0)
    @test r0.max_imW > r2.max_imW + 0.1
end

@testset "BDNK-MHD: Harris-sheet representative state is causal (v_max ≲ 1)" begin
    # paper reports the global spatio-temporal max v_max = 0.9990187 (fully causal).
    # We validate v_max ≲ 1 on a representative Harris state (b ~ O(1)).
    for b in (0.3, 0.5, 0.8)
        c = bdnk_coeffs_from_Dmaps(Du=1e-2, Dε=5e-3, rb=1e-2, τu=1e-1, τX=5e-2,
                                   τb=1e-1, ε=1.0, b2=b^2)
        r = front_velocity_max(c; nθ=361)
        @test r.vmax ≤ 1.0 + 1e-6                            # subluminal / causal
        @test r.max_imW < 1e-6                               # stable (real fronts)
    end
end
