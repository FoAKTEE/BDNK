using Test
using BDNKStar
using LinearAlgebra: eigvals
using BDNKStar.NonRadialModes: nonradial_cowling_spectrum

# TRUE BDNK SHEAR STRESS (η_sh) ON THE BOUNDARY-CONFORMING GRID, read in the
# FREQUENCY DOMAIN via the ℓ-reduced operator.
#
# WHY THIS PATH. SphBDNK's legacy dissipator is ν_mom∇²δS — a scalar friction on
# the momentum density, NOT ∇_j(−2ησ^{ij}). It gets p1 right (~90% of the
# dissipation integral) but over-estimates the f-mode 3.5×, because the f-mode is
# surface-dominated with no single k. `η_sh` switches on the genuine shear stress
# (η = η_sh·ε₀ — the SAME normalisation as CowlingBDNK3D and as the dissipation
# integral in test_cross_method, so the damping rates are directly comparable).
#
# WHY ℓ-REDUCED. On the full (r,θ) grid the spectrum mixes every ℓ, every radial
# overtone and the frame modes — 574 oscillatory modes at Nr=48, all Q=1–6 — and
# the f-mode cannot be identified: continuation in η_sh returns anything from −69%
# to +632% of the expected damping depending on which eigenvalue is followed.
# Projected onto ℓ=2 the tower is clean, and the frequencies land on the validated
# `nonradial_cowling_spectrum` values.
#
# CONVERGENCE (not asserted here — too slow for the suite — but measured). dγ/dη_sh
# as a percentage of the dissipation-integral values (f 0.08746, p1 0.58510):
#     Nr (Nθ=32)   300   450   600   800  1000  1200      Nθ (Nr=800)  24    32    48
#     f            82.6  94.6 101.1 106.1 109.1 111.1      f          104.2 106.1 107.5
#     p1            ---   ---  47.9  49.0  49.6  50.0      p1          48.5  49.0  49.3
# It does NOT converge to the reference: f keeps climbing (increments +5.0,+3.0,+2.0,
# ratio ≈0.65 ⇒ ~115%) and p1 settles near ~51%. Nr=600 crossing 100% is a CROSSING,
# not a limit — do not quote it as agreement.
# WHY IT DOES NOT AGREE — RESOLVED, and NOT for the reason first suspected. The
# reference is Newtonian-form (σ_rr = ∂_r u_r, … no metric factors, no Φ′, dV with
# no √γ), so the natural guess was that its missing O(2M/R)≈0.3 terms explain the
# f-mode's +15%. They do not. Computing BOTH integrals on the CODE'S OWN
# eigenfunction (which also sidesteps the (W,V) convention) gives, at r<0.95R and
# resolution-independent (0.1181 at Nr=400 vs 0.1177 at Nr=600):
#     γ_rel/γ_newt = 0.95   ⇒ the relativistic correction is only −5%
#     γ_newt/γ_ref = 1.42   ⇒ the code's eigenfunction carries ~40% MORE SHEAR
# So the discrepancy is the EIGENFUNCTION, not the reference. Traced to a spike in
# the outermost cell: |δS_r| falls monotonically outward then jumps back up, because
# its equation differences δp against the frozen outer cell (which is identically
# zero for any oscillatory eigenmode, λx=0 ⇒ x=0) — a Dirichlet-zero wall at the
# surface. `free_surface=true` imposes the Lagrangian Δp=0 there instead and cuts
# the contamination 3.1× while leaving the bulk untouched.
# PRACTICAL RULE: any dissipation integral on this engine must EXCLUDE the outer
# ~2%, or it measures the spike (full-star γ_rel = 0.65 vs bulk 0.118). Eigenvalues
# are unaffected — only integrals of (∇v)² are.
@testset "SphBDNK true shear stress: ℓ-reduced frequency-domain QNMs" begin
    eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2
    M2c = BDNKStar.Units.Msun_to_km * BDNKStar.Units.kHz_to_km
    fkHz(λ) = abs(imag(λ)) / (2π * M2c)
    feig, _, _ = nonradial_cowling_spectrum(eos, εc; l=2, nmodes=3, N=1200, nscan=300)

    s = build_sphstar(eos, εc; Nr=200, Nθ=24)
    tower(ηsh) = begin
        e = setup_sphbdnk(s; η̂=0.0, τ̂=0.06, κ̂=0.03, ν̂=0.0, σ_ko=0.02, η_sh=ηsh)
        ev = eigvals(ell_reduced_operator(e, 2))
        sort([λ for λ in ev if imag(λ) > 1e-6 && 0.4 < fkHz(λ) < 14.0 &&
              abs(imag(λ))/(2*max(-real(λ), 1e-30)) > 2.0], by=fkHz)
    end

    b0 = tower(0.0)
    @test length(b0) ≥ 2
    # the ℓ=2 tower reproduces the validated eigensolver
    @test isapprox(fkHz(b0[1]), feig[1]; rtol=0.02)      # f  (measured +0.35%)
    @test isapprox(fkHz(b0[2]), feig[2]; rtol=0.03)      # p1 (measured +1.19%)
    @test -real(b0[1]) > 0                               # inviscid floor damps, not grows

    # viscous damping: track f and p1 by nearest-eigenvalue continuation in η_sh
    ηs = [0.0, 0.01, 0.02, 0.04]
    towers = [b0, (tower(η) for η in ηs[2:end])...]
    γ = map(1:2) do k
        tr = [towers[1][k]]
        for m in 2:length(ηs)
            c = towers[m]; _, idx = findmin(abs.(c .- tr[end])); push!(tr, c[idx])
        end
        [-real(z) for z in tr]
    end
    _slope(x,y) = (n=length(x); (n*sum(x.*y) - sum(x)*sum(y))/(n*sum(abs2,x) - sum(x)^2))
    sf = _slope(ηs, γ[1]); sp = _slope(ηs, γ[2])

    @test sf > 0                        # shear viscosity DAMPS the f-mode
    @test sp > 0
    @test sp > sf                       # p1 damps more (genuine k²-shear ordering)
    # against the dissipation integral (0.08746 for f). At this deliberately modest
    # Nr the floor still suppresses the slope — 65.5% measured — so the window is
    # wide; the convergence to 102% at Nr=600 is documented in the header.
    @test 0.35 < sf/0.08746 < 1.35
end
