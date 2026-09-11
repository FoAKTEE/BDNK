using Test
using BDNKStar

# Ideal-gas EOS through the full stellar pipeline. p=(Γ-1)ρϵ on an adiabat ⇒ p=Kρ^Γ,
# ε=ρ+p/(Γ-1), tabulated as a barotrope (`isentropic_idealgas`). Demonstrates that the
# TOV background, the f/p eigensolver, and BOTH viscous-QNM solvers (polar matrix +
# axial full-GR) are EOS-general — they run on the ideal gas, not just the Shum polytrope.
@testset "Ideal-gas EOS: TOV + f/p spectrum + polar & axial viscous QNM" begin
    ig = isentropic_idealgas(Γ=2.0, K=100.0)

    # valid causal barotrope
    for e in (1e-4, 1e-3, 5e-3)
        @test 0 < sound_speed2(ig, e) < 1
        @test pressure(ig, e) > 0
    end

    # TOV background (sensible M–R for a Γ=2 ideal-gas star)
    εc = 1.4e-3; st = solve_tov(ig, εc)
    @test 0.5 < mass_solar(st) < 1.5
    @test 10 < st.R * BDNKStar.Units.Msun_to_km < 18

    # polar f/p eigenspectrum — real and ordered
    fr,_,_ = nonradial_cowling_spectrum(ig, εc; l=2, nmodes=3)
    @test issorted(fr)
    @test 1.5 < fr[1] < 2.2                         # f-mode in the expected band

    # polar viscous QNM — viscosity damps the f-mode (γ rises with η̂)
    g0 = qnm_damping(polar_qnm(ig, εc; l=2, η̂=0.0,  Nr=96, nmodes=1, warn=false)[1])
    gv = qnm_damping(polar_qnm(ig, εc; l=2, η̂=0.03, Nr=96, nmodes=1, nstep=6, warn=false)[1])
    @test abs(g0) < 5e-3
    @test gv > g0 + 3e-3

    # axial full-GR w-mode — converges with a finite damping time
    ax = axial_qnm(ig, εc; l=2, ηc_cgs=0.0)
    @test ax.converged
    @test ax.f_kHz > 5 && ax.tau_us > 0
end
