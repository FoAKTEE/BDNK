using Test
using BDNKStar

const G2KM = BDNKStar.Units.gram_per_cm3_to_km_minus2

@testset "Radial Cowling modes: convergence + stability (Caballero-Yunes sector)" begin
    eos = PolytropeEnergy(100.0, 1.0)
    εc = 5.5e15 * G2KM
    f4, ω4, R4 = radial_cowling_spectrum(eos, εc; N=400,  h_tov=5e-5, nmodes=4)
    f8, ω8, R8 = radial_cowling_spectrum(eos, εc; N=800,  h_tov=5e-5, nmodes=4)
    f16, ω16, R16 = radial_cowling_spectrum(eos, εc; N=1600, h_tov=5e-5, nmodes=4)
    @info "radial Cowling spectrum [kHz]" f400=f4 f800=f8 f1600=f16
    # eigenfrequencies converge with grid refinement
    @test all(abs.(f16 .- f8) ./ f16 .< 1e-4)
    @test all(abs.(f8  .- f4) ./ f8  .< 1e-3)
    # all ω² real & positive ⇒ star stable to radial perturbations (CY: stable
    # to bulk/shear); fundamental < overtones; physical kHz frequencies
    @test all(ω16 .> 0)
    @test issorted(f16)
    @test 1.0 < f16[1] < 10.0          # fundamental f-mode in the kHz band
end
