#=
    Units — geometrized units (G = c = 1), lengths in km.

    Ported from the conventions of the NeutronStarOscillations.jl reference
    (Keeble & Redondo-Yuste, Zenodo 19207244, src/Constants.jl) so that
    background TOV stars and BDNK transport coefficients are comparable to the
    benchmark codes. In these units energy density and pressure carry dimension
    km^-2, the solar mass is `Msun_to_km` km, and frequencies in kHz map to
    inverse km via `kHz_to_km`.
=#
module Units

const c_SI       = 299_792_458.0            # m / s
const G_SI       = 6.6743015e-11            # m^3 kg^-1 s^-2
const M_sun_kg   = 1.988416e30              # kg

const sec_to_cm  = c_SI * 1e2
const sec_to_km  = c_SI * 1e-3
const kg_to_m    = G_SI / c_SI^2            # mass -> length (geometrized)
const cm_to_km   = 1e-5
const kg_to_g    = 1e3

# g/cm^3 -> km^-2 (mass-energy density in geometrized units)
const gram_per_cm3_to_km_minus2 = kg_to_m * 1e-3 / kg_to_g / cm_to_km^3
const dyne_per_cm2_to_km_minus2 = gram_per_cm3_to_km_minus2 / sec_to_cm^2
const Msun_to_km = M_sun_kg * kg_to_m * 1e-3
const kHz_to_km  = 1e3 / sec_to_km          # multiply a kHz frequency to get km^-1

export Msun_to_km, kHz_to_km, gram_per_cm3_to_km_minus2, dyne_per_cm2_to_km_minus2

end # module Units
