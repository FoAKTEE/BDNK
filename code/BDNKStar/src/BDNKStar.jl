#=
    BDNKStar — causal first-order (BDNK) viscous neutron-star toolkit (Julia).

    Reproduces, then extends, published BDNK neutron-star results:
      STEP 0  shared EOS + primitive-recovery module      (this trunk)
      1A      radial linear benchmark   (Caballero–Yunes 2506.09149)
      1B      axial linear benchmark    (Redondo-Yuste 2411.16841 / Bussières)
      1C      nonlinear Cowling core    (Shum 2509.15303)
      2       1+1D dynamical GR + realistic EOS + collapse
      3       3+1D Cowling non-radial modes
      4       production

    This file wires the STEP-0 submodules. Layout is deliberately modular
    (one concern per file) so the EOS/recovery trunk stays reusable by every
    later stage, matching the project's "reusable module" code-quality rule.
=#
module BDNKStar

include("Numerics.jl")
include("Units.jl")
include("eos/EquationOfState.jl")
include("tov/TOV.jl")
include("transport/Transport.jl")
include("transport/Causality.jl")
include("recovery/Recovery.jl")
include("conformal/ConformalBDNK.jl")
include("conformal/ConformalEvolution.jl")
include("perturbations/RadialModes.jl")
include("flows/Bjorken.jl")
include("dispersion/Kovtun.jl")

using .Numerics
using .Units
using .EquationOfState
using .TOV
using .Transport
using .Causality
using .Recovery
using .ConformalBDNK
using .ConformalEvolution
using .RadialModes
using .Bjorken
using .Kovtun

# Re-export the STEP-0 public surface.
export Numerics, Units, EquationOfState, Transport, Causality, Recovery
# EOS
export AbstractEOS, BarotropicEOS, GeneralEOS,
       PolytropeEnergy, ShumPolytrope, IdealGas, TabulatedBarotrope, tabulate,
       pressure, sound_speed2, cn2, heat_conduction_stable, energy_from_pressure,
       dpdrho_eps, dpdeps_rho, specific_enthalpy, total_energy_density,
       temperature, is_thermodynamically_valid, apply_floor
# TOV background
export TOVStar, solve_tov, mass_solar
# Transport + causality
export TransportCoefficients, conformal_frame_PMP, shum_frame_speeds, shum_frame_wellposed,
       characteristic_speeds, causality_flag, is_causal
# Recovery
export prim2cons_barotropic, cons2prim_barotropic,
       prim2cons_general, cons2prim_general,
       cons2prim_bdnk_barotropic, lorentz_W
# Conformal BDNK (flat-space reference)
export ConformalFrame, pmp_luminal_frame, rankine_hugoniot,
       recover_time_derivs
# Conformal flat-space evolution (1C engine)
export ConfState, init_gaussian, init_smooth_shock, init_step, evolve!, energy_density
# Radial perturbations (STAGE 1A)
export radial_cowling_spectrum, hz_per_invkm
# Bjorken flow (PMP test)
export bjorken_pressure, bjorken_inviscid_analytic, bjorken_diagnostic, bjorken_evolve_rk4
# Kovtun dispersion relations (1907.08191)
export kovtun_cv, kovtun_shear_modes, kovtun_shear_speed

end # module BDNKStar
