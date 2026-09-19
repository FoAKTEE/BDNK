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
include("eos/PiecewisePolytrope.jl")
include("tov/TOV.jl")
include("rotation/SlowRotation.jl")
include("rotation/RModes.jl")
include("transport/Transport.jl")
include("transport/Causality.jl")
include("recovery/Recovery.jl")
include("conformal/ConformalBDNK.jl")
include("conformal/ConformalEvolution.jl")
include("perturbations/RadialModes.jl")
include("perturbations/NonRadialModes.jl")
include("perturbations/TidalDeformability.jl")
include("flows/Bjorken.jl")
include("dispersion/Kovtun.jl")
include("viscous/IsraelStewart.jl")
include("cowling3d/Background3D.jl")
include("cowling3d/CowlingEvolve3D.jl")
include("cowling3d/CowlingBDNK3D.jl")
include("spherical/SphBackground.jl")
include("spherical/SphEvolve.jl")
include("spherical/SphBDNK.jl")
include("perturbations/PolarViscousModes.jl")
include("perturbations/AxialViscousModes.jl")
include("perturbations/PolarGRModes.jl")
include("perturbations/GravityModes.jl")
include("perturbations/DoublyDiffusive.jl")
include("fvcartesian/FVCommon.jl")
include("fvcartesian/FVRadial.jl")
include("fvcartesian/FVCartesian.jl")
include("dyngr/DynGR1D.jl")
include("dg/DGCommon.jl")
include("dg/DGExactRiemann.jl")
include("dg/DGSRHydro.jl")
include("dg/DGStar.jl")
include("dg/DGSubcell.jl")
include("dg/DGStarHP.jl")
include("dg/DGStarFD.jl")
include("dg/DGBall3D.jl")
include("dg/DGCart2D.jl")
include("dg/DGCart3D.jl")
include("mhd/BDNKMHD.jl")
include("mhd/BDNKMHDConstitutive.jl")
include("mhd/BDNKMHD1D.jl")
include("mhd/BDNKMHD2D.jl")

using .Numerics
using .Units
using .EquationOfState
using .PiecewisePolytrope
using .TOV
using .SlowRotation
using .RModes
using .Transport
using .Causality
using .Recovery
using .ConformalBDNK
using .ConformalEvolution
using .RadialModes
using .NonRadialModes
using .TidalDeformability
using .Bjorken
using .Kovtun
using .IsraelStewart
using .Background3D
using .CowlingEvolve3D
using .CowlingBDNK3D
using .SphBackground
using .SphEvolve
using .SphBDNK
using .PolarViscousModes
using .AxialViscousModes
using .PolarGRModes
using .GravityModes
using .DoublyDiffusive
using .FVCommon
using .FVRadial
using .FVCartesian
using .DynGR1D
using .DGCommon
using .DGExactRiemann
using .DGSRHydro
using .DGStar
using .DGSubcell
using .DGStarHP
using .DGStarFD
using .DGBall3D
using .DGCart2D
using .DGCart3D
using .BDNKMHD
using .BDNKMHD1D
using .BDNKMHD2D

# FV Cartesian/radial reproduction (Stage 1 radial MVP + Stage 2 Cartesian)
export RadialGrid, RadialState, setup_fvradial, evolve_fvradial!,
       seed_radial_velocity!, fmode_radial_reference, radial_periodogram_freqs
export setup_fvcart, evolve_fvcart!, seed_l2_velocity!, fvcart_quadrupole,
       fvcart_central_density, fvcart_periodogram_freqs
# DYNAMICAL-GR 1+1D engine (Stage 2): metric evolved (NOT Cowling) — TOV
# stationarity, full-GR radial mode, max-mass marginal stability, BH collapse
export DynGRGrid, DynGRState, DynGREngine, setup_dyngr, evolve_dyngr!,
       seed_dyngr_velocity!, seed_dyngr_toroidal!, deplete_pressure!, dyngr_metric!,
       dyngr_radial_freq, chandrasekhar_radial_omega2, cowling_radial_omega2,
       dyngr_central_lapse, dyngr_max_2mor, dyngr_central_density
# DG (RKDG) fallback scheme — Stage 1 shock tube, Stage 2 radial star, Stage 3 2D
export build_lgl_basis, lgl_nodes_weights, LGLBasis
export setup_srdg, set_initial!, evolve_srdg!, srdg_primitives, srdg_cell_means,
       shocktube_initial!, blast_initial!
export exact_riemann_sr, sample_riemann_sr, RiemannSol
export setup_dgstar, evolve_dgstar!, seed_dgstar_radial!, dgstar_central_density,
       dgstar_surface_width, dgstar_radial_freq
export setup_dgstarfd, evolve_dgstarfd!, seed_dgstarfd_radial!, dgstarfd_central_density,
       dgstarfd_errD, dgstarfd_baryon_mass, dgstarfd_fd_fraction, dgstarfd_active_map,
       dgstarfd_spectrum
export hp_grid, setup_dgstarhp, evolve_dgstarhp!, seed_dgstarhp_radial!, dgstarhp_central_density,
       dgstarhp_errD, dgstarhp_baryon_mass, dgstarhp_spectrum, dgstarhp_surface_state
export ball_grid, setup_dgball3d, evolve_dgball3d!, seed_dgball3d_radial!, seed_dgball3d_l2!,
       dgball3d_central_density, dgball3d_errD, dgball3d_baryon_mass, dgball3d_moment,
       dgball3d_static_residual, dgball3d_volume, dgball3d_rhs_norm, dgball3d_metric_identity
export setup_dgcart2d, evolve_dgcart2d!, seed_dgcart2d_l2!, dgcart2d_quadrupole,
       dgcart2d_central_density
export setup_dgcart3d, evolve_dgcart3d!, seed_dgcart3d_l2!, seed_dgcart3d_Y22!,
       seed_dgcart3d_Y21!, dgcart3d_central_density, dgcart3d_quadrupole,
       dgcart3d_quadrupole_m2, dgcart3d_shocktube_diagonal!, dgcart3d_prim_minmax, dgcart3d_limiter_census

# Re-export the STEP-0 public surface.
export Numerics, Units, EquationOfState, Transport, Causality, Recovery
# EOS
export AbstractEOS, BarotropicEOS, GeneralEOS,
       PolytropeEnergy, ShumPolytrope, IdealGas, TabulatedBarotrope, tabulate, isentropic_idealgas,
       pressure, sound_speed2, cn2, heat_conduction_stable, energy_from_pressure,
       dpdrho_eps, dpdeps_rho, specific_enthalpy, total_energy_density,
       temperature, is_thermodynamically_valid, apply_floor
# Read et al. (2009) realistic piecewise-polytrope EOS (SLy/APR4/H4/MS1 presets)
export piecewise_polytrope
# TOV background
export TOVStar, solve_tov, mass_solar
# Finite-temperature ideal-gas star (heat-conduction background)
export IdealGasStar, solve_tov_idealgas
# Slow-rotation (Hartle 1967) frame dragging + moment of inertia
export SlowRotResult, moment_of_inertia, IBAR_FROM_C_BREU, MR2_FROM_C_LS, KM3_TO_G_CM2
# r-mode (l=m=2) CFS instability window + BDNK causal ζ_eff extension
export rmode_frequencies, rmode_timescales, rmode_instability_window,
       rmode_structure_constants, zeta_eff_factor, lane_emden_n1_star,
       rmode_validate_lom98, RModeWindow, RMODE_GW_COEFF, kepler_frequency,
       TAU_GW_PUB, TAU_SV_PUB, TAU_BV_PUB, JTILDE_PUB, ITILDE_PUB
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
# Non-radial (polar ℓ≥2) Cowling f/p modes (STAGE 3 eigensolver)
export nonradial_cowling_spectrum, freq_kHz_from_omega2
# ℓ=2 tidal Love number k₂, deformability Λ (static) + Λ_eff(ω) (dynamical tide)
export tidal_love_number, lambda_from_k2, lambda_eff, lambda_eff_modes,
       LambdaEffModes, ibar_yagi_yunes, YY_ILOVE_COEFFS
# Gravity (g) modes: relativistic Cowling polar eigensolver WITH buoyancy/Schwarzschild discriminant
export gmode_spectrum, brunt_vaisala, freq_kHz_from_omega2_g
# Local (WKB) doubly-diffusive dispersion: GSF/rotational + thermohaline/heat-coupled
# g-modes in the BDNK causal framework (known instabilities; new BDNK realization)
export gsf_roots_NS, gsf_roots_BDNK, gsf_max_growth_NS, gsf_max_growth_BDNK,
       gsf_epicyclic, gsf_unstable_cone,
       dd_cubic_roots, dd_bdnk_roots, dd_max_growth, dd_is_overstable,
       gmode_root, gmode_growth_slope, most_unstable, fastest_rate
# 3+1D Cowling evolution (STAGE 3 Phase 2)
export CartesianGrid, Star3D, build_star3d
export EvolState, Evo3D, setup_evo3d, seed_l2!, evolve3d!, l2_quadrupole,
       periodogram, freq_kHz_cyclic, damping_rate,
       ylm_real, seed_ylm!, ylm_moment, evolve3d_moments!,
       spectral_peak, pencil_modes, envelope_ratio, analyze_qnm
# Full-frame (causal/hyperbolic) BDNK 3+1D Cowling (STAGE 3, frame recovery)
export BDNKState, BDNK3D, setup_bdnk3d, seed_bdnk_l2!, evolve_bdnk3d!,
       bdnk_causal_denominator, bdnk_bound_violation,
       seed_bdnk_ylm!, evolve_bdnk3d_moments!, bdnk3d_qnm
# Boundary-conforming (r,θ) spherical BDNK (STAGE 3, surface rewrite)
export SphGrid, SphStar, build_sphstar
export SphState, SphEvo, setup_sphevo, seed_sph_l2!, evolve_sph!, l2_quad_sph
export SphBDNKState, SphBDNKEvo, setup_sphbdnk, seed_sphbdnk_l2!, seed_sphbdnk_n!,
       seed_sphbdnk_noise!, evolve_sphbdnk!, l2_quad_sphbdnk, lℓ_quad_sphbdnk, mode_energy,
       ell_reduced_operator
# Polar (even-parity) VISCOUS quasi-normal modes (STAGE 3, frequency-domain eigenproblem)
export polar_bdnk_operator, polar_qnm, qnm_freq_kHz, qnm_damping
# Axial (odd-parity) VISCOUS quasi-normal modes (STAGE 1B, full-GR shooting + Wronskian matching)
export axial_qnm, Viscosity, inviscid,
       omega_to_ftau, ftau_to_omega, frameA_viscosity, build_axial_star
# Polar (even-parity) FULL-GR quasi-normal modes (STAGE 3, Lindblom–Detweiler interior + Zerilli exterior)
export polar_gr_qnm, polar_omega_to_ftau, polar_ftau_to_omega,
       build_polar_star, zerilli_potential, interior_metric_response
# Bjorken flow (PMP test)
export bjorken_pressure, bjorken_inviscid_analytic, bjorken_diagnostic, bjorken_evolve_rk4
# Kovtun dispersion relations (1907.08191)
export kovtun_cv, kovtun_shear_modes, kovtun_shear_speed
# Israel–Stewart bulk-viscous closure + causal limiter (STAGE 2, 2311.13027)
export bulk_pressure_NS, enthalpy_prime, cs2_viscous, tauPi_causal,
       apply_causality_fix, reynolds_inv, reynolds_inv_min, reynolds_inv_max,
       DPi_conserved
# BDNK viscoresistive relativistic MHD (STAGE 1a, Lier–Armas–Porth 2026, 2606.22691):
# ideal one-form constitutive tensors + the causality/front-velocity analysis
# (sound Eq.15/16, telegrapher Eq.9, Alfvén Eq.B5, magnetosonic Eq.B9, v_max Eq.B13)
# + the boosted-telegrapher analytic benchmark + 2nd-order 1D solver
export MHDState, mhd_lorentz, magnetic_fourvector, ideal_Tmunu, ideal_Jmunu
export bdnk_coeffs_from_Dmaps, BDNKMHDCoeffs, second_law_ok
export sound_front_W2, sound_subluminal_window, sound_is_subluminal
export telegrapher_causal, telegrapher_analytic, solve_telegrapher_1d, telegrapher_convergence
export alfven_x, alfven_W, magnetosonic_quartic, magnetosonic_x, magnetosonic_W
export front_velocity_max, kh_initial_b, ot_b2_at, ot_max_b
# BDNK viscoresistive MHD 1D shock-tube EVOLVER (STAGE 1b, Eq.19/20/22, Fig.3):
# 7-component conservative FV with local-matrix primitive recovery (Eq.20)
export MHD1DGrid, MHD1DState, MHD1DEngine, setup_mhd1d_shocktube,
       evolve_mhd1d!, mhd1d_primitives, mhd1d_Jty, mhd1d_pressure,
       mhd1d_assemble_M, mhd1d_cond_number, mhd1d_conserved_totals,
       mhd1d_vmax_field, ST_PARAMS, st_coeffs,
       setup_mhd1d_sinmode, mhd1d_mode_amplitude, mhd1d_growth_rate
# BDNK viscoresistive MHD 2D EVOLVER (STAGE 1c, Eq.23/27/28-32, Fig.4/5/6/8/9):
# extends the 7-component conservative FV to (x,y) with both fluxes, SSP-RK2,
# Tóth cell-centered div-B=0 projection; reproduces Orszag–Tang (the τ_X-essential
# τ_X=0-unstable vs τ_X=0.2-stable contrast), Kelvin–Helmholtz (viscosity/
# resistivity trend), and the double Harris current sheet (causal v_max<1).
export MHD2DGrid, MHD2DState, MHD2DEngine
export setup_mhd2d_orszagtang, setup_mhd2d_kelvinhelmholtz, setup_mhd2d_harris
export evolve_mhd2d!, mhd2d_primitives, mhd2d_energy_density, mhd2d_Jti,
       mhd2d_divB, mhd2d_divB_max, mhd2d_out_of_plane_current, mhd2d_tracer,
       mhd2d_vmax_field, mhd2d_front_velocity, mhd2d_conserved_totals,
       OT_PARAMS, KH_PARAMS, HARRIS_PARAMS, ot_coeffs, kh_coeffs, harris_coeffs

end # module BDNKStar
