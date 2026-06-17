#=
    ConformalBDNK — conformal (traceless, ζ=0, c_s²=1/3) BDNK in flat space,
    slab symmetry. A verbatim Julia port of the reference C solver
    (ref-code/1D_conformal_bdnk/solver.c, A. Pandya; method arXiv:2201.12317),
    which is the cleanest worked example of (i) the BDNK stress tensor with its
    first-derivative dissipative corrections and (ii) the BDNK *primitive
    recovery* `compute_xiD`/`compute_uxD` (the linear solve for the
    time-derivative primitives ξ̇,u̇ given the conserved T^{tt},T^{tx}).

    Variables: ξ = ln ε (log energy density, evolved for positivity),
    u^x = u (spatial 4-velocity), u^t = √(1+u²) = W. Spatial derivs ξ_x,u_x;
    time derivs ξ̇,u̇. Conformal frame coefficients (η0, λ0, χ0); the
    PMP luminal frame is (χ0,λ0)=(25/4,25/7)η0 and the local viscous scale is
    ∝ ε^{3/4}. Reproduces the Rankine–Hugoniot steady-shock state.
=#
module ConformalBDNK

export ConformalFrame, pmp_luminal_frame, rankine_hugoniot,
       T_tt, T_tx, T_xx, compute_A, compute_Qx, compute_m2sxx,
       recover_time_derivs

struct ConformalFrame
    η0::Float64
    λ0::Float64
    χ0::Float64
end

"""PMP luminal conformal frame: (χ0,λ0)=(25/4,25/7)η0 with η0=ε0^{1/4}/(3π)."""
function pmp_luminal_frame(ε0::Real)
    η0 = ε0^0.25 / (3π)
    return ConformalFrame(η0, (25/7)*η0, (25/4)*η0)
end

# --- dissipative corrections (solver.c compute_A/compute_Qx/compute_m2sxx) ---
function compute_A(fr::ConformalFrame, ξ, u, ξx, ux, ξt, ut)
    e = exp(ξ); W = sqrt(1 + u^2)
    return (fr.χ0 * e^0.75 * (4u*ut + 4W*ux + 3*(1+u^2)*ξt + 3u*W*ξx)) / (4W)
end
function compute_Qx(fr::ConformalFrame, ξ, u, ξx, ux, ξt, ut)
    e = exp(ξ); W = sqrt(1 + u^2)
    return (e^0.75 * fr.λ0 * (4W*ut + 4u*ux + u*W*ξt + (1+u^2)*ξx)) / 4
end
function compute_m2sxx(fr::ConformalFrame, ξ, u, ξx, ux, ξt, ut)
    e = exp(ξ); W = sqrt(1 + u^2)
    return (-4 * e^0.75 * fr.η0 * W * (u*ut + W*ux)) / 3
end

# --- stress tensor components (solver.c T_tt/T_tx/T_xx) ----------------------
function T_tt(fr::ConformalFrame, ξ, u, ξx, ux, ξt, ut)
    A = compute_A(fr,ξ,u,ξx,ux,ξt,ut); Qx = compute_Qx(fr,ξ,u,ξx,ux,ξt,ut)
    m2sxx = compute_m2sxx(fr,ξ,u,ξx,ux,ξt,ut)
    W = sqrt(1+u^2); Δtt = -1 + W^2; Qt = u*Qx/W
    m2etasigmatx = u*m2sxx/W; m2etasigmatt = u*m2etasigmatx/W; e = exp(ξ)
    return e*(W^2 + Δtt/3) + A*(W^2 + Δtt/3) + 2Qt*W + m2etasigmatt
end
function T_tx(fr::ConformalFrame, ξ, u, ξx, ux, ξt, ut)
    A = compute_A(fr,ξ,u,ξx,ux,ξt,ut); Qx = compute_Qx(fr,ξ,u,ξx,ux,ξt,ut)
    m2sxx = compute_m2sxx(fr,ξ,u,ξx,ux,ξt,ut)
    W = sqrt(1+u^2); Δtx = W*u; Qt = u*Qx/W
    m2etasigmatx = u*m2sxx/W; e = exp(ξ)
    return e*(W*u + Δtx/3) + A*(W*u + Δtx/3) + Qt*u + W*Qx + m2etasigmatx
end
function T_xx(fr::ConformalFrame, ξ, u, ξx, ux, ξt, ut)
    A = compute_A(fr,ξ,u,ξx,ux,ξt,ut); Qx = compute_Qx(fr,ξ,u,ξx,ux,ξt,ut)
    m2sxx = compute_m2sxx(fr,ξ,u,ξx,ux,ξt,ut)
    W = sqrt(1+u^2); Δxx = 1 + u^2; e = exp(ξ)
    return e*(u^2 + Δxx/3) + A*(u^2 + Δxx/3) + 2Qx*u + m2sxx
end

# --- BDNK primitive recovery (solver.c compute_xiD/compute_uxD) --------------
"""
    recover_time_derivs(fr, ξ, u, ξx, ux, T00, T01) -> (ξ̇, u̇)

Linear BDNK primitive solve: given the conserved densities (T00,T01) and the
state + FROZEN spatial gradients (ξ,u,ξx,ux), invert for the time-derivative
primitives (ξ̇,u̇). Verbatim port of solver.c. The perfect-fluid part TttPF/TtxPF
uses ξ̇=u̇=0; the deficit (TttPF-T00)/η0 drives the dissipative correction.
"""
function recover_time_derivs(fr::ConformalFrame, ξ, u, ξx, ux, T00, T01)
    e = exp(ξ); U = u
    ch = fr.χ0 / fr.η0; l = fr.λ0 / fr.η0
    TttPF = (9e + 8e*U^4 + 6U^2*(3e - 2*e^0.75*fr.η0*ux) + 3*e^0.75*fr.η0*U^3*ξx) / (9 + 6U^2)
    TtxPF = (U*sqrt(1+U^2)*(8e*U^2 + 12*(e - e^0.75*fr.η0*ux) + 3*e^0.75*fr.η0*U*ξx)) / (9 + 6U^2)
    dtt = (TttPF - T00) / fr.η0
    dtx = (TtxPF - T01) / fr.η0
    DEN = e^0.75 * (9*ch*l + 12*ch*(-1 + l)*U^2 + 4*(ch*(-3 + l) - l)*U^4)
    ξt = (-2*(2ux + U*(1+U^2)*ξx)) / (sqrt(1+U^2)*(3 + 2U^2)) -
         (4*sqrt(1+U^2)*(3l + (-4 + 4ch + 6l)*U^2)*dtt) / DEN +
         (4*(3*(ch + 2l)*U + (-4 + 4ch + 6l)*U^3)*dtx) / DEN
    ut = -((sqrt(1+U^2)*(8U*ux + 3ξx)) / (12 + 8U^2)) +
         (3U*sqrt(1+U^2)*(4ch + l + 2*(2ch + l)*U^2)*dtt) / DEN -
         (3*(1+U^2)*(3ch + 2*(2ch + l)*U^2)*dtx) / DEN
    return ξt, ut
end

# --- Rankine–Hugoniot steady shock (solver.c set_initial_data SMOOTH_SHOCK) ---
"""
    rankine_hugoniot(εL, vL) -> (εR, vR)

Conformal (p=ε/3) perfect-fluid jump conditions giving the right state of a
steady planar shock from the left state (εL,vL). For (εL,vL)=(1,0.8) this is the
PMP/Pandya 2201.12317 benchmark (εR=4.4074, vR=0.41667).
"""
function rankine_hugoniot(εL::Real, vL::Real)
    εR = (εL - 9*vL^2*εL) / (3*(-1 + vL^2))
    vR = 1 / (3*vL)
    return εR, vR
end

end # module ConformalBDNK
