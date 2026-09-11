#=
    BDNKMHD — FOUNDATION of BDNK viscoresistive relativistic MHD
    (Lier, Armas, Porth 2026, arXiv:2606.22691).  Stage 1a of the reproduction:

      (1) the IDEAL one-form MHD constitutive tensors T^{μν}_(0), J^{μν}_(0)
          + the BDNK coefficient maps (Eqs. 10,11,17,18) for ε=3p (w=4/3 ε);
      (2) the CAUSALITY / FRONT-VELOCITY analysis — the paper's key theoretical
          result — in every sector:
            • SOUND  front velocity W² (Eq.15) + the subluminal window (Eq.16);
            • the magnetic TELEGRAPHER causality predicate τ_b > r'_b (Eq.9);
            • ALFVÉN  front velocities x_A^± (App. B, Eq. B5, τ_ε:=2τ_u imposed);
            • MAGNETOSONIC front velocities — roots of det M_ms(x)=0 (Eq. B9),
              the 4×4 matrix transcribed verbatim and expanded to a quartic in x;
            • v_max = max over angle θ and over ALL branches of |Re W| (Eq. B13);
      (3) the BOOSTED-TELEGRAPHER analytic benchmark (Eq.21) + a 2nd-order
          1D FV/FD solver of  τ_b ∂_tt b + ∂_t b − r'_b ∂_xx b = 0  that
          converges to it at 2nd order (Fig. 2).

    Flat space, mostly-plus metric η=diag(−1,1,1,1).  STDLIB-only (LinearAlgebra).

    ── MAGNETOSONIC TRANSCRIPTION NOTE (the main reproduction risk) ────────────
    The 4×4 M_ms(x) (Eq. B9) has entries that are polynomials in s=√x.  For
    det M_ms to be a genuine polynomial in x (only EVEN powers of s), the matrix
    must carry a √x-grading: row/col parities r=(0,1,1,0), c=(0,1,1,0), so every
    entry M_{ij} = √x^{r_i+c_j}·P_{ij}(x).  All sixteen printed entries respect
    this EXCEPT the bare `wD_ε` term in entry (1,3), whose grade demands a √x.
    Reading (1,3) = (wD_ε − 3τ_u(w+b⊥²))√x makes ALL odd powers of s cancel
    exactly (verified numerically), leaving a clean quartic a4 x⁴+…+a0.  We use
    that graded reading; with it the OT τ_X=0 anti-diffusive imaginary front
    velocity reproduces the paper's 0.18834 to 5 sig figs (see test).
=#
module BDNKMHD

using LinearAlgebra

export MHDState, mhd_lorentz, magnetic_fourvector, ideal_Tmunu, ideal_Jmunu
export bdnk_coeffs_from_Dmaps, BDNKMHDCoeffs, second_law_ok
export MHDThermo, ConformalThermo, BarotropeThermo, mhd_pressure, mhd_Tslope
export sound_front_W2, sound_subluminal_window, sound_is_subluminal
export telegrapher_causal, telegrapher_analytic, solve_telegrapher_1d, telegrapher_convergence
export alfven_x, alfven_W, magnetosonic_quartic, magnetosonic_x, magnetosonic_W
export front_velocity_max, kh_initial_b, ot_b2_at, ot_max_b

# ===========================================================================
# 0.  KINEMATICS + IDEAL ONE-FORM MHD CONSTITUTIVE TENSORS
# ===========================================================================

const ETA = Float64[-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]   # mostly-plus η^{μν}

"Lorentz factor Γ = 1/√(1−v²) for a 3-velocity vector v=(vx,vy,vz)."
mhd_lorentz(v::AbstractVector) = 1.0/sqrt(1.0 - dot(v,v))

"""
    MHDState(v, B, E=B×v)

A flat-space ideal-MHD fluid cell: 3-velocity `v`, magnetic field `B`,
electric field `E` (defaults to the ideal-MHD value **E** = **B**×**v**).
"""
struct MHDState
    v::Vector{Float64}   # 3-velocity
    B::Vector{Float64}   # magnetic 3-field
    E::Vector{Float64}   # electric 3-field (ideal MHD: B×v)
end
MHDState(v, B) = MHDState(collect(float.(v)), collect(float.(B)),
                          cross(collect(float.(B)), collect(float.(v))))

"u^μ = Γ(1, **v**)."
function four_velocity(s::MHDState)
    Γ = mhd_lorentz(s.v)
    return Γ .* [1.0, s.v[1], s.v[2], s.v[3]]
end

"""
    magnetic_fourvector(s) -> (bμ, b2)

The one-form-MHD magnetic 4-vector  b^μ = Γ(**B·v**, **B** − **v**×**E**)  and
its norm-squared b² = b^μ b_μ.  Satisfies b^μ u_μ = 0, by construction.
"""
function magnetic_fourvector(s::MHDState)
    Γ = mhd_lorentz(s.v)
    bt = Γ * dot(s.B, s.v)
    bvec = Γ .* (s.B .- cross(s.v, s.E))
    bμ = [bt, bvec[1], bvec[2], bvec[3]]
    b2 = -bμ[1]^2 + bμ[2]^2 + bμ[3]^2 + bμ[4]^2
    return bμ, b2
end

"Projector P^{μν} = η^{μν} + u^μ u^ν (orthogonal to u)."
function projector(uμ::AbstractVector)
    P = copy(ETA)
    for μ in 1:4, ν in 1:4
        P[μ,ν] += uμ[μ]*uμ[ν]
    end
    return P
end

"""
    ideal_Tmunu(ε, s) -> 4×4

Ideal one-form-MHD stress tensor (ε=3p, p=ε/3, w=ε+p=4/3 ε):
  T^{μν}_(0) = (ε+p+b²) u^μu^ν + (p+½b²) g^{μν} − b^μb^ν .
"""
function ideal_Tmunu(ε::Real, s::MHDState)
    p = ε/3
    uμ = four_velocity(s)
    bμ, b2 = magnetic_fourvector(s)
    w = ε + p
    T = zeros(4,4)
    for μ in 1:4, ν in 1:4
        T[μ,ν] = (w + b2)*uμ[μ]*uμ[ν] + (p + 0.5*b2)*ETA[μ,ν] - bμ[μ]*bμ[ν]
    end
    return T
end

"""
    ideal_Jmunu(s) -> 4×4

Ideal one-form magnetisation current  J^{μν}_(0) = u^μ b^ν − u^ν b^μ
(antisymmetric).
"""
function ideal_Jmunu(s::MHDState)
    uμ = four_velocity(s)
    bμ, _ = magnetic_fourvector(s)
    J = zeros(4,4)
    for μ in 1:4, ν in 1:4
        J[μ,ν] = uμ[μ]*bμ[ν] - uμ[ν]*bμ[μ]
    end
    return J
end

# ===========================================================================
# 1.  BDNK COEFFICIENT MAPS  (Eqs. 10, 11, 17, 18)
# ===========================================================================

"""
    BDNKMHDCoeffs

The transport + BDNK-relaxation data of the viscoresistive theory, in the
ε=3p conformal MHD parametrisation actually used by the paper's tests:
the dimensionless control set (D_u, D_ε, r_b, τ_u, τ_X, τ_b) plus the derived
physical coefficients (σ, η, ζ, r∥, r⊥) and τ_ε := 2τ_u.
"""
# ---- EOS thermodynamics for the verbatim constitutive (general ε↔p) --------
# The per-cell verbatim BDNK-MHD constitutive needs, at the LOCAL ε: the pressure
# p(ε) (⇒ enthalpy w=ε+p), and — for the heat-flux temperature gradient
# ∂_ρ lnT — the log-slope s_T ≡ d lnT/d lnε.  Conformal ε=3p ⇒ p=ε/3, T∝ε^{1/4}
# ⇒ s_T=1/4.  A general barotrope supplies p(ε) and s_T(ε).  This is what lets the
# verbatim constitutive couple to a REALISTIC (ε≠3p) neutron-star EOS in the GR
# engines.  The general cs²=dp/dε enters the divergence terms AUTOMATICALLY via the
# finite-difference ideal-current Jacobian (which perturbs ε through p(ε)).
abstract type MHDThermo end

"Conformal ultra-relativistic closure ε=3p (the paper's default)."
struct ConformalThermo <: MHDThermo end
@inline mhd_pressure(::ConformalThermo, ε) = ε/3
@inline mhd_Tslope(::ConformalThermo, ε)   = 0.25

"""
General barotrope closure for the constitutive: pressure law `pfun(ε)` and
temperature log-slope `sTfun(ε)=d lnT/d lnε` (use `ε->0.0` to disable the thermal
gradient, e.g. a cold star).  Stored as type-parameter fields ⇒ type-stable and
allocation-free in the per-cell hot loop.
"""
struct BarotropeThermo{P,S} <: MHDThermo
    pfun::P
    sTfun::S
end
@inline mhd_pressure(t::BarotropeThermo, ε) = t.pfun(ε)
@inline mhd_Tslope(t::BarotropeThermo, ε)   = t.sTfun(ε)

Base.@kwdef struct BDNKMHDCoeffs{T<:MHDThermo}
    Du::Float64       # D_u  (sets η = w D_u, ζ = 2 w D_u/3)
    Dε::Float64       # D_ε  (sets σ = w D_ε)
    rb::Float64       # resistivity control r_b
    τu::Float64       # heat/momentum BDNK relaxation
    τX::Float64       # shear/pressure BDNK relaxation
    τb::Float64       # magnetic BDNK relaxation
    w::Float64 = 4/3  # reference enthalpy w = ε+p (dispersion analyses; =4/3 ε conformal)
    b2::Float64 = 0.0 # local b² (for the isotropic resistivity maps)
    thermo::T = ConformalThermo()   # per-cell ε↔p closure (constitutive)
end

"""
    bdnk_coeffs_from_Dmaps(; Du, Dε, rb, τu, τX, τb, ε=1.0, b2=0.0)

Build `BDNKMHDCoeffs` and the derived physical coefficients from the
dimensionless D-maps (Eqs. 10,11): σ = w D_ε, η = w D_u, ζ = (2/3) w D_u
(⇒ ζ−⅔η = 0); resistivity (Eq.17) (w+b²)/w · r⊥ = r∥ = r_b, isotropic-limit
r⊥=r∥=r'_b; and τ_ε := 2τ_u (Eq.18).
"""
function bdnk_coeffs_from_Dmaps(; Du, Dε, rb, τu, τX, τb, ε=1.0, b2=0.0,
                                  thermo::MHDThermo=ConformalThermo())
    # reference enthalpy w=ε+p(ε) at the analysis ε (conformal default ⇒ 4/3 ε)
    w = ε + mhd_pressure(thermo, ε)
    return BDNKMHDCoeffs(Du=Du, Dε=Dε, rb=rb, τu=τu, τX=τX, τb=τb,
                         w=w, b2=b2, thermo=thermo)
end

"σ = w D_ε."           sigma(c::BDNKMHDCoeffs)      = c.w*c.Dε
"η = w D_u."           shear_eta(c::BDNKMHDCoeffs)  = c.w*c.Du
"ζ = (2/3) w D_u."     bulk_zeta(c::BDNKMHDCoeffs)  = (2/3)*c.w*c.Du
"τ_ε := 2τ_u (Eq.18)." tau_eps(c::BDNKMHDCoeffs)    = 2*c.τu
"r∥ = r_b (Eq.17)."    r_par(c::BDNKMHDCoeffs)      = c.rb
"r⊥ = w r_b/(w+b²)."   r_perp(c::BDNKMHDCoeffs)     = c.rb*c.w/(c.w+c.b2)
"isotropic r'_b = r_b (r⊥=r∥)."   r_iso(c::BDNKMHDCoeffs) = c.rb

"""
    second_law_ok(c) -> Bool

Second-law positivity of the first-order transport (Eq.6 footnote): η,ζ,r∥ ≥ 0
and r⊥ + b²σ/(ε+p+b²)² ≥ 0.  (η=wD_u, ζ=2wD_u/3, σ=wD_ε with w>0 ⇒ these are
positive whenever D_u,D_ε,r_b ≥ 0.)
"""
function second_law_ok(c::BDNKMHDCoeffs; tol=0.0)
    η = shear_eta(c); ζ = bulk_zeta(c); rpar = r_par(c)
    σ = sigma(c); W = c.w + c.b2
    rperp = r_perp(c)
    cond = (η ≥ -tol) && (ζ ≥ -tol) && (rpar ≥ -tol) &&
           (rperp + c.b2*σ/W^2 ≥ -tol)
    return cond
end

# ===========================================================================
# 2.  SOUND-SECTOR FRONT VELOCITY  (Eq.15)  +  SUBLUMINAL WINDOW  (Eq.16)
# ===========================================================================

"""
    sound_front_W2(Du, Dε, τu, τX; τε=2τu, branch=:plus) -> ComplexF64

Sound-sector front velocity-squared W² = lim_{k→∞}(ω/k)² (Eq.15):

  W² = [ 6 D_u τ_ε + (3τ_X+τ_ε)(τ_u−D_ε)
         ± √( 12τ_ε(2D_u−τ_X)(τ_u−D_ε)² + [6D_uτ_ε+(τ_u−D_ε)(3τ_X+τ_ε)]² ) ]
       / [ 6τ_ε(τ_u−D_ε) ] .

Reality of both branches needs τ_u>D_ε, τ_ε>0, τ_X ≥ 2D_u.  `branch=:plus`
returns the larger root.  Returned as `Complex` so the loss-of-reality regime
(τ_X<2D_u) is representable (a non-zero imaginary part flags it).
"""
function sound_front_W2(Du::Real, Dε::Real, τu::Real, τX::Real;
                        τε::Real=2*τu, branch::Symbol=:plus)
    common = 6*Du*τε + (3*τX+τε)*(τu-Dε)
    disc   = 12*τε*(2*Du-τX)*(τu-Dε)^2 + (6*Du*τε + (τu-Dε)*(3*τX+τε))^2
    rt     = sqrt(complex(disc))
    sgn    = branch === :plus ? 1.0 : -1.0
    den    = 6*τε*(τu-Dε)
    return (common + sgn*rt) / den
end

"""
    sound_subluminal_window(Du, Dε, τu; τε=2τu) -> (lo, hi)

The Eq.16 window in which the larger sound front is real AND subluminal:
  2 D_u ≤ τ_X ≤ τ_ε − D_u − 3 D_u τ_ε/(τ_u − D_ε).
"""
function sound_subluminal_window(Du::Real, Dε::Real, τu::Real; τε::Real=2*τu)
    lo = 2*Du
    hi = τε - Du - 3*Du*τε/(τu-Dε)
    return (lo, hi)
end

"sound_is_subluminal(Du,Dε,τu,τX): the Eq.16 predicate 2Du ≤ τX ≤ hi."
function sound_is_subluminal(Du::Real, Dε::Real, τu::Real, τX::Real;
                             τε::Real=2*τu, tol::Real=1e-12)
    lo, hi = sound_subluminal_window(Du, Dε, τu; τε=τε)
    return (τX ≥ lo - tol) && (τX ≤ hi + tol)
end

# ===========================================================================
# 3.  MAGNETIC TELEGRAPHER  (Eq.9)  — causality predicate + analytic + solver
# ===========================================================================

"""
    telegrapher_causal(τb, rpb) -> Bool

The rest-frame magnetic equation  ḃ − r'_b ∆b + τ_b b̈ = 0  is subluminal /
causal iff  τ_b > r'_b  (Eq.9).  (Its high-k front speed is √(r'_b/τ_b).)
"""
telegrapher_causal(τb::Real, rpb::Real) = τb > rpb

"""
    telegrapher_analytic(t, x, k, τb, rpb) -> Float64

Analytic plane-wave solution of the magnetic telegrapher equation (Eq.21a):
  b(t,x) = exp(−t/2τ_b) sin(kx − Θ_b t),   Θ_b² = (r'_b/τ_b)k² − 1/(4τ_b²).
Requires Θ_b²>0 (propagating regime).
"""
function telegrapher_analytic(t::Real, x::Real, k::Real, τb::Real, rpb::Real)
    Θ2 = (rpb/τb)*k^2 - 1/(4*τb^2)
    Θ  = sqrt(Θ2)
    return exp(-t/(2*τb)) * sin(k*x - Θ*t)
end

"Θ_b² of the telegrapher dispersion (>0 ⇒ propagating)."
telegrapher_Theta2(k, τb, rpb) = (rpb/τb)*k^2 - 1/(4*τb^2)

"""
    solve_telegrapher_1d(N, tend, k, τb, rpb; L=2π, cfl=0.4) -> (x, b, L1err)

A 2nd-order 1D solver of  τ_b ∂_tt b + ∂_t b − r'_b ∂_xx b = 0  on a periodic
domain [0,L].  First-order system  ∂_t b = q,  ∂_t q = (r'_b ∂_xx b − q)/τ_b
with a 2nd-order central Laplacian and a 2nd-order TVD (Heun / SSPRK2)
time-stepper.  Initial data is the analytic mode at t=0; returns the grid, the
final field, and the discrete L1 error against `telegrapher_analytic` at `tend`.
"""
function solve_telegrapher_1d(N::Integer, tend::Real, k::Real, τb::Real, rpb::Real;
                              L::Real=2π, cfl::Real=0.4)
    dx = L/N
    x  = [(i-1)*dx for i in 1:N]
    Θ2 = telegrapher_Theta2(k, τb, rpb)
    Θ  = sqrt(Θ2)
    b  = [sin(k*xi) for xi in x]
    q  = [-1/(2*τb)*sin(k*xi) - Θ*cos(k*xi) for xi in x]   # ∂_t b at t=0
    cmax = sqrt(rpb/τb)                                     # front speed
    dt = cfl*dx/cmax
    nst = max(1, ceil(Int, tend/dt)); dt = tend/nst
    lap = similar(b)
    function laplacian!(out, u)
        @inbounds for i in 1:N
            im1 = i==1 ? N : i-1
            ip1 = i==N ? 1 : i+1
            out[i] = (u[im1] - 2*u[i] + u[ip1])/dx^2
        end
    end
    bp = similar(b); qp = similar(q)
    k1b = similar(b); k1q = similar(q); k2b = similar(b); k2q = similar(q)
    for _ in 1:nst
        # stage 1
        laplacian!(lap, b)
        @inbounds for i in 1:N
            k1b[i] = q[i]
            k1q[i] = (rpb*lap[i] - q[i])/τb
            bp[i]  = b[i] + dt*k1b[i]
            qp[i]  = q[i] + dt*k1q[i]
        end
        # stage 2
        laplacian!(lap, bp)
        @inbounds for i in 1:N
            k2b[i] = qp[i]
            k2q[i] = (rpb*lap[i] - qp[i])/τb
            b[i]  += 0.5*dt*(k1b[i] + k2b[i])
            q[i]  += 0.5*dt*(k1q[i] + k2q[i])
        end
    end
    err = 0.0
    @inbounds for i in 1:N
        err += abs(b[i] - telegrapher_analytic(tend, x[i], k, τb, rpb))
    end
    return x, b, err*dx
end

"""
    telegrapher_convergence(Ns, tend, k, τb, rpb; L=2π) -> (errs, orders)

Run `solve_telegrapher_1d` over a list of resolutions and return the L1 errors
and the pairwise observed convergence orders log2(err_{i-1}/err_i).  Should be
≈2 (Fig.2).
"""
function telegrapher_convergence(Ns, tend, k, τb, rpb; L=2π)
    errs = [solve_telegrapher_1d(N, tend, k, τb, rpb; L=L)[3] for N in Ns]
    orders = [log2(errs[i-1]/errs[i]) for i in 2:length(errs)]
    return errs, orders
end

# ===========================================================================
# 4.  ALFVÉN SECTOR  (App. B, Eq. B5;  τ_ε := 2τ_u imposed)
# ===========================================================================

"""
    alfven_x(c::BDNKMHDCoeffs, θ) -> (x_plus, x_minus)

The Alfvén front-velocity² roots x_A^± (Eq.B5), with 𝒲≡w+b², τ_ε:=2τ_u, at
propagation angle θ to **b** (b∥=b cosθ, b⊥=b sinθ):

  α_A = 𝒲 τ_u − w D_ε,   r_A = r_b(w+b⊥²)/𝒲,
  K   = [α_A r_A + w D_u τ_b + τ_u τ_b b∥²] / (2 τ_b α_A),
  x_A^± = K ± √( K² − (w D_u r_A)/(τ_b α_A) ).

Returned as a `(ComplexF64, ComplexF64)`; W_A = ±√x_A.
"""
function alfven_x(c::BDNKMHDCoeffs, θ::Real)
    b  = sqrt(max(c.b2, 0.0))
    bp = b*cos(θ); bo = b*sin(θ)
    W  = c.w + c.b2                       # 𝒲
    αA = W*c.τu - c.w*c.Dε
    rA = c.rb*(c.w + bo^2)/W
    K  = (αA*rA + c.w*c.Du*c.τb + c.τu*c.τb*bp^2)/(2*c.τb*αA)
    disc = K^2 - (c.w*c.Du*rA)/(c.τb*αA)
    rt = sqrt(complex(disc))
    return (K + rt, K - rt)
end

"alfven_W(c,θ): the four Alfvén front velocities ±√x_A^±."
function alfven_W(c::BDNKMHDCoeffs, θ::Real)
    xp, xm = alfven_x(c, θ)
    sp = sqrt(complex(xp)); sm = sqrt(complex(xm))
    return (sp, -sp, sm, -sm)
end

# ===========================================================================
# 5.  MAGNETOSONIC SECTOR  (App. B, Eq. B9)
#     det M_ms(x)=0  →  quartic a4 x⁴ + a3 x³ + a2 x² + a1 x + a0
# ===========================================================================
#
# M_ms is built with entries that are polynomials in s=√x.  We represent each
# entry as a coefficient vector in powers of s (index k ↔ s^{k-1}), expand the
# 4×4 determinant as a polynomial in s by Laplace cofactors, verify the odd
# powers cancel (the √x-grading), and read off the even-s coefficients as the
# quartic-in-x coefficients.

# --- small polynomial-in-s helpers ---
@inline _padd(a, b) = begin
    n = max(length(a), length(b)); c = zeros(n)
    @inbounds for i in eachindex(a); c[i] += a[i]; end
    @inbounds for i in eachindex(b); c[i] += b[i]; end
    c
end
@inline _psub(a, b) = _padd(a, -b)
@inline function _pmul(a, b)
    c = zeros(length(a)+length(b)-1)
    @inbounds for i in eachindex(a), j in eachindex(b)
        c[i+j-1] += a[i]*b[j]
    end
    c
end
_cst(v) = [float(v)]
_smono(v, p) = (c = zeros(p+1); c[p+1] = float(v); c)   # v·s^p

# 3×3 determinant of polynomial-in-s entries
function _det3(A)
    t1 = _pmul(A[1,1], _psub(_pmul(A[2,2],A[3,3]), _pmul(A[2,3],A[3,2])))
    t2 = _pmul(A[1,2], _psub(_pmul(A[2,1],A[3,3]), _pmul(A[2,3],A[3,1])))
    t3 = _pmul(A[1,3], _psub(_pmul(A[2,1],A[3,2]), _pmul(A[2,2],A[3,1])))
    _padd(_psub(t1, t2), t3)
end

# 4×4 determinant (Laplace along row 1) of polynomial-in-s entries
function _det4(M)
    acc = zeros(1)
    @inbounds for j in 1:4
        rows = (2,3,4); cols = ntuple(k -> k < j ? k : k+1, 3)
        minor = Array{Vector{Float64}}(undef, 3, 3)
        for ri in 1:3, ci in 1:3
            minor[ri,ci] = M[rows[ri], cols[ci]]
        end
        sgn = isodd(j) ? 1.0 : -1.0
        acc = _padd(acc, _pmul(_cst(sgn), _pmul(M[1,j], _det3(minor))))
    end
    acc
end

# Build M_ms(s) (Eq.B9).  Entry (1,3) is taken in its √x-graded reading
# (wD_ε − 3τ_u(w+b⊥²))√x so the determinant is a genuine quartic in x.
function _build_Mms(τu, τX, τb, De, Du, w, bp, bo, rb)
    Wc = w + bp^2 + bo^2
    bp2 = bp^2; bo2 = bo^2
    M = Array{Vector{Float64}}(undef, 4, 4)
    # Row 1
    M[1,1] = _padd(_cst((τu-De)/3), _smono(2*τu, 2))           # (τu−De)/3 + 2τu x
    M[1,2] = _smono(3*τu*bp*bo, 1)                             # 3τu b∥b⊥ √x
    M[1,3] = _smono(w*De - 3*τu*(w+bo2), 1)                    # (wDe − 3τu(w+b⊥²))√x
    M[1,4] = _padd(_cst(τu*bo), _smono(2*τu*bo, 2))            # τu b⊥ (1+2x)
    # Row 2
    M[2,1] = _cst(0.0)
    M[2,2] = _padd(_cst(-w*Du), _smono(τu*(w+bp2) - w*De, 2))  # −wDu + (τu(w+b∥²)−wDe) x
    M[2,3] = _smono(-τu*bp*bo, 2)                             # −τu b∥b⊥ x
    M[2,4] = _smono(τu*bp, 1)                                 # τu b∥ √x
    # Row 3
    M[3,1] = _smono((De - τu - 3*τX)/3, 1)                    # (De−τu−3τX)/3 √x
    M[3,2] = _padd(_cst(-bp*bo*τX), _smono(-τu*bp*bo, 2))     # −b∥b⊥τX − τu b∥b⊥ x
    M[3,3] = _padd(_cst(τX*(w+bo2) - 2*w*Du),
                   _smono(τu*(w+bo2) - w*De, 2))              # τX(w+b⊥²)−2wDu + (…)x
    M[3,4] = _smono(-bo*(τu+τX), 1)                          # −b⊥(τu+τX)√x
    # Row 4
    M[4,1] = _cst(rb*bo/(3*Wc))                              # rb b⊥/(3𝒲)
    M[4,2] = _smono(τb*bp, 1)                                # τb b∥ √x
    M[4,3] = _smono(-τb*bo, 1)                               # −τb b⊥ √x
    M[4,4] = _padd(_cst(-w*rb/Wc), _smono(τb, 2))            # −w rb/𝒲 + τb x
    return M
end

"""
    magnetosonic_quartic(c::BDNKMHDCoeffs, θ) -> (a0,a1,a2,a3,a4), oddmax

Quartic coefficients of det M_ms(x)=0 (Eq.B9) at propagation angle θ:
returns `(a0,a1,a2,a3,a4)` with poly a4 x⁴+a3 x³+a2 x²+a1 x+a0, and the max
absolute odd-power-of-√x coefficient `oddmax` (a numerical witness that the
√x-grading holds — should be ≈0).
"""
function magnetosonic_quartic(c::BDNKMHDCoeffs, θ::Real)
    b  = sqrt(max(c.b2, 0.0))
    bp = b*cos(θ); bo = b*sin(θ)
    d  = _det4(_build_Mms(c.τu, c.τX, c.τb, c.Dε, c.Du, c.w, bp, bo, c.rb))
    getc(k) = (k+1 ≤ length(d)) ? d[k+1] : 0.0
    oddmax = 0.0
    for k in 1:2:length(d)-1
        oddmax = max(oddmax, abs(d[k+1]))
    end
    coeffs = (getc(0), getc(2), getc(4), getc(6), getc(8))   # a0..a4
    return coeffs, oddmax
end

# roots of a4 x⁴ + a3 x³ + a2 x² + a1 x + a0 via the companion matrix
function _quartic_roots(a0, a1, a2, a3, a4)
    if abs(a4) < 1e-300
        return ComplexF64[]
    end
    c = (a0/a4, a1/a4, a2/a4, a3/a4)   # monic x⁴ + c3 x³ + c2 x² + c1 x + c0
    C = [0.0 0.0 0.0 -c[1];
         1.0 0.0 0.0 -c[2];
         0.0 1.0 0.0 -c[3];
         0.0 0.0 1.0 -c[4]]
    return ComplexF64.(eigvals(C))   # always complex: front velocity² may be <0 or complex
end

"""
    magnetosonic_x(c::BDNKMHDCoeffs, θ) -> Vector{ComplexF64}

The four magnetosonic front-velocity² roots x_ms (roots of det M_ms(x)=0).
"""
function magnetosonic_x(c::BDNKMHDCoeffs, θ::Real)
    coeffs, _ = magnetosonic_quartic(c, θ)
    return _quartic_roots(coeffs...)
end

"magnetosonic_W(c,θ): the magnetosonic front velocities ±√x_ms."
function magnetosonic_W(c::BDNKMHDCoeffs, θ::Real)
    Ws = ComplexF64[]
    for x in magnetosonic_x(c, θ)
        s = sqrt(complex(x)); push!(Ws, s); push!(Ws, -s)
    end
    return Ws
end

# ===========================================================================
# 6.  v_max  (Eq. B13):  max over θ and over ALL branches of |Re W|
# ===========================================================================

"""
    front_velocity_max(c::BDNKMHDCoeffs; nθ=361, sectors=(:alfven,:magnetosonic))
        -> NamedTuple

Scan propagation angle θ∈[0,π/2] over all front-velocity branches (Alfvén +
magnetosonic) on the background encoded by `c` (its b²) and return
  • `vmax`      = max over θ, branches of |Re W|         (Eq. B13);
  • `max_imW`   = max over θ, branches of |Im W|         (anti-diffusive growth);
  • `complex_front` = vmax + i·max_imW                   (the paper's combined
        diagnostic — a non-zero imaginary part signals an UNSTABLE state);
  • `W_at_vmax`, `θ_at_vmax`, `θ_at_imax`.
A non-zero `max_imW` is the signature the paper uses to declare a state
acausal/anti-diffusive (e.g. Orszag–Tang with τ_X=0).
"""
function front_velocity_max(c::BDNKMHDCoeffs; nθ::Integer=361,
                            sectors=(:alfven, :magnetosonic))
    vmax = 0.0;  Wbest = 0.0 + 0.0im;  θbest = 0.0
    imax = 0.0;  θimax = 0.0
    for θ in range(0, π/2; length=nθ)
        Ws = ComplexF64[]
        if :alfven in sectors
            xp, xm = alfven_x(c, θ)
            for x in (xp, xm); push!(Ws, sqrt(complex(x))); end
        end
        if :magnetosonic in sectors
            for x in magnetosonic_x(c, θ); push!(Ws, sqrt(complex(x))); end
        end
        for W in Ws
            rW = abs(real(W)); iW = abs(imag(W))
            if rW > vmax; vmax = rW; Wbest = W; θbest = θ; end
            if iW > imax; imax = iW; θimax = θ; end
        end
    end
    return (vmax=vmax, max_imW=imax, complex_front=complex(vmax, imax),
            W_at_vmax=Wbest, θ_at_vmax=θbest, θ_at_imax=θimax)
end

# ===========================================================================
# 7.  TEST-BACKGROUND HELPERS  (KH Eq.23, Orszag–Tang Eq.27)
# ===========================================================================

"""
    kh_initial_b(; Jtx=0.08, Jty=0.0, Jtz=0.8, ux=0.15) -> b

The magnitude b=√b² of the KH initial magnetic field (Eq.23) at the shear-layer
bulk (where u^x≈0.15, u^y≈0): from J^{ti} and u^i, b^t=J^{ti}u_i,
b^i=(J^{ti}+b^t u^i)/Γ, b²=−(b^t)²+b·b.
"""
function kh_initial_b(; Jtx=0.08, Jty=0.0, Jtz=0.8, ux=0.15)
    u = [ux, 0.0, 0.0]; Γ = sqrt(1 + dot(u,u))   # u^i given directly ⇒ Γ=√(1+|u|²)
    J = [Jtx, Jty, Jtz]
    bt = dot(J, u)                                # b^t = J^{ti} u_i
    bvec = (J .+ bt .* u) ./ Γ
    b2 = -bt^2 + dot(bvec, bvec)
    return sqrt(max(b2, 0.0))
end

"""
    ot_b2_at(x, y; Γvinit=0.8) -> b2

Orszag–Tang b² at point (x,y) (Eq.27): (u^x,u^y)=(−Γv sin y, Γv sin x),
(J^{tx},J^{ty})=(−sin y, sin 2x), b^t=J^{ti}u_i, b^i=(J^{ti}+b^t u^i)/Γ.
"""
function ot_b2_at(x::Real, y::Real; Γvinit::Real=0.8)
    u = [-Γvinit*sin(y), Γvinit*sin(x), 0.0]; Γ = sqrt(1 + dot(u,u))
    J = [-sin(y), sin(2x), 0.0]
    bt = dot(J, u)
    bvec = (J .+ bt .* u) ./ Γ
    return -bt^2 + dot(bvec, bvec)
end

"ot_max_b(; n=200): the max b over the OT periodic [0,2π]² grid."
function ot_max_b(; n::Integer=200, Γvinit::Real=0.8)
    b2max = 0.0
    for x in range(0, 2π; length=n), y in range(0, 2π; length=n)
        b2max = max(b2max, ot_b2_at(x, y; Γvinit=Γvinit))
    end
    return sqrt(max(b2max, 0.0))
end

end # module BDNKMHD
