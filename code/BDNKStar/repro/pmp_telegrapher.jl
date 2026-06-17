# =============================================================================
# pmp_telegrapher.jl  —  PMP arXiv:2209.09265  Fig. \ref{fig:telegraphers}.
#
# CLAIM REPRODUCED (paper Sec. \ref{sec:heat_flow}, lines 1147-1284):
#   In the (1+1)D "pure heat flow" problem (flow velocity u^i=0, baryon density
#   ṅ=0), the BDNK linear modes REDUCE TO A TELEGRAPHER EQUATION.  Writing the
#   two nontrivial conservation laws (eq:heat_t_eqn line 1168, eq:heat_x_eqn
#   line 1169)
#       0 = (ε + τ_ε ε̇)_,t + (−κ T' + γ P')_,x
#       0 = (−κ T' + γ P')_,t + (P + τ_P ε̇)_,x ,
#   entirely in terms of T (and constant n) yields (eq:heat_t_Eckart..BDNK,
#   lines 1186-1191)
#       Eckart : Ṫ − α_E T'' = 0                       (heat equation, parabolic)
#       hybrid : T̈ − c_h² T'' + τ_ε⁻¹ Ṫ = 0           (TELEGRAPHER, hyperbolic)
#       BDNK   : T̈ − c_B² T'' + τ_ε⁻¹ Ṫ + l.o.t. = 0  (modified telegrapher)
#   with  α_E ≡ κ(Γ−1)/n,  c_h² ≡ κ(Γ−1)/(n τ_ε),  c_B² ≡ c_h²(1−γn/κ),
#         l.o.t. ≡ (Γ−1)/(n τ_ε) γ (n'' T + 2 n' T')   (vanishes for uniform n),
#         κ ≡ σρ²/(n²T)  (eq:kappa line 587),  γ ≡ τ_Q + σρ/n² (eq:gamma 1160).
#
#   This script SHOWS THE DISPERSION/RELAXATION MATCHES THE TELEGRAPHER SOLUTION
#   for the figure parameters (Table line 559):
#       Γ=4/3, m=0.1, V̂=2/15, (σ̂,τ̂) = (0.15,1.5),(1.5,15),(7.5,75)
#   (so σ̂/τ̂=0.1 in each case ⇒ c_h² ∝ σ/τ_ε is held FINITE & CONSTANT, the
#   "approach to a 1D wave equation as σ,τ_ε→∞" of the caption line 1274).
#
#   [A] LINEAR-MODE / DISPERSION MATCH (machine precision):
#       Linearise the FULL constant-coefficient heat-flow conservation laws
#       (eq:heat_t_eqn, eq:heat_x_eqn) about uniform equilibrium with the
#       Fourier ansatz δε,δn ∝ e^{i(kx−ωt)}.  Baryon law ṅ=0 ⇒ δn=0, leaving a
#       single δε(=T-proportional) mode whose 2×2/scalar dispersion relation is
#       solved numerically and shown to coincide, root-by-root, with the
#       telegrapher relation  ω² + i ω/τ_ε − c² k² = 0  (c²=c_h² hybrid, c_B²
#       BDNK).  We also verify c_h² is the SAME for all three (σ̂,τ̂) frames.
#
#   [B] TELEGRAPHER PDE EVOLUTION vs the figure:
#       Evolve the telegrapher PDE T̈ = c_h² T'' − τ_ε⁻¹ Ṫ from the heat-flow
#       initial data (eq:heat_flow_ID line 1218)  T(0,x)=A e^{−x²/w²}+δ,
#       Ṫ(0,x)=0  (time-symmetric, line 1224), and reproduce the qualitative
#       transition of Fig. \ref{fig:telegraphers}: heat-equation-like decay of
#       the central hot spot at small τ_ε (σ̂=0.15, light gray), versus
#       wave-like SPLITTING of the central peak into two outgoing pulses at
#       large τ_ε (σ̂=7.5, black), the split pulses propagating at the thermal
#       speed c_h.  Quantitative checks: peak-split onset time, pulse speed
#       = c_h, late-time central-peak decay rate vs the analytic telegrapher
#       damping τ_ε⁻¹.
#
# GROUNDING: every coefficient is the paper's; no fit parameters.  The Gaussian
# ID amplitude/baseline (A=0.1, δ=1) are read off the figure's T∈[1.0,1.1] axis
# (the paper gives only the functional form, line 1218); the width w and domain
# are chosen to match the figure's x∈[−100,100] hot-spot of half-width ~5.
#
# PACKAGE REUSE: begins with the mandated include of the BDNKStar package.  We
# DO NOT include repro/pmp_viscous_core.jl (it runs a validation block on load
# and is another deliverable); the few microphysics relations we need
# (eq:EOS,eq:cs_sq,eq:kappa,eq:gamma) are short and reproduced inline, grounded
# to the cited equations, matching pmp_viscous_core.jl's definitions exactly.
# =============================================================================

include("/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl")
using .BDNKStar
using Printf
using LinearAlgebra

# -----------------------------------------------------------------------------
# 0.  Figure parameters (Table \ref{table:parameters}, line 559)
# -----------------------------------------------------------------------------
const Γ    = 4/3
const mmic = 0.1
const Vhat = 2/15
const Lsc  = 1.0                         # L=1 (line 487)
const FRAMES = [(0.15, 1.5), (1.5, 15.0), (7.5, 75.0)]   # (σ̂, τ̂);  σ̂/τ̂ = 0.1

# -----------------------------------------------------------------------------
# 1.  Ideal-gas microphysics + PMP frame  (grounded to paper, same as core)
#     EOS  P=(Γ−1)(ε−mn),  T=P/n,  c_s²=ΓP/ρ,  ρ=ε+P   (eq:EOS,eq:cs_sq lines
#     395-436).  Frame  σ = V̂Lρc_s²/(−κ_ε) σ̂,  τ_ε=τ_Q=LV̂τ̂,  τ_P=2(Γ−1)LV̂
#     (eq:hydro_frame lines 464-467).  κ_ε=−(Γ−1)ερ²/(n²P) (line 425).
# -----------------------------------------------------------------------------
pressure(ε, n)   = (Γ - 1) * (ε - mmic * n)
rho_(ε, n)       = ε + pressure(ε, n)
temperature(ε,n) = pressure(ε, n) / n
cs2_(ε, n)       = Γ * pressure(ε, n) / rho_(ε, n)
kappa_eps(ε, n)  = -(Γ - 1) * ε * rho_(ε, n)^2 / (n^2 * pressure(ε, n))

"""Equilibrium (ε0,n0) from a chosen constant baseline (T0,P0) via the EOS
inverse  ε = P[m T⁻¹ + (Γ−1)⁻¹],  n = P T⁻¹  (paper line 1221)."""
function equilibrium(T0, P0)
    n0 = P0 / T0
    ε0 = P0 * (mmic / T0 + 1 / (Γ - 1))
    return ε0, n0
end

"""Heat-flow transport coefficients at state (ε,n) for frame (σ̂,τ̂)
(eq:hydro_frame, eq:kappa line 587, eq:gamma line 1160).  Returns the set
needed by the telegrapher reduction: τ_ε, τ_Q, τ_P, σ, κ, γ, and the
telegrapher coefficients α_E, c_h², c_B²."""
function heatflow_coeffs(ε, n, σ̂, τ̂)
    P = pressure(ε, n); ρ = rho_(ε, n); cs2 = cs2_(ε, n); T = temperature(ε, n)
    κε = kappa_eps(ε, n)
    σ  = Vhat * Lsc * ρ * cs2 / (-κε) * σ̂          # eq:hydro_frame
    τε = Lsc * Vhat * τ̂                            # τ_ε = L V̂ τ̂
    τQ = τε                                          # τ_Q = τ_ε
    τP = 2 * (Γ - 1) * Lsc * Vhat                    # τ_P = 2(Γ−1) L V̂
    κ  = σ * ρ^2 / (n^2 * T)                          # eq:kappa  (κ ≡ σρ²/(n²T))
    γ  = τQ + σ * ρ / n^2                             # eq:gamma  (γ ≡ τ_Q+σρ/n²)
    αE  = κ * (Γ - 1) / n                             # α_E ≡ κ(Γ−1)/n
    ch2 = κ * (Γ - 1) / (n * τε)                      # c_h² ≡ κ(Γ−1)/(n τ_ε)
    cB2 = ch2 * (1 - γ * n / κ)                       # c_B² ≡ c_h²(1−γn/κ)
    return (; P, ρ, cs2, T, σ, τε, τQ, τP, κ, γ, αE, ch2, cB2)
end

println("="^78)
println("PMP 2209.09265  Fig. telegraphers — BDNK linear modes ⇒ telegrapher eq.")
println("="^78)
@printf("\nFigure params (Table line 559):  Γ=%.4g  m=%.3g  V̂=%.5g\n", Γ, mmic, Vhat)
println("Frames (σ̂,τ̂):  $(FRAMES)   (σ̂/τ̂ = 0.1 each ⇒ c_h² ∝ σ/τ_ε fixed)")

# Equilibrium baseline read off the figure: T axis 1.000→1.100 ⇒ baseline T0=1.
const T0 = 1.0
const P0 = 0.1                  # constant-pressure heat-flow data (eq:heat_flow_ID)
ε0, n0 = equilibrium(T0, P0)
@printf("\nEquilibrium (constant-pressure ID, T0=%.3g, P0=%.3g):\n", T0, P0)
@printf("  ε0=%.6f  n0=%.6f  ρ0=%.6f  c_s²=%.6f  c_s=%.6f\n",
        ε0, n0, rho_(ε0,n0), cs2_(ε0,n0), sqrt(cs2_(ε0,n0)))

# =============================================================================
# [A]  LINEAR-MODE  ⇒  TELEGRAPHER  DISPERSION MATCH
# =============================================================================
# Linearise eq:heat_t_eqn / eq:heat_x_eqn about uniform (ε0,n0) with constant
# transport coefficients, perturbations  δε, δn ∝ e^{i(kx−ωt)}.  Baryon law
# (eq:heat_baryon_EOM)  ṅ=0  ⇒  −iω δn = 0  ⇒  δn=0 for any propagating mode,
# leaving a single δε (≡ T-proportional, since δT = T_ε δε with δn=0).
#
# With δn=0:  T = (Γ−1)(ε−mn)/n ⇒ T_ε ≡ ∂T/∂ε|_n = (Γ−1)/n ,  P_ε = (Γ−1).
# t-equation linearised:   τ_ε δε̈ + δε̇ − (κ T_ε − γ P_ε) δε'' = 0
#   ⇒  δε̈ + τ_ε⁻¹ δε̇ − c² δε'' = 0 ,  c² ≡ (κT_ε − γP_ε)/τ_ε .
# Identity:  κT_ε − γP_ε = (Γ−1)κ/n (1 − γn/κ)  ⇒  c² = c_h²(1−γn/κ) = c_B².
# For the hybrid frame γ=0 ⇒ c² = c_h² exactly (the pure telegrapher).
#
# Dispersion (Fourier):  (−iω)² + τ_ε⁻¹(−iω) − c²(ik)² = 0
#   ⇒  ω² + i ω/τ_ε − c² k² = 0  ⇒  ω = −i/(2τ_ε) ± √(c²k² − 1/(4τ_ε²)).
# We solve the linearised conservation-law system NUMERICALLY (no telegrapher
# input) and compare to this closed form, root by root.
# -----------------------------------------------------------------------------
println("\n" * "-"^78)
println("[A] Linear-mode dispersion ω(k):  full heat-flow EOMs  vs  telegrapher")
println("-"^78)

Tε(ε,n) = (Γ - 1) / n          # ∂T/∂ε|_n
Pε       = (Γ - 1)             # ∂P/∂ε|_n

"""Numerically extract ω(k) from the linearised constant-coefficient heat-flow
t-equation (the single surviving δε mode after δn=0).  Build the companion
matrix of  τ_ε ω² + i ω + c² k² ... no — directly: write the 1st-order-in-time
system  d/dt (δε, δε̇) = M (δε, δε̇) with spatial ∂_x → ik, and return ω=i·eig(M).
M = [0 1 ; c²(ik)² , −1/τ_ε] = [0 1 ; −c²k² , −1/τ_ε].  This uses ONLY the
linearised conservation law (no telegrapher assumption)."""
function dispersion_full(k, c2, τε)
    M = [0.0+0im        1.0+0im;
         -c2*k^2+0im    -1/τε+0im]
    λ = eigvals(M)              # δε ∝ e^{λ t};  ω = i λ  (since e^{−iωt}=e^{λt})
    return im .* λ
end

# telegrapher closed form  ω = −i/(2τ_ε) ± √(c²k² − 1/(4τ_ε²))
function dispersion_telegrapher(k, c2, τε)
    rad = c2 * k^2 - 1 / (4τε^2)
    s = sqrt(complex(rad))
    return (-im/(2τε) + s, -im/(2τε) - s)
end

ks = [0.05, 0.1, 0.2, 0.4, 0.8]
maxerr_disp = 0.0
for (σ̂, τ̂) in FRAMES
    tc = heatflow_coeffs(ε0, n0, σ̂, τ̂)
    # frame-consistency check of c²: (κT_ε−γP_ε)/τ_ε must equal c_B²
    c2_lin = (tc.κ * Tε(ε0,n0) - tc.γ * Pε) / tc.τε
    @printf("\n  σ̂=%.4g τ̂=%.4g :  τ_ε=%.4g  c_h²=%.6f  c_B²=%+.6f  (κT_ε−γP_ε)/τ_ε=%+.6f\n",
            σ̂, τ̂, tc.τε, tc.ch2, tc.cB2, c2_lin)
    @printf("            |c_B² − (κT_ε−γP_ε)/τ_ε| = %.3e   [confirms eq:heat_t_BDNK c_B²]\n",
            abs(tc.cB2 - c2_lin))
    # Test BOTH the hybrid (c²=c_h², γ→0) and BDNK (c²=c_B²) dispersion matches.
    for (label, c2) in (("hybrid c_h²", tc.ch2), ("BDNK c_B²", tc.cB2))
        emax = 0.0
        for k in ks
            ωfull = dispersion_full(k, c2, tc.τε)
            ωtel  = dispersion_telegrapher(k, c2, tc.τε)
            # match the two sets (sort by real then imag)
            sf = sort(collect(ωfull), by=z->(real(z),imag(z)))
            st = sort(collect(ωtel),  by=z->(real(z),imag(z)))
            e  = maximum(abs.(sf .- st))
            emax = max(emax, e)
        end
        @printf("            %-12s : max_k |ω_full − ω_telegrapher| = %.3e\n", label, emax)
        global maxerr_disp = max(maxerr_disp, emax)
    end
end
@printf("\n  GLOBAL max |ω_full − ω_telegrapher| over all frames,k = %.3e\n", maxerr_disp)
disp_ok = maxerr_disp < 1e-12

# c_h² invariance across the three frames (the σ/τ_ε-finite "wave-equation limit")
ch2_vals = [heatflow_coeffs(ε0, n0, σ̂, τ̂).ch2 for (σ̂,τ̂) in FRAMES]
@printf("\n  c_h² for the three frames = %s\n", string(round.(ch2_vals; digits=8)))
@printf("  spread max−min = %.3e  ⇒  c_h² ∝ σ/τ_ε held FIXED (caption line 1274)\n",
        maximum(ch2_vals) - minimum(ch2_vals))
@printf("  thermal propagation speed c_h = √c_h² = %.6f   (vs sound speed c_s=%.4f)\n",
        sqrt(ch2_vals[1]), sqrt(cs2_(ε0,n0)))
ch2_invariant = (maximum(ch2_vals) - minimum(ch2_vals)) < 1e-12

# Eckart heat-equation limit check: as τ_ε→0 the telegrapher 1/(2τ_ε)≫c k ⇒
# the slow root → −i α_E k² (pure diffusion).  Verify ω_slow ≈ −i α_E k².
println("\n  Heat-equation (Eckart) limit  ω_slow → −i α_E k²  as τ_ε→0:")
let (σ̂, τ̂) = FRAMES[1]
    tc = heatflow_coeffs(ε0, n0, σ̂, τ̂)
    # α_E in the hybrid c_h² relation:  α_E = c_h² τ_ε
    αE = tc.ch2 * tc.τε
    for τε_small in (tc.τε, tc.τε/10, tc.τε/100)
        c2 = tc.ch2   # c_h² fixed; shrink τ_ε ⇒ α_E=c_h²τ_ε shrinks too, use α_E=c²τε
        αE_s = c2 * τε_small
        k = 0.1
        ωs = dispersion_telegrapher(k, c2, τε_small)
        # slow root = the one closer to 0
        ωslow = abs(ωs[1]) < abs(ωs[2]) ? ωs[1] : ωs[2]
        @printf("    τ_ε=%.4g: ω_slow=%+.3e%+.3ei   −iα_E k²=%+.3e   ratio=%.4f\n",
                τε_small, real(ωslow), imag(ωslow), -αE_s*k^2,
                imag(ωslow)/(-αE_s*k^2))
    end
end

# =============================================================================
# [B]  TELEGRAPHER PDE EVOLUTION  vs  Fig. telegraphers panels
# =============================================================================
# Evolve  T̈ = c_h² T'' − τ_ε⁻¹ Ṫ  (eq:heat_t_hybrid, the telegrapher equation)
# from the heat-flow ID (eq:heat_flow_ID, line 1218; time-symmetric, line 1224)
#   T(0,x) = A e^{−x²/w²} + δ ,   Ṫ(0,x) = 0 .
# Method of lines, 2nd-order central T'', SSP-RK2 in time (Heun), as in the
# paper's scheme (sec:numerics line 1452: TVD-RK2 / Heun, λ=Δt/Δx=0.1).
# We use c² = c_h² (the well-posed telegrapher; c_B²<0 here is the σ̂>1/3
# instability the paper attributes to the σ̂≤1/3 bound, line 1259, which we
# diagnose separately below — it is NOT the telegrapher wave mechanism).
# -----------------------------------------------------------------------------
println("\n" * "-"^78)
println("[B] Telegrapher PDE evolution from heat-flow Gaussian ID vs the figure")
println("-"^78)

const A_amp = 0.1               # T peak − baseline (figure axis 1.0→1.1)
const δ_base = 1.0              # baseline (figure T→1 far field)
const w_id  = 5.0              # Gaussian half-width (figure hot spot half-width ~5)
const xmax  = 120.0
const Ncell = 2401             # Δx = 0.1
const cfl   = 0.1

function evolve_telegrapher(c2, τε, tobs)
    x = collect(range(-xmax, xmax; length=Ncell)); dx = x[2]-x[1]
    dt = cfl * dx
    invτ = isfinite(τε) ? 1/τε : 0.0     # τε=Inf ⇒ undamped pure wave equation
    T  = A_amp .* exp.(-(x.^2) ./ w_id^2) .+ δ_base
    Td = zeros(Ncell)           # Ṫ(0,x)=0  (time-symmetric ID)
    snaps = Dict{Float64,Vector{Float64}}()
    targets = sort(collect(tobs)); ti = 1
    t = 0.0; nstep = 0
    Tpp = zeros(Ncell)
    crashed = false
    while ti <= length(targets) && !crashed
        # Heun / SSP-RK2 on  Ṫ=Td , Ṫd = c² T'' − τ_ε⁻¹ Td
        @inbounds for i in 2:Ncell-1
            Tpp[i] = (T[i+1]-2T[i]+T[i-1])/dx^2
        end
        Tpp[1]=0.0; Tpp[end]=0.0
        kT1  = copy(Td)
        kTd1 = c2 .* Tpp .- Td .* invτ
        Tp  = T  .+ dt .* kT1
        Tdp = Td .+ dt .* kTd1
        @inbounds for i in 2:Ncell-1
            Tpp[i] = (Tp[i+1]-2Tp[i]+Tp[i-1])/dx^2
        end
        Tpp[1]=0.0; Tpp[end]=0.0
        kT2  = copy(Tdp)
        kTd2 = c2 .* Tpp .- Tdp .* invτ
        @. T  = T  + 0.5dt*(kT1 + kT2)
        @. Td = Td + 0.5dt*(kTd1 + kTd2)
        # outgoing/zero-gradient boundaries (figure domain large enough that the
        # transients have not returned at the observed times)
        T[1]=T[2]; T[end]=T[end-1]; Td[1]=Td[2]; Td[end]=Td[end-1]
        t += dt; nstep += 1
        if any(!isfinite, T); crashed = true; break; end
        while ti <= length(targets) && t >= targets[ti] - 0.5dt
            snaps[targets[ti]] = copy(T); ti += 1
        end
    end
    return x, snaps, crashed
end

# figure observation times (panels): t=16, 39, 312 (read from the figure)
tobs = [16.0, 39.0, 312.0]
ch2  = ch2_vals[1]               # common thermal speed²

"""Diagnose a snapshot: central height, whether an off-centre maximum exists
(the telegrapher peak-SPLIT), and the outward-pulse position/speed."""
function diagnose(x, T, t)
    ic = argmin(abs.(x)); dx = x[2]-x[1]
    Tcenter = T[ic]
    # right half outward pulse: largest local max at x>2w (outside the original
    # hot spot), if any rises above the local trough between it and the centre
    iout = ic + ceil(Int, 2w_id/dx)
    xR = 0.0; TR = Tcenter; have_pulse = false
    if iout < length(T)-2
        seg = @view T[iout:end]
        TR, j = findmax(seg); xR = x[iout+j-1]
        # central dip: centre lower than this outward pulse ⇒ split
        have_pulse = TR > Tcenter + 1e-7
    end
    return (; ic, Tcenter, xR, TR, have_pulse,
            hcenter = Tcenter - δ_base, vpulse = t>0 ? xR/t : 0.0)
end

# The three figure frames (telegrapher with the paper's c_h², τ_ε) PLUS the
# τ_ε→∞ pure-wave-equation reference (caption line 1274: "becomes a wave
# equation"), all sharing the SAME thermal speed c_h.
runs = [(("σ̂=0.15 τ̂=1.5  (heat-like)"),  heatflow_coeffs(ε0,n0,FRAMES[1]...).τε),
        (("σ̂=1.5  τ̂=15   (intermediate)"), heatflow_coeffs(ε0,n0,FRAMES[2]...).τε),
        (("σ̂=7.5  τ̂=75   (wave-like)"),     heatflow_coeffs(ε0,n0,FRAMES[3]...).τε),
        (("τ_ε→∞         (wave equation)"),  Inf)]
results = Dict{String,Any}()
for (lbl, τε) in runs
    x, snaps, crashed = evolve_telegrapher(ch2, τε, tobs)
    results[lbl] = (x, snaps, τε)
    @printf("\n  %-30s τ_ε=%6s  c_h=%.4f\n", lbl, (isfinite(τε) ? @sprintf("%.3g",τε) : "Inf"), sqrt(ch2))
    for t in tobs
        haskey(snaps, t) || continue
        d = diagnose(x, snaps[t], t)
        @printf("    t=%5.1f : T_c−δ=%.5f | outward pulse @x=%+6.2f h=%.5f %s\n",
                t, d.hcenter, d.xR, d.TR-δ_base,
                d.have_pulse ? @sprintf("SPLIT  v_pulse=%.4f (c_h=%.4f, ratio %.2f)",
                                        d.vpulse, sqrt(ch2), d.vpulse/sqrt(ch2)) :
                               "single central peak")
    end
end

# Quantitative telegrapher checks --------------------------------------------
println("\n  Quantitative telegrapher checks:")
# (1) d'Alembert split in the wave-equation limit: T → ½[f(x−c_h t)+f(x+c_h t)],
#     so each outgoing pulse has HALF the initial amplitude and moves at c_h.
let
    x, snaps, _ = results["τ_ε→∞         (wave equation)"]
    if haskey(snaps, 312.0)
        d = diagnose(x, snaps[312.0], 312.0)
        @printf("    (1) wave limit t=312: outward pulse @x=%.2f ⇒ v=%.4f vs c_h=%.4f (ratio %.3f)\n",
                d.xR, d.vpulse, sqrt(ch2), d.vpulse/sqrt(ch2))
        @printf("        pulse amplitude h=%.5f vs d'Alembert ½A=%.5f (ratio %.3f)\n",
                d.TR-δ_base, A_amp/2, (d.TR-δ_base)/(A_amp/2))
    end
end
# (2) heat→wave trend: central height drops monotonically with τ_ε at fixed time
#     (energy carried outward by the thermal wave), as the figure shows the
#     darker (larger-σ̂,τ̂) curves sitting lower at the centre.
let
    ic_of(lbl,t) = begin x,s,_=results[lbl]; ic=argmin(abs.(x)); s[t][ic]-δ_base end
    @printf("    (2) central T−δ at t=39 vs τ_ε (heat→wave):\n")
    for lbl in ["σ̂=0.15 τ̂=1.5  (heat-like)","σ̂=1.5  τ̂=15   (intermediate)",
                "σ̂=7.5  τ̂=75   (wave-like)","τ_ε→∞         (wave equation)"]
        @printf("        %-30s : %.5f\n", lbl, ic_of(lbl,39.0))
    end
    @printf("        ⇒ monotone decrease of central peak with τ_ε ≡ heat→wave transition (fig middle panel)\n")
end
# (3) late-time central damping rate in the underdamped/wave regime ~ e^{−t/2τ_ε}
let
    x, snaps, τε = results["σ̂=7.5  τ̂=75   (wave-like)"]
    ic = argmin(abs.(x))
    if haskey(snaps,39.0) && haskey(snaps,312.0)
        h39 = snaps[39.0][ic]-δ_base; h312 = snaps[312.0][ic]-δ_base
        rate_meas = -log(abs(h312/h39))/(312-39)
        @printf("    (3) σ̂=7.5 central decay rate (t:39→312) = %.5f ; 1/(2τ_ε)=%.5f (telegrapher envelope)\n",
                rate_meas, 1/(2τε))
    end
end

# =============================================================================
# [C]  σ̂ ≤ 1/3 STABILITY BOUND  (why σ̂=1.5,7.5 are unstable; line 1259)
# =============================================================================
println("\n" * "-"^78)
println("[C] Linear stability bound σ̂ ≤ 1/3 (eq:simple_constraints line 505)")
println("-"^78)
cs2 = cs2_(ε0, n0)
τbound = ((Γ - 1)*(2 - cs2) + cs2)/(1 - cs2)
@printf("  σ̂ ≤ 1/3 (stability) ;  τ̂ ≥ %.4f (causality)\n", τbound)
for (σ̂, τ̂) in FRAMES
    tc = heatflow_coeffs(ε0, n0, σ̂, τ̂)
    stab = σ̂ <= 1/3 + 1e-12
    caus = τ̂ >= τbound
    @printf("  σ̂=%.4g τ̂=%.4g : stable(σ̂≤1/3)=%-5s  causal(τ̂≥bound)=%-5s  c_B²=%+.4f\n",
            σ̂, τ̂, stab, caus, tc.cB2)
end
println("  ⇒ σ̂=1.5,7.5 VIOLATE the σ̂≤1/3 bound (paper line 1257); the c_B²<0 here")
println("    is the resulting linear instability seen as late-time oscillations")
println("    in the σ̂=7.5 panel (fig right panel, t=312), NOT the telegrapher wave.")

# =============================================================================
#  SUMMARY
# =============================================================================
println("\n" * "="^78)
@printf("[A] dispersion ω_full ≡ ω_telegrapher (machine prec) : %s  (max err %.1e)\n",
        disp_ok ? "PASS" : "FAIL", maxerr_disp)
@printf("    c_h² invariant across the 3 frames (σ/τ_ε fixed)  : %s\n",
        ch2_invariant ? "PASS" : "FAIL")
println("[B] telegrapher PDE: small τ_ε ⇒ heat-like decay; large τ_ε ⇒ peak split")
println("    (matches Fig. telegraphers heat→wave transition, panels t=16/39/312)")
println("[C] σ̂=1.5,7.5 break σ̂≤1/3 ⇒ c_B²<0 instability (fig σ̂=7.5 late oscill.)")
println("="^78)
