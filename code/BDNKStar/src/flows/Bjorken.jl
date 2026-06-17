#=
    Bjorken — boost-invariant (0+1D) flow in Milne coordinates, the PMP
    (2209.09265 Sec. Bjorken_flow) test problem. Ideal-gas microphysics
    P = (Γ-1) m n e with ε = m n (1+e), n(τ)=n_0/τ, so P = (Γ-1)(ε − m n_0/τ).

    Inviscid baseline (the reference underlying fig:bjorken): the τ-component of
    ∇_a T^{ab}=0 gives  dε/dτ = −(ε+P)/τ,  with exact solution
        ε(τ) = m n_0 τ^{-1} [ 1 + e_0 τ^{-(Γ-1)} ]   (eq:inviscid_bjorken),
    and the frame-independent diagnostic  ε' + Γ ε/τ = (Γ-1) m n_0 / τ².
    Integrated with classical RK4 ⇒ 4th-order self-convergence (Q_N → 16).
=#
module Bjorken

export bjorken_pressure, bjorken_inviscid_analytic, bjorken_diagnostic,
       bjorken_evolve_rk4

@inline bjorken_pressure(ε, τ; Γ=4/3, m=1.0, n0=1.0) = (Γ-1)*(ε - m*n0/τ)

bjorken_inviscid_analytic(τ; Γ=4/3, m=1.0, n0=1.0, e0=1.0) =
    m*n0/τ * (1 + e0*τ^(-(Γ-1)))

"""diagnostic ε' + Γε/τ; analytically = (Γ-1) m n0 / τ² on the inviscid flow."""
bjorken_diagnostic_exact(τ; Γ=4/3, m=1.0, n0=1.0) = (Γ-1)*m*n0/τ^2

"""
    bjorken_evolve_rk4(τ0, τf, ε0; N, Γ, m, n0) -> (τs, εs)

Classical RK4 of dε/dτ = −(ε+P)/τ from τ0 to τf in N steps.
"""
function bjorken_evolve_rk4(τ0, τf, ε0; N=2000, Γ=4/3, m=1.0, n0=1.0)
    h = (τf - τ0)/N
    τs = Vector{Float64}(undef, N+1); εs = similar(τs)
    τs[1] = τ0; εs[1] = ε0
    f(τ, ε) = -(ε + bjorken_pressure(ε, τ; Γ=Γ, m=m, n0=n0))/τ
    τ = τ0; ε = ε0
    for i in 1:N
        k1 = f(τ, ε)
        k2 = f(τ+h/2, ε+h/2*k1)
        k3 = f(τ+h/2, ε+h/2*k2)
        k4 = f(τ+h,   ε+h*k3)
        ε += h/6*(k1+2k2+2k3+k4); τ += h
        τs[i+1] = τ; εs[i+1] = ε
    end
    return τs, εs
end

# convenience: diagnostic ε'+Γε/τ from a numeric (τ,ε) using the ODE rhs for ε'
bjorken_diagnostic(τ, ε; Γ=4/3, m=1.0, n0=1.0) =
    -(ε + bjorken_pressure(ε, τ; Γ=Γ, m=m, n0=n0))/τ + Γ*ε/τ

end # module Bjorken
