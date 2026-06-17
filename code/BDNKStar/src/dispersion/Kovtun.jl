#=
    Kovtun — linearized dispersion relations of first-order (general-frame)
    relativistic hydrodynamics (Kovtun 1907.08191, "First-order relativistic
    hydrodynamics is stable"). Reproduces the paper's analytic mode figures.

    Phase velocity of a linearly-dispersing wave in a moving fluid (eq:cvphi):
      c_v(φ) = v0(1-c0²)/(1-c0²v0²) cosφ ± c0/(1-c0²v0²) √{(1-v0²)[1-v0²c0²-v0²(1-c0²)cos²φ]}

    Shear-channel eigenfrequencies (eq. line 394), units w0=η=1 (ω,k in (ε0+p0)/η):
      (θ-v0²η)ω² + (i w0/γ0 - 2(θ-η)k·v0)ω - (i w0/γ0)(k·v0) - k²η/γ0² + (k·v0)²(θ-η) = 0
    with k·v0 = k v0 cosφ, γ0 = 1/√(1-v0²). Stable iff Im ω ≤ 0 for all k (true
    for θ/η>0): the central claim of the paper.
=#
module Kovtun

export kovtun_cv, kovtun_shear_modes, kovtun_shear_speed

"""c_v(φ) phase velocity (± branches) of a linear mode (c0) in a fluid moving at v0."""
function kovtun_cv(v0::Real, φ::Real; c0::Real=0.5)
    den = 1 - c0^2*v0^2
    pre = v0*(1-c0^2)/den * cos(φ)
    arg = (1-v0^2)*(1 - v0^2*c0^2 - v0^2*(1-c0^2)*cos(φ)^2)
    rad = c0/den * sqrt(max(arg, 0.0))
    return pre+rad, pre-rad
end

"""Shear-channel eigenfrequencies ω±(k,φ) as complex roots of the F_shear quadratic."""
function kovtun_shear_modes(k::Real, φ::Real; v0::Real=0.9, θη::Real=2.0)
    η = 1.0; w0 = 1.0; θ = θη*η; γ0 = 1/sqrt(1-v0^2)
    kv = k*v0*cos(φ)
    a = θ - v0^2*η
    b = im*w0/γ0 - 2*(θ-η)*kv
    c = -(im*w0/γ0)*kv - k^2*η/γ0^2 + kv^2*(θ-η)
    s = sqrt(b^2 - 4*a*c + 0im)
    return (-b + s)/(2a), (-b - s)/(2a)
end

"""Large-k shear phase velocity c_shear(φ) (eq:cs-shear), the ± roots of the quadratic
   (θ-v0²η)c² - 2v0cosφ(θ-η)c + v0²(θcos²φ+ηsin²φ) - η = 0."""
function kovtun_shear_speed(φ::Real; v0::Real=0.9, θη::Real=2.0)
    η = 1.0; θ = θη*η
    a = θ - v0^2*η
    b = -2*v0*cos(φ)*(θ-η)
    c = v0^2*(θ*cos(φ)^2 + η*sin(φ)^2) - η
    s = sqrt(b^2 - 4a*c + 0im)
    return (-b+s)/(2a), (-b-s)/(2a)
end

end # module Kovtun
