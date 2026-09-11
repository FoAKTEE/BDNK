using Test
using BDNKStar
using LinearAlgebra

# FLAT-SPACE PLANE-WAVE VALIDATION of the 3+1D full-frame BDNK engine, against the
# CLOSED-FORM BDNK dispersion relations. This is the only check on CowlingBDNK3D
# that is independent of the derivations that produced the code: the prediction
# comes from the BDNK constitutive relations, the measurement from calling the
# shipped `_rhs!` on an exact eigenvector.
#
# SETUP. The Cowling background is FROZEN and arbitrary, so it may be flattened
# (α=1, γ_ij=δ_ij, Φ′=0) with uniform w₀,c_s²,η,ζ,τ and every regularizer off
# (σ_ko=0, excise_frac=0, den_frac→0). Then `_rhs!` is a REAL linear operator L
# with exact plane-wave eigenmodes, tested POINTWISE in the deep interior — no
# boundary conditions enter at all.
#
# WHY IT IS EXACT AND NOT O(h²): the code composes two CENTRAL differences, so the
# exact discrete symbol is ∂_x → i·κ with κ = sin(kh)/h. Using κ (not k) in the
# dispersion relations makes these identities hold to MACHINE PRECISION.
#
#   SHEAR (transverse, δv_y ∝ e^{ikx}):   D s² + w₀ s + η κ² = 0,  D = τ_Q w₀ − η
#   SOUND (longitudinal, 4×4 in δE,δS,δε,δv):
#       s δE = −iκ δS
#       s δS = −iκ c_s² δE − κ²(ζ + 4η/3) δv
#       s δε = (δE−δε)/τ − iκ w₀ δv
#       s δv = (δS − w₀δv − iκ τ c_s² δε)/D
#   The δS equation contains NO δε — that is the Π-closure τ_P = c_s²τ_ε, which
#   makes the isotropic stress depend on the CONSERVED energy only (iso = c_s²δE).
#
# COVERAGE: shear tensor construction, momentum stress divergence, bulk ζ, heat
# flux τ_Q c_s²∂δε, frame times τ_ε/τ_Q, the Π closure, and both primitive
# recovery equations. NOT covered: the curved-space O(2M/R) terms, which vanish
# identically in flat space (Φ′ = gΓ1 = gΓ2 = qinv = A²−1 = 0).

const CB = BDNKStar.CowlingBDNK3D

# uniform state of the flat box; τ_Q w₀ − η > 0 is the BDN causality condition (a)
const W0 = 1.0; const CS2 = 0.25; const ETA = 0.05; const ZETA = 0.02; const TAU = 0.30
const DEN = TAU * W0 - ETA
const NPW = 32

"""Flatten the frozen background and impose the uniform BDNK coefficients."""
function _flat_box(; α::Float64=1.0)
    s = build_star3d(ShumPolytrope(100.0), 0.00128 + 100*0.00128^2; N=NPW, Lfac=1.2)
    s.α .= α; s.e2λ .= 1.0; s.sqrtγ .= 1.0; s.interior .= true
    e = setup_bdnk3d(s; η̂=0.0, ζ̂=0.0, τ̂=TAU, σ_ko=0.0,
                     den_frac=1e-14, excise_frac=0.0)
    e.Φp .= 0.0; e.qinv .= 0.0; e.gΓ1 .= 0.0; e.gΓ2 .= 0.0
    e.w0 .= W0; e.cs2 .= CS2; e.η0 .= ETA; e.ζ0 .= ZETA
    e.τε .= TAU; e.τP .= CS2*TAU; e.τQ .= TAU
    s, e
end

_flds(x) = (x.δE, x.δSx, x.δSy, x.δSz, x.δε, x.δvx, x.δvy, x.δvz)

"""max pointwise |L q − s q| / max|s q| over the deep interior, for q = amp·e^{ikx}.
`_rhs!` is real-linear, so for a complex mode we form L(Re q) + i·L(Im q)."""
function _residual(s, e, amp::Vector{ComplexF64}, k::Float64, sval::ComplexF64)
    N = NPW; xs = s.grid.x
    st = BDNKState(N); d = BDNKState(N); scr = CB.Scratch(N)
    out = map((:re, :im)) do part
        for A in _flds(st); fill!(A, 0.0); end
        for kk in 1:N, jj in 1:N, ii in 1:N
            ph = cis(k*xs[ii])
            for (fi, A) in enumerate(_flds(st))
                A[ii,jj,kk] = part === :re ? real(amp[fi]*ph) : imag(amp[fi]*ph)
            end
        end
        CB._rhs!(d, st, e, scr)
        [copy(A) for A in _flds(d)]
    end
    R, I = out
    num = 0.0; den = 0.0
    for kk in 6:N-5, jj in 6:N-5, ii in 6:N-5, fi in 1:8
        Lq  = R[fi][ii,jj,kk] + im*I[fi][ii,jj,kk]
        tgt = sval * amp[fi] * cis(k*xs[ii])
        num = max(num, abs(Lq - tgt)); den = max(den, abs(tgt))
    end
    num / den
end

"""the two shear branches at wavenumber k, and their eigenvectors δS=(sD+w₀)δv"""
function _shear_modes(h, k)
    κ = sin(k*h)/h
    rt = sqrt(complex(W0^2 - 4*DEN*ETA*κ^2))
    map(((-W0+rt)/(2DEN), (-W0-rt)/(2DEN))) do sv
        amp = zeros(ComplexF64, 8)
        amp[7] = 1.0 + 0im          # δv_y
        amp[3] = sv*DEN + W0        # δS_y
        (sv, amp)
    end
end

"""the 4×4 longitudinal system; returns (s, amp) for all four branches"""
function _sound_modes(h, k)
    κ = sin(k*h)/h; ik = im*κ
    M = ComplexF64[  0        -ik          0             0
                  -ik*CS2      0           0     -κ^2*(ZETA + 4ETA/3)
                   1/TAU       0        -1/TAU        -ik*W0
                     0       1/DEN  -ik*TAU*CS2/DEN    -W0/DEN ]
    F = eigen(M)
    map(1:4) do m
        ev = F.vectors[:, m]; amp = zeros(ComplexF64, 8)
        amp[1] = ev[1]; amp[2] = ev[2]; amp[5] = ev[3]; amp[6] = ev[4]
        (F.values[m], amp)
    end
end

@testset "3+1D BDNK: flat-space plane waves reproduce the BDNK dispersion relations" begin
    s, e = _flat_box()
    h = s.grid.dx

    # POSITIVE CONTROL on the setup: if any regularizer were live the modes below
    # would not be eigenmodes at all, and the test would be measuring something else.
    @test e.wfloor == 0.0                       # nothing excised
    @test e.σ_ko == 0.0                         # no Kreiss–Oliger damping
    @test e.den_floor < 1e-12                   # recovery denominator unclamped
    @test DEN > 0                               # BDN causality condition (a)
    @test all(==(1.0), s.sqrtγ) && all(==(1.0), s.α)

    @testset "shear channel (transverse)" begin
        for k in (0.15, 0.35, 0.7, 1.2), (sv, amp) in _shear_modes(h, k)
            @test _residual(s, e, amp, k, sv) < 1e-9
        end
        # hydrodynamic branch → the Navier–Stokes shear rate −ηk²/w₀ as k→0, with
        # the causal BDNK correction growing with k (0.2% at k=0.15, 14% at k=1.2)
        s_hydro(k) = real(_shear_modes(h, k)[1][1])
        @test isapprox(s_hydro(0.15), -ETA*0.15^2/W0; rtol=5e-3)
        @test s_hydro(0.15) > -ETA*0.15^2/W0     # BDNK damps LESS than NS
        # non-hydrodynamic branch is the frame mode at −w₀/D
        @test isapprox(real(_shear_modes(h, 0.15)[2][1]), -W0/DEN; rtol=1e-3)
    end

    @testset "sound channel (longitudinal)" begin
        for k in (0.35, 0.7), (sv, amp) in _sound_modes(h, k)
            @test _residual(s, e, amp, k, sv) < 1e-9
        end
        # the acoustic pair: phase speed → c_s, damping → the NS value
        k = 0.35; κ = sin(k*h)/h
        ac = [v for (v, _) in _sound_modes(h, k) if abs(imag(v)) > 1e-8 && real(v) > -1.0]
        @test length(ac) == 2
        @test isapprox(abs(imag(ac[1]))/κ, sqrt(CS2); rtol=2e-3)
        @test isapprox(-real(ac[1]), (ZETA + 4ETA/3)*κ^2/(2W0); rtol=1e-2)
    end

    # A UNIFORM lapse is still flat (Φ′=0), but the BDNK constitutive relations are
    # covariant — D = u^μ∂_μ = (1/α)∂_t — so the ENTIRE rhs must scale by exactly α
    # and every eigenvalue by α. This discriminates the recovery lapse specifically:
    # putting α on the conserved equations but not on the recovery ones leaves a
    # residual of exactly 1 − 1/α.
    @testset "uniform lapse: every eigenvalue scales by exactly α" begin
        α = 1.7
        sα, eα = _flat_box(; α=α)
        for k in (0.35, 0.7)
            for (sv, amp) in _shear_modes(h, k)
                @test _residual(sα, eα, amp, k, α*sv) < 1e-9
            end
            for (sv, amp) in _sound_modes(h, k)
                @test _residual(sα, eα, amp, k, α*sv) < 1e-9
            end
        end
    end
end
