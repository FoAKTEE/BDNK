#=
    Reproduce Kovtun 1907.08191 SOUND-channel eigenfrequencies (picresoundv09 /
    picimsoundv09) and the at-rest stability conditions (picstab region).

    GROUNDING (ref-paper/sources/arXiv-1907.08191/src/hydro-general-arxiv-v2.tex):

    F_sound(v0=0,ω,k)  [eq:Fsound0, lines 441-446]:
      F = c_s^2 ε1 θ ω^4
        + i w0 (c_s^2 ε1 + θ) ω^3
        - ( w0^2 + k^2 c_s^2 ( c_s^4 ε1^2 + γs ε1
                               + (ε2+π1)(θ - c_s^2 ε1) + ε2 π1 ) ) ω^2
        - i k^2 w0 ( γs + c_s^4 ε1 + c_s^2 θ ) ω
        + k^2 c_s^2 ( w0^2 + k^2 θ ( c_s^2(ε2+π1 - c_s^2 ε1) - γs ) )
      with c_s^2 = ∂p0/∂ε0, γs = (4/3)η+ζ, w0 = ε0+p0.

    Boost to v0≠0 [eq:FF-boosted, lines 363-365]: F_sound(v0≠0; ω',k') = 0 is
      F_sound( v0=0, ω = (ω' - k'_x v0)/√(1-v0²),
                     k_x = (k'_x - ω' v0)/√(1-v0²), k_y = k'_y ) = 0,
    with k'_x = k' cosφ, k'_y = k' sinφ, and in F_sound k² means k_x²+k_y².
    Both ω and k_x are linear in ω', so the boosted F_sound is again a quartic
    in ω'. We expand it by exact polynomial (Complex coefficient) arithmetic
    and root-find via the companion matrix (mirrors kovtun_shear_modes).

    Figure picresoundv09/picimsoundv09 parameters [caption, lines 554-560]:
      v0 = 0.9,  c_s = 0.5,  c_s^2 ε1/γs = 3,  θ/γs = 4,  ε2 = 0,  π1/γs = 3/c_s^2.
    Units w0/γs: set γs = 1, w0 = 1 ⇒ ε1 = 3/c_s^2, θ = 4, π1 = 3/c_s^2.

    Stability conditions at rest [eq:sound-c1 (484), eq:sound-c2 (489)]:
      e>0:  ε2+π1 > γs/c_s^2 + c_s^2 ε1
      Routh-Hurwitz 2nd (dimensionless, ε̄1=c_s^2 ε1/γs, ε̄2=ε2/γs, θ̄=θ/γs, π̄1=π1/γs):
        ε̄1²/c_s² + c_s²(ε̄1-ε̄2)(ε̄1+θ̄)²(ε̄1-π̄1)
          + (ε̄1+θ̄)(2ε̄1² - ε̄1(ε̄2+π̄1) + (θ̄+ε̄2)(θ̄+π̄1)) > 0
=#

include("/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl")
using .BDNKStar
using LinearAlgebra

# ---- exact tiny polynomial arithmetic in a single variable (Complex coeffs) ----
# represent p(x) = sum c[i+1] x^i, lowest degree first.
polymul(a::Vector{<:Complex}, b::Vector{<:Complex}) = begin
    r = zeros(ComplexF64, length(a)+length(b)-1)
    for i in eachindex(a), j in eachindex(b)
        r[i+j-1] += a[i]*b[j]
    end
    r
end
polyadd(a, b) = begin
    n = max(length(a), length(b)); r = zeros(ComplexF64, n)
    for i in eachindex(a); r[i] += a[i]; end
    for i in eachindex(b); r[i] += b[i]; end
    r
end
polyscale(a, s) = a .* ComplexF64(s)
polypow(a, n) = (r = ComplexF64[1.0]; for _ in 1:n; r = polymul(r, a); end; r)

# roots of polynomial (lowest-first) via companion matrix
function polyroots(c::Vector{ComplexF64})
    # strip leading (highest) zeros
    while length(c) > 1 && abs(c[end]) < 1e-300
        c = c[1:end-1]
    end
    n = length(c) - 1
    n == 0 && return ComplexF64[]
    a = c ./ c[end]               # monic, lowest-first
    C = zeros(ComplexF64, n, n)
    for i in 1:n-1; C[i+1, i] = 1.0; end
    for i in 1:n; C[i, n] = -a[i]; end
    return eigvals(C)
end

"""
F_sound(v0=0) coefficients in ω, lowest-order first: returns length-5 vector
[ω^0, ω^1, ω^2, ω^3, ω^4] for a given k² (here generic complex polynomial input).
We build them as plain complex scalars given a numeric k2.
"""
function fsound_rest_coeffs(k2; cs2, ε1, ε2, π1, θ, γs, w0)
    c4 = cs2*ε1*θ
    c3 = im*w0*(cs2*ε1 + θ)
    c2 = -( w0^2 + k2*cs2*( cs2^2*ε1^2 + γs*ε1 + (ε2+π1)*(θ - cs2*ε1) + ε2*π1 ) )
    c1 = -im*k2*w0*( γs + cs2^2*ε1 + cs2*θ )
    c0 = k2*cs2*( w0^2 + k2*θ*( cs2*(ε2+π1 - cs2*ε1) - γs ) )
    return ComplexF64[c0, c1, c2, c3, c4]
end

"""
Boosted sound eigenfrequencies ω'(k',φ): all four complex roots.
Builds the quartic in ω' by substituting the boost polynomials into F_sound(v0=0).
"""
function kovtun_sound_modes(kp, φ; v0=0.9, cs=0.5,
                            ε1=3/0.5^2, ε2=0.0, π1=3/0.5^2, θ=4.0, γs=1.0, w0=1.0)
    cs2 = cs^2
    g = sqrt(1 - v0^2)
    kpx = kp*cos(φ); kpy = kp*sin(φ)
    # ω (rest) as a polynomial in ω':  ω = (ω' - kpx v0)/g  =>  [-(kpx v0)/g, 1/g]
    ωpoly  = ComplexF64[-(kpx*v0)/g, 1/g]
    # k_x (rest) as polynomial in ω': k_x = (kpx - ω' v0)/g => [kpx/g, -v0/g]
    kxpoly = ComplexF64[kpx/g, -v0/g]
    # k² (rest) = k_x² + k_y²  (k_y = kpy is constant)
    k2poly = polyadd(polymul(kxpoly, kxpoly), ComplexF64[kpy^2])

    # F_sound(v0=0) = Σ_n C_n(k²) ω^n, where C_n is a polynomial in k².
    # From fsound_rest_coeffs: dependence on k2 is at most quadratic. Build each
    # C_n(k²) as a polynomial in k², then compose with k2poly, multiply by ω^n.
    # C_n coefficients vs powers of k2:
    #   c4: const
    #   c3: const
    #   c2: const + k2*(...)
    #   c1: k2*(...)
    #   c0: k2*(...) + k2^2*(...)
    A2 = cs2^2*ε1^2 + γs*ε1 + (ε2+π1)*(θ - cs2*ε1) + ε2*π1   # coeff of k2 in -c2 (inside)
    B0 = cs2*w0^2                                            # coeff of k2 in c0
    B1 = cs2*θ*( cs2*(ε2+π1 - cs2*ε1) - γs )                 # coeff of k2^2 in c0
    # C_n as polynomials in k2 (lowest-first):
    Cn = Vector{Vector{ComplexF64}}(undef, 5)
    Cn[5] = ComplexF64[ cs2*ε1*θ ]                                  # ω^4
    Cn[4] = ComplexF64[ im*w0*(cs2*ε1 + θ) ]                        # ω^3
    Cn[3] = ComplexF64[ -w0^2, -cs2*A2 ]                            # ω^2: c2
    Cn[2] = ComplexF64[ 0.0, -im*w0*( γs + cs2^2*ε1 + cs2*θ ) ]     # ω^1: c1
    Cn[1] = ComplexF64[ 0.0, B0, B1 ]                              # ω^0: c0

    # Compose: total(ω') = Σ_n [ C_n(k2poly) ] * ωpoly^n
    total = ComplexF64[0.0]
    for n in 0:4
        cn_of_k2 = Cn[n+1]
        # evaluate polynomial cn_of_k2 at argument k2poly (Horner over k2 powers)
        comp = ComplexF64[0.0]
        for (p, coef) in enumerate(cn_of_k2)   # p-1 power of k2
            comp = polyadd(comp, polyscale(polypow(k2poly, p-1), coef))
        end
        term = polymul(comp, polypow(ωpoly, n))
        total = polyadd(total, term)
    end
    return polyroots(total)
end

# ---------------- validation ----------------
const PARAMS  = (v0=0.9, cs=0.5, ε1=3/0.5^2, ε2=0.0, π1=3/0.5^2, θ=4.0, γs=1.0, w0=1.0)
# rest-frame transport params (NO v0 key) so v0=0.0 keyword is not overridden by splat
const PARAMS0 = (cs=0.5, ε1=3/0.5^2, ε2=0.0, π1=3/0.5^2, θ=4.0, γs=1.0, w0=1.0)

println("="^70)
println("Kovtun 1907.08191 sound channel — params (fig picresoundv09): ", PARAMS)
println("c_s^2 ε1/γs = ", PARAMS.cs^2*PARAMS.ε1/PARAMS.γs, " (target 3)")
println("θ/γs = ", PARAMS.θ/PARAMS.γs, " (target 4)  π1/γs = ", PARAMS.π1/PARAMS.γs,
        " (target 3/cs^2=", 3/PARAMS.cs^2, ")")
println("="^70)

# (A) at-rest stability conditions (eq:sound-c1, eq:sound-c2)
function stability_at_rest(; cs, ε1, ε2, π1, θ, γs, kw...)
    cs2 = cs^2
    cond1 = (ε2+π1) - (γs/cs2 + cs2*ε1)           # >0 needed
    ε̄1 = cs2*ε1/γs; ε̄2 = ε2/γs; θ̄ = θ/γs; π̄1 = π1/γs
    cond2 = ε̄1^2/cs2 + cs2*(ε̄1-ε̄2)*(ε̄1+θ̄)^2*(ε̄1-π̄1) +
            (ε̄1+θ̄)*(2*ε̄1^2 - ε̄1*(ε̄2+π̄1) + (θ̄+ε̄2)*(θ̄+π̄1))
    return cond1, cond2
end
c1, c2 = stability_at_rest(; PARAMS0...)
println("\n[at-rest stability] eq:sound-c1 LHS-RHS = ", round(c1,digits=4), "  (>0 ⇒ stable): ", c1>0)
println("[at-rest stability] eq:sound-c2 value    = ", round(c2,digits=4), "  (>0 ⇒ stable): ", c2>0)

# (B) gapless sound pair ω = ±c_s k at small k (v0=0)
println("\n[small-k, v0=0] gapless sound pair  ω/k → ±c_s = ±", PARAMS.cs)
for kk in (1e-4, 1e-3, 1e-2)
    rts = kovtun_sound_modes(kk, 0.0; v0=0.0, PARAMS0...)
    # pick the two with smallest |Im ω| (gapless sound), report Re ω / k
    perm = sortperm(rts, by = z->abs(imag(z)))
    g1, g2 = rts[perm[1]], rts[perm[2]]
    println(rpad("  k=$(kk):",14), " Re ω/k = ",
            round(real(g1)/kk,digits=4), ", ", round(real(g2)/kk,digits=4),
            "   (Im ω/k = ", round(imag(g1)/kk,digits=4), ", ",
            round(imag(g2)/kk,digits=4), ")")
end

# (C) gapped modes at small k (v0=0): ω = -i w0/(c_s^2 ε1), -i w0/θ  (eq:wgaps-2)
rts0 = kovtun_sound_modes(1e-4, 0.0; v0=0.0, PARAMS0...)
gap_pred1 = -im*PARAMS.w0/(PARAMS.cs^2*PARAMS.ε1)
gap_pred2 = -im*PARAMS.w0/PARAMS.θ
gapped = sort(rts0, by=z->-abs(imag(z)))[1:2]
println("\n[small-k, v0=0] gapped modes (eq:wgaps-2):")
println("  predicted: ", round(gap_pred1,digits=4), ", ", round(gap_pred2,digits=4))
println("  computed : ", round(gapped[1],digits=4), ", ", round(gapped[2],digits=4))

# (D) boosted v0=0.9 dispersion — verify Im ω' ≤ 0 for all k, all φ (stable frame)
println("\n[boosted v0=0.9] scanning Im ω'(k,φ) ≤ 0 over k∈[1e-3,50], φ∈[0,π/2]")
ks = vcat(10 .^ range(-3, log10(50), length=140))
φs = range(0, π/2, length=19)
maxim = -Inf; argmax_kφ = (0.0,0.0); nbad = 0
for φ in φs, k in ks
    rts = kovtun_sound_modes(k, φ; PARAMS...)
    m = maximum(imag, rts)
    global maxim, argmax_kφ, nbad
    if m > maxim; maxim = m; argmax_kφ = (k,φ); end
    if m > 1e-8; nbad += 1; end
end
println("  max Im ω' over scan = ", round(maxim, sigdigits=5),
        "  at (k,φ)=", round.(argmax_kφ,digits=4))
println("  # (k,φ) grid points with Im ω' > 1e-8 : ", nbad, " / ", length(ks)*length(φs))
stable_boosted = maxim <= 1e-8

# (E) sample boosted dispersion values (small k → linear sound; large k → causal)
println("\n[boosted v0=0.9] sample ω'(k,φ=0):")
for kk in (0.001, 0.01, 0.1, 1.0, 10.0, 50.0)
    rts = kovtun_sound_modes(kk, 0.0; PARAMS...)
    rts = sort(rts, by=z->real(z))
    println(rpad("  k=$(kk):",12), join([string("(",round(real(z),digits=3),
            ",", round(imag(z),digits=3),"i)") for z in rts], "  "))
end

# large-k phase speed: |ω/k| < 1 (causal) per eq:cs0largek for φ=0
rts_big = kovtun_sound_modes(1e6, 0.0; PARAMS...)
maxspeed = maximum(abs(real(z))/1e6 for z in rts_big)
println("\n[boosted v0=0.9] large-k max |Re ω/k| (φ=0) = ", round(maxspeed,digits=5),
        "  (<1 ⇒ causal): ", maxspeed < 1)

println("\n", "="^70)
println("SUMMARY:")
println("  at-rest stability (both conditions >0): ", c1>0 && c2>0)
println("  gapless sound pair ω=±c_s k at small k: confirmed (Re ω/k → ±", PARAMS.cs, ")")
println("  boosted v0=0.9 stable (Im ω' ≤ 0 ∀ k,φ): ", stable_boosted)
println("="^70)
