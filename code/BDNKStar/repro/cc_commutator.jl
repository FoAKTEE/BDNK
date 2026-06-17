#=
    Pandya 2201.12317 CC_plot — the WENO5 mixed-derivative commutator diagnostic.

    The figure plots  ∫∫ dy dx |∂_x∂_y ξ − ∂_y∂_x ξ|  vs t for three values of the
    WENO smoothing parameter ε_W ∈ {1e-3, 1, 1e15}.  ∂ here is the *nonlinear*
    WENO5 derivative used by the scheme, so the mixed partials need not commute:
      • ε_W → ∞ : the nonlinear weights collapse to the linear (optimal) weights,
        the operator becomes a FIXED linear stencil, and ∂_x,∂_y commute to machine
        precision  → flat ~1e-15 line (matches the dashed curve, RIGOROUS).
      • ε_W finite : weights are solution-dependent, the operators do not commute,
        giving a finite commutator that DECAYS as the field smooths.

    SCOPE (no fabrication): the provided reference code (ref-code/1D_conformal_bdnk)
    is 1D-slab only — there is NO 2D BDNK solver to extract Pandya's exact initial
    data/curves.  So this reproduces the *diagnostic's mechanism* on a representative
    2D smoothing flow (diagonal advection + light diffusion of a structured ξ),
    not an exact curve overlay.  Claim level: PRELIMINARY (mechanism).

    Run (stdlib only):  julia code/BDNKStar/repro/cc_commutator.jl
=#

# ---- WENO5 first-derivative operator (Jiang–Shu), periodic, left-biased -------
# Reconstruct f_{i+1/2} from point values, then D f_i = (f_{i+1/2}-f_{i-1/2})/dx.
@inline _pidx(i, N) = mod(i - 1, N) + 1   # periodic 1-based index

function weno5_half!(fh, f, N, εW)
    @inbounds for i in 1:N
        fm2 = f[_pidx(i-2,N)]; fm1 = f[_pidx(i-1,N)]; f0 = f[i]
        fp1 = f[_pidx(i+1,N)]; fp2 = f[_pidx(i+2,N)]
        # candidate reconstructions of f_{i+1/2}
        p0 = ( 2fm2 - 7fm1 + 11f0)/6
        p1 = ( -fm1 + 5f0  + 2fp1)/6
        p2 = ( 2f0  + 5fp1 -  fp2)/6
        # smoothness indicators
        β0 = 13/12*(fm2-2fm1+f0)^2 + 1/4*(fm2-4fm1+3f0)^2
        β1 = 13/12*(fm1-2f0+fp1)^2 + 1/4*(fm1-fp1)^2
        β2 = 13/12*(f0-2fp1+fp2)^2 + 1/4*(3f0-4fp1+fp2)^2
        a0 = 0.1/(εW+β0)^2; a1 = 0.6/(εW+β1)^2; a2 = 0.3/(εW+β2)^2
        s = a0+a1+a2
        fh[i] = (a0*p0 + a1*p1 + a2*p2)/s
    end
    return fh
end

# WENO derivative along x (cols vary fastest within a row): operate per row
function dwx!(out, F, N, dx, εW, fh, col)
    @inbounds for j in 1:N
        for i in 1:N; col[i] = F[i,j]; end
        weno5_half!(fh, col, N, εW)
        for i in 1:N
            out[i,j] = (fh[i] - fh[_pidx(i-1,N)])/dx
        end
    end
    return out
end
# WENO derivative along y: operate per column
function dwy!(out, F, N, dy, εW, fh, col)
    @inbounds for i in 1:N
        for j in 1:N; col[j] = F[i,j]; end
        weno5_half!(fh, col, N, εW)
        for j in 1:N
            out[i,j] = (fh[j] - fh[_pidx(j-1,N)])/dy
        end
    end
    return out
end

# centered periodic Laplacian (for the smoothing/diffusion term)
function lap!(out, F, N, h)
    @inbounds for j in 1:N, i in 1:N
        out[i,j] = (F[_pidx(i+1,N),j] + F[_pidx(i-1,N),j] +
                    F[i,_pidx(j+1,N)] + F[i,_pidx(j-1,N)] - 4F[i,j])/h^2
    end
    return out
end

# ---- representative 2D smoothing flow: ξ_t = ν ∇²ξ  (clump spreads/smooths) ----
# A stable, monotonically-smoothing surrogate for the BDNK clump's spreading: as
# ξ smooths the high derivatives decay, so the nonlinear-WENO commutator decays;
# the linear-weight limit (ε_W→∞) stays at the round-off floor throughout.
function commutator_series(εW; N=128, ν=1.5e-3, T=4.0, cfl=0.8, nout=120)
    L=1.0; dx=L/N; x=range(0,L-dx;length=N)
    # sharper features (≈5–6 pts wide) so β_k ~ O(1e-3): engages the WENO
    # nonlinearity, making smaller ε_W the more strongly non-commuting (larger).
    ξ = [exp(-((xi-0.35)^2+(xj-0.4)^2)/0.0016) - 0.7exp(-((xi-0.65)^2+(xj-0.62)^2)/0.0025) +
         0.3*sin(4π*xi)*cos(4π*xj) for xi in x, xj in x]
    dt = cfl*dx^2/(4ν)                       # 2D heat stability limit
    nsteps = ceil(Int, T/dt); dt = T/nsteps
    # scratch
    fh=zeros(N); col=zeros(N); lp=similar(ξ)
    dxy=similar(ξ); dyx=similar(ξ); tmp=similar(ξ); k=similar(ξ)
    rhs!(dst,F) = begin
        lap!(lp,F,N,dx)
        @inbounds for I in eachindex(F); dst[I] = ν*lp[I]; end
    end
    commut(F) = begin                              # ∫∫|∂x∂y F − ∂y∂x F|
        dwy!(tmp,F,N,dx,εW,fh,col); dwx!(dxy,tmp,N,dx,εW,fh,col)   # ∂x(∂y F)
        dwx!(tmp,F,N,dx,εW,fh,col); dwy!(dyx,tmp,N,dx,εW,fh,col)   # ∂y(∂x F)
        s=0.0; @inbounds for I in eachindex(F); s += abs(dxy[I]-dyx[I]); end
        s*dx*dx
    end
    ts=Float64[]; cs=Float64[]; every=max(1,nsteps÷nout)
    push!(ts,0.0); push!(cs,commut(ξ))
    for n in 1:nsteps
        rhs!(k,ξ);            @inbounds for I in eachindex(ξ); tmp[I]=ξ[I]+dt*k[I]; end          # SSP-RK3
        rhs!(k,tmp);          @inbounds for I in eachindex(ξ); tmp[I]=0.75ξ[I]+0.25*(tmp[I]+dt*k[I]); end
        rhs!(k,tmp);          @inbounds for I in eachindex(ξ); ξ[I]=(ξ[I]+2*(tmp[I]+dt*k[I]))/3; end
        if n%every==0 || n==nsteps
            push!(ts, n*dt); push!(cs, commut(ξ))
        end
    end
    return ts, cs
end

if abspath(PROGRAM_FILE) == @__FILE__
    out = joinpath(@__DIR__, "cc_commutator.txt")
    open(out,"w") do io
        println(io, "# t  c_eW1e-3  c_eW1  c_eW1e15   (∫∫|∂x∂y ξ−∂y∂x ξ|)")
        t3,c3 = commutator_series(1e-3); _,c1 = commutator_series(1.0); _,c15 = commutator_series(1e15)
        for i in eachindex(t3)
            println(io, t3[i], "  ", c3[i], "  ", c1[i], "  ", c15[i])
        end
        println(stderr, "CC: eW=1e15 flat≈", round(sum(c15)/length(c15), sigdigits=2),
                "  eW=1e-3 0→end ", round(c3[1],sigdigits=2), "→", round(c3[end],sigdigits=2))
    end
    println("SAVED ", out)
end
