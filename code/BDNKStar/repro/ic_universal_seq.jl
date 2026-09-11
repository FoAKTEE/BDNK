using BDNKStar
using BDNKStar.Units: gram_per_cm3_to_km_minus2
using Printf

# energy-density range (geometric, km^-2). Central ENERGY densities (not rest-mass)
# scan in g/cm^3 equivalent of ε for typical NS cores.
gcgs = gram_per_cm3_to_km_minus2

eos_names = [:SLy, :APR4, :H4, :MS1]

# Build M(εc) along the sequence, find M_max, then build the requested sequence
# from M ≈ 1.0 Msun up to near M_max.
function build_sequence(name::Symbol)
    eos = piecewise_polytrope(name)
    # scan central energy density in g/cm^3 (cgs) over a broad NS core range
    εc_cgs = exp.(range(log(3e14), log(4e16); length=400))
    Ms = Float64[]
    εs = Float64[]
    for ec in εc_cgs
        εc = ec * gcgs
        try
            res = moment_of_inertia(eos, εc)
            push!(Ms, res.M_Msun); push!(εs, εc)
        catch
        end
    end
    # locate maximum mass index
    imax = argmax(Ms)
    Mmax = Ms[imax]
    εc_max = εs[imax]
    # build a sequence: stable branch only (εc up to εc_max), masses from ~1.0 to near Mmax
    out = NamedTuple[]
    # target masses spanning 1.0 .. ~0.99*Mmax
    Mtargets = collect(range(1.0, 0.99*Mmax; length=8))
    for Mt in Mtargets
        # find stable-branch εc giving mass closest to Mt (εc ≤ εc_max)
        best = Inf; bi = 0
        for i in eachindex(Ms)
            εs[i] ≤ εc_max || continue
            d = abs(Ms[i] - Mt)
            if d < best; best = d; bi = i; end
        end
        bi == 0 && continue
        res = moment_of_inertia(eos, εs[bi])
        push!(out, (M=res.M_Msun, R=res.R_km, C=res.C, Ibar=res.Ibar))
    end
    # also append the actual maximum-mass config
    resmax = moment_of_inertia(eos, εc_max)
    push!(out, (M=resmax.M_Msun, R=resmax.R_km, C=resmax.C, Ibar=resmax.Ibar))
    return Mmax, out
end

function main()
    println("EOS    M[Msun]   R[km]    C        Ibar_comp   Ibar_Breu   dev%")
    allrows = NamedTuple[]
    maxdev = 0.0
    for name in eos_names
        Mmax, seq = build_sequence(name)
        for s in seq
            pred = IBAR_FROM_C_BREU(s.C)
            dev = 100*(s.Ibar - pred)/pred
            maxdev = max(maxdev, abs(dev))
            push!(allrows, (eos=String(name), M=s.M, R=s.R, C=s.C, Ibar=s.Ibar, pred=pred, dev=dev))
            @printf("%-5s  %6.4f  %7.3f  %6.4f   %8.5f   %8.5f  %+6.3f\n",
                    String(name), s.M, s.R, s.C, s.Ibar, pred, dev)
        end
        @printf("# %-5s M_max = %.4f Msun\n", String(name), Mmax)
    end
    @printf("\nMAX |dev| across all EOS and all masses = %.4f %%\n", maxdev)

    # universality check: spread in Ibar at fixed C across the 4 EOS.
    println("\n# Universality: Ibar spread at matched compactness")
    Cgrid = 0.12:0.02:0.30
    for Cc in Cgrid
        vals = Float64[]
        for name in eos_names
            rows = filter(r -> r.eos == String(name), allrows)
            Cs = [r.C for r in rows]; Is = [r.Ibar for r in rows]
            (minimum(Cs) ≤ Cc ≤ maximum(Cs)) || continue
            p = sortperm(Cs); Cs = Cs[p]; Is = Is[p]
            k = searchsortedfirst(Cs, Cc)
            k = clamp(k, 2, length(Cs))
            t = (Cc - Cs[k-1])/(Cs[k]-Cs[k-1])
            push!(vals, Is[k-1] + t*(Is[k]-Is[k-1]))
        end
        length(vals) < 2 && continue
        spread = 100*(maximum(vals)-minimum(vals))/((maximum(vals)+minimum(vals))/2)
        @printf("C=%.2f  n=%d  Ibar in [%.4f, %.4f]  spread=%.3f%%\n",
                Cc, length(vals), minimum(vals), maximum(vals), spread)
    end
    return allrows, maxdev
end

main()
