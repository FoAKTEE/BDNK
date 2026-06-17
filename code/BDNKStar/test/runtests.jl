using Test

@testset "BDNKStar — STEP 0 (EOS + primitive recovery + causality)" begin
    include("test_eos.jl")
    include("test_recovery.jl")
    include("test_causality.jl")
    include("test_tov.jl")
    include("test_conformal.jl")
    include("test_radial.jl")
end
