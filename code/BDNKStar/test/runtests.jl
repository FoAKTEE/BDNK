using Test

@testset "BDNKStar — STEP 0 (EOS + primitive recovery + causality)" begin
    include("test_eos.jl")
    include("test_recovery.jl")
    include("test_causality.jl")
end
