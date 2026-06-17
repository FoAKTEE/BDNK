using Test

@testset "BDNKStar — STEP 0 (EOS + primitive recovery + causality)" begin
    include("test_eos.jl")
    include("test_recovery.jl")
    include("test_causality.jl")
    include("test_tov.jl")
    include("test_conformal.jl")
    include("test_conformal_evolution.jl")
    include("test_conformal_convergence.jl")
    include("test_radial.jl")
    include("test_heat_criterion.jl")
    include("test_bjorken.jl")
    include("test_kovtun.jl")
end
