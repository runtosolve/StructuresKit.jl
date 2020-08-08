using Test

@testset "InternalForces" begin

    include(string(@__DIR__) * "/InternalForces/" * "InternalForcesTest1.jl")
    include(string(@__DIR__) * "/InternalForces/" * "InternalForcesTest2.jl")
    include(string(@__DIR__) * "/InternalForces/" * "InternalForcesTest3.jl")
    include(string(@__DIR__) * "/InternalForces/" * "InternalForcesTest4.jl")

end
