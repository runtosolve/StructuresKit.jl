using Test

@testset "Beam" begin

    include(string(@__DIR__) * "/Beam/" * "BeamTest1.jl")
    include(string(@__DIR__) * "/Beam/" * "BeamTest2.jl")
    include(string(@__DIR__) * "/Beam/" * "BeamTest3.jl")

end
