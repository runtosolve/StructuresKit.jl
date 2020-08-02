using StructuresKit
using Test

@testset "PlautBeam" begin

    include(string(@__DIR__) * "/PlautBeam/" * "PlautBeamTest1.jl")
    include(string(@__DIR__) * "/PlautBeam/" * "PlautBeamTest2.jl")
    include(string(@__DIR__) * "/PlautBeam/" * "PlautBeamTest3.jl")

end
