using Test

@testset "AISIS10016" begin

    include(string(@__DIR__) * "/AISIS10016/" * "AISIS10016Test1.jl")
    include(string(@__DIR__) * "/AISIS10016/" * "AISIS10016Test2.jl")

end
