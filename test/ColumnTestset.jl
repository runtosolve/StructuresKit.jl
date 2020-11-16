using Test

@testset "Column" begin

    include(string(@__DIR__) * "/Column/" * "ColumnTest1.jl")
    include(string(@__DIR__) * "/Column/" * "ColumnTest2.jl")

end
