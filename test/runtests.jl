using SafeTestsets

@safetestset "Beam" begin include("BeamTestset.jl") end
@safetestset "Column" begin include("ColumnTestset.jl") end
@safetestset "InternalForces" begin include("InternalForcesTestset.jl") end
@safetestset "AISIS10016" begin include("AISIS10016Testset.jl") end
# @safetestset "AISIS10024" begin include("AISIS10024Testset.jl") end
