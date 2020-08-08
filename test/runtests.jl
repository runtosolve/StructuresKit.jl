using SafeTestsets

@safetestset "PlautBeam" begin include("PlautBeamTestset.jl") end
@safetestset "InternalForces" begin include("InternalForcesTestset.jl") end
@safetestset "AISIS10016" begin include("AISIS10016Testset.jl") end
