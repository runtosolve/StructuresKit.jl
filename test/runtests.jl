using SafeTestsets

@safetestset "PlautBeam" begin include("PlautBeamTestset.jl") end
@safetestset "InternalForces" begin include("InternalForcesTestset.jl") end
