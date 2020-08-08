module StructuresKit

export PlautBeam
include("PlautBeam.jl")
using .PlautBeam

export InternalForces
include("InternalForces.jl")
using .InternalForces

export AISIS10016
include("AISIS10016.jl")
using .AISIS10016

export AISIS10024
include("AISIS10024.jl")
using .AISIS10024

end
