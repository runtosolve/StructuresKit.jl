using StructuresKit
using Test

#Verify Julia solution against Mathematica solution from
#Example 2 from Plaut, R.H., Moen, C.D.(2020). "Lateral-Torsional Deformations of C-Section and Z-Section Beams with Continuous Bracing".  Proceedings of the Structural Stability Research Council Annual Conference, Atlanta, Georgia.

#Z-section, no slope, simple span
#kϕ=300 N*mm/rad/mm, kx=0, gravity load
#This solution was calculated with Mathematica
#  max ϕ    q
PlautSolution=[0.0 0.0
0.0221 0.2
0.0476 0.4
0.0771 0.6
0.1116 0.8
0.1313 0.9
0.1527 1
0.1763 1.1
0.2022 1.2
0.2310 1.3
0.2630 1.4
0.2989 1.5
0.3393 1.6
0.3615 1.65
0.3756 1.68
0.3853 1.7
0.4003 1.73]

#*********** PlautBeam solution

#inputs

#Ix Iy Ixy J Cw
SectionProperties = [(3.230E6,449530,-865760, 397.09, 3.4104E9)]


#E  ν
MaterialProperties = [(200,0.30)]

#ax ay
LoadLocation = [(27.826,101.6)]


#kx kϕ
SpringStiffness = [(0.0,300/1000)]

#roof slope
RoofSlope = [0.0]   #degrees

#member information
#L dL SectionProperties MaterialProperties LoadLocation SpringStiffness RoofSlope
MemberDefinitions = [(7620,7620/12,1,1,1,1,1)]

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
EndBoundaryConditions = [1 1]

#supports
#location where u=v=ϕ=0
Supports = [0.0 7620]

#initialize twist vector
ϕmax=zeros(length(PlautSolution[:,1]))

#calculate max twist with PlautBeam
for i=1:length(PlautSolution[:,1])

    UniformLoad = (0.0, PlautSolution[i,2]/1000)  #kN/mm

    z, u, v, ϕ, BeamProperties = PlautBeam.solve(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, SpringStiffness, EndBoundaryConditions, Supports, UniformLoad)

    ϕmax[i]=maximum(ϕ)

end

#calculate percent error at every load compared to Mathematica solution
Error= abs.((ϕmax .- PlautSolution[:,1]))./ PlautSolution[:,1]

MaxError=maximum(Error[2:end])

#accept 1% error from numerical solution
@test MaxError <= 0.01
