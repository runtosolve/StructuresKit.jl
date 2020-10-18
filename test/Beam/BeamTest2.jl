using StructuresKit

#Verify Julia solution against Mathematica solution from
#Example 2 from Plaut, R.H., Moen, C.D.(2020). "Lateral-Torsional Deformations of Single-Span and Two-Span Thin-Walled Beams with Continuous Bracing".
#Journal of Constructional Steel Research.

#Z-section, no slope, single span, fixed-fixed
#kϕ=1500 N*mm/rad/mm, kx=0.1 N/mm^2, gravity load
#q = 1 kN/m
#Check again Figure 14 in the Plaut and Moen (2020) manuscript.
#  max ϕ
PlautSolution=0.0164

#*********** PlautBeam solution

#inputs

#Ix Iy Ixy J Cw
SectionProperties = [(3.230E6,449530,-865760, 397.09, 3.4104E9)]


#E  ν
MaterialProperties = [(200,0.30)]

#ax ay
LoadLocation = [(27.826,101.6)]


#kx kϕ
SpringStiffness = [(0.1/1000, 1500/1000)]

#roof slope
RoofSlope = [0.0]   #degrees

#member information
#L dL SectionProperties MaterialProperties LoadLocation SpringStiffness RoofSlope
MemberDefinitions = [(7620,7620/48,1,1,1,1,1)]

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
EndBoundaryConditions = [2 2]

#supports
#location where u=v=ϕ=0
Supports = [0.0 7620]

UniformLoad = (0.0, 5.0/1000)  #kN/mm

z, u, v, ϕ, BeamProperties = Beam.solve(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, SpringStiffness, EndBoundaryConditions, Supports, UniformLoad)

ϕmax=maximum(ϕ)

#calculate percent error as compared to Mathematica solution
Error= abs.((ϕmax .- PlautSolution))./ PlautSolution

#accept 1% error from numerical solution
@test Error <= 0.01
