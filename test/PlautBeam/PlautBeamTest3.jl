using StructuresKit
using Test

#Compare Julia solution with LS-DYNA thin-shell FEA
#Example 2 from Plaut, R.H., Moen, C.D.(2020). "Lateral-Torsional Deformations of C-Section and Z-Section Beams with Continuous Bracing".  Proceedings of the Structural Stability Research Council Annual Conference, Atlanta, Georgia.

#Z-section, no slope, simple span
#ax=0, ay=0
#kϕ=0, kx=0, gravity load

#inputs

#Ix Iy Ixy J Cw
SectionProperties = [(3.230E6,449530,-865760, 397.09, 3.4104E9)]


#E  ν
MaterialProperties = [(200,0.30)]

#ax ay
LoadLocation = [(0.0, 0.0)]


#kx kϕ
SpringStiffness = [(0.0, 0.0)]

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

UniformLoad = (0.0, 0.008/1000)  #kN/mm

z, u, v, ϕ, BeamProperties = PlautBeam.solve(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, SpringStiffness, EndBoundaryConditions, Supports, UniformLoad)

umax=maximum(u)

vmax=maximum(v)

ϕmax=maximum(ϕ)

#results from LS-DYNA
umax_LSDYNA = 2.19 #at Node 16556 in LS-DYNA model
vmax_LSDYNA = 1.14 #at Node 16556 in LS-DYNA model

#(y disp. @ Node 21372 - y disp. @ Node 16556)/(z loc @ Node 21372 - zloc @ Node 11740)
ϕmax_LSDYNA = atan((2.26 - 2.19)/(195.936 - 100.851))

#calculate percent error
uError = abs.((umax .- umax_LSDYNA))./ umax_LSDYNA
vError = abs.((vmax .- vmax_LSDYNA))./ vmax_LSDYNA


#accept 1.5% error from numerical solution
@test uError <= 0.015
@test vError <= 0.015

#PlautBeam predicts ϕ=0, judgement call here that 0.0007 rad from LS-DYNA is small
@test ϕmax ≈ ϕmax_LSDYNA atol=ϕmax_LSDYNA
