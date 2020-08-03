using StructuresKit
using Test


#Test torsion calculation for a simply-supported beam with
#twist fixed warping free ends.
#ax=1.0 in., offset uniform load

#Ix Iy Ixy J Cw
sectionProperties = [(100.0, 100.0, 0.0, 10.0, 10.0)]

#E  ν
materialProperties = [(29500, 0.30)]

#ax ay
loadLocation = [(1.0, 0.0)]

#kx kϕ
springStiffness = [(0.0, 0.0)]

#member information
#L dL SectionProperties MaterialProperties LoadLocation SpringStiffness
memberDefinitions = [(30*12,3,1,1,1,1)]

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
endBoundaryConditions = [1 1]

#supports
#location where u=v=ϕ=0
supports = [0.0 30*12]

#load  qx   qy
uniformLoad = (0.0, 5.0/1000)

z, u, v, ϕ, BeamProperties = PlautBeam.solve(memberDefinitions, sectionProperties, materialProperties, loadLocation, springStiffness, endBoundaryConditions, supports, uniformLoad)

T = InternalForces.torsion(ϕ, z, BeamProperties.E, BeamProperties.G, BeamProperties.J, BeamProperties.Cw)

using Plots

plot(z, T)

maxT = maximum(T)
maxT_theory = (uniformLoad[2] *BeamProperties.ax[1] *z[end])/2

Error= abs((maxT - maxT_theory))/ maxT_theory


#accept 1% error from numerical solution
@test Error <= 0.01
