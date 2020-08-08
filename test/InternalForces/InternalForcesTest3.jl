using StructuresKit

#Test shear and moment calculations for a cantilever beam.


#Ix Iy Ixy J Cw
sectionProperties = [(100.0, 100.0, 0.0, 10.0, 10.0)]

#E  ν
materialProperties = [(29500, 0.30)]

#ax ay
loadLocation = [(0.0, 0.0)]

#kx kϕ
springStiffness = [(0.0, 0.0)]

#member information
#L dL SectionProperties MaterialProperties LoadLocation SpringStiffness
memberDefinitions = [(30*12,3,1,1,1,1)]

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
endBoundaryConditions = [2 3]

#supports
#location where u=v=ϕ=0
supports = [0.0]

#load  qx   qy
uniformLoad = (0.0, 5.0/1000)

z, u, v, ϕ, beamProperties = PlautBeam.solve(memberDefinitions, sectionProperties, materialProperties, loadLocation, springStiffness, endBoundaryConditions, supports, uniformLoad)


Mxx = -InternalForces.moment(z, v, beamProperties.E, beamProperties.Ix)
Vyy = InternalForces.shear(z, v, beamProperties.E, beamProperties.Ix)


maxMxx = minimum(Mxx)
maxMxx_AISC = -uniformLoad[2]*z[end]^2/2

maxVyy = maximum(Vyy)
maxVyy_AISC = uniformLoad[2]*z[end]

MxxError= abs((maxMxx - maxMxx_AISC))/ maxMxx_AISC
VyyError= abs((maxVyy - maxVyy_AISC))/ maxVyy_AISC


#accept 1% error from numerical solution
@test MxxError <= 0.01
@test VyyError <= 0.01
