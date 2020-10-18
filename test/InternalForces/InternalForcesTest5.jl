using StructuresKit

#Test bimoment calculation for a single span beam.
#Example II-11, AISI D100-17, 9CS2.5x059, 25 ft. span, load is 10 lbs/ft
#AISC assumes vertical location of load is at center of web, so do that here also


#Ix Iy Ixy J Cw
sectionProperties = [(10.2911, 0.69678, 0.0, 0.001022, 11.9)]

#E  ν
materialProperties = [(29500, 0.30)]

#ax ay
loadLocation = [(-1.05, 0.0)]

#kx kϕ
springStiffness = [(0.0, 0.0)]

#member information
#L dL SectionProperties MaterialProperties LoadLocation SpringStiffness
memberDefinitions = [(25*12,3.0,1,1,1,1)]

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
endBoundaryConditions = [1 1]

#supports
#location where u=v=ϕ=0
supports = [0.0 25.0*12]

#load  qx   qy
uniformLoad = (0.0, 10/12/1000)  #10 lbs/ft to kips/in.

z, u, v, ϕ, beamProperties = Beam.solve(memberDefinitions, sectionProperties, materialProperties, loadLocation, springStiffness, endBoundaryConditions, supports, uniformLoad)

B = InternalForces.bimoment(z, beamProperties.dm, ϕ, beamProperties.E, beamProperties.Cw)

Bmax = maximum(B)

B_D100 = 7.44   #from Bob Glauz's ballot example, see /testfiles/InternalForcesTest5/495A_Examples_v2.pdf

Error= abs((Bmax - B_D100)/ B_D100)

#accept 1% error from numerical solution
@test Error <= 0.01
