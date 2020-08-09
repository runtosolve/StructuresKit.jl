using StructuresKit

#Test shear and moment calculations for a three span continuous beam.


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
memberDefinitions = [(30*12,3.0,1,1,1,1)]

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
endBoundaryConditions = [1 1]

#supports
#location where u=v=ϕ=0
supports = [0.0 30.0*12/3 30.0*12*(2/3) 30.0*12]

#load  qx   qy
uniformLoad = (0.0, 5.0/1000)

z, u, v, ϕ, beamProperties = PlautBeam.solve(memberDefinitions, sectionProperties, materialProperties, loadLocation, springStiffness, endBoundaryConditions, supports, uniformLoad)

Mxx = InternalForces.moment(z, -v, beamProperties.E, beamProperties.Ix)

Vyy = InternalForces.shear(z, -v, beamProperties.E, beamProperties.Ix)

#check interior moment over a support

intIndex = findall(x-> x==30.0*12/3, z)

intMxx = Mxx[intIndex[1]]
intVyy = Vyy[intIndex[1]-1]


#http://faculty-legacy.arch.tamu.edu/anichols/index_files/courses/arch331/NS8-2beamdiagrams.pdf
#Figure 36
intMxx_AISC = -0.100*uniformLoad[2]*z[intIndex[1]].^2
intVyy_AISC = -0.600*uniformLoad[2]*z[intIndex[1]] + uniformLoad[2]*3.0  #shear is checked 3 in. away from support


MxxError= abs((intMxx - intMxx_AISC)/ intMxx_AISC)
VyyError= abs((intVyy - intVyy_AISC)/ intVyy_AISC)


#accept 1% error from numerical solution
@test MxxError <= 0.01
@test VyyError <= 0.01
