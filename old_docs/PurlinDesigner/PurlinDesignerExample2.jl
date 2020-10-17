using StructuresKit

#RunToSolve solution for the AISI 4 span example.


#inputs

#Ix Iy Ixy J Cw Wn Mcrlxx Mcrlyy
sectionProperties = [(9.18,1.28,-2.47,0.00159,15.1, 151.86, 66.33),
                     (7.76,1.08,-2.08,0.000954,12.7, 91.48, 40.06),
                     (9.18+7.76,1.28+1.08,-(2.47+2.08),0.00159+0.000954, 15.1+12.7, 151.86+91.48, 66.33+40.06),
                     (7.76*2,1.08*2,-2.08*2,0.000954*2,12.7*2, 2*91.48, 2*40.06)]

#t, ho, b, d, θ, CorZ, h     h is flat portion of web, for shear
crossSectionDimensions = [(0.070, 8.0, 2.25, 0.930, 50, 1, (7.7450-0.1850)),
                          (0.059, 8.0, 2.25, 0.910, 50, 1, (7.7615-0.1795)),
                          (0.070+0.059, 8.0, 2.25, 0.910, 50, 1, (7.7450-0.1850)),
                          (0.059*2, 8.0, 2.25, 0.910, 50, 1, (7.7615-0.1795))]



#E  ν Fy
materialProperties = [(29500,0.30, 55)]

#ax ay
loadLocation = [((2.250-0.070/2)/2,4),
                ((2.250-0.059/2)/2,4)]

#change to BracingProperties
#kx kϕ  Lm a      a is shear stiffener spacing, make equal to span since there are none provided
bracingProperties = [(0.000,0.000, 25.0*12, 25.0*12)]


#roof slope
roofSlope = 4.76   #degrees


purlinSpacing=5*12

#member information
#L dL SectionProperties MaterialProperties LoadLocation BracingProperties CrossSectionDimensions
 memberDefinitions = [(25*12/3, 5.0,      1,1,1,1,1),
                      (25*12/3, 5.0,      1,1,1,1,1)]

# memberDefinitions = [(25*12, 5.0,      1,1,1,1,1)]

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
endBoundaryConditions = [1 1]

#supports
#location where u=v=ϕ=0
supports = [0.0 200.00]


pressure = 10/1000/144  #kips/in.^2


q = pressure * purlinSpacing

qx = -q*sin(deg2rad(roofSlope))
qy = q*cos(deg2rad(roofSlope))
uniformLoad=(qx,qy)



z, u, v, ϕ, beamProperties = PlautBeam.solve(memberDefinitions, sectionProperties, materialProperties, loadLocation, bracingProperties, endBoundaryConditions, supports, uniformLoad)


Mxx = InternalForces.moment(z, beamProperties.dm, -v, beamProperties.E, beamProperties.Ix)
Myy = InternalForces.moment(z, beamProperties.dm, -u, beamProperties.E, beamProperties.Iy)
Vyy = InternalForces.shear(z, beamProperties.dm, -v, beamProperties.E, beamProperties.Ix)
Vxx = InternalForces.shear(z, beamProperties.dm, -u, beamProperties.E, beamProperties.Iy)
T = InternalForces.torsion(z, beamProperties.dm, ϕ, beamProperties.E, beamProperties.G, beamProperties.J, beamProperties.Cw)
B = InternalForces.bimoment(z, beamProperties.dm, ϕ, beamProperties.E, beamProperties.Cw)

using Plots

plot(z, v, markershape = :circle)
plot(z, Mxx, markershape = :circle)
plot(z, T, markershape = :circle)


diffdm = abs.(diff(beamProperties.dm))
  jumps = findall(x-> x > 0.0, diffdm)


plot!(z, beamProperties.dm, markershape = :circle)
plot!(z, beamProperties.Ix/4, markershape = :circle)


diffdm = diff(beamProperties.dm)
jumps = findall(x-> x > 0.0, diffdm)

  # #shift first jump
  # dm[jumps[1]] = dm[jumps[1]]+1
  #
  # #calculate jumps again
  # diffdm = diff(dm)
  # jumps = findall(x-> x > 0.0, diffdm)


plot(z, ϕ, markershape = :circle)
plot!(z, beamProperties.Ix/4000, markershape = :circle)

)
plot(z,ϕ, markershape = :circle)
plot(z, Mxx, markershape = :circle)
plot!(z, beamProperties.dm*2000, markershape = :circle)
plot!(z, beamProperties.Ix*1000, markershape = :circle)
plot(z, Myy, markershape = :circle)
plot(z, Vyy, markershape = :circle)
plot!(z, beamProperties.Ix*10, markershape = :circle)



plot(z, Vxx, markershape = :circle)
plot(z, T, markershape = :circle)
plot(z, B, markershape = :circle)
