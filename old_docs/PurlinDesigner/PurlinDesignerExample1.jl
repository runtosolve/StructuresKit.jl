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
bracingProperties = [(0.100,0.100, 25.0*12, 25.0*12)]


#roof slope
roofSlope = 4.76   #degrees


purlinSpacing=5*12

#member information
#L dL SectionProperties MaterialProperties LoadLocation BracingProperties CrossSectionDimensions
 memberDefinitions = [(1*12,1.0,      1,1,1,1,1),
                      (0.5*12,1.0,    1,1,1,1,1),
                      (22.0*12,3.0,  1,1,1,1,1),
                      (2.5*12,3.0,    3,1,1,1,3),
                      (3.5*12,3.0,    3,1,1,1,3),
                      (20.0*12,3.0,  2,1,1,1,2),
                      (0.5*12,1.0,    2,1,1,1,2),
                      (2*12,3.0,      4,1,1,1,4),
                      (0.5*12,1.0,    2,1,1,1,2),
                      (20.0*12,3.0,  2,1,1,1,2),
                      (3.5*12,3.0,    3,1,1,1,3),
                      (2.5*12,3.0,    3,1,1,1,3),
                      (22.0*12,3.0,  1,1,1,1,1),
                      (0.5*12,1.0,    1,1,1,1,1),
                      (1*12,1.0,      1,1,1,1,1)]





#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
endBoundaryConditions = [3 3]

#supports
#location where u=v=ϕ=0
supports = [1.0*12 26.0*12 51.0*12 76.0*12 101.0*12]


pressure = 10/1000/144  #kips/in.^2


q = pressure * purlinSpacing

qx = -q*sin(deg2rad(roofSlope))
qy = q*cos(deg2rad(roofSlope))
uniformLoad=(qx,qy)



z, u, v, ϕ, beamProperties = PlautBeam.solve(memberDefinitions, sectionProperties, materialProperties, loadLocation, bracingProperties, endBoundaryConditions, supports, uniformLoad)


plot(z, v, markershape = :circle)
plot!(z, beamProperties.dm, markershape = :circle)
Mxx = InternalForces.moment(z, beamProperties.dm, -v, beamProperties.E, beamProperties.Ix)
Myy = InternalForces.moment(z, beamProperties.dm, -u, beamProperties.E, beamProperties.Iy)
Vyy = InternalForces.shear(z, beamProperties.dm, -v, beamProperties.E, beamProperties.Ix)
Vxx = InternalForces.shear(z, beamProperties.dm, -u, beamProperties.E, beamProperties.Iy)
T = InternalForces.torsion(z, beamProperties.dm, ϕ, beamProperties.E, beamProperties.G, beamProperties.J, beamProperties.Cw)
B = InternalForces.bimoment(z, beamProperties.dm, ϕ, beamProperties.E, beamProperties.Cw)

plot(z, B)
plot!(z, beamProperties.dm, markershape = :circle)
xlims!(270, 310)

using Plotsz

plot(z, beamProperties.Ix, markershape = :circle)


using Dierckx

spl=Spline1D(z, beamProperties.Ix, s=0.1)

plot(z, spl(z))
plot!(z, beamProperties.Ix)




plot(z, Mxx)


#define code approach
ASDorLRFD=0
GravityOrUplift=0


kx=[0.0:0.02:0.1;0.15:0.05:0.3;0.4:0.2:1.0]
kϕ=0.0:0.1:0.3


ExpectedPurlinLineStrength=zeros(length(kx),length(kϕ))
for i in eachindex(kx)
    for j in eachindex(kϕ)

        show(i)
        show(j)

        BracingProperties=[(kx[i],kϕ[j],25.0*12, 25.0*12)]

    ExpectedPurlinLineStrength[i,j], z, eMnℓxx, eMnℓyy, eBn, eMnd, eVn, Mxx, Myy, Vxx, Vyy, M1, M2, T, B, BTActionM1, BTActionM2, BTActionB, BTTotalInteraction, BTDemandToCapacity, DistDemandToCapacity, MVDemandToCapacity, BBActionP, BBActionM1, BBActionM2, BBTotalInteraction, BBDemandToCapacity, DemandToCapacity = FindPurlinLineStrength(ASDorLRFD, GravityOrUplift, MemberDefinitions, SectionProperties, CrossSectionDimensions, MaterialProperties, LoadLocation, BracingProperties, RoofSlope, EndBoundaryConditions, Supports)


    end
end

FailurePressure=ExpectedPurlinLineStrength./(PurlinSpacing.*cos.(deg2rad(RoofSlope))) .*1000 .*144;

using Plots
plot(kx, FailurePressure[:,1], markershape=:circle)
plot!(kx, FailurePressure[:,2], markershape=:circle)
plot!(kx, FailurePressure[:,3], markershape=:circle)
plot!(kx, FailurePressure[:,4], markershape=:circle)
#
# i=2
# j=1
# BracingProperties=[(kx[i],kϕ[j],25.0*12, 25.0*12)]
#
# ExpectedPurlinLineStrength, z, eMnℓxx, eMnℓyy, eBn, eMnd, eVn, Mxx, Myy, Vxx, Vyy, M1, M2, T, B, BTActionM1, BTActionM2, BTActionB, BTTotalInteraction, BTDemandToCapacity, DistDemandToCapacity, MVDemandToCapacity, BBActionP, BBActionM1, BBActionM2, BBTotalInteraction, BBDemandToCapacity, DemandToCapacity = FindPurlinLineStrength(ASDorLRFD, GravityOrUplift, MemberDefinitions, SectionProperties, CrossSectionDimensions, MaterialProperties, LoadLocation, BracingProperties, RoofSlope, EndBoundaryConditions, Supports)



dz, dm = DefineMesh(MemberDefinitions)

DemandToCapacityRange=zeros(4)

using AISIS10024
#calculate expected capacities
eMnℓxx, eMnℓyy, eBn, eMnd,  eVn  = CalculateExpectedStrengths(ASDorLRFD, dz, dm, MemberDefinitions, SectionProperties, CrossSectionDimensions, MaterialProperties, LoadLocation, BracingProperties, RoofSlope)



#pick a small load
if GravityOrUplift==0
  q=0.0001  #initial guess
  LoadAngle=RoofSlope
elseif GravityOrUplift==1
  q=-0.0001
  LoadAngle=0.0   #for uplift load is always perpendicular to purlin flange
end

qx = -qStart*sin(deg2rad(LoadAngle))
qy = qStart*cos(deg2rad(LoadAngle))
UniformLoad=(qx,qy)
DemandToCapacity=FindPurlinDemandToCapacity(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn)

#initialize residual
Residual=abs(1.0-DemandToCapacity)
eps=0.0001

qfinal = RootFinder(q, eps, Residual, DemandToCapacity, MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn, LoadAngle)

qx = -qfinal*sin(deg2rad(LoadAngle))
qy = qfinal*cos(deg2rad(LoadAngle))
UniformLoad=(qx,qy)
DemandToCapacity=FindPurlinDemandToCapacity(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn)



function FindExpectedPurlinStrength(q, eps, Residual, DemandToCapacity, MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn, LoadAngle)

    #use bisection method of rootfinding to solve for purlin strength

    while Residual>eps

    newq=q/DemandToCapacity
    q=q + (newq-q)/2

    qx = -q*sin(deg2rad(LoadAngle))
    qy = q*cos(deg2rad(LoadAngle))
    UniformLoad=(qx,qy)

    DemandToCapacity=CalculatePurlinDemandToCapacity(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn)

    Residual=abs(1.0-DemandToCapacity)

    end

    return q

end

q=q+ (newq-q)/2

qx = -q*sin(deg2rad(LoadAngle))
qy = q*cos(deg2rad(LoadAngle))
UniformLoad=(qx,qy)

DemandToCapacity=FindPurlinDemandToCapacity(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn)

newq=q/DemandToCapacity

q=q+ (newq-q)/2

qx = -q*sin(deg2rad(LoadAngle))
qy = q*cos(deg2rad(LoadAngle))
UniformLoad=(qx,qy)

DemandToCapacity=FindPurlinDemandToCapacity(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn)

newq=q/DemandToCapacity

q=q+ (newq-q)/2

qx = -q*sin(deg2rad(LoadAngle))
qy = q*cos(deg2rad(LoadAngle))
UniformLoad=(qx,qy)

DemandToCapacity=FindPurlinDemandToCapacity(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn)

newq=q/DemandToCapacity

q=q+ (newq-q)/2

qx = -q*sin(deg2rad(LoadAngle))
qy = q*cos(deg2rad(LoadAngle))
UniformLoad=(qx,qy)

DemandToCapacity=FindPurlinDemandToCapacity(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn)


#####
#scale load linearly so that DtoC is approximately equal to 1
qFail=qStart/DemandToCapacityStart

#now calculate what the load is at qFail

qx = -qFail*sin(deg2rad(LoadAngle))
qy = qFail*cos(deg2rad(LoadAngle))

UniformLoad=(qx,qy)

u, v, ϕ, z, Mxx, Myy, Vxx, Vyy, M1, M2, T, B = SolvePlautBeam(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad)
BTActionM1, BTActionM2, BTActionB, BTTotalInteraction, BTDemandToCapacity = BendingTorsionDemandToCapacity(Mxx, Myy, B, eMnℓxx, eMnℓyy, eBn)
DistDemandToCapacity = DistortionalDemandToCapacity(Mxx,eMnd)
MVDemandToCapacity = BendingShearDemandToCapacity(Vyy, Mxx, eMnℓxx, eVn)
BBActionP, BBActionM1, BBActionM2, BBTotalInteraction, BBDemandToCapacity = BiaxialBendingDemandToCapacity(Mxx, Myy, eMnℓxx, eMnℓyy)
DemandToCapacityRange[2]=maximum([BTDemandToCapacity; DistDemandToCapacity; MVDemandToCapacity; BBDemandToCapacity])

qRange=[qFail*0.90; qFail; qFail*1.05; qFail*1.10]

qx = -qRange[1]*sin(deg2rad(LoadAngle))
qy = qRange[1]*cos(deg2rad(LoadAngle))

UniformLoad=(qx,qy)

u, v, ϕ, z, Mxx, Myy, Vxx, Vyy, M1, M2, T, B = SolvePlautBeam(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad)
BTActionM1, BTActionM2, BTActionB, BTTotalInteraction, BTDemandToCapacity = BendingTorsionDemandToCapacity(Mxx, Myy, B, eMnℓxx, eMnℓyy, eBn)
DistDemandToCapacity = DistortionalDemandToCapacity(Mxx,eMnd)
MVDemandToCapacity = BendingShearDemandToCapacity(Vyy, Mxx, eMnℓxx, eVn)
BBActionP, BBActionM1, BBActionM2, BBTotalInteraction, BBDemandToCapacity = BiaxialBendingDemandToCapacity(Mxx, Myy, eMnℓxx, eMnℓyy)
DemandToCapacityRange[1]=maximum([BTDemandToCapacity; DistDemandToCapacity; MVDemandToCapacity; BBDemandToCapacity])

qx = -qRange[3]*sin(deg2rad(LoadAngle))
qy = qRange[3]*cos(deg2rad(LoadAngle))

UniformLoad=(qx,qy)

u, v, ϕ, z, Mxx, Myy, Vxx, Vyy, M1, M2, T, B = SolvePlautBeam(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad)
BTActionM1, BTActionM2, BTActionB, BTTotalInteraction, BTDemandToCapacity = BendingTorsionDemandToCapacity(Mxx, Myy, B, eMnℓxx, eMnℓyy, eBn)
DistDemandToCapacity = DistortionalDemandToCapacity(Mxx,eMnd)
MVDemandToCapacity = BendingShearDemandToCapacity(Vyy, Mxx, eMnℓxx, eVn)
BBActionP, BBActionM1, BBActionM2, BBTotalInteraction, BBDemandToCapacity = BiaxialBendingDemandToCapacity(Mxx, Myy, eMnℓxx, eMnℓyy)
DemandToCapacityRange[3]=maximum([BTDemandToCapacity; DistDemandToCapacity; MVDemandToCapacity; BBDemandToCapacity])

qx = -qRange[4]*sin(deg2rad(LoadAngle))
qy = qRange[4]*cos(deg2rad(LoadAngle))

UniformLoad=(qx,qy)

u, v, ϕ, z, Mxx, Myy, Vxx, Vyy, M1, M2, T, B = SolvePlautBeam(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad)
BTActionM1, BTActionM2, BTActionB, BTTotalInteraction, BTDemandToCapacity = BendingTorsionDemandToCapacity(Mxx, Myy, B, eMnℓxx, eMnℓyy, eBn)
DistDemandToCapacity = DistortionalDemandToCapacity(Mxx,eMnd)
MVDemandToCapacity = BendingShearDemandToCapacity(Vyy, Mxx, eMnℓxx, eVn)
BBActionP, BBActionM1, BBActionM2, BBTotalInteraction, BBDemandToCapacity = BiaxialBendingDemandToCapacity(Mxx, Myy, eMnℓxx, eMnℓyy)
DemandToCapacityRange[4]=maximum([BTDemandToCapacity; DistDemandToCapacity; MVDemandToCapacity; BBDemandToCapacity])

show(DemandToCapacityRange)
show(qRange)

DtoCSpline=Spline1D(DemandToCapacityRange, qRange)

ExpectedPurlinLineStrength=DtoCSpline(1.0)

qx = -ExpectedPurlinLineStrength*sin(deg2rad(LoadAngle))
qy = ExpectedPurlinLineStrength*cos(deg2rad(LoadAngle))

UniformLoad=(qx,qy)

u, v, ϕ, z, Mxx, Myy, Vxx, Vyy, M1, M2, T, B = SolvePlautBeam(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad)
BTActionM1, BTActionM2, BTActionB, BTTotalInteraction, BTDemandToCapacity = BendingTorsionDemandToCapacity(Mxx, Myy, B, eMnℓxx, eMnℓyy, eBn)
DistDemandToCapacity = DistortionalDemandToCapacity(Mxx,eMnd)
MVDemandToCapacity = BendingShearDemandToCapacity(Vyy, Mxx, eMnℓxx, eVn)
BBActionP, BBActionM1, BBActionM2, BBTotalInteraction, BBDemandToCapacity = BiaxialBendingDemandToCapacity(Mxx, Myy, eMnℓxx, eMnℓyy)
DemandToCapacity=maximum([BTDemandToCapacity; DistDemandToCapacity; MVDemandToCapacity; BBDemandToCapacity])
