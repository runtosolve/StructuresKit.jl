function DefineMesh(MemberDefinitions)

   for i in eachindex(MemberDefinitions)

      L = MemberDefinitions[i][1]
      dL = MemberDefinitions[i][2]
      NumSegments = floor(Int64, L/dL)

      if i == 1
         dz = ones(NumSegments)*dL   #member discretization
         dm = ones(Int8,NumSegments+1)*i    #member properties at each node
      else
         dz = [dz; ones(NumSegments)*dL]
         dm = [dm; ones(Int8, NumSegments)*i]
      end

   end

   dz = dz   #need this to get dz out of function, not sure why
   dm = dm

   return dz, dm

end


function BuildPropertyVector(MemberDefinitions, dm, dz, Property, PropertyOrder, PropertyType)


   z = [0; cumsum(dz)]

   NumberOfNodes = length(dm)

   A = zeros(NumberOfNodes)

   for i=1:NumberOfNodes

      PropertyIndex = MemberDefinitions[dm[i]][PropertyOrder]  #maps properties to each member along beam
      A[i] = Property[PropertyIndex][PropertyType]
   end

   return A

end

function CalculateExpectedStrengths(ASDorLRFD, dz, dm, MemberDefinitions, SectionProperties, CrossSectionDimensions, MaterialProperties, LoadLocation, BracingProperties, RoofSlope)

    NumberOfNodes=length(dz)+1

    #calculate strength limit state capacities

    #strong axis flexure, local-global interaction
    Fy=BuildPropertyVector(MemberDefinitions, dm, dz, MaterialProperties, 4, 3)
    Ixx=BuildPropertyVector(MemberDefinitions, dm, dz, SectionProperties, 3, 1)
    ho=BuildPropertyVector(MemberDefinitions, dm, dz, CrossSectionDimensions, 7, 2)
    ycy=ho./2  #distance from neutral axis to outer fiber
    Sxx=Ixx./ycy
    Myxx=Fy.*Sxx
    Mnexx=Myxx
    Mcrℓxx=BuildPropertyVector(MemberDefinitions, dm, dz, SectionProperties, 3, 6)
    eMnℓxx=S10016F321.(Mnexx, Mcrℓxx, ASDorLRFD)

    #weak axis flexure, local-global interaction
    Iyy=BuildPropertyVector(MemberDefinitions, dm, dz, SectionProperties, 3, 2)

    t = BuildPropertyVector(MemberDefinitions, dm, dz, CrossSectionDimensions, 7, 1)
    b = BuildPropertyVector(MemberDefinitions, dm, dz, CrossSectionDimensions, 3, 3)
    d = BuildPropertyVector(MemberDefinitions, dm, dz, CrossSectionDimensions, 3, 4)
    θc = BuildPropertyVector(MemberDefinitions, dm, dz, CrossSectionDimensions, 3, 5)

    #distance from neutral axis to outer fiber
    ycx = b.+d.*cos.(deg2rad.(θc)) .-t./2
    Syy = Iyy./ycx
    Myyy=Fy.*Syy
    Mneyy=Myyy
    Mcrℓyy=BuildPropertyVector(MemberDefinitions, dm, dz, SectionProperties, 3, 7)
    eMnℓyy=S10016F321.(Mneyy, Mcrℓyy, ASDorLRFD)

    #torsion
    Cw=BuildPropertyVector(MemberDefinitions, dm, dz, SectionProperties, 3, 5)
    Fy=BuildPropertyVector(MemberDefinitions, dm, dz, MaterialProperties, 4, 3)
    Wn=BuildPropertyVector(MemberDefinitions, dm, dz, SectionProperties, 3, 6)
    eBn = S10024H411.(Cw, Fy, Wn, ASDorLRFD)

    #distortional buckling
    CorZ = BuildPropertyVector(MemberDefinitions, dm, dz, CrossSectionDimensions, 3, 6)
    CorZ=trunc.(Int, CorZ)
    E = BuildPropertyVector(MemberDefinitions, dm, dz, MaterialProperties, 4, 1)
    μ = BuildPropertyVector(MemberDefinitions, dm, dz, MaterialProperties, 4, 2)
    G= E./(2 .*(1 .+ μ))
    f1=1.0*ones(length(NumberOfNodes))
    f2=-1.0*ones(length(NumberOfNodes))
    Lm = BuildPropertyVector(MemberDefinitions, dm, dz, BracingProperties, 6, 3)

    #define moment gradient factor
    #assume it is 1
    M1=1.0*ones(length(NumberOfNodes))
    M2=1.0*ones(length(NumberOfNodes))
    CurvatureSign = -1.0*ones(length(NumberOfNodes))

    #assume clip stiffness does not help distortional buckling for a standing seam roof
    #clips are spaced too far apart
    kϕ=0.0*zeros(length(NumberOfNodes))

    Sf=Sxx

    Mcrd=S1001623331.(CorZ, t, ho, b, d, θc, E, μ, G, f1, f2, M1, M2, CurvatureSign, Lm, kϕ, Sf)

    My=Fy.*Sf

    eMnd = S10016F411.(My, Mcrd, ASDorLRFD)


    #Web shear strength
    a = BuildPropertyVector(MemberDefinitions, dm, dz, BracingProperties, 6, 4)
    h = BuildPropertyVector(MemberDefinitions, dm, dz, CrossSectionDimensions, 7, 7)
    kv  = S10016G233.(a, h)
    Fcr = S10016G232.(E, μ, kv, h, t)
    Vcr = S10016G231.(h, t, Fcr)
    eVn = S10016G21.(E, h, t, Fy, Vcr, ASDorLRFD)


    #note nominal capacities are divided by ASD safety factor or
    #multiplied by LRFD resistance factor here
    return eMnℓxx, eMnℓyy, eBn, eMnd, eVn

end


function BendingTorsionDemandToCapacity(Mxx, Myy, B, eMnℓxx, eMnℓyy, eBn)

    #check bending + torsion interaction
    #ActionM1, ActionM2, ActionB, Interaction
    InteractionCheck=S10024H42.(Mxx,Myy,B,eMnℓxx,eMnℓyy,eBn)

    ActionMxx = [x[1] for x in InteractionCheck]
    ActionMyy = [x[2] for x in InteractionCheck]
    ActionB = [x[3] for x in InteractionCheck]
    TotalInteraction = [x[4] for x in InteractionCheck]

    DemandToCapacity = TotalInteraction./1.15

    return ActionMxx, ActionMyy, ActionB, TotalInteraction, DemandToCapacity

end

function DistortionalDemandToCapacity(Mxx,eMnd)

    #check distortional buckling
    DemandToCapacity=abs.(Mxx./eMnd)

end

function BendingShearDemandToCapacity(Vyy, Mxx, eMnℓxx, eVn)

    Interaction = S10016H21.(Mxx, Vyy, eMnℓxx, eVn)

    DemandToCapacity = Interaction

    return DemandToCapacity

end

function BiaxialBendingDemandToCapacity(Mxx, Myy, eMnℓxx, eMnℓyy)

    #no axial force for now
    Pbar=zeros(length(Mxx))
    Pa=ones(length(Mxx))

    InteractionCheck = S10016H121.(Pbar, Mxx, Myy, Pa, eMnℓxx, eMnℓyy)

    ActionP = [x[1] for x in InteractionCheck]
    ActionMxx = [x[2] for x in InteractionCheck]
    ActionMyy = [x[3] for x in InteractionCheck]
    TotalInteraction = [x[4] for x in InteractionCheck]

    DemandToCapacity = TotalInteraction

    return ActionP, ActionMxx, ActionMyy, TotalInteraction, DemandToCapacity

end


function FindPurlinLineStrength(ASDorLRFD, GravityOrUplift, MemberDefinitions, SectionProperties, CrossSectionDimensions, MaterialProperties, LoadLocation, BracingProperties, RoofSlope, EndBoundaryConditions, Supports)

    dz, dm = DefineMesh(MemberDefinitions)

    DemandToCapacityRange=zeros(4)

    #pick a small load
    if GravityOrUplift==0
      qStart=0.00001
      LoadAngle=RoofSlope
    elseif GravityOrUplift==1
      qStart=-0.00001
      LoadAngle=0.0   #for uplift load is always perpendicular to purlin flange
    end

    qx = -qStart*sin(deg2rad(LoadAngle))
    qy = qStart*cos(deg2rad(LoadAngle))

    UniformLoad=(qx,qy)

    #calculate expected capacities
    eMnℓxx, eMnℓyy, eBn, eMnd,  eVn  = CalculateExpectedStrengths(ASDorLRFD, dz, dm, MemberDefinitions, SectionProperties, CrossSectionDimensions, MaterialProperties, LoadLocation, BracingProperties, RoofSlope)


    #calculate DtoC at this small load
    u, v, ϕ, z, Mxx, Myy, Vxx, Vyy, M1, M2, T, B = SolvePlautBeam(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad)
    BTActionM1, BTActionM2, BTActionB, BTTotalInteraction, BTDemandToCapacity = BendingTorsionDemandToCapacity(Mxx, Myy, B, eMnℓxx, eMnℓyy, eBn)
    DistDemandToCapacity = DistortionalDemandToCapacity(Mxx,eMnd)
    MVDemandToCapacity = BendingShearDemandToCapacity(Vyy, Mxx, eMnℓxx, eVn)
    BBActionP, BBActionM1, BBActionM2, BBTotalInteraction, BBDemandToCapacity = BiaxialBendingDemandToCapacity(Mxx, Myy, eMnℓxx, eMnℓyy)
    DemandToCapacityStart=maximum([BTDemandToCapacity; DistDemandToCapacity; MVDemandToCapacity; BBDemandToCapacity])

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


    return ExpectedPurlinLineStrength, z, eMnℓxx, eMnℓyy, eBn, eMnd, eVn, Mxx, Myy, Vxx, Vyy, M1, M2, T, B, BTActionM1, BTActionM2, BTActionB, BTTotalInteraction, BTDemandToCapacity, DistDemandToCapacity, MVDemandToCapacity, BBActionP, BBActionM1, BBActionM2, BBTotalInteraction, BBDemandToCapacity, DemandToCapacity

end

function CalculatePurlinDemandToCapacity(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad, eMnℓxx, eMnℓyy, eBn, eMnd, eVn)

    u, v, ϕ, z, Mxx, Myy, Vxx, Vyy, M1, M2, T, B = SolvePlautBeam(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad)
    BTActionM1, BTActionM2, BTActionB, BTTotalInteraction, BTDemandToCapacity = BendingTorsionDemandToCapacity(Mxx, Myy, B, eMnℓxx, eMnℓyy, eBn)
    DistDemandToCapacity = DistortionalDemandToCapacity(Mxx,eMnd)
    MVDemandToCapacity = BendingShearDemandToCapacity(Vyy, Mxx, eMnℓxx, eVn)
    BBActionP, BBActionM1, BBActionM2, BBTotalInteraction, BBDemandToCapacity = BiaxialBendingDemandToCapacity(Mxx, Myy, eMnℓxx, eMnℓyy)
    PurlinDemandToCapacity=maximum([BTDemandToCapacity; DistDemandToCapacity; MVDemandToCapacity; BBDemandToCapacity])

    return PurlinDemandToCapacity

end


using AISIS10016
using PlautBeam


using MetalBuildingRoofDesigner


#RunToSolve solution for the AISI 4 span example.


#inputs

#Ix Iy Ixy J Cw Wn Mcrlxx Mcrlyy
SectionProperties = [(9.18,1.28,-2.47,0.00159,15.1, 151.86, 66.33),
                     (7.76,1.08,-2.08,0.000954,12.7, 91.48, 40.06),
                     (9.18+7.76,1.28+1.08,-(2.47+2.08),0.00159+0.000954, 15.1+12.7, 151.86+91.48, 66.33+40.06),
                     (7.76*2,1.08*2,-2.08*2,0.000954*2,12.7*2, 2*91.48, 2*40.06)]

#t, ho, b, d, θ, CorZ, h     h is flat portion of web, for shear
CrossSectionDimensions = [(0.070, 8.0, 2.25, 0.930, 50, 1, (7.7450-0.1850)),
                          (0.059, 8.0, 2.25, 0.910, 50, 1, (7.7615-0.1795)),
                          (0.070+0.059, 8.0, 2.25, 0.910, 50, 1, (7.7450-0.1850)),
                          (0.059*2, 8.0, 2.25, 0.910, 50, 1, (7.7615-0.1795))]



#E  ν Fy
MaterialProperties = [(29500,0.30, 55)]

#ax ay
LoadLocation = [((2.250-0.070/2)/2,4),
                ((2.250-0.059/2)/2,4)]

#change to BracingProperties
#kx kϕ  Lm a      a is shear stiffener spacing, make equal to span since there are none provided
BracingProperties = [(0.100,0.100, 25.0*12, 25.0*12)]


#roof slope
RoofSlope = 4.76   #degrees


PurlinSpacing=5*12

#member information
#L dL SectionProperties MaterialProperties LoadLocation BracingProperties CrossSectionDimensions
 MemberDefinitions = [(1*12,1.0,      1,1,1,1,1),
                      (0.5*12,6.0,    1,1,1,1,1),
                      (22.0*12,12.0,  1,1,1,1,1),
                      (2.5*12,6.0,    3,1,1,1,3),
                      (3.5*12,6.0,    3,1,1,1,3),
                      (20.0*12,12.0,  2,1,2,1,2),
                      (0.5*12,6.0,    2,1,2,1,2),
                      (2*12,3.0,      4,1,2,1,4),
                      (0.5*12,6.0,    2,1,2,1,2),
                      (20.0*12,12.0,  2,1,2,1,2),
                      (3.5*12,6.0,    3,1,1,1,3),
                      (2.5*12,6.0,    3,1,1,1,3),
                      (22.0*12,12.0,  1,1,1,1,1),
                      (0.5*12,6.0,    1,1,1,1,1),
                      (1*12,1.0,      1,1,1,1,1)]

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
EndBoundaryConditions = [3 3]

#supports
#location where u=v=ϕ=0
Supports = [0.5*12 1.5*12 25.5*12 26.5*12 50.5*12 51.5*12 75.5*12 76.5*12 100.5*12 101.5*12]


Pressure=10/1000/144  #kips/in.^2
PurlinSpacing=5*12   #in.

q=Pressure*PurlinSpacing

qx = -q*sin(deg2rad(RoofSlope))
qy = q*cos(deg2rad(RoofSlope))
UniformLoad=(qx,qy)


using PlautBeam
u, v, ϕ, z, Mxx, Myy, Vxx, Vyy, M1, M2, T, B = SolvePlautBeam(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad)

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
