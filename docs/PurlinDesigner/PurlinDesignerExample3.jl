using StructuresKit

#RunToSolve solution for the AISI 4 span example.


ASDorLRFD=0;    #ASDorLRFD=0 for ASD, =1 for LRFD


    MemberDefinitions = [(1*12,1.0,      1,1,1,1,1),
                         (1.5*12,1.0,    1,1,1,1,1),
                         (21.0*12,6.0,  1,1,1,1,1),
                         (2.5*12,3.0,    3,1,1,1,3),
                         (3.5*12,3.0,    3,1,1,1,3),
                         (20.0*12,12.0,  2,1,2,1,2),
                         (0.5*12,1.0,    2,1,2,1,2),
                         (2*12,1.0,      4,1,2,1,4),
                         (0.5*12,1.0,    2,1,2,1,2),
                         (20.0*12,12.0,  2,1,2,1,2),
                         (3.5*12,3.0,    3,1,1,1,3),
                         (2.5*12,3.0,    3,1,1,1,3),
                         (21.0*12,6.0,  1,1,1,1,1),
                         (1.5*12,1.0,    1,1,1,1,1),
                         (1*12,1.0,      1,1,1,1,1)];


PurlinSpacing=5*12;  #in.


#θ
RoofSlope = 4.76;   #degrees


#location where u=v=ϕ=0
Supports = [1.0*12 26.0*12 51.0*12 76.0*12 101.0*12]

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
EndBoundaryConditions = [3 3];


#Each row below defines a set of section properties.  There are 4 rows because there are 4 cross-sections in this example -  8ZS2.25x070, 8ZS2.25x059, 8ZS2.25x070+ 8ZS2.25x059 at a splice, and 8ZS2.25x059+8ZS2.25x059 at a splice.

#Ixy is negative here because the AISI D100 coordinate system has y pointing up however the PlautBeam.jl coordinate system has y pointing down.

#Ix Iy Ixy J Cw Wn Mcrℓx Mcrℓy
SectionProperties = [
  (9.18,1.28,-2.47,0.00159,15.1, 151.86, 66.33),
  (7.76,1.08,-2.08,0.000954,12.7, 91.48, 40.06),
  (9.18+7.76,1.28+1.08,-(2.47+2.08),0.00159+0.000954, 15.1+12.7, 151.86+91.48, 66.33+40.06),
  (7.76*2,1.08*2,-2.08*2,0.000954*2,12.7*2, 2*91.48, 2*40.06)];


  #t is base metal thickness
  #ho, b, d, and h are outside purlin depth, flange width, lip length, and web flat height
  #θ is lip angle from the horizontal
  #CorZ=0 for C, CorZ=1 for Z
  #This nomenclature is consistent with AISI S100-16.

  #t, ho, b, d, θ, CorZ, h
  CrossSectionDimensions =
  [(0.070, 8.0, 2.25, 0.930, 50, 1, 7.560),
   (0.059, 8.0, 2.25, 0.910, 50, 1, 7.582),
   (0.070+0.059, 8.0, 2.25, 0.910, 50, 1, 7.56),
   (0.059*2, 8.0, 2.25, 0.910, 50, 1, 7.582)];


   #ax ay
LoadLocation = [((2.250-0.070/2)/2,4), ((2.250-0.059/2)/2,4)];


#E  ν  Fy
MaterialProperties = [(29500,0.30, 55)];


#kx  kϕ  Lm  a
#kx has units of kips/in./in., and kϕ of kips-in./rad/in.
#Lm is the spacing between bracing that restrains distortional buckling
#a is the web shear stiffener spacing, assumed equal to the span length here since none are provided
BracingProperties = [(0.100,0.100, 25.0*12, 25.0*12)];


#Solve for the purlin line deformations
#10 psf pressure
Pressure=10/1000/144  #kips/in.^2

q=Pressure*PurlinSpacing*cos(deg2rad(RoofSlope))  #cosine here for slope
qx = -q*sin(deg2rad(RoofSlope))
qy = q*cos(deg2rad(RoofSlope))
UniformLoad=(qx,qy)


z, u, v, ϕ, beamProperties = PlautBeam.solve(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, BracingProperties, EndBoundaryConditions, Supports, UniformLoad);

using Plots
plot(z/12,u, legend=false, markershape= :circle)
ylabel!("u [in.]")
xlabel!("z [ft.]")

plot(z/12,v, legend=false)
ylabel!("v [in.]")
xlabel!("z [ft.]")

plot(z/12,rad2deg.(ϕ), legend=false)
ylabel!("phi [degrees]")
xlabel!("z [ft.]")

Mxx = InternalForces.moment(z, beamProperties.dm, -v, beamProperties.E, beamProperties.Ix)
Myy = InternalForces.moment(z, beamProperties.dm, -u, beamProperties.E, beamProperties.Iy)
Vyy = InternalForces.shear(z, beamProperties.dm, -v, beamProperties.E, beamProperties.Ix)
Vxx = InternalForces.shear(z, beamProperties.dm, -u, beamProperties.E, beamProperties.Iy)
T = InternalForces.torsion(z, beamProperties.dm, ϕ, beamProperties.E, beamProperties.G, beamProperties.J, beamProperties.Cw)
B = InternalForces.bimoment(z, beamProperties.dm, ϕ, beamProperties.E, beamProperties.Cw)


plot(z/12,Mxx, legend=false)
ylabel!("Mxx [kip-in.]")
xlabel!("z [ft.]")

plot(z/12,Vyy, legend=false)
ylabel!("Vyy [kips]")
xlabel!("z [ft.]")

plot(z/12,Myy, legend=false)
ylabel!("Myy [kip-in.]")
xlabel!("z [ft.]")

plot(z/12,Vxx, legend=false)
ylabel!("Vxx [kips]")
xlabel!("z [ft.]")

plot(z/12,T, legend=false)
ylabel!("T [kip-in.]")
xlabel!("z [ft.]")

plot(z/12,B, legend=false)
ylabel!("B [kip-in.^2]")
xlabel!("z [ft.]")


GravityOrUplift=0   #GravityOrUplift=0 for gravity loading

expectedPurlinLineStrength, z, eMnℓxx, eMnℓyy, eBn, eMnd, eVn, Mxx, Myy, Vyy, T, B, BTActionMxx, BTActionMyy, BTActionB, BTTotalInteraction, BTDemandToCapacity, distDemandToCapacity, MVDemandToCapacity, BBActionP, BBActionMxx, BBActionMyy, BBTotalInteraction, BBDemandToCapacity, demandToCapacity = PurlinDesigner.lineStrength(ASDorLRFD, GravityOrUplift, MemberDefinitions, SectionProperties, CrossSectionDimensions, MaterialProperties, LoadLocation, BracingProperties, RoofSlope, EndBoundaryConditions, Supports)

FailurePressure=expectedPurlinLineStrength/(PurlinSpacing*cos(deg2rad(RoofSlope)))*1000*144

println("ASD expected gravity roof capacity = ",round(FailurePressure,digits=1), " psf")

plot(z/12,distDemandToCapacity, legend=false)
plot!([0, 102],[1.0, 1.0], linecolor=:red, linewidth=2)
ylabel!("Distortional buckling demand/capacity ratio")
xlabel!("z [ft.]")
ylims!((0,1.41))
xlims!((0,102))


plot(z/12,BTTotalInteraction, label="Interaction")
plot!(z/12, BTActionMxx, label="ActionMxx")
plot!(z/12, BTActionMyy, label="ActionMyy")
plot!(z/12, BTActionB, label="ActionB")
plot!([0, 102],[1.15, 1.15], label="Limit", linecolor=:red, linewidth=2)
ylabel!("Flexure+torsion interaction")
xlabel!("z [ft.]")
ylims!((0,1.41))
xlims!((0,102))

plot(z/12,BBTotalInteraction, label="Interaction")
plot!(z/12, BBActionMxx, label="ActionMxx")
plot!(z/12, BBActionMyy, label="ActionMyy")
plot!([0, 102],[1.0, 1.0], label="Limit", linecolor=:red, linewidth=2)
ylabel!("Biaxial bending")
xlabel!("z [ft.]")
ylims!((0,1.41))
xlims!((0,102))

plot(z/12,MVDemandToCapacity, legend=false)
plot!([0, 102],[1.0, 1.0], linecolor=:red, linewidth=2)
ylabel!("Shear+flexure demand/capacity ratio")
xlabel!("z [ft.]")
ylims!((0,1.41))
xlims!((0,102))


kx=[0.0:0.02:0.1;0.15:0.05:0.3;0.4:0.2:1.0]
kϕ=0.0:0.1:0.3


ExpectedPurlinLineStrength=zeros(length(kx),length(kϕ))
for i in eachindex(kx)
    for j in eachindex(kϕ)

        BracingProperties=[(kx[i],kϕ[j],25.0*12, 25.0*12)]

    ExpectedPurlinLineStrength[i,j], z, eMnℓxx, eMnℓyy, eBn, eMnd, eVn, Mxx, Myy, Vyy, T, B, BTActionMxx, BTActionMyy, BTActionB, BTTotalInteraction, BTDemandToCapacity, DistDemandToCapacity, MVDemandToCapacity, BBActionP, BBActionMxx, BBActionMyy, BBTotalInteraction, BBDemandToCapacity, DemandToCapacity = PurlinDesigner.lineStrength(ASDorLRFD, GravityOrUplift, MemberDefinitions, SectionProperties, CrossSectionDimensions, MaterialProperties, LoadLocation, BracingProperties, RoofSlope, EndBoundaryConditions, Supports)

    end

end

FailurePressure=ExpectedPurlinLineStrength./(PurlinSpacing.*cos.(deg2rad(RoofSlope))) .*1000 .*144;

plot(kx,FailurePressure[:,1],  markershape=:circle, label="kphi=0.0")
plot!(kx,FailurePressure[:,2], markershape=:square, label="kphi=0.1")
plot!(kx,FailurePressure[:,3], markershape=:diamond, label="kphi=0.2")
plot!(kx,FailurePressure[:,4], markershape=:star, label="kphi=0.3")

ylabel!("ASD expected gravity roof capacity [psf]")
xlabel!("kx [kips/in./in.]")
