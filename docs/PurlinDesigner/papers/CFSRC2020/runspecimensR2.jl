using StructuresKit
using CSV, DataFrames

#Import Gao and Moen (2013) specimen section properties and Mcrl.

path = "/Users/crismoen/.julia/dev/StructuresKit/docs/PurlinDesigner/papers/CFSRC2020/inputs/"

filename = string(path, "julia_sectionprops.csv")
sectprops = CSV.read(filename, DataFrame, header=false)

filename = string(path, "julia_Mcrl_xx.csv")
Mcrlxxall = CSV.read(filename, DataFrame, header=false)

filename = string(path, "Wn.csv")
Wnall = CSV.read(filename, DataFrame, header=false)

filename = string(path, "/Gao_Moen_2013_Table1", ".csv")
dims = CSV.read(filename, DataFrame, header=false)


i  = 1

ASDorLRFD=2    #ASDorLRFD=0 for ASD, =1 for LRFD, 2 for nominal

MemberDefinitions = [(7468.0,7468.0/100,    1,1,1,1,1)]

PurlinSpacing=(4140.0-2286.0)/2 + 2286.0/2;  #in.

RoofSlope = 0.0;   #degrees

#location where u=v=ϕ=0
Supports = [0.0 7468.0]

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
EndBoundaryConditions = [1 1];

#Each row below defines a set of section properties.  There are 4 rows because there are 4 cross-sections in this example -  8ZS2.25x070, 8ZS2.25x059, 8ZS2.25x070+ 8ZS2.25x059 at a splice, and 8ZS2.25x059+8ZS2.25x059 at a splice.

#Ixy is negative here because the AISI D100 coordinate system has y pointing up however the PlautBeam.jl coordinate system has y pointing down.

#Ix Iy Ixy J Cw Wn Mcrℓx Mcrℓy Bcrℓ
Ixx = sectprops[i,1]
Iyy = sectprops[i,2]
Ixy = -sectprops[i,3]
J = sectprops[i,4]
Cw = sectprops[i,5]
Wn = Wnall[i,1]
Mcrlxx = Mcrlxxall[i,1]
Mcrlyy = 99999999999.0
Bcrl = 99999999999.0

SectionProperties = [(Ixx, Iyy, Ixy, J, Cw, Wn, Mcrlxx, Mcrlyy, Bcrl)];


#t is base metal thickness
#ho, b, d, and h are outside purlin depth, flange width, lip length, and web flat height
#θ is lip angle from the horizontal
#CorZ=0 for C, CorZ=1 for Z
#This nomenclature is consistent with AISI S100-16.

#t, ho, b, d, θ, CorZ, h
b = dims[i,1]
d = dims[i,2]
θ = dims[i,4]
ho = dims[i,9]
r = dims[i,10]
t = dims[i,11]
h = ho - 2*(r+t)
CorZ = dims[i, 13]


CrossSectionDimensions =
[(t, ho, b, d, θ , CorZ, h)];


#ax ay


#location of centroid
Xc = sectprops[i,6]   #from centerline web
Yc = sectprops[i,7]   #from bottom fiber

#distance from centroid to shear center
Xs = sectprops[i,8]
Ys = sectprops[i,9]

#location of shear center
Xs_loc = Xc+Xs
Ys_loc = Yc+Ys

c = dims[i, 14]

#for a C
if CorZ == 0
        ax = abs(Xs_loc) - t/2 + b - c
elseif CorZ == 1  #for a Z
        ax = abs(Xs_loc) - t/2 + c
end

ay = ho - abs(Ys_loc)  #top of top flange

LoadLocation = [(ax, ay)];


#E  ν  Fy
Fy = dims[i, 12]
MaterialProperties = [(203395,0.30, Fy)];



#Bracing stiffness
#Lm is the spacing between bracing that restrains distortional buckling
#a is the web shear stiffener spacing, assumed equal to the span length here since none are provided
#kx  kϕ  Lm  a

kϕ = dims[i,15]*1.0
kx = dims[i,16]*1.0


BracingProperties = [(kx,kϕ, 7468.0, 7468.0)];



GravityOrUplift=1   #GravityOrUplift=0 for gravity loading

using Plots
plot(z, ϕ)



eqn, z, beamProperties, deformation, strengths, forces, interactions, demand_to_capacity = PurlinDesigner.lineStrength(ASDorLRFD, GravityOrUplift, MemberDefinitions, SectionProperties, CrossSectionDimensions, MaterialProperties, LoadLocation, BracingProperties, RoofSlope, EndBoundaryConditions, Supports)

FailurePressure=eqn/PurlinSpacing*10^6   #N/m^2


plot(z, forces.Mxx)

using Plots
plot(node[:, 2], node[:, 3])


#A Ix Iy J Cw xo yo

FlangeProperties, Springs, Loads = PurlinDesigner.free_flange_define(MemberDefinitions, MaterialProperties, CrossSectionDimensions, forces.Mxx)


u, v, ϕ, properties = BeamColumn.solve(MemberDefinitions, SectionProperties, MaterialProperties, Loads, Springs, EndBoundaryConditions, Supports)
