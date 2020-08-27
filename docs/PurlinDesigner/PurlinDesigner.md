# PurlinDesigner

Determine the expected strength of a single or multi-span purlin line under gravity or uplift loading.  Demand moments, shear, and torsion are calculated considering load eccentricity and roof slope.  Strength limit states are defined by AISI S100-16 and include: flexure+shear, flexure+torsion, biaxial bending, and distortional buckling.  Continuous lateral and rotational springs are available to simulate through fastened or standing seam panel bracing.  

## Nomenclature

![Module nomenclature](/docs/PurlinDesigner/figures/gravity.png)

## Example
Calculate the deflection of a 4 span Z-section purlin line supporting a standing seam roof.  Consistent units of kips and inches are used.  The purlin is loaded at the center of the top flange with a uniform downward gravity load. Continuous bracing from roof sheathing is provided as `kx=0.100 kips/in.\in.` and `kϕ=0.100 kip-in./rad/in.`.   More commentary on this example is available in a [Nextjournal notebook](https://nextjournal.com/runtosolve/metal-building-standing-seam-roof-design-example/) and on YouTube.

```julia
using StructuresKit

ASDorLRFD=0;    #ASDorLRFD=0 for ASD, =1 for LRFD

MemberDefinitions = [(1.0*12,2.0,    1,1,1,1,1),
                      (22.5*12,6.0,   1,1,1,1,1),
                      (6.0*12,6.0,    3,1,1,1,3),
                      (20.5*12,6.0,   2,1,2,1,2),
                      (2*12,3.0,      4,1,2,1,4),
                      (20.5*12,6.0,   2,1,2,1,2),
                      (6.0*12,6.0,    3,1,1,1,3),
                      (22.5*12,6.0,   1,1,1,1,1),
                      (1.0*12,2.0,    1,1,1,1,1)];


PurlinSpacing=5*12;  #in.

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
(9.18+7.76,1.28+1.08,-(2.47+2.08),0.00159+0.000954, (15.1+12.7)*2, 151.86+91.48, 66.33+40.06),
(7.76*2,1.08*2,-2.08*2,0.000954*2,12.7*4, 2*91.48, 2*40.06)];


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


GravityOrUplift=0   #GravityOrUplift=0 for gravity loading

eqn, z, strengths, forces, interactions, dc = PurlinDesigner.lineStrength(ASDorLRFD, GravityOrUplift, MemberDefinitions, SectionProperties, CrossSectionDimensions, MaterialProperties, LoadLocation, BracingProperties, RoofSlope, EndBoundaryConditions, Supports)

FailurePressure=eqn/(PurlinSpacing*cos(deg2rad(RoofSlope)))*1000*144

println("ASD expected gravity roof capacity = ",round(FailurePressure,digits=1), " psf")



```
## Background
The PurlinDesigner module uses PlautBeam to perform second-order analysis of the purlin line and AISIS10016 and AISIS10024 to calculate strength limit states.   A [bisection-based root-finding method](https://en.wikipedia.org/wiki/Bisection_method) determines the expected strength by using the demand-to-capacity envelope along the purlin line.

## Verification and testing log
Comparison of PurlinDesigner predicted strengths to physical experiments is ongoing.  The first study will be completed for the 2020 Cold-Formed Steel Research Consortium Colloquium.
