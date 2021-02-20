using StructuresKit

#Verify Julia solution against Mathematica solution from
#Example 2 from Plaut, R.H., Moen, C.D.(2020). "Lateral-Torsional Deformations of Single-Span and Two-Span Thin-Walled Beams with Continuous Bracing".
#Journal of Constructional Steel Research.

#Z-section, no slope, single span, fixed-fixed
#kϕ=1500 N*mm/rad/mm, kx=0.1 N/mm^2, gravity load
#q = 1 kN/m
#Check again Figure 14 in the Plaut and Moen (2020) manuscript.
#  max ϕ

plaut_solution = 0.0164

#*********** Julia solution

#inputs

#Ix Iy Ixy J Cw
section_properties = [(3.230E6,449530,-865760, 397.09, 3.4104E9)]

#E  ν
material_properties = [(200,0.30)]

#ax ay
load_location = [(27.826,101.6)]

#kx kϕ
spring_stiffness = [(0.1/1000, 1500/1000)]

#ay_kx
spring_location = [(101.6)]

#member information
#L dL SectionProperties MaterialProperties LoadLocation SpringStiffness spring_location
member_definitions = [(7620,7620/48,1,1,1,1,1)]

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
end_boundary_conditions = [2 2]

#supports
#location where u=v=ϕ=0
supports = [0.0 7620]

load = (0.0, 5.0/1000)  #kN/mm

z, u, v, ϕ, beam_properties = ThinWalledBeam.solve(member_definitions, section_properties, material_properties, spring_stiffness, spring_location, supports, load, load_location, end_boundary_conditions)

ϕmax = maximum(ϕ)

#calculate percent error as compared to Mathematica solution
error= abs.((ϕmax .- plaut_solution))./ plaut_solution

#accept 1% error from numerical solution
@test error <= 0.01
