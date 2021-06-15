using StructuresKit

#Verify Julia solution against Mathematica solution from
#Example 2 from Plaut, R.H., Moen, C.D.(2020). "Lateral-Torsional Deformations of C-Section and Z-Section Beams with Continuous Bracing".  Proceedings of the Structural Stability Research Council Annual Conference, Atlanta, Georgia.

#Z-section, no slope, simple span
#kϕ=300 N*mm/rad/mm, kx=0, gravity load
#This solution was calculated with Mathematica
#  max ϕ    q
plaut_solution=[0.0 0.0
0.0221 0.2
0.0476 0.4
0.0771 0.6
0.1116 0.8
0.1313 0.9
0.1527 1
0.1763 1.1
0.2022 1.2
0.2310 1.3
0.2630 1.4
0.2989 1.5
0.3393 1.6
0.3615 1.65
0.3756 1.68
0.3853 1.7
0.4003 1.73]

#*********** Julia solution

#inputs

#Ix Iy Ixy J Cw
section_properties = [(3.230E6,449530,-865760, 397.09, 3.4104E9)]

#E  ν
material_properties = [(200,0.30)]

#ax ay
load_location = [(27.826,101.6)]

#kx kϕ
spring_stiffness = [(0.0,300/1000)]

#ay_kx
spring_location = [(101.6)]

#member information
#L dL section_properties material_properties load_location spring_stiffness spring_location
member_definitions = [(7620, 7620/12 ,1, 1, 1, 1, 1)]

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
end_boundary_conditions = [1 1]

#supports
#location where u=v=ϕ=0
supports = [0.0 7620]

#initialize twist vector
ϕmax=zeros(length(plaut_solution[:,1]))

#calculate max twist with PlautBeam
for i=1:length(plaut_solution[:,1])

    load = (0.0, plaut_solution[i,2]/1000)  #kN/mm

    z, u, v, ϕ, beam_properties = ThinWalledBeam.solve(member_definitions, section_properties, material_properties, spring_stiffness, spring_location, supports, load, load_location, end_boundary_conditions)

    ϕmax[i] = maximum(ϕ)

end

#calculate percent error at every load compared to Mathematica solution
error= abs.((ϕmax .- plaut_solution[:,1]))./ plaut_solution[:,1]

max_error=maximum(error[2:end])

#accept 1% error from numerical solution
@test max_error <= 0.01




#update with this format!!

using StructuresKit

#Ix Iy Ixy J Cw
section_properties = [(3.230E6,449530,-865760, 397.09, 3.4104E9)]

#E  ν
material_properties = [(200,0.30)]

#kx kϕ
spring_stiffness = [(0.0,300/1000)]

#ay_kx
spring_location = [(101.6)]

#qx qy
loads = [(0.0001, 0.00002)]

#ax ay
load_locations = [(27.826,101.6)]

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
end_boundary_conditions = [1 1]

#supports
#location where u=v=ϕ=0
supports = [0.0 7620]

#member information
#L(1) dL(2) section_properties(3) material_properties(4) spring_stiffness(5) spring_location(6) load(7) load_location(8) 
member_definitions = [(7620, 7620/12, 1, 1, 1, 1, 1, 1)]

#Define model inputs.
model = ThinWalledBeam.user_interface(member_definitions, section_properties, material_properties, spring_stiffness, spring_location, loads, load_locations, end_boundary_conditions, supports)

#Calculate model stiffness and external force.
model = ThinWalledBeam.define(model)

#Solve model.
model = ThinWalledBeam.solve(model)