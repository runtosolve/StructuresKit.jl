using StructuresKit

shape_name = "W14X90"
#define discretization along member
zcoords = 0:1:20
n=(4, 4, 4, 4, 4)

shape_info = StructuresKit.CrossSection.AISC(shape_name)

#A Ix Iy J Cw xc yc xs ys
section_properties = [(shape_info.A, shape_info.Ix, shape_info.Iy, shape_info.J, shape_info.Cw, shape_info.d/2, 0.0, shape_info.d/2, 0.0)]

#E  ν
material_properties = [(29000,0.30)]

#kx ky kϕ hx hy
springs = [0.0, 0.0, 0.0, 0.0, 0.0]

#member information
#L dL SectionProperties MaterialProperties Springs
member_definitions = [(13.0*12, (13*12)/12,1,1,1)]

num_nodes = Int(member_definitions[1][1]/member_definitions[1][2]) + 1

#P qx qy ax ay
loads = [(1000.0 * ones(num_nodes)),(0.0 * ones(num_nodes)), (0.0 * ones(num_nodes)),(0.0 * ones(num_nodes)),(0.0 * ones(num_nodes))]

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
end_boundary_conditions = [1 1]

#supports
#location where u=v=ϕ=0
supports = [0.0 13.0*12]

#imperfections
L = member_definitions[1][1]
dL = member_definitions[1][2]
z = 0: dL :L
Δo = L / 1000 * sin.(π*z/L) 
uo = Δo
vo = zeros(Float64, num_nodes)
ϕo = zeros(Float64, num_nodes)

imperfections = [uo, vo, ϕo]

u, v, ϕ, properties = BeamColumn.solve(member_definitions, section_properties, material_properties, loads, springs, end_boundary_conditions, supports)

using Plots
plot(properties.z, u)
plot!(properties.z, Δo) 