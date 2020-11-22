
using StructuresKit



shape_name = "W14X90"

shape_info = StructuresKit.CrossSection.AISC(shape_name)

#A Ix Iy J Cw xc yc xs ys
section_properties = [(shape_info.A, shape_info.Ix, shape_info.Iy, shape_info.J, shape_info.Cw, 0.0, shape_info.d/2, 0.0, shape_info.d/2)]

#E  ν
material_properties = [(29000,0.30)]

#kx ky kϕ hx hy
springs = [0.0, 0.0, 0.0, 0.0, 0.0]

#member information
#L dL SectionProperties MaterialProperties Springs
member_definitions = [(13.0*12, 1.0,1,1,1)]

num_nodes = Int(member_definitions[1][1]/member_definitions[1][2]) + 1

#P 
loads = 10000.0 * ones(num_nodes)

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
end_boundary_conditions = [1 1]

#supports
#location where u=v=ϕ=0
supports = [0.0 13.0*12]

#imperfection
L = member_definitions[1][1]
dL = member_definitions[1][2]
z = 0: dL :L
Δo = L / 1000 * sin.(π*z/L) 

uo = Δo  #in the weak axis
vo = zeros(Float64, num_nodes)
ϕo = zeros(Float64, num_nodes)

imperfections = [uo, vo, ϕo]

#Calculate weak axis flexural buckling load
E = material_properties[1][1]
Iy = shape_info.Iy
Pcry = π^2*E*Iy/L^2   #weak axis buckling load

#calculate yield load
Py = 50*shape_info.A

#define load range for a VR test
P = 0:0.99*Pcry/10:0.99*Pcry

for i = 1:length(P)

    loads = P[i] * ones(num_nodes)

    #solve for column deformation
    u, v, ϕ, properties = Column.solve(member_definitions, section_properties, material_properties, loads, springs, end_boundary_conditions, supports, imperfections)

    #add initial imperfection to solution
    u_total = uo .+ u
    v_total = vo .+ v
    ϕ_total = ϕo .+ ϕ


    #create surface mesh of undeformed column

    #define discretization along member
    zcoords = 0:dL:L
    n=(4, 4, 4, 4, 4)

    xcoords, ycoords = StructuresKit.CrossSection.wshape_nodes(shape_info, n)

    #undeformed shape
    coordinates, connectivity = StructuresKit.Mesh.surface(xcoords, ycoords, zcoords)

    #deformed shape

    σ_total_normal = Column.normal_stresses(xcoords, ycoords, zcoords, u, v, properties)
    ϵ_total_normal = Column.normal_stains(properties, σ_total_normal)
    w = Column.warping_displacements(zcoords, ϵ_total_normal)

    scale_u = 1.0
    scale_v = 1.0
    scale_ϕ = 1.0
    scale_w = 1.0

    x_deformed, y_deformed, z_deformed = Column.deformed_shape(xcoords, ycoords, zcoords, u, v, ϕ, w, scale_u, scale_v, scale_ϕ, scale_w)

    filename = string(@__DIR__, "/", "W14x90", "_SS_", string(i-1), ".ply")

    def_coordinates = [x_deformed y_deformed z_deformed]

    Visualize.create_ply_file(def_coordinates, connectivity, filename)

end


using Makie
scene = poly(def_coordinates, connectivity, color=:lightgray, shading=:true, show_axis=:true, overdraw=:false, strokecolor = (:black, 0.6), strokewidth = 1)




 








