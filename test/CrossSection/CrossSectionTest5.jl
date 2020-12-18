
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
supports = [0.0, 13.0*12]

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

#inelasticity
inelasticity_flag = 0


# mutable struct Properties

#     A::Array{Float64,1}
#     Ix::Array{Float64,1}
#     Iy::Array{Float64,1}
#     J::Array{Float64,1}
#     Cw::Array{Float64,1}
#     xc::Array{Float64,1}
#     yc::Array{Float64,1}
#     xs::Array{Float64,1}
#     ys::Array{Float64,1}
#     xo::Array{Float64,1}
#     yo::Array{Float64,1}
#     Io::Array{Float64,1}
#     E::Array{Float64,1}
#     ν::Array{Float64,1}
#     G::Array{Float64,1}
#     kx::Float64
#     ky::Float64
#     kϕ::Float64
#     hx::Float64
#     hy::Float64
#     uo::Array{Float64,1}
#     vo::Array{Float64,1}
#     ϕo::Array{Float64,1}
#     P::Array{Float64,1}
#     dz::Array{Float64,1}
#     z::Array{Float64,1}
#     dm::Array{Float64,1}
#     Azz::Array{Float64,2}
#     Azzzz::Array{Float64,2}
#     supports::Array{Float64,1}
    
# end


# dz, z, dm = Mesh.define_line_element(member_definitions)

#    NumberOfNodes=length(dz)+1


#    #define property vectors
#    A =  Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 1)
#    Ix = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 2)
#    Iy = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 3)
#    J = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 4)
#    Cw = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 5)
#    xc = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 6)
#    yc = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 7)
#    xs = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 8)
#    ys = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 9)

#    xo = -(xc .- xs)
#    yo = yc .- ys

#    Io = Ix .+ Iy .+ A .* (xo.^2 + yo.^2)

#    E = Mesh.create_line_element_property_array(member_definitions, dm, dz, material_properties, 4, 1)
#    ν = Mesh.create_line_element_property_array(member_definitions, dm, dz, material_properties, 4, 2)
#    G = E./(2 .*(1 .+ ν))

#    kx = springs[1]
#    ky = springs[2]
#    kϕ = springs[3]

#    hx = springs[4]
#    hy = springs[5]

#    #define imperfections 
#    uo = imperfections[1]
#    vo = imperfections[2]
#    ϕo = imperfections[3] 

#    #define load
#    P =  loads

#    #calculate derivative operators
#    Azzzz,Azz = Column.calculate_derivative_operators(dz) 

#    #apply left and right end boundary condition stencils to derivative operators
#    NthDerivative = 4
#    Azzzz = Column.apply_end_boundary_conditions(Azzzz, end_boundary_conditions, NthDerivative, dz)
    
#    NthDerivative = 2
#    Azz = Column.apply_end_boundary_conditions(Azz, end_boundary_conditions, NthDerivative, dz)

#    column_properties = Properties(A, Ix, Iy, J, Cw, xc, yc, xs, ys, xo, yo, Io, E, ν, G, kx, ky, kϕ, hx, hy, uo, vo, ϕo, P, dz, z, dm, Azz, Azzzz, supports)











#********************************************************************************************

# dz, z, dm = StructuresKit.Mesh.define_line_element(member_definitions)
# A =  StructuresKit.Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 1)

def_coordinates = []
connectivity = []

for i = 1:length(P)

    loads = P[i] * ones(num_nodes)

    #define column information
    column_properties = Column.initialize(member_definitions, section_properties, material_properties, loads, springs, end_boundary_conditions, supports, imperfections)

    #solve for column deformation
    u, v, ϕ = Column.solve(column_properties, inelasticity_flag)
  

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

    σ_total_normal = Column.normal_stresses(xcoords, ycoords, zcoords, u, v, column_properties)
    ϵ_total_normal = Column.normal_stains(column_properties, σ_total_normal)
    w = Column.warping_displacements(zcoords, ϵ_total_normal)

    scale_u = 1.0
    scale_v = 1.0
    scale_ϕ = 1.0
    scale_w = 1.0

    x_deformed, y_deformed, z_deformed = Column.deformed_shape(xcoords, ycoords, zcoords, u, v, ϕ, w, scale_u, scale_v, scale_ϕ, scale_w)

    filename = string(@__DIR__, "/", "W14x90", "_SS_", string(i-1), ".ply")

    def_coordinates = [x_deformed y_deformed z_deformed]

    # Visualize.create_ply_file(def_coordinates, connectivity, filename)

end


using Makie
scene = poly(def_coordinates, connectivity, color=:lightgray, shading=:true, show_axis=:true, overdraw=:false, strokecolor = (:black, 0.6), strokewidth = 1)




 








