
using StructuresKit



shape_name = "W14X90"

shape_info = StructuresKit.CrossSection.AISC(shape_name)

#A Ix Iy J Cw xc yc xs ys
section_properties = [(shape_info.A, shape_info.Ix, shape_info.Iy, shape_info.J, shape_info.Cw, 0.0, shape_info.d/2, 0.0, shape_info.d/2)]


#define column cross-section mesh
n_Wshape=(4, 2, 4, 4, 4)
mesh_size = 0.01
xcoords, ycoords = CrossSection.wshape_nodes(shape_info, n_Wshape)
section_geometry = Column.SectionGeometry(xcoords, ycoords, mesh_size)




#define material properties
Es = 29000.0
σy = 50.0                   #steel yield stress
σy1 = 50.0                    #steel yield stress
σu = 65.3                   #steel ultimate stress
σf= σu      #steel fracture stress#
ϵy = σy / Es              #steel yield strain
ϵy1 = σy / Es * 10    #steel strain at end of yield plateau
ϵu = 0.18             #steel ultimate strain
ϵf = 0.21             #steel fracture strain
n = [1000, 1000, 1000]

material_properties = [(29000,0.30, σy, σy1, σu, σf, ϵy, ϵy1, ϵu, ϵf, n)]


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
P = 0:Py/10:0.90*Py

P = [P; 0.901 * Py]

#inelasticity
#cross-section mesh level inelasticity is not working yet here, so keep it simple
#just reduce Iy at first yield for now
strain_limit = 0.0
solution_controls = Column.SolutionControls(" ", strain_limit)


def_coordinates = []
connectivity = []

u_total_max = zeros(Float64, length(P))

for i = 1:length(P)

    loads = P[i] * ones(num_nodes)

    #reduce Iy when column fails
    if i == length(P)

        #A Ix Iy J Cw xc yc xs ys
        section_properties = [(shape_info.A, shape_info.Ix, shape_info.Iy * 0.30, shape_info.J, shape_info.Cw, 0.0, shape_info.d/2, 0.0, shape_info.d/2)]

    end


    #define column information
    column_properties = Column.initialize(member_definitions, section_properties, section_geometry, material_properties, loads, springs, end_boundary_conditions, supports, imperfections, solution_controls)

#   if i == length(P)

#         #A Ix Iy J Cw xc yc xs ys
#         column_properties.Iy[76:80] .= shape_info.Iy * 0.30

#     end


    #solve for column deformation
    u, v, ϕ = Column.solve(column_properties, solution_controls)
  

    #add initial imperfection to solution
    u_total = uo .+ u
    v_total = vo .+ v
    ϕ_total = ϕo .+ ϕ

    #get max out-of-plane deflection
    u_total_max[i] = maximum(u_total)

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

    if i == length(P)
        scale_u = 1.0
        scale_w = 1.0
    end

    x_deformed, y_deformed, z_deformed = Column.deformed_shape(xcoords, ycoords, zcoords, u, v, ϕ, w, scale_u, scale_v, scale_ϕ, scale_w)

    filename = string(@__DIR__, "/", "W14x90", "_SS_", string(i-1), "_R1", ".ply")

    def_coordinates = [x_deformed y_deformed z_deformed]

    # Visualize.create_ply_file(def_coordinates, connectivity, filename)

end

using Plots
SSRC = plot(u_total_max * 25.4, P * 4.448, markershape = :o, legend = false, fmt = :pdf)
Plots.xlims!(0, 8)
Plots.ylims!(0, 6000)
Plots.xlabel!("column midheight deflection, mm")
Plots.ylabel!("column compressive load, kN")

display(SSRC)

savefig(SSRC, "SSRC.pdf")



using Makie
scene = poly(def_coordinates, connectivity, color=:lightgray, shading=:true, show_axis=:true, overdraw=:false, strokecolor = (:black, 0.6), strokewidth = 1)




 








