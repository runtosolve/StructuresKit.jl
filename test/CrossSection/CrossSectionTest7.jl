
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


loads = P[end] * ones(num_nodes)


   #calculate moment from deformation
   #get cross-section cell discretization
   #calculate strains from P+M along the column
   #assign strain to each cell
   #reduce area of cell based on tangent elastic modulus,  Acell * Etangent/Ei  
   #calculate updated I


#define column information
properties = Column.initialize(member_definitions, section_properties, material_properties, loads, springs, end_boundary_conditions, supports, imperfections)

#solve for column deformation
u, v, ϕ = Column.solve(properties, inelasticity_flag)

#define column cross-section mesh

n_Wshape=(4, 2, 4, 4, 4)
xcoords, ycoords = CrossSection.wshape_nodes(shape_info, n_Wshape)

mesh_size = 0.01
mesh = CrossSection.triangular_mesh(xcoords, ycoords, mesh_size)

Ai, cxi, cyi = CrossSection.triangular_mesh_properties(mesh)

cx, cy = CrossSection.centroid_from_cells(Ai, cxi, cyi)

# Ixx = CrossSection.moment_of_inertia_from_cells(Ai, cyi, cy)
# Iyy = CrossSection.moment_of_inertia_from_cells(Ai, cxi, cx)

#calculate moment
Myy = InternalForces.moment(properties.z, properties.dm, -u, properties.E, properties.Iyy)

#calculate strain
ϵ = zeros(Float64, mesh.n_cell)


Es = 29000.0
σy = 50.0                   #steel yield stress
σy1 = 50.0                    #steel yield stress
σu = 65.3                   #steel ultimate stress
σf= σu      #steel fracture stress#
ϵy = σy / Es              #steel yield strain
ϵy1 = σy / Es * 10    #steel strain at end of yield plateau
ϵu = 0.18             #steel ultimate strain
ϵf = 0.21             #steel fracture strain
n = [10, 10, 10]



# dx = ϵy1 - ϵy
# dy = σy1 - σy

# normAB = norm([dx, dy])

# unitnormals_test = [-dy, dx] / normAB



closed_or_open = 1
unit_normals = CrossSection.surface_normals([0.0, ϵy, ϵy1], [0.0, σy, σy1], closed_or_open)
node_normal = CrossSection.avg_node_normals(unit_normals, closed_or_open)

interior_angle = 90


radius = 2.0
xy_PI = [ϵy, σy]
PI_unit_normal = node_normal[2,:]
n_radius = 4
xy_curve, Δ, E, T, xy_o, xy_BC, xy_EC, BC_unit_tangent, EC_unit_tangent, radius_unit_vector_BC = Geometry.circular_curve(radius, interior_angle, xy_PI, PI_unit_normal, n_radius)






ϵ_s, σ_s  = MaterialModels.steel(σy, σy1, σu, σf, ϵy, ϵy1, ϵu, ϵf, n)








PI_unit_normal
    






using Dierckx
spl = Spline1D(ϵ_s, σ_s)

ϵ_s_check = 0:0.001:0.15

 steel = LinearInterpolation(ϵ_s, σ_s)

#number of cross-sections 
num_sections = length(z)

k = 11  #load
# for i = 1:num_sections

    # for j = 1:1:mesh.n_cells

    i=78
    j=1

        ϵ_axial = -P[k]/(properties.A[i] * properties.E[i]) 
        ϵ_flexure = Myy[i]/(properties.E[i]*properties.Iy[i]) * (cxi[j] - cx)

        ϵ_total[j] = ϵ_axial + ϵ_flexure

        #calculate tangent elastic modulus



    # end

# end






 








