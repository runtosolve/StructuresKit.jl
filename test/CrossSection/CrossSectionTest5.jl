
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
member_definitions = [(13.0*12, (13*12)/12,1,1,1)]

num_nodes = Int(member_definitions[1][1]/member_definitions[1][2]) + 1

#P 
loads = 4000.0 * ones(num_nodes)

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


#solve for column deformation
u, v, ϕ, properties = Column.solve(member_definitions, section_properties, material_properties, loads, springs, end_boundary_conditions, supports, imperfections)



#create surface mesh of undeformed column

#define discretization along member
zcoords = 0:dL:L
n=(4, 4, 4, 4, 4)

xcoords, ycoords = StructuresKit.CrossSection.wshape_nodes(shape_info, n)

coordinates, connectivity = StructuresKit.Mesh.surface(xcoords, ycoords, zcoords)

#now work on deformed shape

function deformed_shape(x, y, z, u, v, ϕ)

    x_deformed = []
    y_deformed = []
    num_sections = length(z)


    for i = 1:num_sections

        if i == 1
            x_deformed = x.+u[i]+x.*cos.(ϕ[i]).-y.*sin.(ϕ[i])
            y_deformed = y.+v[i]+x.*sin.(ϕ[i]).+y.*cos.(ϕ[i])
        else
            x_deformed = [x_deformed; x.+u[i]+x.*cos.(ϕ[i]).-y.*sin.(ϕ[i])]
            y_deformed = [y_deformed; y.+v[i]+x.*sin.(ϕ[i]).+y.*cos.(ϕ[i])]
        end

    end

    return x_deformed, y_deformed

end

Az, Azz, Azzz = InternalForces.calculateDerivativeOperators(properties.z, properties.dm)

#warping displacements

function axial_stress(P, A)

    σ_axial = P / A

end

function flexural_stress(M, y, I)

    σ_flexural = M * y / I

end


num_sections = length(zcoords)
num_cross_section_nodes = length(xcoords)

Mxx = InternalForces.moment(z, properties.dm, -v, properties.E, properties.Ix)
Myy = InternalForces.moment(z, properties.dm, -u, properties.E, properties.Iy)

σ_axial = zeros(Float64, (num_sections, num_cross_section_nodes))
σ_flexural_xx = zeros(Float64, (num_sections, num_cross_section_nodes))
σ_flexural_yy = zeros(Float64, (num_sections, num_cross_section_nodes))
σ_total_normal = zeros(Float64, (num_sections, num_cross_section_nodes))

for i=1:num_sections

    A=properties.A[i]
    P=properties.P[i]
    Ix = properties.Ix[i]
    Iy = properties.Iy[i]

    σ_axial[i, :] = axial_stress(-properties.P[i], properties.A[i]) * ones(Float64, num_cross_section_nodes)
    σ_flexural_xx[i, :] = flexural_stress.(Mxx[i], -(ycoords .- properties.yc[i]), properties.Ix[i]) 
    σ_flexural_yy[i, :] = flexural_stress.(Myy[i], -(xcoords .- properties.xc[i]), properties.Iy[i]) 
    #add warping stresses here one day...
    σ_total_normal[i,:] = σ_axial[i, :] + σ_flexural_xx[i, :] + σ_flexural_yy[i, :]

end

#calculate strains
ϵ_total_normal = σ_total_normal ./ properties.E

# #calculate change in strain
# dϵ_total_normal = zeros(Float64, (num_sections, num_cross_section_nodes))
# strain_trend = zeros(Float64, (num_sections, num_cross_section_nodes))

# for i = 1:num_cross_section_nodes
#     dϵ_total_normal[:, i] = Az * ϵ_total_normal[:,i]
#     strain_trend[:, i] = sign.(dϵ_total_normal[:, i])
# end


using NumericalIntegration

#integrate strains to calculate warping displacements

w = zeros(Float64, (num_sections, num_cross_section_nodes))
w_total = zeros(Float64, (num_sections, num_cross_section_nodes))


for i=1:num_sections

    for j = 1:num_cross_section_nodes

        w_total[i,j] = integrate(zcoords[1:end], ϵ_total_normal[1:end,j])
        
        w[i, j] = -w_total[i,j]/2 + integrate(zcoords[1:i], ϵ_total_normal[1:i,j])  #not sure if this is general or not

    end

end

z_deformed = []
for i=1:num_sections

    if i == 1
        z_deformed = zcoords[i] .+ w[i, :]
    else
        z_deformed = [z_deformed; zcoords[i] .+ w[i, :]]
    end

end




using Plots
plot(xcoords, ycoords, σ_flexural_yy[1,:])


# function normal_warping_stress(E, ϕzz, Wns)

#     σ_normal_warping_stress = E * Wns * ϕzz

# end
# i=5
# integrate(zcoords[i:7], strain_trend[i:7, 4] .* ϵ_total_normal[i:7,4])

#calculate cross-sectional stresses
    #axial
    #flexural
    #warping torsion

#convert stresses to strains
#integrate strains along the column


# #calculate derivatives to calculate warping displacements
# uz = Az*u
# vz = Az*v
# ϕz = Az*ϕ

# num_sections = length(zcoords)

# dwarp = zeros(Float64, (num_sections, length(xcoords)))

# xs = section_properties[1][8]
# ys = section_properties[1][9]

# for i=1:num_sections

#     dwarp[i, :] = uz[i].^2 + vz[i].^2 .+ ((xcoords .-xs).^2 .+(ycoords.-ys).^2).*ϕz[i].^2 .- 2 .*(ycoords .-ys).*uz[i].*ϕz[i] .+ 2 .*(xcoords .-xs).*vz[i].*ϕz[i]
    
# end


# x_deformed, y_deformed = deformed_shape(xcoords, ycoords, zcoords, u, v, ϕ)

# def_coordinates = [x_deformed y_deformed z_deformed] .* 1000

using Makie
scene = poly(def_coordinates, connectivity, color=:lightgray, shading=:true, show_axis=:false, overdraw=:false, strokecolor = (:black, 0.6), strokewidth = 1)
# poly!(scene, coordinates, connectivity, color=:lightgray, shading=:true, show_axis=:false, overdraw=:false, strokecolor = (:black, 0.6), strokewidth = 1)
# scene

# #  using GeometryBasics

 

#  filename = string(@__DIR__, "/", "column.ply")

#  Visualize.create_ply_file(coordinates, connectivity, filename)






