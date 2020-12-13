using StructuresKit

E = 29000.0
σy = 50.0                   #steel yield stress
σy1 = 50.0                  #steel yield stress
σu = 65.3                   #steel ultimate stress
σf=σy * (1 + 0.10)      #steel fracture stress#
ϵy = σy / E              #steel yield strain
ϵy1 = σy / E * 10    #steel strain at end of yield plateau
ϵu = 0.18             #steel ultimate strain
ϵf = 0.21             #steel fracture strain
n = 11



σ, ϵ  = MaterialModels.steel(σy, σy1, σu, σf, ϵy, ϵy1, ϵu, ϵf, n)

using Plots
plot(σ, ϵ, markershape = :o)


shape_name = "W14X90"

shape_info = CrossSection.AISC(shape_name)

n=(4, 4, 4, 4, 4)

xcoords, ycoords = CrossSection.wshape_nodes(shape_info, n)


shape_mesh, cross_section_edges = Mesh.open_cross_section_tessellation(xcoords, ycoords)

connectivity = zeros(Int, length(shape_mesh), 3)
for i=1:length(shape_mesh)
    connectivity[i, :] = shape_mesh[i]
end

zcoords = zeros(Float64, length(xcoords))
coordinates = [xcoords ycoords zcoords]


using Makie
scene = poly(coordinates, connectivity, color=:lightgray, shading=:true, show_axis=:true, overdraw=:false, strokecolor = (:black, 0.6), strokewidth = 1)
   