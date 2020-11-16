
using StructuresKit



shape_name = "W14X90"
#define discretization along member
zcoords = 0:1:20
n=(4, 4, 4, 4, 4)

shape_info = StructuresKit.CrossSection.AISC(shape_name)
xcoords, ycoords = StructuresKit.CrossSection.wshape_nodes(shape_info, n)


coordinates, connectivity = StructuresKit.Mesh.surface(xcoords, ycoords, zcoords)

# using Makie
# scene = poly(coordinates, connectivity, color=:lightgray, shading=:true, show_axis=:false, overdraw=:false, strokecolor = (:black, 0.6), strokewidth = 1)


#  using GeometryBasics

 

 filename = string(@__DIR__, "/", "column.ply")

 Visualize.create_ply_file(coordinates, connectivity, filename)






