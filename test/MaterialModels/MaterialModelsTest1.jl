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

n_Wshape=(4, 2, 4, 4, 4)

xcoords, ycoords = CrossSection.wshape_nodes(shape_info, n_Wshape)

plot(xcoords, ycoords, markershape = :o)


# using TriangleMesh

mesh_size = 0.01
mesh = CrossSection.triangular_mesh(xcoords, ycoords, mesh_size)

Ai, cxi, cyi = CrossSection.triangulation_properties(mesh)

cx, cy = CrossSection.centroid_from_cells(Ai, cxi, cyi)


# function moment_of_inertia_from_cells(Ai, ci, c)

#     I = sum(((c .- ci) .^2 .* Ai))

#     return I

# end

Ix = moment_of_inertia_from_cells(Ai, cyi, cy)
Iy = moment_of_inertia_from_cells(Ai, cxi, cx)




# #discretize cross-section with a triangular mesh
# function triangular_mesh(xcoords, ycoords, mesh_size)

#     num_nodes = length(xcoords)
#     num_segments = num_nodes

#     # n_point, n_point_marker, n_point_attribute, n_segment, n_holes
#     poly = TriangleMesh.Polygon_pslg(num_nodes, 1, 0, num_segments, 0)

#     node = [xcoords ycoords]
#     set_polygon_point!(poly, node)

#     node_marker = ones(Int, num_nodes, 1)
#     set_polygon_point_marker!(poly, node_marker)

#     segments = zeros(Int, num_segments, 2)
#     for i=1:num_segments

#         if i == num_segments
#             segments[i, 1:2] = [i, 1]
#         else
#             segments[i, 1:2] = [i, i+1]
#         end

#     end

#     set_polygon_segment!(poly, segments)

#     segment_markers = ones(Int, num_segments)
#     set_polygon_segment_marker!(poly, segment_markers)

#     #switches from https://www.cs.cmu.edu/~quake/triangle.html
#     switches = "penvVa" * string(mesh_size) * "D"

#     mesh = create_mesh(poly, switches)

#     return mesh

# end


#calculate triangle areas


# # https://keisan.casio.com/has10/SpecExec.cgi?path=05000000.Mathematics%252F01000500.Plane%2520geometry%252F10010300.Area%2520of%2520a%2520triangle%2520with%2520three%2520points%252Fdefault.xml&charset=utf-8
# function triangle_area(x1, y1, x2, y2, x3, y3)

#     A = abs((x1*y2 + x2*y3 + x3*y1 - y1*x2 - y2*x3 - y3*x1)/2)

# end

# #https://www.mathopenref.com/coordcentroid.html
# function triangle_centroid(x1, y1, x2, y2, x3, y3)

#     cx = (x1 + x2 + x3)/3
#     cy = (y1 + y2 + y3)/3

#     return cx, cy

# end


# function triangulation_properties(mesh)

#     #calculate cell area and centroid
#     Ai = zeros(Float64, mesh.n_cell)
#     cxi = zeros(Float64, mesh.n_cell)
#     cyi = zeros(Float64, mesh.n_cell)

#     for i = 1:mesh.n_cell

#         p1 = mesh.cell[1, i]
#         p2 = mesh.cell[2, i]
#         p3 = mesh.cell[3, i]

#         x1 = mesh.point[1, p1]
#         y1 = mesh.point[2,p1]
#         x2 = mesh.point[1, p2]
#         y2 = mesh.point[2,p2]
#         x3 = mesh.point[1, p3]
#         y3 = mesh.point[2,p3]

#     Ai[i] = triangle_area(x1, y1, x2, y2, x3, y3)

#     cxi[i], cyi[i] = triangle_centroid(x1, y1, x2, y2, x3, y3)

#     end

#     return Ai, cxi, cyi

# end


# #calculate cross-sectional area from cells

# # function area_from_cells(Ai)

# #     A = sum(Ai)

# # end

# # #calculate cross-section centroid from cells

# # function centroid_from_cells(Ai, cxi, cyi)

# #     A = area_from_cells(Ai)

# #     cx = sum(cxi .* Ai) / A
# #     cy = sum(cyi .* Ai) / A

# #     return cx, cy

# # end


# cx, cy = CrossSection.centroid_from_cells(Ai, cxi, cyi)


# # function moment_of_inertia_from_cells(Ai, ci, c)

# #     I = sum(((c .- ci) .^2 .* Ai))

# #     return I

# # end

# Ix = moment_of_inertia_from_cells(Ai, cyi, cy)
# Iy = moment_of_inertia_from_cells(Ai, cxi, cx)





# #calculate triangle centroids




# # for i = 1:mesh.n_cell

# #     p1 = mesh.cell[1, i]
# #     p2 = mesh.cell[2, i]
# #     p3 = mesh.cell[3, i]

# #     x1 = mesh.point[1, p1]
# #     y1 = mesh.point[2,p1]
# #     x2 = mesh.point[1, p2]
# #     y2 = mesh.point[2,p2]
# #     x3 = mesh.point[1, p3]
# #     y3 = mesh.point[2,p3]

# #    cx[i] = triangle_area(x1, y1, x2, y2, x3, y3)

# # end

# zcoords = zeros(Float64, mesh.n_point)
# coordinates = [mesh.point[1,:] mesh.point[2,:] zcoords]

# connectivity = mesh.cell'

# using Makie
# scene = Makie.poly(coordinates, connectivity, color=:lightgray, shading=:true, show_axis=:true, overdraw=:false, strokecolor = (:black, 0.6), strokewidth = 1)
   



# poly = polygon_Lshape()
# mesh = create_mesh(poly, info_str="my mesh", voronoi=true, delaunay=true, set_area_max=true)


# node = [1.0 0.0 ; 0.0 1.0 ; -1.0 0.0 ; 0.0 -1.0 ;
#         0.25 0.25 ; -0.25 0.25 ; -0.25 -0.25 ; 0.25 -0.25] 

#         # size is number_segments x 2
# seg = [1 2 ; 2 3 ; 3 4 ; 4 1 ; 5 6 ; 6 7 ; 7 8 ; 8 5] 


# # all points get marker 1
# node_marker = [ones(Int,4,1) ; 2*ones(Int,4,1)]
# # last segment gets a different marker
# seg_marker = [ones(Int,4) ; 2*ones(Int,4)]


# # size is number_points x number_attr
# node_attr = rand(8,2) 


# # size is number_holes x 2
# hole = [0.5 0.5] 









# polygon_struct_from_points


# shape_mesh, cross_section_edges = Mesh.open_cross_section_tessellation(xcoords, ycoords)

# connectivity = zeros(Int, length(shape_mesh), 3)
# for i=1:length(shape_mesh)
#     connectivity[i, :] = shape_mesh[i]
# end

