using TriangleMesh


# size is number_points x 2
node = [1.0 0.0 ; 0.0 1.0 ; -1.0 0.0 ; 0.0 -1.0 ;
0.25 0.25 ; -0.25 0.25 ; -0.25 -0.25 ; 0.25 -0.25]


# size is number_segments x 2
seg = [1 2 ; 2 3 ; 3 4 ; 4 1 ; 5 6 ; 6 7 ; 7 8 ; 8 5]


# all points get marker 1
node_marker = [ones(Int,4,1) ; 2*ones(Int,4,1)]
# last segment gets a different marker
seg_marker = [ones(Int,4) ; 2*ones(Int,4)]

# size is number_points x number_attr
node_attr = rand(8,2)


# size is number_holes x 2
hole = [0.5 0.5]


poly = Polygon_pslg(8, 1, 2, 8, 1)

set_polygon_point!(poly, node)
set_polygon_point_marker!(poly, node_marker)
set_polygon_point_attribute!(poly, node_attr)
set_polygon_segment!(poly, seg)
set_polygon_segment_marker!(poly, seg_marker)
set_polygon_hole!(poly, hole)

# mesh = create_mesh(poly, info_str="my mesh", voronoi=true, delaunay=true, set_area_max=true)



m = create_mesh(poly, info_str = "Mesh test",
		verbose = false,
		check_triangulation = true,
		voronoi = true,
		delaunay = true,
		output_edges = true,
		output_cell_neighbors = true,
		quality_meshing = true,
		prevent_steiner_points_boundary  = false,
		prevent_steiner_points = false,
		set_max_steiner_points = false,
		set_area_max = false,
		set_angle_min = false,
		add_switches = "")


m.point

# using Makie
#
# mesh(m)
#
# # p = polygon_unitSimplex()
# # 	    m = create_mesh(p, info_str = "Mesh test",
# #                             verbose = false,
#                             check_triangulation = true,
#                             voronoi = true,
#                             delaunay = true,
#                             output_edges = true,
#                             output_cell_neighbors = true,
#                             quality_meshing = true,
#                             prevent_steiner_points_boundary  = false,
#                             prevent_steiner_points = false,
#                             set_max_steiner_points = false,
#                             set_area_max = false,
#                             set_angle_min = false,
#                             add_switches = "")
