

using TetGen
using TetGen: TetgenIO
using GeometryBasics
using GeometryBasics: Mesh, QuadFace
# Construct a cube out of Quads
points = Point{3, Float64}[
    (0.0, 0.0, 0.0), (2.0, 0.0, 0.0),
    (2.0, 2.0, 0.0), (0.0, 2.0, 0.0),
    (0.0, 0.0, 12.0), (2.0, 0.0, 12.0),
    (2.0, 2.0, 12.0), (0.0, 2.0, 12.0)
]

# using Plots

# plot(points, markershape = :o)


facets = QuadFace{Cint}[
    1:4,
    5:8,
    [1,5,6,2],
    [2,6,7,3],
    [3, 7, 8, 4],
    [4, 8, 5, 1]
]
# markers = Cint[-1, -2, 0, 0, 0, 0]
# attach some additional information to our faces!
mesh = Mesh(points, facets)
# result = tetrahedralize(mesh, "vpq1.414a0.1")
# using FileIO, MeshIO
# # and save it in your favorite format
# FileIO.save("column.ply", mesh)
# # using GLMakie, AbstractPlotting
# # GLMakie.mesh(normal_mesh(result), color=(:blue, 0.1), transparency=true)
# # GLMakie.wireframe!(result)



#try a triangular mesh

points = Point{3, Float64}[
    (0.0, 0.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0)]

# using Plots

# plot(points, markershape = :o)


facets = TriangleFace{Cint}[
    [1, 2, 3]]
# markers = Cint[-1, -2, 0, 0, 0, 0]
# attach some additional information to our faces!
mesh = GeometryBasics.Mesh(points, facets)