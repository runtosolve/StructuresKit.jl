module Visualize

using GeometryBasics
using FileIO
using Plots

export Mesh3D, create_ply_file, show_closed_cross_section, show_multi_branch_cross_section


function Mesh3D(CrossSectionNodes, z, u, v, ϕ)

    #define planar cross-section geometry
    x=view(CrossSectionNodes,:,1)
    y=view(CrossSectionNodes,:,2)

    # write deformed cross-section in 3D coordinates, section by section along the z axis
    for i in 1:length(z)
        CrossSectionShape=hcat(x.+u[i]+x.*cos.(ϕ[i]).-y.*sin.(ϕ[i]),y.+v[i]+x.*sin.(ϕ[i]).+y.*cos.(ϕ[i]),ones(length(x))*z[i])
        if i==1
            global coordinates=CrossSectionShape   #start with first cross-section
        else
            coordinates=vcat(coordinates, CrossSectionShape)  #add following cross-sections along z
        end
    end

    # apply triangular mesh to 3D coordinates
    NumCrossSectionNodes=length(x)
    NumZNodes=length(z)

    #uppoly are triangles pointing towards +z, downpoly are triangles pointing towards -z
    uppoly=[[i i+1 i+1+NumCrossSectionNodes].+NumCrossSectionNodes*(j-1) for i=1:NumCrossSectionNodes-1 for j=1:NumZNodes-1]
    downpoly=[[i i+NumCrossSectionNodes i+NumCrossSectionNodes+1].+NumCrossSectionNodes*(j-1) for i=1:NumCrossSectionNodes-1 for j=1:NumZNodes-1]

    connectivity=vcat(uppoly,downpoly)  #combine all the triangles
    connectivity=permutedims(reshape(hcat(connectivity...), (length(connectivity[1]), length(connectivity))))  #convert into format to be read by Makie

    return coordinates, connectivity   #these are ready for Makie scenes, e.g., poly()

end


"""
    create_ply_file(coordinates, connectivity, filename)

Accepts triangular surface mesh `coordinates` and `connectivity` along with a `filename' and saves a .ply file.

https://en.wikipedia.org/wiki/PLY_(file_format)

The coordinates and connectivity are converted to GeometryBasics.jl primitives before saving with FileIO.jl.  

"""



function create_ply_file(coordinates, connectivity, filename)

    points = [Point{3, Float64}(coordinates[i,1], coordinates[i,2], coordinates[i,3]) for i = 1:size(coordinates)[1]]

    facets = [TriangleFace{Cint}(connectivity[i, 1:3]) for i in eachindex(connectivity[:,1])]

    mesh = GeometryBasics.Mesh(points, facets)

    # using FileIO, MeshIO
    FileIO.save(filename, mesh)

end



function show_closed_cross_section(xcoords, ycoords, markershape, xlims, ylims)

    num_nodes = length(xcoords)

    plt = []

    for i = 1:num_nodes

        if i == 1

            plt = plot([xcoords[i], xcoords[i+1]], [ycoords[i], ycoords[i+1]], size = (600, 600), legend = false, linecolor = :black, markercolor = :black, markershape = markershape, seriestype = :line, xlims = xlims, ylims = ylims)

        elseif i<num_nodes
            
            plot!(plt, [xcoords[i], xcoords[i+1]], [ycoords[i], ycoords[i+1]], size = (600, 600), legend = false, linecolor = :black, markercolor = :black, markershape = markershape, seriestype = :line, xlims = xlims, ylims = ylims)
        
        elseif i==num_nodes
            
            plot!(plt, [xcoords[i], xcoords[1]], [ycoords[i], ycoords[1]], size = (600, 600), legend = false, linecolor = :black, markercolor = :black, markershape = markershape, seriestype = :line, xlims = xlims, ylims = ylims)
        

        end

    end

    return plt

end


function show_multi_branch_cross_section(xcoords, ycoords, element_connectivity, markershape, xlims, ylims)

    num_elem = size(element_connectivity)[1]

    plt = []

    for i = 1:num_elem

        node_i = Int(element_connectivity[i, 1])
        node_j = Int(element_connectivity[i, 2])

        if i == 1

            plt = plot([xcoords[node_i], xcoords[node_j]], [ycoords[node_i], ycoords[node_j]], size = (600, 600), legend = false, linecolor = :black, markercolor = :black, markershape = markershape, seriestype = :line, xlims = xlims, ylims = ylims)

        else
            
            plot!(plt, [xcoords[node_i], xcoords[node_j]], [ycoords[node_i], ycoords[node_j]], size = (600, 600), legend = false, linecolor = :black, markercolor = :black, markershape = markershape, seriestype = :line, xlims = xlims, ylims = ylims)
        
        end

    end

    return plt

end


end #module

