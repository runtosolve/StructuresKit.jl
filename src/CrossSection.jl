module CrossSection

using Statistics
using LinearAlgebra
using CSV
using DataFrames
using TriangleMesh

using ..Geometry

export AISC, wshape_nodes,
       assemble, Feature, Deck, surface_normals, avg_node_normals, xycoords_along_normal, create_CUFSM_node_elem,
       discretize_feature, feature_geometry, get_xy_coordinates, area_from_cells, centroid_from_cells, moment_of_inertia_from_cells,
       triangular_mesh_properties, mesh, SectionProperties, rectangular_tube_geometry


struct WShape

    d::Float64
    tw::Float64
    bf::Float64
    tf::Float64
    kdes::Float64
    k1::Float64
    A::Float64
    Ix::Float64
    Iy::Float64
    J::Float64
    Cw::Float64
    Zx::Float64
    Zy::Float64
    Wno::Float64


end


#primitive, line element defined as vector, with n segments
struct Feature

    ΔL::Array{Float64,1}
    θ::Array{Float64,1}
    n::Array{Int,1}
    radius::Array{Float64,1}
    n_radius::Array{Int,1}
    closed_or_open::Int

end

struct Point

    x::Float64
    y::Float64

end

#deck cross-section definition
#this could be a model for other cross-section types
struct Deck

    features::Tuple{Feature,Feature,Feature}
    feature_map::Array{Int,1}

end


struct Open

    features::Feature
    feature_map::Array{Int,1}

end


#Get the node and element properties just for the flange and lip of a C or Z section.
#This code grabs the bottom flange and lip.
function CZflange_template(CorZ,H,Bc,Bt,Dc,Dt,r1,r2,r3,r4,θc,θt,t,nh,nb1,nb2,nd1,nd2,nr1,nr2,nr3,nr4,kipin,center)

    prop,node,elem,lengths,springs,constraints,geom,cz = CrossSection.CUFSMtemplate(CorZ,H,Bc,Bt,Dc,Dt,r1,r2,r3,r4,θc,θt,t,nh,nb1,nb2,nd1,nd2,nr1,nr2,nr3,nr4,kipin,center)

    if CorZ == 2
        node[:, 2] = -node[:, 2]
    end

    index_xo = findall(x-> x==0.0, node[:, 2])
    index_yo = findall(y-> y==0.0, node[:, 3])
    index_o = intersect(index_xo, index_yo)
    index = 1:(index_o[1]+1)

	nodeflange = node[index,:]
    elemflange = elem[index[1:end-1],:]

    return nodeflange, elemflange

end


"""
    AISC(shape_name)

Accepts `shape_name` as a String and returns cross-section dimensions and section properties in the Struct `shape_info`.

For now only W-shapes have been tested, eg., `AISC("W14x90")` where the struct contains `d, tw, bf, tf, kdes, k1, A, Ix, Iy, J, Cw, Zx, Wno`.

"""


function AISC(shape_name)

    filename = string(@__DIR__, "/assets/aisc-shapes-database-v15.0.csv")

    data = CSV.read(filename, DataFrame, header=true)

    shape_row = findfirst(==(shape_name), data.AISC_Manual_Label)

    section_type = data[shape_row, :Type]

    if section_type == "W"

        d = parse(Float64, data[shape_row, :d])
        tw = parse(Float64, data[shape_row, :tw])
        bf = parse(Float64, data[shape_row, :bf])
        tf = parse(Float64, data[shape_row, :tf])
        kdes = parse(Float64, data[shape_row, :kdes])

        #get k1 from AISC table, it is in fraction format
        k1 = data[shape_row, :k1]
        index = findfirst("/", k1)

        whole = Int(data[shape_row, :k1][1] - '0')

        if isempty(index) == false
            top_fraction = parse(Float64, data[shape_row, :k1][index[1]-2:index[1]-1])
            bottom_fraction = parse(Float64, data[shape_row, :k1][index[1]+1:index[1]+2])
        else
            top_fraction = 0.0
            bottom_fraction = 0.0
        end

        k1 = whole + top_fraction/bottom_fraction

        A = data[shape_row, :A]
        Ix = data[shape_row, :Ix]
        Iy = data[shape_row, :Iy]
        J = parse(Float64, data[shape_row, :J])
        Cw = parse(Float64, data[shape_row, :Cw])
        Zx = data[shape_row, :Zx]
        Zy = data[shape_row, :Zy]
        Wno = parse(Float64, data[shape_row, :Wno])

        shape_info = WShape(d, tw, bf, tf, kdes,k1, A, Ix, Iy, J, Cw, Zx, Zy, Wno)

    end

    return shape_info

end

"""
    wshape_nodes(shape_info, n)

Accepts the Struct `shape_info` generated using CrossSection.AISC and the discretization Vector `n` and outputs the outline x-y coordinates of a W shape 'xcoords' and 'ycoords'.

The Vector 'n' describes the number of segments in a quarter cross-section, i.e., `n = [half of outside flange face, flange thickness, half of inside flange face, flange-web radius, half of web]`.

"""


function wshape_nodes(shape_info, n)

    #from bottom of bottom flange, web centerline to left edge
    xcoords = zeros(n[1]+1)
    ycoords = zeros(n[1]+1)

    flange_range = 0.0 : -shape_info.bf / 2 / n[1] : -shape_info.bf / 2
    [xcoords[i] =  flange_range[i] for i in eachindex(flange_range)]
    ycoords .= 0.0

    #up along bottom flange thickness
    flange_thickness_range = shape_info.tf/n[2]:shape_info.tf/n[2]:shape_info.tf
    xcoords = [xcoords; ones(n[2])*xcoords[end]]
    ycoords = [ycoords; flange_thickness_range]

    #over to fillet radius at bottom flange - web intersection

    # flange_flat = shape_info.bf/2 - shape_info.k1
    flange_flat = shape_info.bf/2 - shape_info.tw/2 - (shape_info.kdes - shape_info.tf)

    inside_flange_range = (xcoords[end] + flange_flat/n[3]) : flange_flat/n[3] : (xcoords[end] + flange_flat)

    xcoords = [xcoords; inside_flange_range]
    ycoords = [ycoords; ones(n[3])*ycoords[end]]

    #go around the fillet
    radius = -xcoords[end] - shape_info.tw/2
    θ = (-π/2 + π/2/n[4]):π/2/n[4]: 0.0

    xo = xcoords[end]
    yo = ycoords[end] + radius

    x_radius = xo .+ radius .* cos.(θ)
    y_radius = yo .+ radius .* sin.(θ)

    # plot(x_radius, y_radius, markershape = :o, linetype = :scatter)

    xcoords = [xcoords; x_radius]
    ycoords = [ycoords; y_radius]

    #add web flat
    web_flat = shape_info.d/2 - shape_info.tf - radius

    web_flat_range = (ycoords[end] + web_flat/n[5]): web_flat/n[5]: (ycoords[end] + web_flat)
    xcoords = [xcoords; ones(n[5])*xcoords[end]]
    ycoords = [ycoords; web_flat_range]

    #mirror about horizontal axis
    ycoords_horz_flip = ycoords .- ycoords[end]
    ycoords_horz_flip = -ycoords_horz_flip
    ycoords_horz_flip = ycoords_horz_flip .+ ycoords[end]

    xcoords = [xcoords; reverse(xcoords)[2:end]]
    ycoords = [ycoords; reverse(ycoords_horz_flip)[2:end]]

    #mirror about vertical axis
    xcoords_vert_flip = reverse(-xcoords)[2:end-1]

    xcoords = [xcoords; xcoords_vert_flip]
    ycoords = [ycoords; reverse(ycoords)[2:end-1]]


    return xcoords, ycoords

end


#define nodal coordinates within a feature
function discretize_feature(feature)

    dx = feature.Δx ./ feature.n
    dy = feature.Δy ./ feature.n

    return dx, dy

end

#define feature geometry, typically a feature is repetitive
function feature_geometry(feature, dx, dy)

    xcoords = []
    ycoords = []

    num_lines = length(feature.Δx)

    for i = 1:num_lines

        if i==1
            xcoords = range(0.0, feature.Δx[i], length = feature.n[i] + 1)
            ycoords = range(0.0, feature.Δy[i], length = feature.n[i] + 1)
        else
            xcoords = [xcoords; xcoords[end] .+ range(dx[i], feature.Δx[i], length = feature.n[i])]
            ycoords = [ycoords; ycoords[end] .+ range(dy[i], feature.Δy[i], length = feature.n[i])]
        end

    end

    return xcoords, ycoords

end


#assemble an open cross-section as a series of line element features
function assemble(OpenSection)

    xcoords = []
    ycoords = []
    num_features = length(OpenSection.feature_map)

    for i =1:num_features

        feature_type = OpenSection.feature_map[i]

        # dx, dy = discretize_feature(OpenSection.features[feature_type])

        if i==1
            xcoords, ycoords = get_xy_coordinates(OpenSection.features[feature_type])
        else
            feature_xcoords, feature_ycoords = get_xy_coordinates(OpenSection.features[feature_type])
            xcoords = [xcoords; xcoords[end] .+ feature_xcoords[2:end]]
            ycoords = [ycoords; ycoords[end] .+ feature_ycoords[2:end]]
        end

    end

    return xcoords, ycoords

end








# calculate surface normals for each line segment in a 2D cross-section
function surface_normals(xcoords, ycoords, closed_or_open)


    if closed_or_open == 0
        numel = length(xcoords)   #closed
    elseif closed_or_open ==1
        numel = length(xcoords) - 1  #open
    end

    unitnormals = zeros(Float64, (numel, 2))

    for i=1:numel

        if (i == numel) & (closed_or_open == 0)   #for tubes
            pointA = Point(xcoords[i], ycoords[i])
            pointB = Point(xcoords[1], ycoords[1])
        else
            pointA = Point(xcoords[i], ycoords[i])
            pointB = Point(xcoords[i + 1], ycoords[i + 1])
        end

        dx = pointB.x - pointA.x
        dy = pointB.y - pointA.y

        normAB = norm([dx, dy])

        unitnormals[i, :] = [-dy, dx] / normAB

        if unitnormals[i,1] == -0.0
            unitnormals[i,1]= 0.0
        end

        if unitnormals[i,2] == -0.0
            unitnormals[i,2]= 0.0
        end

    end

    return unitnormals

end


# calculate average unit normals at each node in a 2D cross-section from element unit normals
function avg_node_normals(unitnormals, closed_or_open)

    if closed_or_open == 0
        numnodes = size(unitnormals)[1]
        numel = size(unitnormals)[1]
    elseif closed_or_open ==1
        numnodes = size(unitnormals)[1]+1
        numel = size(unitnormals)[1]
    end

    nodenormals = zeros(Float64, (numnodes, 2))

    for i=1:numnodes

        if (i == 1) & (closed_or_open == 0)  # where nodes meet in the tube
            nodenormals[i, :] = mean(unitnormals[[numel, 1], :], dims=1)
        elseif (i != 1) & (closed_or_open == 0)  #tube
            nodenormals[i, :] = mean(unitnormals[i-1:i, :], dims=1)
        elseif (i == 1) & (closed_or_open == 1)  #open, first node is element norm
            nodenormals[i, :] = unitnormals[i, :]
        elseif (i != 1) & (i != numnodes) & (closed_or_open == 1)  #open
            nodenormals[i, :] = mean(unitnormals[i-1:i, :], dims=1)
        elseif (i == numnodes) & (closed_or_open == 1)  #open, last node is element norm
            nodenormals[i, :] = unitnormals[i-1, :]

        end

        #make sure unit normal always = 1.0
        unitnorm = norm(nodenormals[i, :])
        if unitnorm < 0.99
            scale = 1.0/unitnorm
            nodenormals[i,:] = scale .* nodenormals[i,:]
        end

    end

    return nodenormals

end


function xycoords_along_normal(xcoords, ycoords, nodenormals, Δ)

    numnodes = size(xcoords)[1]
    xcoords_normal = zeros(Float64, numnodes)
    ycoords_normal = zeros(Float64, numnodes)

    for i=1:numnodes
        xcoords_normal[i] = xcoords[i] + nodenormals[i, 1] * Δ
        ycoords_normal[i] = ycoords[i] + nodenormals[i, 2] * Δ
    end

    return xcoords_normal, ycoords_normal

end

function create_CUFSM_node_elem(xcoords, ycoords, connectivity, t)

    num_elem = size(connectivity)[1]
    num_nodes = length(xcoords)

    node = zeros((num_nodes, 8))
    elem = zeros((num_elem, 5))

    for i=1:num_nodes
        node[i,:] = [i, xcoords[i], ycoords[i], 1, 1, 1, 1, 1.0]
    end

    for i=1:num_elem
        elem[i, :] = [i, connectivity[i,1], connectivity[i,2], t[i], 100]
    end

    return node, elem

end


function get_xy_coordinates(feature)

    xcoords = []
    ycoords = []

    #convert feature vectors to xy components
    Δxy = Geometry.vector_components(feature.ΔL, feature.θ)

    #number of straight line segments in the feature
    num_lines = size(Δxy)[1]

    #get xy line coordinates of feature
    xcoords_straight, ycoords_straight = Geometry.line_coordinates(Δxy, feature.closed_or_open)

    #calculate feature surface normals
    unitnormals = CrossSection.surface_normals(xcoords_straight, ycoords_straight, feature.closed_or_open)

    #calculate feature node normals
    nodenormals = CrossSection.avg_node_normals(unitnormals, feature.closed_or_open)

    #calculate interior corner angle magnitudes
    num_radius = length(feature.radius)
    interior_radius_angle = zeros(Float64, num_radius)
    for i=1:num_radius

        θ1 = feature.θ[i]
        θ2 = feature.θ[i+1]


        if sign(θ1) != sign(θ2)   #don't think this is fully general
            interior_radius_angle[i] = 180 - (abs(θ1) + abs(θ2))
        end

    end

    #calculate distance from curve PI to start and end of curve, along tangents
    Tc = zeros(Float64, length(feature.radius))

    for i = 1:length(feature.radius)

        radius = feature.radius[i]

        if radius > 0.0

            xy_PI = [xcoords_straight[i+1], ycoords_straight[i+1]]
            n = feature.n_radius[i]
            γ = interior_radius_angle[i]
            θ1 = feature.θ[i]
            θ2 = feature.θ[i+1]

            if (sign(θ1) <=0) & (sign(θ2) >= 0)
                PI_unit_normal = -nodenormals[i+1,:]
            else
                PI_unit_normal = nodenormals[i+1,:]
            end

            xy_curve, Δ, E, Tc[i], xy_o, xy_BC, xy_EC, BC_unit_tangent, EC_unit_tangent, radius_unit_vector_BC = Geometry.circular_curve(radius, γ, xy_PI, PI_unit_normal, n)

        elseif radius == 0.0
            Tc[i] = 0.0
        end

    end

    #shorten vector lengths to include curves
    ΔLc = zeros(Float64, num_lines)
    for i = 1:length(feature.ΔL)

        if i == 1  #first
            ΔLc[i] = feature.ΔL[i] - Tc[i]   #j end
        elseif i == length(feature.ΔL)  #last
            ΔLc[i] = feature.ΔL[i] - Tc[end]  #j end
        else  #others
            ΔLc[i] = feature.ΔL[i] - Tc[i-1] - Tc[i] #i and j ends
        end

    end

    Δxyc = Geometry.vector_components(ΔLc, feature.θ)

    #discretize feature
    dx = Δxyc[:,1] ./ feature.n
    dy = Δxyc[:,2] ./ feature.n

    #assemble feature
    for i = 1:num_lines

        if i==1   #first vector
            xcoords = range(0.0, Δxyc[i,1], length = feature.n[i] + 1)
            ycoords = range(0.0, Δxyc[i,2], length = feature.n[i] + 1)
        elseif i != num_lines
            xcoords = [xcoords; xcoords[end] .+ range(dx[i], Δxyc[i,1], length = feature.n[i])]
            ycoords = [ycoords; ycoords[end] .+ range(dy[i], Δxyc[i,2], length = feature.n[i])]
        end

        #add radius to end of vector
        if i < num_lines
            if feature.radius[i]>0.0
                xy_PI = [xcoords_straight[i+1], ycoords_straight[i+1]]
                n = feature.n_radius[i]
                radius = feature.radius[i]
                γ = interior_radius_angle[i]
                θ1 = feature.θ[i]
                θ2 = feature.θ[i+1]

                if (sign(θ1) <=0) & (sign(θ2) >= 0)
                    PI_unit_normal = -nodenormals[i+1,:]  #concave up
                else
                    PI_unit_normal = nodenormals[i+1,:]  #concave down
                end

                xy_curve, Δ, E, T, xy_o, xy_BC, xy_EC, BC_unit_tangent, EC_unit_tangent, radius_unit_vector_BC = Geometry.circular_curve(radius, γ, xy_PI, PI_unit_normal, n)

                if (sign(θ1) <=0) & (sign(θ2) >= 0)
                    xy_curve = reverse(xy_curve, dims=1)    #concave up
                end

                xcoords = [xcoords; xy_curve[2:end,1]]
                ycoords = [ycoords; xy_curve[2:end,2]]

            end
        end

        if i== num_lines  #last vector
            xcoords = [xcoords; xcoords[end] .+ range(dx[i], Δxyc[i,1], length = feature.n[i])]
            ycoords = [ycoords; ycoords[end] .+ range(dy[i], Δxyc[i,2], length = feature.n[i])]
        end

    end

    return xcoords, ycoords

end



#calculate cross-sectional area from cells

function area_from_cells(Ai)

    A = sum(Ai)

end

#calculate cross-section centroid from cells

function centroid_from_cells(Ai, cxi, cyi)

    A = area_from_cells(Ai)

    cx = sum(cxi .* Ai) / A
    cy = sum(cyi .* Ai) / A

    return cx, cy

end


function moment_of_inertia_from_cells(Ai, ci, c)

    I = sum(((c .- ci) .^2 .* Ai))

    return I

end

#discretize cross-section with a triangular mesh
function triangular_mesh(xcoords, ycoords, mesh_size)

    num_nodes = length(xcoords)
    num_segments = num_nodes

    # n_point, n_point_marker, n_point_attribute, n_segment, n_holes
    poly = TriangleMesh.Polygon_pslg(num_nodes, 1, 0, num_segments, 0)

    node = [xcoords ycoords]
    set_polygon_point!(poly, node)

    node_marker = ones(Int, num_nodes, 1)
    set_polygon_point_marker!(poly, node_marker)

    segments = zeros(Int, num_segments, 2)
    for i=1:num_segments

        if i == num_segments
            segments[i, 1:2] = [i, 1]
        else
            segments[i, 1:2] = [i, i+1]
        end

    end

    set_polygon_segment!(poly, segments)

    segment_markers = ones(Int, num_segments)
    set_polygon_segment_marker!(poly, segment_markers)

    #switches from https://www.cs.cmu.edu/~quake/triangle.html
    switches = "penvVa" * string(mesh_size) * "D"

    mesh = create_mesh(poly, switches)

    return mesh

end


function triangular_mesh_properties(mesh)

    #calculate cell area and centroid
    Ai = zeros(Float64, mesh.n_cell)
    cxi = zeros(Float64, mesh.n_cell)
    cyi = zeros(Float64, mesh.n_cell)

    for i = 1:mesh.n_cell

        p1 = mesh.cell[1, i]
        p2 = mesh.cell[2, i]
        p3 = mesh.cell[3, i]

        x1 = mesh.point[1, p1]
        y1 = mesh.point[2,p1]
        x2 = mesh.point[1, p2]
        y2 = mesh.point[2,p2]
        x3 = mesh.point[1, p3]
        y3 = mesh.point[2,p3]

    Ai[i] = Geometry.triangle_area(x1, y1, x2, y2, x3, y3)

    cxi[i], cyi[i] = Geometry.triangle_centroid(x1, y1, x2, y2, x3, y3)

    end

    return Ai, cxi, cyi

end



function mesh(xcoords, ycoords, mesh_size)

    section_mesh = CrossSection.triangular_mesh(xcoords, ycoords, mesh_size)
    Ai, cxi, cyi = CrossSection.triangular_mesh_properties(section_mesh)

    return Ai, cxi, cyi

 end


 function rectangular_tube_geometry(feature)


    xcoords = []
    ycoords = []

    #convert feature vectors to xy components
    Δxy = Geometry.vector_components(feature.ΔL, feature.θ)

    #number of straight line segments in the feature
    num_lines = size(Δxy)[1]

    #get xy line coordinates of feature
    xcoords_straight, ycoords_straight = Geometry.line_coordinates(Δxy, feature.closed_or_open)

    #calculate feature surface normals
    unitnormals = CrossSection.surface_normals(xcoords_straight, ycoords_straight, feature.closed_or_open)

    #calculate feature node normals
    nodenormals = CrossSection.avg_node_normals(unitnormals, feature.closed_or_open)

    #calculate interior corner angle magnitudes
    num_radius = length(feature.radius)
    interior_radius_angle = zeros(Float64, num_radius)

    for i=1:num_radius

        θ1 = feature.θ[i] - 180

        if (feature.closed_or_open == 0) & (i==num_radius)  #for closed section, return to beginning node
            θ2 = feature.θ[1]
        else
            θ2 = feature.θ[i+1]
        end

        unit_vector_i = [cos(deg2rad(θ1)), sin(deg2rad(θ1))]
        unit_vector_j = [cos(deg2rad(θ2)), sin(deg2rad(θ2))]

        interior_radius_angle[i] = rad2deg(acos(dot(unit_vector_i, unit_vector_j)))

    end

    #calculate distance from curve PI to start and end of curve, along tangents
    Tc = zeros(Float64, length(feature.radius))

    for i = 1:length(feature.radius)

        radius = feature.radius[i]

        if radius > 0.0


            if feature.closed_or_open == 0
                xy_PI = [xcoords_straight[i], ycoords_straight[i]]
            elseif feature.closed_or_open == 1
                xy_PI = [xcoords_straight[i+1], ycoords_straight[i+1]]
            end

            n = feature.n_radius[i]
            γ = interior_radius_angle[i]


            θ1 = feature.θ[i]


            if (feature.closed_or_open == 0) & (i == length(feature.radius))  #for closed section, return to beginning node
                θ2 = feature.θ[1]
            else
                θ2 = feature.θ[i+1]
            end


            if (feature.closed_or_open == 0) & (i == length(feature.radius))

                if (sign(θ1) <=0) & (sign(θ2) >= 0)
                    PI_unit_normal = -nodenormals[1,:]
                else
                    PI_unit_normal = nodenormals[1,:]
                end

            else

                if (sign(θ1) <=0) & (sign(θ2) >= 0)
                    PI_unit_normal = -nodenormals[i+1,:]
                else
                    PI_unit_normal = nodenormals[i+1,:]
                end

            end


            xy_curve, Δ, E, Tc[i], xy_o, xy_BC, xy_EC, BC_unit_tangent, EC_unit_tangent, radius_unit_vector_BC = Geometry.circular_curve(radius, γ, xy_PI, PI_unit_normal, n)

        elseif radius == 0.0
            Tc[i] = 0.0
        end

    end

    #Shorten vector lengths to include curves.
    ΔLc = zeros(Float64, num_lines)

    for i = 1:length(feature.ΔL)

        if i == 1  #first
            if feature.closed_or_open == 1
                ΔLc[i] = feature.ΔL[i] - Tc[i]   #j end
            elseif feature.closed_or_open == 0
                ΔLc[i] = feature.ΔL[i] - Tc[i] - Tc[end]
            end
        elseif i == length(feature.ΔL)  #last
            if feature.closed_or_open == 1
                ΔLc[i] = feature.ΔL[i] - Tc[end]  #j end
            elseif feature.closed_or_open == 0
                ΔLc[i] = feature.ΔL[i] - Tc[end] - Tc[i - 1]
            end
        else  #others
            ΔLc[i] = feature.ΔL[i] - Tc[i-1] - Tc[i] #i and j ends
        end

    end

    #Get xy components of flat lengths around tube.
    Δxyc = Geometry.vector_components(ΔLc, feature.θ)


    #Discretize feature.
    dx = Δxyc[:,1] ./ feature.n
    dy = Δxyc[:,2] ./ feature.n

    #Define vertical tube segment at y=0.

    ycoords = Tc[1]:dy[1]:(Δxyc[1,2] + Tc[1])
	xcoords = zeros(Float64, length(ycoords))

    #Define top left corner.

    n = feature.n_radius[1]
    radius = feature.radius[1]
    γ = interior_radius_angle[1]
    PI_unit_normal = nodenormals[2,:]
    xy_PI = [xcoords_straight[2], ycoords_straight[2]]

    xy_curve, Δ, E, T, xy_o, xy_BC, xy_EC, BC_unit_tangent, EC_unit_tangent, radius_unit_vector_BC = Geometry.circular_curve(radius, γ, xy_PI, PI_unit_normal, n)

    xcoords = [xcoords; xy_curve[2:end, 1]]
    ycoords = [ycoords; xy_curve[2:end, 2]]

    #Define top flat.

    x_range = xcoords[end]:dx[2]:xcoords[end]+Δxyc[2,1]
    y_range = ycoords[end] * ones(length(x_range))
    xcoords = [xcoords; x_range[2:end]]
    ycoords = [ycoords; y_range[2:end]]

    #Define top right corner.

    n = feature.n_radius[2]
    radius = feature.radius[2]
    γ = interior_radius_angle[2]
    PI_unit_normal = nodenormals[3,:]
    xy_PI = [xcoords_straight[3], ycoords_straight[3]]

    xy_curve, Δ, E, T, xy_o, xy_BC, xy_EC, BC_unit_tangent, EC_unit_tangent, radius_unit_vector_BC = Geometry.circular_curve(radius, γ, xy_PI, PI_unit_normal, n)

    xcoords = [xcoords; xy_curve[2:end, 1]]
    ycoords = [ycoords; xy_curve[2:end, 2]]

    #Define right vertical flat.

    y_range = ycoords[end]:dy[3]:(ycoords[end] + Δxyc[3,2])
    x_range = xcoords[end] * ones(length(y_range))
    xcoords = [xcoords; x_range[2:end]]
    ycoords = [ycoords; y_range[2:end]]

    #Define bottom right corner.

    n = feature.n_radius[3]
    radius = feature.radius[3]
    γ = interior_radius_angle[3]
    PI_unit_normal = nodenormals[4,:]
    xy_PI = [xcoords_straight[4], ycoords_straight[4]]

    xy_curve, Δ, E, T, xy_o, xy_BC, xy_EC, BC_unit_tangent, EC_unit_tangent, radius_unit_vector_BC = Geometry.circular_curve(radius, γ, xy_PI, PI_unit_normal, n)

    xcoords = [xcoords; xy_curve[2:end, 1]]
    ycoords = [ycoords; xy_curve[2:end, 2]]

    #Define bottom horizontal flat.

    x_range = xcoords[end]:dx[4]:xcoords[end]+Δxyc[4,1]
    y_range = ycoords[end] * ones(length(x_range))
    xcoords = [xcoords; x_range[2:end]]
    ycoords = [ycoords; y_range[2:end]]

    #Define lower left corner.

    n = feature.n_radius[4]
    radius = feature.radius[4]
    γ = interior_radius_angle[4]
    PI_unit_normal = nodenormals[1,:]
    xy_PI = [xcoords_straight[1], ycoords_straight[1]]

    xy_curve, Δ, E, T, xy_o, xy_BC, xy_EC, BC_unit_tangent, EC_unit_tangent, radius_unit_vector_BC = Geometry.circular_curve(radius, γ, xy_PI, PI_unit_normal, n)

    xcoords = [xcoords; xy_curve[2:(end-1), 1]]
    ycoords = [ycoords; xy_curve[2:(end-1), 2]]

    return xcoords, ycoords

end


end #module
