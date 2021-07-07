module Geometry

using StaticArrays
using LinearAlgebra

export rotation_matrix, vector_components, line_coordinates, circular_curve, triangle_area, triangle_centroid

#counterclockwise, right hand rule
function rotation_matrix(θ)

    R = [cos(θ) -sin(θ)
         sin(θ) cos(θ)]

end

function vector_components(ΔL, θ)

    num_lines = length(ΔL)

    Δxy = zeros(Float64, (num_lines, 2))

    for i = 1:num_lines

        Δxy[i,:] = [ΔL[i] * cos(deg2rad(θ[i]))  ΔL[i] * sin(deg2rad(θ[i]))]

    end

    return Δxy

end



function line_coordinates(Δxy, closed_or_open)

    xcoords = []
    ycoords = []

    num_lines = size(Δxy)[1]

    for i = 1:num_lines

        if i==1
            xcoords = [0.0, Δxy[i,1]]
            ycoords = [0.0, Δxy[i,2]]
        else

            if ((closed_or_open == 0) & (i == num_lines))
                #avoid last node for closed section
            else
                xcoords = [xcoords; xcoords[end] .+ Δxy[i,1]]
                ycoords = [ycoords; ycoords[end] .+ Δxy[i,2]]
            end
        end

    end

    return xcoords, ycoords

end

function circular_curve(radius, interior_angle, xy_PI, PI_unit_normal, n)

    #follow this nomenclature for a horizontal curve
    #https://www.cpp.edu/~hturner/ce220/circular_curves.pdf

    #interior angle
    θ = deg2rad(interior_angle)

    #angle over which the arc sweeps
    Δ = π - θ

    #define orthogonal distance from curve to PI
    E = radius*(1/(cos(Δ/2))-1)

    #define distance along tangent from PI to BC and EC
    T = radius*tan(Δ/2)

    #calculate vector that runs tangent to curve through BC
    #rotate PI_normal vector

    # if concave_up_or_down == 1
        BC_rotation = Δ + (π/2 - Δ/2)
    # elseif concave_up_or_down == 0
        # BC_rotation = -(Δ + (π/2 - Δ/2))
    # end
    R_BC = rotation_matrix(BC_rotation)
    BC_unit_tangent = -R_BC*PI_unit_normal  #- sign is to flip direction of vector, face towards direction of curve calculation

    #calculate vector that runs tangent to curve through EC
    #rotate PI_normal vector
    #if concave_up_or_down == 1
        EC_rotation = -(Δ + (π/2 - Δ/2))
    #elseif concave_up_or_down == 0
    #    EC_rotation = Δ + (π/2 - Δ/2)
    #end
    R_EC = rotation_matrix(EC_rotation)
    EC_unit_tangent = R_EC*PI_unit_normal

    #define x-y coordinates of BC
    xy_BC = xy_PI - T * BC_unit_tangent

    #define x-y coordinates of EC
    xy_EC = xy_PI + T * EC_unit_tangent

    #define x-y coordinates at center of circle
    xy_o = xy_PI - (E+radius) * PI_unit_normal

    #define unit radius vector going through BC
    #this is the unit vector normal to the tangent vector at BC
    radius_unit_vector_BC = [-BC_unit_tangent[2], BC_unit_tangent[1]]

    #rotate unit radius vector
    Δ_sweep = range(0.0, -Δ, length=n+1)
    xy_curve = zeros(Float64, (n+1, 2))
    xy_curve[1,:] = xy_BC
    xy_curve[end,:] = xy_EC

    for i = 2:n
        R = rotation_matrix(Δ_sweep[i])
        xy_curve[i,:] = xy_o + radius*(R*radius_unit_vector_BC)
    end

    return xy_curve, Δ, E, T, xy_o, xy_BC, xy_EC, BC_unit_tangent, EC_unit_tangent, radius_unit_vector_BC

end


# https://keisan.casio.com/has10/SpecExec.cgi?path=05000000.Mathematics%252F01000500.Plane%2520geometry%252F10010300.Area%2520of%2520a%2520triangle%2520with%2520three%2520points%252Fdefault.xml&charset=utf-8
function triangle_area(x1, y1, x2, y2, x3, y3)

    A = abs((x1*y2 + x2*y3 + x3*y1 - y1*x2 - y2*x3 - y3*x1)/2)

end

#https://www.mathopenref.com/coordcentroid.html
function triangle_centroid(x1, y1, x2, y2, x3, y3)

    cx = (x1 + x2 + x3)/3
    cy = (y1 + y2 + y3)/3

    return cx, cy

end


#Added from https://discourse.julialang.org/t/point-to-line-distance-in-geometrybasics-or-another-geometry-package/42501/5

"Point in N dimensions"
struct Point{N,T}
    x::SVector{N,T}
end

"Line in N dimensions. `p` is a point on the line and `u` is the direction vector
(not necessarily normalized). Parametrised as \$p + ut\$"
struct Line{N,T}
    p::SVector{N,T}
    u::SVector{N,T}
end

distance_between_two_points(p1::Point{N}, p2::Point{N}) where {N} = norm(p1.x - p2.x)

function distance_between_point_and_line(y::Point{N}, l::Line{N}) where {N}
    p, u = l.p, l.u
    
    t = (y.x - p) ⋅ u / (u ⋅ u) 
    x = Point(p + t*u)
    
    return distance_between_two_points(x, y)
end

#Added from https://discourse.julialang.org/t/point-to-line-distance-in-geometrybasics-or-another-geometry-package/42501/5


end #module
