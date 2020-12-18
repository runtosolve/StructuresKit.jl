module Geometry

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



function line_coordinates(Δxy)

    xcoords = []
    ycoords = []

    num_lines = size(Δxy)[1]

    for i = 1:num_lines

        if i==1
            xcoords = [0.0, Δxy[i,1]]    
            ycoords = [0.0, Δxy[i,2]]
        else
            xcoords = [xcoords; xcoords[end] .+ Δxy[i,1]]
            ycoords = [ycoords; ycoords[end] .+ Δxy[i,2]]
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


end #module