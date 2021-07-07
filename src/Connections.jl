module Connections

using Interpolations
using CSV, DataFrames

export cfs_rot_screwfastened_k, cfs_trans_screwfastened_k, cfs_pull_through_plate_stiffness, cfs_load_deformation_interpolation_model

#Cold-formed steel
#Rotational stiffness for a screw-fastened connection
#Distributed spring!
#https://doi.org/10.1016/j.tws.2012.06.005
function cfs_rot_screwfastened_k(b, c, S, t, kp, E, CorZ)

    I = 1/12*S*t^3

    if CorZ == 0   #C

        kϕ = (S/(c^2*kp)+(((b^2*c/2)+c^2*b+(c^3/3))*S/(c^2*E*I)))^-1

    elseif CorZ == 1  #Z

        kϕ = (S/(c^2*kp) + c*S/(3*E*I))^-1

    end

    return kϕ

end


#Cold-formed steel
#Translational stiffness for a screw-fastened connection
#Discrete spring!
#https://www.buildusingsteel.org/aisi-design-resources/-/media/doc/buildusingsteel/research-reports/CFSD%20-%20Report%20-%20RP17-2.pdf
function cfs_trans_screwfastened_k(t1, t2, E1, E2, Fss, Fu1, Fu2, D)

    Ka = (1/(E1*t1)+1/(E2*t2))^-1

    ψ=(Fss/(t1*D*Fu1))*(Fss/(t2*D*Fu2))

    #monotonic elastic, steel to steel
    α = 0.27
    β = - 0.69

    Ke = α*ψ^β*Ka

    return Ka, ψ, α, β, Ke

end


"""
        cfs_pull_through_plate_stiffness(x, c, tw)

`x` is the fastener distance from an R-panel major deck rib, in mm.
`c` is the fastener distance from the purlin rotation point, in mm.
'tw' is the R-panel deck base metal thickness, in mm.

The fastener pull-through elastic stiffness 'kp' is interpolated and output in N/mm.

These results come from elastic finite element studies in
https://www.buildusingsteel.org/aisi-design-resources/-/media/doc/buildusingsteel/research-reports/CFSD%20-%20Report%20-%20RP17-2.pdf

"""

function cfs_pull_through_plate_stiffness(x, c, tw)

    #Define source file for pullover stiffness values.
    filename = string(@__DIR__, "/assets/kp.csv")

    #Import the parameters and stiffness values.

    data = CSV.File(filename)
    data = DataFrame(x = Vector(data.x), c = Vector(data.c), tw = Vector(data.tw), kp=Vector(data.kp))
    # data = CSV.read(filename, DataFrame, header=true)

    #Define the range of screw locations relative to an R-panel deck rib.
    x_range = sort(unique(data.x))

    #Define the range of fastener locations from the purlin pivot point.
    c_range = sort(unique(data.c))

    #Define the range of deck panel thicknesses considered.
    t_range = sort(unique(data.tw))

    #Check if the screw location in the R-panel is within the tested range.

    # if (tw > maximum(t_range)) | (tw < minimum(t_range))
    #     print("Beware!  The deck base metal thickness, tw, is outside the tested range.")
    # end

    # if (c > maximum(c_range)) | (c < minimum(c_range))
    #     print("Beware!  The fastener distance from the purlin pivot point, c, is outside the tested range.")
    # end

    # if (x > maximum(x_range)) | (x < minimum(x_range))
    #     print("Beware!  The fastener distance from the R-panel major rib, x, is outside the tested range.")
    # end

    #Define the number of parameters that will be used to predict kp.

    x_dim = length(x_range)
    c_dim = length(c_range)
    t_dim = length(t_range)

    #Initialize the parameter array.
    A = zeros(Float64, (x_dim,c_dim, t_dim))

    #Fill the parameter array with kp.
    for i = 1:x_dim
        for j = 1:c_dim
            for k = 1:t_dim
                x_entry = x_range[i]
                c_entry = c_range[j]
                t_entry = t_range[k]

                index_x = findall(==(x_entry), data.x)
                index_c = findall(==(c_entry), data.c)
                index_t = findall(==(t_entry), data.tw)

                index = intersect(index_x, index_c)
                index = intersect(index, index_t)

                A[i, j, k] = data.kp[index[1]]

            end
        end
    end

    #Define the interpolation operator.  Only linear interpolation is available for multi-dimensional parameter sets in Interpolations.jl.
    #Allow extrapolation with warnings above.
    stiffness_model = LinearInterpolation((x_range, c_range, t_range), A, extrapolation_bc = Line())

    kp = stiffness_model(x, c, tw)

    return kp

end


function cfs_load_deformation_interpolation_model(t1, t2, fy1, fy2, fu1, fu2, screw_diameter)

    #Define source file for experimental data used to define the interpolation model.
    filename = string(@__DIR__, "/assets/Connections/cfs_fastener_test_data.csv")

    #Load experimental data from the CSV file.
    data = CSV.File(filename)

    #Define sorted, unique ranges of each parameter.
    t1_range = sort(unique(data.t1))
    t2_range = sort(unique(data.t2))
    fy1_range = sort(unique(data.fy1))
    fy2_range = sort(unique(data.fy2))
    fu1_range = sort(unique(data.fu1))
    fu2_range = sort(unique(data.fu2))
    screw_diameter_range = sort(unique(data.screw_diameter))

    #Define interpolation operators for all the deformation points (D1, D2, D3, D4) and the peak load (P3).
    D1_operator = define_interpolation_operator(data, t1_range, t2_range, fy1_range, fy2_range, fu1_range, fu2_range, screw_diameter_range, data.D1)
    D2_operator = define_interpolation_operator(data, t1_range, t2_range, fy1_range, fy2_range, fu1_range, fu2_range, screw_diameter_range, data.D2)
    D3_operator = define_interpolation_operator(data, t1_range, t2_range, fy1_range, fy2_range, fu1_range, fu2_range, screw_diameter_range, data.D3)
    D4_operator = define_interpolation_operator(data, t1_range, t2_range, fy1_range, fy2_range, fu1_range, fu2_range, screw_diameter_range, data.D4)
    P3_operator = define_interpolation_operator(data, t1_range, t2_range, fy1_range, fy2_range, fu1_range, fu2_range, screw_diameter_range, data.P3)
    
    #Define interpolation models for all the deformation points (D1, D2, D3, D4) and the peak load (P3).
    D1_model =  LinearInterpolation((t1_range, t2_range, fy1_range, fy2_range, fu1_range, fu2_range, screw_diameter_range), D1_operator, extrapolation_bc = Line())
    D2_model =  LinearInterpolation((t1_range, t2_range, fy1_range, fy2_range, fu1_range, fu2_range, screw_diameter_range), D2_operator, extrapolation_bc = Line())
    D3_model =  LinearInterpolation((t1_range, t2_range, fy1_range, fy2_range, fu1_range, fu2_range, screw_diameter_range), D3_operator, extrapolation_bc = Line())
    D4_model =  LinearInterpolation((t1_range, t2_range, fy1_range, fy2_range, fu1_range, fu2_range, screw_diameter_range), D4_operator, extrapolation_bc = Line())
    P3_model =  LinearInterpolation((t1_range, t2_range, fy1_range, fy2_range, fu1_range, fu2_range, screw_diameter_range), P3_operator, extrapolation_bc = Line())
    
    #Perform interpolation predictions based on user-defined values.
    D1 = D1_model(t1, t2, fy1, fy2, fu1, fu2, screw_diameter)
    D2 = D2_model(t1, t2, fy1, fy2, fu1, fu2, screw_diameter)
    D3 = D3_model(t1, t2, fy1, fy2, fu1, fu2, screw_diameter)
    D4 = D4_model(t1, t2, fy1, fy2, fu1, fu2, screw_diameter)
    P3 = P3_model(t1, t2, fy1, fy2, fu1, fu2, screw_diameter)

    #Calculate P1, P2, and P4. 
    P1 = 0.4 * P3
    P2 = 0.8 * P3
    P4 = 0.1 * P3

    #Warn user if the interpolation is out of range.
    interpolation_warning(t1, t2, fy1, fy2, fu1, fu2, screw_diameter, t1_range, t2_range, fy1_range, fy2_range, fu1_range, fu2_range, screw_diameter_range)

    return D1, D2, D3, D4, P1, P2, P3, P4

end


function define_interpolation_operator(data, t1_range, t2_range, fy1_range, fy2_range, fu1_range, fu2_range, screw_diameter_range, D)
    
    t1_dim = length(t1_range)
    t2_dim = length(t2_range)
    fy1_dim = length(fy1_range)
    fy2_dim = length(fy2_range)
    fu1_dim = length(fu1_range)
    fu2_dim = length(fu2_range)
    screw_diameter_dim = length(screw_diameter_range)

    #Initialize the parameter array.
    A = zeros(Float64, (t1_dim, t2_dim, fy1_dim, fy2_dim, fu1_dim, fu2_dim, screw_diameter_dim))

    #Fill the parameter array with kp.
    for i = 1:t1_dim
        for j = 1:t2_dim
            for k = 1:fy1_dim
                for n = 1:fy2_dim
                    for m = 1:fu1_dim
                        for h = 1:fu2_dim
                            for s  = 1:screw_diameter_dim
                                t1_entry = t1_range[i]
                                t2_entry = t2_range[j]
                                fy1_entry = fy1_range[k]
                                fy2_entry = fy2_range[n]
                                fu1_entry = fu1_range[m]
                                fu2_entry = fu2_range[h]
                                screw_diameter_entry = screw_diameter_range[s]

                                index_t1 = findall(==(t1_entry), data.t1)
                                index_t2 = findall(==(t2_entry), data.t2)
                                index_fy1 = findall(==(fy1_entry), data.fy1)
                                index_fy2 = findall(==(fy2_entry), data.fy2)
                                index_fu1 = findall(==(fu1_entry), data.fu1)
                                index_fu2 = findall(==(fu2_entry), data.fu2)
                                index_screw_diameter = findall(==(screw_diameter_entry), data.screw_diameter)

                                index = intersect(index_t1, index_t2)
                                index = intersect(index, index_fy1)
                                index = intersect(index, index_fy2)
                                index = intersect(index, index_fu1)
                                index = intersect(index, index_fu2)
                                index = intersect(index, index_screw_diameter)
                                

                                if size(index,1) == 0
                                    A[i, j, k, n, m, h, s] = 0
                                else
                                    A[i, j, k, n, m, h, s] = D[index[1]]
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return A
end


function interpolation_warning(t1, t2, fy1, fy2, fu1, fu2, screw_diameter, t1_range, t2_range, fy1_range, fy2_range, fu1_range, fu2_range, screw_diameter_range)

    if (t1 > maximum(t1_range))
        println("Beware!  The thickness of ply 1, t1, is outside the tested range, too large.")
    elseif (t1 < minimum(t1_range))
        println("Beware!  The thickness of ply 1, t1, is outside the tested range, too small.")
    end

    if (t2 > maximum(t2_range))
        println("Beware!  The thickness of ply 2, t2, is outside the tested range, too large.")
    elseif (t2 < minimum(t2_range))
        println("Beware!  The thickness of ply 2, t2, is outside the tested range. too small.")
    end

    if (fy1 > maximum(fy1_range))
        println("Beware!  The yielding strength of ply 1, fy1, is outside the tested range, too large.")
    elseif (fy1 < minimum(fy1_range))
        println("Beware!  The yielding strength of ply 1, fy1, is outside the tested range, too small.")
    end

    if (fy2 > maximum(fy2_range))
        println("Beware!  The yielding strength of ply 2, fy2, is outside the tested range, too large.")
    elseif (fy2 < minimum(fy2_range))
        println("Beware!  The yielding strength of ply 2, fy2, is outside the tested range, too small.")
    end

    if (fu1 > maximum(fu1_range))
        println("Beware!  The ultimate strength of ply 1, fu1, is outside the tested range, too large.")
    elseif (fu1 < minimum(fu1_range))
        println("Beware!  The ultimate strength of ply 1, fu1, is outside the tested range, too small.")
    end

    if (fu2 > maximum(fu2_range)) 
        println("Beware!  The ultimate strength of ply 2, fu2, is outside the tested range, too large.")
    elseif (fu2 < minimum(fu2_range))
        println("Beware!  The ultimate strength of ply 2, fu2, is outside the tested range, too small.")
    end

    if (screw_diameter > maximum(screw_diameter_range))
        println("Beware!  The screw diameter, screw_diameter, is outside the tested range, too large.")
    elseif (screw_diameter < minimum(screw_diameter_range))
        println("Beware!  The screw diameter, screw_diameter, is outside the tested range, too small.")
    end

end


end #module
