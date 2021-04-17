module Connections

using Interpolations
using CSV, DataFrames

export cfs_rot_screwfastened_k, cfs_trans_screwfastened_k, cfs_pull_through_plate_stiffness

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

    if (tw > maximum(t_range)) | (tw < minimum(t_range))
        print("Beware!  The deck base metal thickness, tw, is outside the tested range.")
    end

    if (c > maximum(c_range)) | (c < minimum(c_range))
        print("Beware!  The fastener distance from the purlin pivot point, c, is outside the tested range.")
    end

    if (x > maximum(x_range)) | (x < minimum(x_range))
        print("Beware!  The fastener distance from the R-panel major rib, x, is outside the tested range.")
    end

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


end #module
