
using CSV
using Plots
using Interpolations

# Load the experimental data from the CSV file
el = DataFrame(CSV.File("D:/JHU Spring 2021/Fastener Testing Data_Zhidong/Trail1/Mon (2).csv"))
df_matrix = Matrix(el)

# preparting the data information for each dimensions in the following interpolation
t1 = df_matrix[1:20, 10]
t2 = df_matrix[1:20, 11]
fy1 = df_matrix[1:20, 12]
fy2 = df_matrix[1:20, 13]
fu1 = df_matrix[1:20, 14]
fu2 = df_matrix[1:20, 15]
screw_diameter = df_matrix[1:20, 16]
D1_range = df_matrix[1:20, 2]
D2_range = df_matrix[1:20, 3]
D3_range = df_matrix[1:20, 4]
D4_range = df_matrix[1:20, 5]
F_range = df_matrix[1:20,8]


function define_interpolation_operator(t1, t2, fy1, fy2, fu1, fu2, screw_diameter,D)

    t1_range = sort(unique(t1))
    t2_range = sort(unique(t2))
    fy1_range = sort(unique(fy1))
    fy2_range = sort(unique(fy2))
    fu1_range = sort(unique(fu1))
    fu2_range = sort(unique(fu2))
    screw_diameter_range = sort(unique(screw_diameter))
    
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

                                index_t1 = findall(==(t1_entry), t1)
                                index_t2 = findall(==(t2_entry), t2)
                                index_fy1 = findall(==(fy1_entry), fy1)
                                index_fy2 = findall(==(fy2_entry), fy2)
                                index_fu1 = findall(==(fu1_entry), fu1)
                                index_fu2 = findall(==(fu2_entry), fu2)
                                index_screw_diameter = findall(==(screw_diameter_entry), screw_diameter)

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

function target_stiffness_operator(t1, t2, fy1, fy2, fu1, fu2, screw_diameter, D1_range, D2_range, D3_range, D4_range, F_range, t1_target, t2_target, fy1_target, fy2_target, fu1_target, fu2_target, screw_diameter_target)

    t1_range = sort(unique(t1))
    t2_range = sort(unique(t2))
    fy1_range = sort(unique(fy1))
    fy2_range = sort(unique(fy2))
    fu1_range = sort(unique(fu1))
    fu2_range = sort(unique(fu2))
    screw_diameter_range = sort(unique(screw_diameter))

    if (t1_target > maximum(t1_range))
        println("Beware!  The thickness of ply 1, t1, is outside the tested range, too large.")
    elseif (t1_target < minimum(t1_range))
        println("Beware!  The thickness of ply 1, t1, is outside the tested range, too small.")
    end

    if (t2_target > maximum(t2_range))
        println("Beware!  The thickness of ply 2, t2, is outside the tested range, too large.")
    elseif (t2_target < minimum(t2_range))
        println("Beware!  The thickness of ply 2, t2, is outside the tested range. too small.")
    end

    if (fy1_target > maximum(fy1_range))
        println("Beware!  The yielding strength of ply 1, fy1, is outside the tested range, too large.")
    elseif (fy1_target < minimum(fy1_range))
        println("Beware!  The yielding strength of ply 1, fy1, is outside the tested range, too small.")
    end

    if (fy2_target > maximum(fy2_range))
        println("Beware!  The yielding strength of ply 2, fy2, is outside the tested range, too large.")
    elseif (fy2_target < minimum(fy2_range))
        println("Beware!  The yielding strength of ply 2, fy2, is outside the tested range, too small.")
    end

    if (fu1_target > maximum(fu1_range))
        println("Beware!  The ultimate strength of ply 1, fu1, is outside the tested range, too large.")
    elseif (fu1_target < minimum(fu1_range))
        println("Beware!  The ultimate strength of ply 1, fu1, is outside the tested range, too small.")
    end

    if (fu2_target > maximum(fu2_range)) 
        println("Beware!  The ultimate strength of ply 2, fu2, is outside the tested range, too large.")
    elseif (fu2_target < minimum(fu2_range))
        println("Beware!  The ultimate strength of ply 2, fu2, is outside the tested range, too small.")
    end

    if (screw_diameter_target > maximum(screw_diameter_range))
        println("Beware!  The screw diameter, screw_diameter, is outside the tested range, too large.")
    elseif (screw_diameter_target < minimum(screw_diameter_range))
        println("Beware!  The screw diameter, screw_diameter, is outside the tested range, too small.")
    end

    D1_operator = define_interpolation_operator(t1, t2, fy1, fy2, fu1, fu2, screw_diameter,D1_range)
    D1_model =  LinearInterpolation((t1_range, t2_range, fy1_range, fy2_range, fu1_range, fu2_range, screw_diameter_range), D1_operator, extrapolation_bc = Line())
    D1_emperical = D1_model(t1_target, t2_target, fy1_target, fy2_target, fu1_target, fu2_target, screw_diameter_target)

    D2_operator = define_interpolation_operator(t1, t2, fy1, fy2, fu1, fu2, screw_diameter,D2_range)
    D2_model =  LinearInterpolation((t1_range, t2_range, fy1_range, fy2_range, fu1_range, fu2_range, screw_diameter_range), D2_operator, extrapolation_bc = Line())
    D2_emperical = D2_model(t1_target, t2_target, fy1_target, fy2_target, fu1_target, fu2_target, screw_diameter_target)

    D3_operator = define_interpolation_operator(t1, t2, fy1, fy2, fu1, fu2, screw_diameter,D3_range)
    D3_model =  LinearInterpolation((t1_range, t2_range, fy1_range, fy2_range, fu1_range, fu2_range, screw_diameter_range), D3_operator, extrapolation_bc = Line())
    D3_emperical = D3_model(t1_target, t2_target, fy1_target, fy2_target, fu1_target, fu2_target, screw_diameter_target)

    D4_operator = define_interpolation_operator(t1, t2, fy1, fy2, fu1, fu2, screw_diameter,D4_range)
    D4_model =  LinearInterpolation((t1_range, t2_range, fy1_range, fy2_range, fu1_range, fu2_range, screw_diameter_range), D4_operator, extrapolation_bc = Line())
    D4_emperical = D4_model(t1_target, t2_target, fy1_target, fy2_target, fu1_target, fu2_target, screw_diameter_target)

    F_operator = define_interpolation_operator(t1, t2, fy1, fy2, fu1, fu2, screw_diameter,F_range)
    F_model =  LinearInterpolation((t1_range, t2_range, fy1_range, fy2_range, fu1_range, fu2_range, screw_diameter_range), F_operator, extrapolation_bc = Line())
    F_emperical = F_model(t1_target, t2_target, fy1_target, fy2_target, fu1_target, fu2_target, screw_diameter_target)

    F1 = 0.4*F_emperical
    F2 = 0.8*F_emperical
    F4 = 0.1*F_emperical

    return D1_emperical, D2_emperical, D3_emperical, D4_emperical, F1, F2, F_emperical, F4

end

t1_target1 = 0.78
t1_target2 = 3.6
t2_target = 2.54
fy1_target = 150.4
fy2_target = 422.4
fu1_target = 312.1
fu2_target = 534.2
screw_diameter_target = 4.75

prediction = target_stiffness_operator(t1, t2, fy1, fy2, fu1, fu2, screw_diameter, D1_range, D2_range, D3_range, D4_range,F_range, t1_target1, t2_target, fy1_target, fy2_target, fu1_target, fu2_target, screw_diameter_target)
