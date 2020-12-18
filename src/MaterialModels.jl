module MaterialModels

using Polynomials

using ..Geometry

export steel

function steel(σ_prop_limit, σy, σy1, σu, σf, ϵ_prop_limit, ϵy, ϵy1, ϵu, ϵf, n)

    # define elastic part of steel stress-strain curve
    σ = 0.0:σy/n[1]:σ_prop_limit
    ϵ = 0.0:ϵy/n[1]:ϵ_prop_limit

    #connection linear part of curve to yield plateau

    #work on proportional limit!


    closed_or_open = 1
    unit_normals = CrossSection.surface_normals([0.0, ϵy, ϵy1], [0.0, σy, σy1], closed_or_open)
    node_normal = CrossSection.avg_node_normals(unit_normals, closed_or_open)

    interior_angle = 90


    radius = 2.0
    xy_PI = [ϵy, σy]
    PI_unit_normal = node_normal[2,:]
    n_radius = 4
    xy_curve, Δ, E, T, xy_o, xy_BC, xy_EC, BC_unit_tangent, EC_unit_tangent, radius_unit_vector_BC = Geometry.circular_curve(radius, interior_angle, xy_PI, PI_unit_normal, n_radius)







    #define plateau region
    plateau_slope = (σy1 - σy) / (ϵy1 - ϵy)
    ϵ_plateau = (ϵy + (ϵy1 - ϵy)/n[2]):(ϵy1 - ϵy)/n[2]:(ϵy1 - (ϵy1 - ϵy)/n[2])
    σ_plateau = plateau_slope .* (ϵ_plateau .- ϵy) .+ σy

    ϵ = [ϵ; ϵ_plateau]
    σ = [σ; σ_plateau]

    #define plateau and hardening part of steel stress-strain curve
    p = fit([ϵy1, ϵu, ϵf],[σy1, σu, σf], 2)

    ϵ_range = ϵy1:(ϵf-ϵy1)./n[3]:ϵf

    σ_fit = zeros(length(ϵ_range))

    for i=1:length(ϵ_range)
        σ_fit[i] = p(ϵ_range[i])
    end

    ϵ = [ϵ; ϵ_range]
    σ = [σ; σ_fit]

    # add compression part of steel stress-strain curve
    # σ = vcat(-reverse(σ[2:end]), σ)
    # ϵ = vcat(-reverse(ϵ[2:end]), ϵ)

    return ϵ, σ

end





end

