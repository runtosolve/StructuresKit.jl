module MaterialModels

using Polynomials
using LinearAlgebra
using Dierckx

using ..Geometry
using ..CrossSection
using ..Mesh


export steel, ConstitutiveLaw


struct ConstitutiveLaw

    stress_strain::Dierckx.Spline1D
    E_tan_strain::Dierckx.Spline1D

end


function steel(σy, σy1, σu, σf, ϵy, ϵy1, ϵu, ϵf, n)

    # define elastic and yielding portion of the steel stress-strain curve

    seg_1 = norm([σy, ϵy])
    seg_2 = norm([(σy-σy1), (ϵy-ϵy1)])
    
    θ1 = rad2deg(atan(σy/ϵy))
    θ2 = rad2deg(atan((σy-σy1) / (ϵy-ϵy1)))
    
    ΔL = [seg_1, seg_2]
    θ = [θ1, θ2]  #degrees
    n_ep = [n[1], n[2]]          #hard code this for now
    radius = [0.001]  #hard code this for now
    n_radius = [400]     #hard code this for now
    closed_or_open = 1
    
    
   elastic_plastic = CrossSection.Feature(ΔL, θ, n_ep, radius, n_radius, closed_or_open)
   ϵ, σ = CrossSection.get_xy_coordinates(elastic_plastic)


    # #define plateau region
    # plateau_slope = (σy1 - σy) / (ϵy1 - ϵy)
    # ϵ_plateau = (ϵy + (ϵy1 - ϵy)/n[2]):(ϵy1 - ϵy)/n[2]:(ϵy1 - (ϵy1 - ϵy)/n[2])
    # σ_plateau = plateau_slope .* (ϵ_plateau .- ϵy) .+ σy

    # ϵ = [ϵ; ϵ_plateau]
    # σ = [σ; σ_plateau]

    #define plateau and hardening part of steel stress-strain curve
    p = fit([ϵy1, ϵu, ϵf],[σy1, σu, σf], 2)

    ϵ_range = ϵy1:(ϵf-ϵy1)./n[3]:ϵf

    σ_fit = zeros(length(ϵ_range))

    for i=1:length(ϵ_range)
        σ_fit[i] = p(ϵ_range[i])
    end

    ϵ = [ϵ; ϵ_range[2:end]]
    σ = [σ; σ_fit[2:end]]

    #calculate tangent elastic modulus
    num_points = length(ϵ)
    E_tan = zeros(Float64,  num_points)
    order = 1
    for i = 1:num_points
 
        if i < num_points - 2

            x = ϵ[i:i+2]
            x0 = x[1]
    
            stencil = Mesh.calculate_weights(order, x0, x)

            E_tan[i] = sum(stencil .* σ[i:i+2])

        elseif i > num_points - 2

            x = ϵ[i-2:i]
            x0 = x[3]
    
            stencil = Mesh.calculate_weights(order, x0, x)

            E_tan[i] = sum(stencil .* σ[i-2:i])

        end

        #tough to round the corner at yielding, revisit this with better model of proportional limit
        if E_tan[i] > E_tan[1]

            E_tan[i] = E_tan[i-1]

        end

    end

    stress_strain = Spline1D(ϵ, σ)
    E_tan_strain = Spline1D(ϵ, E_tan)

    steel_law = ConstitutiveLaw(stress_strain, E_tan_strain)

    return steel_law

end





end

