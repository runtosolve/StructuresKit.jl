module MaterialModels

using Polynomials

export steel

function steel(σy, σy1, σu, σf, ϵy, ϵy1, ϵu, ϵf, n)

    # define elastic part of steel stress-strain curve
    σ=[0 σy]
    ϵ=[0 ϵy]

    #define plateau and hardening part of steel stress-strain curve
    p = fit([ϵy1, ϵu, ϵf],[σy1, σu, σf], 2)

    σ_fit = zeros(n)
    ϵ_range = ϵy1:(ϵf-ϵy1)./(n-1):ϵf

    for i=1:n 
        σ_fit[i] = p(ϵ_range[i])
    end

    ϵ = vcat(vec(ϵ),collect(range(ϵy1,length=n,stop=ϵf)))
    σ = vcat(vec(σ), σ_fit)

    # add compression part of steel stress-strain curve
    σ = vcat(-reverse(σ[2:end]), σ)
    ϵ = vcat(-reverse(ϵ[2:end]), ϵ)

    return ϵ, σ

end





end

