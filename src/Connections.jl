module Connections

export cfs_rot_screwfastened_k, cfs_trans_screwfastened_k

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




end #module
