module InternalForces

using DiffEqOperators: CenteredDifference

export moment, shear, torsion, bimoment, calculateDerivativeOperators


#############################################################
# Fornberg algorithm

# This implements the Fornberg (1988) algorithm (https://doi.org/10.1090/S0025-5718-1988-0935077-0)
# to obtain Finite Difference weights over arbitrary points to arbitrary order

function calculate_weights(order::Int, x0::T, x::AbstractVector) where T<:Real
    #=
        order: The derivative order for which we need the coefficients
        x0   : The point in the array 'x' for which we need the coefficients
        x    : A dummy array with relative coordinates, eg. central differences
               need coordinates centred at 0 while those at boundaries need
               coordinates starting from 0 to the end point
        The approximation order of the stencil is automatically determined from
        the number of requested stencil points.
    =#
    N = length(x)
    @assert order < N "Not enough points for the requested order."
    M = order
    c1 = one(T)
    c4 = x[1] - x0
    C = zeros(T, N, M+1)
    C[1,1] = 1
    @inbounds for i in 1 : N-1
        i1 = i + 1
        mn = min(i, M)
        c2 = one(T)
        c5 = c4
        c4 = x[i1] - x0
        for j in 0 : i-1
            j1 = j + 1
            c3 = x[i1] - x[j1]
            c2 *= c3
            if j == i-1
                for s in mn : -1 : 1
                    s1 = s + 1
                    C[i1,s1] = c1*(s*C[i,s] - c5*C[i,s1]) / c2
                end
                C[i1,1] = -c1*c5*C[i,1] / c2
           end
            for s in mn : -1 : 1
                s1 = s + 1
                C[j1,s1] = (c4*C[j1,s1] - s*C[j1,s]) / c3
            end
            C[j1,1] = c4 * C[j1,1] / c3
        end
        c1 = c2
    end
    #=
        This is to fix the problem of numerical instability which occurs when the sum of the stencil_coefficients is not
        exactly 0.
        https://scicomp.stackexchange.com/questions/11249/numerical-derivative-and-finite-difference-coefficients-any-update-of-the-fornb
        Stack Overflow answer on this issue.
        http://epubs.siam.org/doi/pdf/10.1137/S0036144596322507 - Modified Fornberg Algorithm
    =#
    _C = C[:,end]
    _C[div(N,2)+1] -= sum(_C)
    return _C
end



function operatorjumps(z, dm, A, order)

    diffdm = diff(dm)
    jumps = findall(x-> x > 0.0, diffdm)


    for i in eachindex(jumps)

        #left pointing single sided stencil, left of jump
        x = z[jumps[i]-4-1:jumps[i]-1]
        x0 = z[jumps[i]-1]
        coeffs=calculate_weights(order, x0, x)
        A[jumps[i]-1, 1:end] .= 0.0
        A[jumps[i]-1, jumps[i]-4-1:jumps[i]-1] = coeffs

        #left pointing single sided stencil, at jump
        x = z[jumps[i]-4:jumps[i]]
        x0 = z[jumps[i]]
        coeffs=calculate_weights(order, x0, x)
        A[jumps[i], 1:end] .= 0.0
        A[jumps[i], jumps[i]-4:jumps[i]] = coeffs

        # right pointing single sided stencil, right of jump
        x = z[jumps[i]+1:jumps[i]+1 + 4]
        x0 = z[jumps[i]+1]
        coeffs=calculate_weights(order, x0, x)
        A[jumps[i]+1, 1:end] .= 0.0
        A[jumps[i]+1, jumps[i]+1:jumps[i]+4+1] = coeffs

        # right pointing single sided stencil, two right of jump
        x = z[jumps[i]+2:jumps[i]+2 + 4]
        x0 = z[jumps[i]+2]
        coeffs=calculate_weights(order, x0, x)
        A[jumps[i]+2, 1:end] .= 0.0
        A[jumps[i]+2, jumps[i]+2:jumps[i]+4+2] = coeffs


    end

    return A

end



function calculateDerivativeOperators(z, dm)

    NumberOfNodes=length(z)
    dz=diff(z)
    dz=[dz[1]; dz; dz[end]]
    dL=dz

    #first derivative

    NthDerivative = 1
    DerivativeOrder = 2
    Az = CenteredDifference(NthDerivative, DerivativeOrder, dL, NumberOfNodes)
    Az=Array(Az)
    Az=Az[:,2:end-1]  #remove ghost nodes
    Az[1,:]=zeros(NumberOfNodes)
    Az[2,:]=zeros(NumberOfNodes)
    Az[(end-1),:]=zeros(NumberOfNodes)
    Az[end,:]=zeros(NumberOfNodes)

    #update boundary condition with single sided stencils
    order=1

    x=z[1:5]
    x0=x[1]
    coeffs=calculate_weights(order, x0, x)
    Az[1,1:5]=coeffs

    x=z[1:5]
    x0=x[2]
    coeffs=calculate_weights(order, x0, x)
    Az[2,1:5]=coeffs

    x=z[end-4:end]
    x0=x[end-1]
    coeffs=calculate_weights(order, x0, x)
    Az[end-1,end-4:end]=coeffs

    x=z[end-4:end]
    x0=x[end]
    coeffs=calculate_weights(order, x0, x)
    Az[end,end-4:end]=coeffs

    #place singled sided stencils at discontinuities
    Az = operatorjumps(z, dm, Az, order)

    #second derivative

    NthDerivative = 2
    DerivativeOrder = 2
    Azz = CenteredDifference(NthDerivative, DerivativeOrder, dL, NumberOfNodes)
    Azz=Array(Azz)
    Azz=Azz[:,2:end-1]
    Azz[1,:]=zeros(NumberOfNodes)
    Azz[2,:]=zeros(NumberOfNodes)
    Azz[(end-1),:]=zeros(NumberOfNodes)
    Azz[end,:]=zeros(NumberOfNodes)

    order=2   #update stencil to consider boundary conditions without ghost nodes

    x=z[1:5]
    x0=x[1]
    coeffs=calculate_weights(order, x0, x)
    Azz[1,1:5]=coeffs

    x=z[1:5]
    x0=x[2]
    coeffs=calculate_weights(order, x0, x)
    Azz[2,1:5]=coeffs

    x=z[end-4:end]
    x0=x[end-1]
    coeffs=calculate_weights(order, x0, x)
    Azz[end-1,end-4:end]=coeffs

    x=z[end-4:end]
    x0=x[end]
    coeffs=calculate_weights(order, x0, x)
    Azz[end,end-4:end]=coeffs

    #place singled sided stencils at discontinuities
    Azz = operatorjumps(z, dm, Azz, order)

    #third derivative

    NthDerivative = 3
    DerivativeOrder = 2
    Azzz = CenteredDifference(NthDerivative, DerivativeOrder, dL, NumberOfNodes)
    Azzz=Array(Azzz)
    Azzz=Azzz[:,2:end-1]
    Azzz[1,:]=zeros(NumberOfNodes)
    Azzz[2,:]=zeros(NumberOfNodes)
    Azzz[(end-1),:]=zeros(NumberOfNodes)
    Azzz[end,:]=zeros(NumberOfNodes)

    order=3

    x=z[1:5]
    x0=x[1]
    coeffs=calculate_weights(order, x0, x)
    Azzz[1,1:5]=coeffs

    x=z[1:5]
    x0=x[2]
    coeffs=calculate_weights(order, x0, x)
    Azzz[2,1:5]=coeffs

    x=z[end-4:end]
    x0=x[end-1]
    coeffs=calculate_weights(order, x0, x)
    Azzz[end-1,end-4:end]=coeffs

    x=z[end-4:end]
    x0=x[end]
    coeffs=calculate_weights(order, x0, x)
    Azzz[end,end-4:end]=coeffs

    #place singled sided stencils at discontinuities
    Azzz = operatorjumps(z, dm, Azzz, order)

    return Az, Azz, Azzz

end


function moment(z, dm, Δ, E, I)

    Az, Azz, Azzz = calculateDerivativeOperators(z, dm)

    M = E .* I .* Azz * Δ

    return M

end

function shear(z, dm, Δ, E, I)

    Az, Azz, Azzz = calculateDerivativeOperators(z, dm)

    V = E .* I .* Azzz * Δ

    return V

end

function torsion(z, dm, ϕ, E, G, J, Cw)

    Az, Azz, Azzz = calculateDerivativeOperators(z, dm)

    #T=GJ*dϕ/dz - ECw*dϕ3/dz3
    T=G.*J.*Az*ϕ .- E.*Cw.*Azzz*ϕ

    return T

end

function bimoment(z, dm, ϕ, E, Cw)

    Az, Azz, Azzz = calculateDerivativeOperators(z, dm)

    #B=ECwϕ''
    B=E.*Cw.*Azz*ϕ

    return B

end


end #module
