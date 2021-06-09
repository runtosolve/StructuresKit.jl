module BeamColumn

using DiffEqOperators: CenteredDifference
using LinearAlgebra
using RecursiveArrayTools: VectorOfArray
using NLsolve

using ..Mesh


export define, solve, Model 



mutable struct Model

   member_definitions::Vector{Tuple{Float64, Float64, Int64, Int64, Int64, Int64, Int64}}
   section_properties::Vector{Tuple{Float64, Float64, Float64, Float64, Float64}}
   material_properties::Vector{Tuple{Float64, Float64}}
   spring_stiffness::Vector{Tuple{Float64, Float64}}
   spring_location::Vector{Tuple{Float64}}
   supports::Array{Float64}
   loads::Tuple{Float64, Float64, }
   load_location::Vector{Tuple{Float64, Float64}}
   end_boundary_conditions::Array{Int64}

   qx::Array{Float64}
   qy::Array{Float64}
   P::Array{Float64}

   K::Matrix{Float64}
   F::Array{Float64}

   free_dof::Array{Int64}
  
   A::Array{Float64}
   Ix::Array{Float64}
   Iy::Array{Float64}
   Io::Array{Float64}
   J::Array{Float64}
   Cw::Array{Float64}
   xc::Array{Float64}
   yc::Array{Float64}
   xs::Array{Float64}
   ys::Array{Float64}
   xo::Array{Float64}  #distance from centroid to shear center
   yo::Array{Float64}  #distance from centroid to shear center

   E::Array{Float64}
   ν::Array{Float64}
   G::Array{Float64}

   ax::Array{Float64} #distance from shear center to qy applied load location
   ay::Array{Float64} #distance from shear center to qx applied load location

   kx::Array{Float64} 
   ky::Array{Float64}
   kϕ::Array{Float64}

   hx::Array{Float64} #distance from centroid to ky spring location
   hy::Array{Float64} #distance from centroid to kx spring location

   z::Array{Float64}
   dm::Array{Int64}
   
   u::Array{Float64}
   v::Array{Float64}
   ϕ::Array{Float64}

end

function calculate_derivative_operators(dz)

   num_nodes=length(dz)+1

   # add extra dz on each end for padding nodes used in CenteredDifference
   dz = [dz[1]; dz; dz[end]]

   nth_derivative = 4
   derivative_order = 2
   Azzzz = CenteredDifference(nth_derivative, derivative_order, dz, num_nodes)

   Azzzz = Array(Azzzz)   #concretization of derivative-free operator
   Azzzz = Azzzz[:,2:end-1]  #trim off ghost nodes, not needed

   nth_derivative = 2
   derivative_order = 2
   Azz = CenteredDifference(nth_derivative, derivative_order, dz, num_nodes)
   Azz = Array(Azz)   #concretization of derivative-free operator
   Azz = Azz[:,2:end-1]  #trim off ghost nodes, not needed

   return Azzzz, Azz

end

function calculate_boundary_stencils(bc_flag, h, nth_derivative)

   #Calculate boundary conditions stencils without ghost nodes using
   #Jorge M. Souza, "Boundary Conditions in the Finite Difference Method"
   #https://cimec.org.ar/ojs/index.php/mc/article/download/2662/2607


   taylor_coeffs =  [1; h; h^2/2; h^3/6; h^4/24]

   rhs = zeros(5)
   #row1*[A;B;C;D;E]*u2
   #row2*[A;B;C;D;E]*u2'
   #row3*[A;B;C;D;E]*u2''
   #row4*[A;B;C;D;E]*u2'''
   #row5*[A;B;C;D;E]*u2''''

   if bc_flag == 1 #simply supported end

      lhs = [1 1 1 1 0
         -1 0 1 2 0
         1 0 1 4 2
         -1 0 1 8 -6
         1 0 1 16 12]

      rhs[nth_derivative+1] = (1/TaylorCoeffs[nth_derivative+1])
      boundary_stencil = lhs\rhs
      boundary_stencil = ((boundary_stencil[1:4]),(zeros(4)))  #since u''=0



   elseif bc_flag == 2 #fixed end

      lhs = [1 1 1 1 0
         -1 0 1 2 1
         1 0 1 4 -2
         -1 0 1 8 3
         1 0 1 16 -4]

      rhs[nth_derivative+1] = (1/TaylorCoeffs[nth_derivative+1])
      boundary_stencil = lhs\rhs
      boundary_stencil = ((boundary_stencil[1:4]),(zeros(4)))  #since u'=0

   elseif bc_flag == 3 #free end
                #u'' u'''
      lhs = [1 1 1  0   0
           0 1 2  0   0
           0 1 4  0.5 0
           0 1 8  0   6
           0 1 16 0  0]

      rhs[nth_derivative+1] = (1/TaylorCoeffs[nth_derivative+1])
      boundary_stencil1 = lhs\rhs  #at free end
      boundary_stencil1 = boundary_stencil1[1:3]   #u''=u'''=0

      # use simply supported BC to find stencil at one node in from free end
      lhs = [1 1 1 1 0
         -1 0 1 2 0
         1 0 1 4 2
         -1 0 1 8 -6
         1 0 1 16 12]

      rhs[nth_derivative+1] = (1/TaylorCoeffs[nth_derivative+1])
      boundary_stencil2 = lhs\rhs #at one node over from free end
      boundary_stencil2 = boundary_stencil2[1:4]
      boundary_stencil = ((boundary_stencil1), (boundary_stencil2))  #two stencils are calculated

   end

   return boundary_stencil

end

function apply_end_boundary_conditions(A, end_boundary_conditions, nth_derivative, dz)

   #left end
   h = dz[1]
   bc_flag = end_boundary_conditions[1]
   boundary_stencil = calculate_boundary_stencils(bc_flag,h,nth_derivative)

   A[1,:] .= 0.0
   A[2,:] .= 0.0

   if (bc_flag == 1) | (bc_flag == 2)   #make this cleaner, combine
      A[2,1:length(boundary_stencil[1])] = boundary_stencil[1]
   else
      A[1,1:length(boundary_stencil[1])] = boundary_stencil[1]
      A[2,1:length(boundary_stencil[2])] = boundary_stencil[2]
   end

   #right end
   h = dz[end]
   bc_flag = end_boundary_conditions[2]
   boundary_stencil = calculateboundary_stencils(bc_flag,h,nth_derivative)

   A[end,:] .= 0.0
   A[end-1,:] .= 0.0

   if (bc_flag == 1) | (bc_flag == 2)
      A[end-1,(end-(length(boundary_stencil[1])-1)):end] = reverse(boundary_stencil[1])
   else
      A[end,end-(length(boundary_stencil[1])-1):end] = reverse(boundary_stencil[1])
      A[end-1,end-(length(boundary_stencil[2])-1):end] = reverse(boundary_stencil[2])
   end

   return A

end


function define(member_definitions, section_properties, material_properties, loads, springs, end_boundary_conditions, supports)

   dz, z, dm = Mesh.define_line_element(member_definitions)

   num_nodes=length(dz)+1

   #define property vectors
   A =  Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 1)
   Ix = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 2)
   Iy = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 3)
   J = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 4)
   Cw = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 5)
   xc = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 6)
   yc = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 7)
   xs = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 8)
   ys = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 9)

   xo = -(xc .- xs)
   yo = yc .- ys

   Io = Ix .+ Iy .+ A .* (xo.^2 + yo.^2)

   E = Mesh.create_line_element_property_array(member_definitions, dm, dz, material_properties, 4, 1)
   ν = Mesh.create_line_element_property_array(member_definitions, dm, dz, material_properties, 4, 2)
   G = E./(2 .*(1 .+ ν))


   kx = springs[1]
   ky = springs[2]
   kϕ = springs[3]

   hx = springs[4]
   hy = springs[5]

   Azzzz,Azz = calculate_derivative_operators(dz) #calculate derivative operators

   #apply left and right end boundary condition stencils to derivative operators
   nth_derivative = 4
   Azzzz = apply_end_boundary_conditions(Azzzz,end_boundary_conditions,nth_derivative,dz)

   nth_derivative = 2
   Azz = apply_end_boundary_conditions(Azz,end_boundary_conditions,nth_derivative,dz)

   #define load and load locations
   P =  loads[1]
   qx = loads[2]
   qy = loads[3]
   ax = loads[4]
   ay = loads[5]

   #build identity matrix for ODE operations
   AI = Matrix(1.0I,num_nodes,num_nodes)

   #build operator matrix that doesn't update with load

   #start building operator matrix, LHS first
   A11 = zeros(num_nodes,num_nodes)
   A12 = zeros(num_nodes,num_nodes)
   A13 = zeros(num_nodes,num_nodes)
   A21 = zeros(num_nodes,num_nodes)
   A22 = zeros(num_nodes,num_nodes)
   A23 = zeros(num_nodes,num_nodes)
   A31 = zeros(num_nodes,num_nodes)
   A32 = zeros(num_nodes,num_nodes)
   A33 = zeros(num_nodes,num_nodes)

   #calculate operator quantities on LHS  AU=B
   for i = 1:num_nodes
      A11[i,:] = E.*Iy.*Azzzz[i,:] .+ P .* Azz[i,:] .+kx.*AI[i,:]
      A13[i,:] = -kx.*hy.*AI[i,:]
      A22[i,:] = E.*Ix.*Azzzz[i,:] .+ P.*Azz[i,:] .+ ky.*AI[i,:]
      A23[i,:] = ky.*hx.*AI[i,:]
      A31[i,:] = -kx.*hy.*AI[i,:]
      A32[i,:] = ky.*hx.*AI[i,:]
      A33[i,:] = E.*Cw.*Azzzz[i,:] .-(G.*J .- (P.*Io./A)).*Azz[i,:] .+kx.*hy.^2 .*AI[i,:] .+ky.*hx.^2 .*AI[i,:] .+kϕ.*AI[i,:] .+qx.*ax.*AI[i,:] - .+qy.*ay.*AI[i,:]
   end

   #calculate RHS of AU=B
   B1 = qx
   B2 = qy
   B3 = qx.*ay .+qy.*ax

   #reduce problem to free dof
   fixed_dof = [findall(x->abs(x-supports[i])<=10e-6,z) for i=1:length(supports)]
   fixed_dof = VectorOfArray(fixed_dof)
   fixed_dof  = convert(Array,fixed_dof)
   free_dof = setdiff(1:num_nodes,fixed_dof)

   K = [A11[free_dof,free_dof] A12[free_dof,free_dof] A13[free_dof,free_dof];
      A21[free_dof,free_dof] A22[free_dof,free_dof] A23[free_dof,free_dof];
      A31[free_dof,free_dof] A32[free_dof,free_dof] A33[free_dof,free_dof]]

   F = [B1[free_dof]; B2[free_dof]; B3[free_dof]]

   properties = NamedTuple{(:dm, :z, :P, :A, :Ix, :Iy, :J, :Cw, :xc, :yc, :xs, :ys, :xo, :yo, :Io, :E, :ν, :G, :ax, :ay, :kx, :ky, :kϕ, :hx, :hy, :qx, :qy)}((dm, z, P, A, Ix, Iy, J, Cw, xc, yc, xs, ys, xo, yo, Io, E, ν, G, ax, ay, kx, ky, kϕ, hx, hy, qx, qy))

   model = Model(member_definitions, section_properties, material_properties, spring_stiffness, spring_location, supports, loads, load_location, end_boundary_conditions, qx, qy, P, K, F, free_dof, A, Ix, Iy, Io, J, Cw, xc, yc, xs, ys, xo, yo, E, ν, G, ax, ay, kx, ky, kϕ, hx, hy, z, dm, u, v, ϕ)

   return model

end

function residual!(R, U, K, F)

   for i=1:length(F)
    R[i] = transpose(K[i,:]) * (U) - F[i]
   end

   return R

end


function solve(model)

   # K, F, free_dof, properties = define(member_definitions, section_properties, material_properties, loads, springs, end_boundary_conditions, supports)

   num_nodes=length(model.z)

   model.u=zeros(num_nodes)
   model.v=zeros(num_nodes)
   model.ϕ=zeros(num_nodes)

   deformation_guess = model.K \ model.F    #consider revising this for large systems, it might be slow...

   solution = nlsolve((R, U) ->residual!(R, U, model.K, model.F), deformation_guess)

   model.u[model.free_dof] = solution.zero[1:length(model.free_dof)]
   model.v[model.free_dof] = solution.zero[length(model.free_dof)+1:2*length(model.free_dof)]
   model.ϕ[model.free_dof] = solution.zero[2*length(model.free_dof)+1:3*length(model.free_dof)]

   return model

end

end #module
