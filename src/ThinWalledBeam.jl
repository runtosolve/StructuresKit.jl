module ThinWalledBeam

using DiffEqOperators: CenteredDifference
using LinearAlgebra
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
   loads::Vector{Vector{Float64}}
   load_location::Vector{Vector{Float64}}
   end_boundary_conditions::Array{Int64}

   qx::Array{Float64}
   qy::Array{Float64}

   K::Matrix{Float64}
   F::Array{Float64}

   free_dof::Array{Int64}
   
   Ix::Array{Float64}
   Iy::Array{Float64}
   Ixy::Array{Float64}
   J::Array{Float64}
   Cw::Array{Float64}
   E::Array{Float64}
   ν::Array{Float64}
   G::Array{Float64}
   ax::Array{Float64}
   ay::Array{Float64}
   ay_kx::Array{Float64}
   kx::Array{Float64}
   kϕ::Array{Float64}

   z::Array{Float64}
   dz::Array{Float64}
   dm::Array{Int64}
   
   u::Array{Float64}
   v::Array{Float64}
   ϕ::Array{Float64}

   Model() = new()

end


function calculate_derivative_operators(dz)

   #Define the number of nodes.
   num_nodes = length(dz) + 1

   #Add extra dz on each end for padding nodes used in the CenteredDifference function.
   dz = [dz[1]; dz; dz[end]]

   #Calculate the 4th derivative operator.
   nth_derivative = 4
   derivative_order = 2
   Azzzz = CenteredDifference(nth_derivative, derivative_order, dz, num_nodes)

   #Convert the operator to a matrix.
   Azzzz = Array(Azzzz)   

   #Trim off the ghost nodes.
   Azzzz = Azzzz[:,2:end-1]  

   #Calculate the 2nd derivative operator.
   nth_derivative = 2
   derivative_order = 2
   Azz = CenteredDifference(nth_derivative, derivative_order, dz, num_nodes)
   Azz = Array(Azz)   
   Azz = Azz[:,2:end-1]  

   return Azzzz, Azz

end

function calculate_boundary_stencils(bc_flag, h, nth_derivative)

   #Calculate boundary conditions stencils without ghost nodes using
   #Jorge M. Souza, "Boundary Conditions in the Finite Difference Method"
   #https://cimec.org.ar/ojs/index.php/mc/article/download/2662/2607


   taylor_coeffs =  [1; h; h^2/2; h^3/6; h^4/24]

   RHS = zeros(Float64, 5)
   #row1*[A;B;C;D;E]*u2
   #row2*[A;B;C;D;E]*u2'
   #row3*[A;B;C;D;E]*u2''
   #row4*[A;B;C;D;E]*u2'''
   #row5*[A;B;C;D;E]*u2''''

   if bc_flag == 1 #simply supported end

      LHS = [1 1 1 1 0
         -1 0 1 2 0
         1 0 1 4 2
         -1 0 1 8 -6
         1 0 1 16 12]

      RHS[nth_derivative + 1] = (1/taylor_coeffs[nth_derivative + 1])
      boundary_stencil = LHS \ RHS
      boundary_stencil = ((boundary_stencil[1:4]),(zeros(4)))  #since u''=0

   elseif bc_flag == 2 #fixed end

      LHS = [1 1 1 1 0
         -1 0 1 2 1
         1 0 1 4 -2
         -1 0 1 8 3
         1 0 1 16 -4]

      RHS[nth_derivative+1] = (1/taylor_coeffs[nth_derivative + 1])
      boundary_stencil = LHS\RHS
      boundary_stencil = ((boundary_stencil[1:4]),(zeros(4)))  #since u'=0

   elseif bc_flag == 3 #free end
                #u'' u'''
      LHS = [1 1 1  0   0
           0 1 2  0   0
           0 1 4  0.5 0
           0 1 8  0   6
           0 1 16 0  0]

      RHS[nth_derivative+1] = (1/taylor_coeffs[nth_derivative+1])
      boundary_stencil1 = LHS\RHS  #at free end
      boundary_stencil1 = boundary_stencil1[1:3]   #u''=u'''=0

      # use simply supported BC to find stencil at one node in from free end
      LHS = [1 1 1 1 0
         -1 0 1 2 0
         1 0 1 4 2
         -1 0 1 8 -6
         1 0 1 16 12]

      RHS[nth_derivative+1] = (1/taylor_coeffs[nth_derivative+1])
      boundary_stencil2 = LHS\RHS #at one node over from free end
      boundary_stencil2 = boundary_stencil2[1:4]
      boundary_stencil = ((boundary_stencil1), (boundary_stencil2))  #two stencils are calculated

   end

   return boundary_stencil

end

function apply_end_boundary_conditions(A, end_boundary_conditions, nth_derivative, dz)

   #Consider the left end boundary condition.  Swap out interior stencil for boundary stencil.
   h = dz[1]
   bc_flag = end_boundary_conditions[1]
   boundary_stencil = calculate_boundary_stencils(bc_flag, h, nth_derivative)

   A[1,:] .= 0.0
   A[2,:] .= 0.0

   if (bc_flag == 1) | (bc_flag == 2)   #make this cleaner, combine
      A[2,1:length(boundary_stencil[1])] = boundary_stencil[1]
   else
      A[1,1:length(boundary_stencil[1])] = boundary_stencil[1]
      A[2,1:length(boundary_stencil[2])] = boundary_stencil[2]
   end

   #Consider the right end boundary condition.
   h = dz[end]
   bc_flag = end_boundary_conditions[2]
   boundary_stencil = calculate_boundary_stencils(bc_flag,h,nth_derivative)

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


function define(member_definitions, section_properties, material_properties, spring_stiffness, spring_location, supports, loads, load_location, end_boundary_conditions)

   #Define discretization along beam and assign a member type to each segment, i.e., dm.
   dz, z, dm = Mesh.define_line_element(member_definitions)

   #Define the number of nodes along the beam.
   num_nodes=length(dz)+1

   #Define properties at each node.
   Ix = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 1)
   Iy = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 2)
   Ixy = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 3)
   J = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 4)
   Cw = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 5)

   E = Mesh.create_line_element_property_array(member_definitions, dm, dz, material_properties, 4, 1)
   ν = Mesh.create_line_element_property_array(member_definitions, dm, dz, material_properties, 4, 2)
   G = E./(2 .*(1 .+ ν))

   # kx = Mesh.create_line_element_property_array(member_definitions, dm, dz, spring_stiffness, 6, 1)
   # kϕ = Mesh.create_line_element_property_array(member_definitions, dm, dz, spring_stiffness, 6, 2)
   # ay_kx = Mesh.create_line_element_property_array(member_definitions, dm, dz, spring_location, 7, 1)

   kx = spring_stiffness[1]
   kϕ = spring_stiffness[2]

   ay_kx = spring_location

   #Calculate the derivative operators. 
   Azzzz,Azz = calculate_derivative_operators(dz) 

   #Apply left and right end boundary condition stencils to derivative operators.
   nth_derivative = 4
   Azzzz = apply_end_boundary_conditions(Azzzz,end_boundary_conditions, nth_derivative, dz)

   nth_derivative = 2
   Azz = apply_end_boundary_conditions(Azz, end_boundary_conditions, nth_derivative, dz)

   #Assign a uniform load magnitude at each node.  
   qx = loads[1]
   qy = loads[2]

   #Assign load location on the cross-section at each node.
   ax = load_location[1]
   ay = load_location[2]

   #Build identity matrix for ODE operations.
   AI = Matrix(1.0I, num_nodes, num_nodes)

   #Build operator matrix that doesn't update with load.

   #LHS first.
   A11 = zeros(Float64, num_nodes,num_nodes)
   A12 = zeros(Float64, num_nodes,num_nodes)
   A13 = zeros(Float64, num_nodes,num_nodes)
   A21 = zeros(Float64, num_nodes,num_nodes)
   A22 = zeros(Float64, num_nodes,num_nodes)
   A23 = zeros(Float64, num_nodes,num_nodes)
   A31 = zeros(Float64, num_nodes,num_nodes)
   A32 = zeros(Float64, num_nodes,num_nodes)
   A33 = zeros(Float64, num_nodes,num_nodes)

   #Calculate operator quantities on LHS  AU=B.
   for i = 1:num_nodes
      A11[i,:] = E.*Iy.*Azzzz[i,:] .+kx.*AI[i,:]
      A12[i,:] = E.*Ixy.*Azzzz[i,:]
      A13[i,:] = kx.*ay_kx.*AI[i,:]
      A21[i,:] = E.*Ixy.*Azzzz[i,:]
      A22[i,:] = E.*Ix.*Azzzz[i,:]
      A31[i,:] = kx.*ay_kx.*AI[i,:]
      A33[i,:] = E.*Cw.*Azzzz[i,:] .-G.*J.*Azz[i,:] .+kx.*ay_kx.*ay_kx.*AI[i,:] .+kϕ.*AI[i,:] .+qx.*ax.*AI[i,:] .-qy.*ay.*AI[i,:]
   end

   #Calculate RHS of AU=B.
   B1 = qx
   B2 = qy
   B3 = qx.*ay .+ qy.*ax

   #Reduce problem to free dof.

   #Find fixed dof.
   fixed_dof = []
   for i = 1:length(supports)

      if i == 1

            fixed_dof = findall(x->abs(x-supports[i])<=10e-6, z)

      elseif i != 1

            new_fixed_dof = findall(x->abs(x-supports[i])<=10e-6, z)
            fixed_dof = [fixed_dof; new_fixed_dof]

      end

   end

   #Find free dof. 
   free_dof = setdiff(1:num_nodes,fixed_dof)

   #Define stiffness matrix.
   K = [A11[free_dof,free_dof] A12[free_dof,free_dof] A13[free_dof,free_dof];
      A21[free_dof,free_dof] A22[free_dof,free_dof] A23[free_dof,free_dof];
      A31[free_dof,free_dof] A32[free_dof,free_dof] A33[free_dof,free_dof]]

   #Define external force vector.
   F = [B1[free_dof]; B2[free_dof]; B3[free_dof]]

   #Add definitions to data structure.
   model = Model(member_definitions, section_properties, material_properties, spring_stiffness, spring_location, supports, loads, load_location, end_boundary_conditions, qx, qy, K,  F, free_dof, Ix, Iy, Ixy, J, Cw, E, ν, G, ax, ay, ay_kx, kx, kϕ, z, dm, [], [], [])
  
   return model

end

function residual!(R, U, K, F)

   for i=1:length(F)
    R[i] = transpose(K[i,:])*U-F[i]
   end

   return R

end


function solve(model)

   # #Set up the beam problem, including the stiffness matrix and external force vector.
   # K, F, free_dof, Ix, Iy, Ixy, J, Cw, E, ν, G, ax, ay, ay_kx, kx, kϕ, z, dm = define(member_definitions, section_properties, material_properties, spring_stiffness, spring_location, supports, load, load_location, end_boundary_conditions)

   #Define the number of nodes along the beam.
   num_nodes = length(model.z)

   #Define the deformation vectors.
   u = zeros(num_nodes)
   v = zeros(num_nodes)
   ϕ = zeros(num_nodes)

   #Define the deformation initial guess for the nonlinear solver.
   u_guess = model.K \ model.F

   #Solve for the beam deformations.
   solution = nlsolve((R,U) ->residual!(R, U, model.K, model.F), u_guess)

   #Pull the displacements and twist from the solution results.
   u[model.free_dof] = solution.zero[1:length(model.free_dof)]
   v[model.free_dof] = solution.zero[length(model.free_dof)+1:2*length(model.free_dof)]
   ϕ[model.free_dof] = solution.zero[2*length(model.free_dof)+1:3*length(model.free_dof)]

   #Add deformations to data structure.
   model.u = u
   model.v = v
   model.ϕ = ϕ

   return model

end

end #module
