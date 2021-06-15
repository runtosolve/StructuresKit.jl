module ThinWalledBeam

using DiffEqOperators: CenteredDifference
using LinearAlgebra
using NLsolve

using ..Mesh

export user_interface, define, solve, Model


Base.@kwdef mutable struct Model   #using Parameters.jl macro here to help assign some constants, leave others for later.

   Ix::Union{Array{Float64}, Nothing} = nothing
   Iy::Union{Array{Float64}, Nothing} = nothing
   Ixy::Union{Array{Float64}, Nothing} = nothing
   J::Union{Array{Float64}, Nothing} = nothing
   Cw::Union{Array{Float64}, Nothing} = nothing
   E::Union{Array{Float64}, Nothing} = nothing
   ν::Union{Array{Float64}, Nothing} = nothing
   G::Union{Array{Float64}, Nothing} = nothing
   ax::Union{Array{Float64}, Nothing} = nothing
   ay::Union{Array{Float64}, Nothing} = nothing
   ay_kx::Union{Array{Float64}, Nothing} = nothing
   kx::Union{Array{Float64}, Nothing} = nothing
   kϕ::Union{Array{Float64}, Nothing} = nothing

   qx::Union{Array{Float64}, Nothing} = nothing
   qy::Union{Array{Float64}, Nothing} = nothing

   end_boundary_conditions::Union{Array{Int64}, Nothing} = nothing
   supports::Union{Array{Float64}, Nothing} = nothing

   z::Union{Array{Float64}, Nothing} = nothing
   dz::Union{Array{Float64}, Nothing} = nothing
   dm::Union{Array{Int64}, Nothing} = nothing

   K::Union{Matrix{Float64}, Nothing} = nothing  
   F::Union{Array{Float64}, Nothing} = nothing

   free_dof::Union{Array{Int64}, Nothing} = nothing
   
   u::Union{Array{Float64}, Nothing} = nothing
   v::Union{Array{Float64}, Nothing} = nothing
   ϕ::Union{Array{Float64}, Nothing} = nothing

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

function user_interface(member_definitions, section_properties, material_properties, spring_stiffness, spring_location, loads, load_locations, end_boundary_conditions, supports)

   #Define discretization along beam and assign a member type to each segment, i.e., dm.
   dz, z, dm = Mesh.define_line_element(member_definitions)

   #Define the number of nodes along the beam.
   num_nodes=length(dz)+1

   #This is the user interface mapping.
   #L(1) dL(2) section_properties(3) material_properties(4) spring_stiffness(5) spring_location(6) load(7) load_location(8) 

   #Define properties at each node.
   Ix = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 1)
   Iy = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 2)
   Ixy = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 3)
   J = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 4)
   Cw = Mesh.create_line_element_property_array(member_definitions, dm, dz, section_properties, 3, 5)

   E = Mesh.create_line_element_property_array(member_definitions, dm, dz, material_properties, 4, 1)
   ν = Mesh.create_line_element_property_array(member_definitions, dm, dz, material_properties, 4, 2)
   G = E./(2 .*(1 .+ ν))

   kx = Mesh.create_line_element_property_array(member_definitions, dm, dz, spring_stiffness, 5, 1)
   kϕ = Mesh.create_line_element_property_array(member_definitions, dm, dz, spring_stiffness, 5, 2)

   ay_kx = Mesh.create_line_element_property_array(member_definitions, dm, dz, spring_location, 6, 1)

   #Assign a uniform load magnitude at each node.  
   qx = Mesh.create_line_element_property_array(member_definitions, dm, dz, loads, 7, 1)
   qy = Mesh.create_line_element_property_array(member_definitions, dm, dz, loads, 7, 2)

   #Assign load location on the cross-section at each node.
   ax = Mesh.create_line_element_property_array(member_definitions, dm, dz, load_locations, 8, 1)
   ay = Mesh.create_line_element_property_array(member_definitions, dm, dz, load_locations, 8, 2)

   #Store model inputs in data structure.
   model = Model(Ix, Iy, Ixy, J, Cw, E, ν, G, ax, ay, ay_kx, kx, kϕ, qx, qy, end_boundary_conditions, supports, z, dz, dm, nothing, nothing, nothing, nothing, nothing, nothing)
  
   return model

end




function define(model)

   #Define the number of nodes.
   num_nodes = length(model.z)

   #Calculate the derivative operators. 
   Azzzz,Azz = calculate_derivative_operators(model.dz) 

   #Apply left and right end boundary condition stencils to derivative operators.
   nth_derivative = 4
   Azzzz = apply_end_boundary_conditions(Azzzz, model.end_boundary_conditions, nth_derivative, model.dz)

   nth_derivative = 2
   Azz = apply_end_boundary_conditions(Azz, model.end_boundary_conditions, nth_derivative, model.dz)

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
      A11[i,:] = model.E .* model.Iy .* Azzzz[i,:] .+ model.kx .* AI[i,:]
      A12[i,:] = model.E .* model.Ixy .* Azzzz[i,:]
      A13[i,:] = model.kx .* model.ay_kx .*AI[i,:]
      A21[i,:] = model.E .* model.Ixy .* Azzzz[i,:]
      A22[i,:] = model.E .* model.Ix .* Azzzz[i,:]
      A31[i,:] = model.kx .* model.ay_kx .* AI[i,:]
      A33[i,:] = model.E .* model.Cw .* Azzzz[i,:] .- model.G .* model.J .* Azz[i,:] .+ model.kx .* model.ay_kx .* model.ay_kx .* AI[i,:] .+ model.kϕ .* AI[i,:] .+ model.qx .* model.ax .* AI[i,:] .- model.qy .* model.ay .* AI[i,:]
   end

   #Calculate RHS of AU=B.
   B1 = model.qx
   B2 = model.qy
   B3 = model.qx .* model.ay .+ model.qy .* model.ax

   #Reduce problem to free dof.

   #Find fixed dof.
   fixed_dof = []
   for i = 1:length(model.supports)

      if i == 1

            fixed_dof = findall(x->abs(x-model.supports[i])<=10e-6, model.z)

      elseif i != 1

            new_fixed_dof = findall(x->abs(x-model.supports[i])<=10e-6, model.z)
            fixed_dof = [fixed_dof; new_fixed_dof]

      end

   end

   #Find free dof. 
   model.free_dof = setdiff(1:num_nodes,fixed_dof)

   #Define stiffness matrix.
   model.K = [A11[model.free_dof,model.free_dof] A12[model.free_dof,model.free_dof] A13[model.free_dof,model.free_dof];
      A21[model.free_dof,model.free_dof] A22[model.free_dof,model.free_dof] A23[model.free_dof,model.free_dof];
      A31[model.free_dof,model.free_dof] A32[model.free_dof,model.free_dof] A33[model.free_dof,model.free_dof]]

   #Define external force vector.
   model.F = [B1[model.free_dof]; B2[model.free_dof]; B3[model.free_dof]]

   return model

end

function residual!(R, U, K, F)

   for i=1:length(F)
    R[i] = transpose(K[i,:])*U-F[i]
   end

   return R

end


function solve(model)

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
