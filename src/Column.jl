module Column

using DiffEqOperators: CenteredDifference
using LinearAlgebra
using RecursiveArrayTools: VectorOfArray
using NLsolve

using ..Mesh
using ..InternalForces



export Properties, initialize, solve, normal_stresses, normal_strains, warping_displacements, deformed_shape, calculate_derivative_operators, apply_end_boundary_conditions



function calculate_derivative_operators(dz)

   NumberOfNodes=length(dz)+1

   # add extra dz on each end for padding nodes used in CenteredDifference
   dz = [dz[1]; dz; dz[end]]

   NthDerivative = 4
   DerivativeOrder = 2
   Azzzz = CenteredDifference(NthDerivative, DerivativeOrder, dz, NumberOfNodes)

   Azzzz = Array(Azzzz)   #concretization of derivative-free operator
   Azzzz = Azzzz[:,2:end-1]  #trim off ghost nodes, not needed

   NthDerivative = 2
   DerivativeOrder = 2
   Azz = CenteredDifference(NthDerivative, DerivativeOrder, dz, NumberOfNodes)
   Azz = Array(Azz)   #concretization of derivative-free operator
   Azz = Azz[:,2:end-1]  #trim off ghost nodes, not needed

   return Azzzz, Azz

end

function calculateBoundaryStencils(BCFlag, h, NthDerivative)

   #Calculate boundary conditions stencils without ghost nodes using
   #Jorge M. Souza, "Boundary Conditions in the Finite Difference Method"
   #https://cimec.org.ar/ojs/index.php/mc/article/download/2662/2607


   TaylorCoeffs =  [1; h; h^2/2; h^3/6; h^4/24]

   RHS = zeros(5)
   #row1*[A;B;C;D;E]*u2
   #row2*[A;B;C;D;E]*u2'
   #row3*[A;B;C;D;E]*u2''
   #row4*[A;B;C;D;E]*u2'''
   #row5*[A;B;C;D;E]*u2''''

   if BCFlag == 1 #simply supported end

      LHS = [1 1 1 1 0
         -1 0 1 2 0
         1 0 1 4 2
         -1 0 1 8 -6
         1 0 1 16 12]

      RHS[NthDerivative+1] = (1/TaylorCoeffs[NthDerivative+1])
      BoundaryStencil = LHS\RHS
      BoundaryStencil = ((BoundaryStencil[1:4]),(zeros(4)))  #since u''=0



   elseif BCFlag == 2 #fixed end

      LHS = [1 1 1 1 0
         -1 0 1 2 1
         1 0 1 4 -2
         -1 0 1 8 3
         1 0 1 16 -4]

      RHS[NthDerivative+1] = (1/TaylorCoeffs[NthDerivative+1])
      BoundaryStencil = LHS\RHS
      BoundaryStencil = ((BoundaryStencil[1:4]),(zeros(4)))  #since u'=0

   elseif BCFlag == 3 #free end
                #u'' u'''
      LHS = [1 1 1  0   0
           0 1 2  0   0
           0 1 4  0.5 0
           0 1 8  0   6
           0 1 16 0  0]

      RHS[NthDerivative+1] = (1/TaylorCoeffs[NthDerivative+1])
      BoundaryStencil1 = LHS\RHS  #at free end
      BoundaryStencil1 = BoundaryStencil1[1:3]   #u''=u'''=0

      # use simply supported BC to find stencil at one node in from free end
      LHS = [1 1 1 1 0
         -1 0 1 2 0
         1 0 1 4 2
         -1 0 1 8 -6
         1 0 1 16 12]

      RHS[NthDerivative+1] = (1/TaylorCoeffs[NthDerivative+1])
      BoundaryStencil2 = LHS\RHS #at one node over from free end
      BoundaryStencil2 = BoundaryStencil2[1:4]
      BoundaryStencil = ((BoundaryStencil1), (BoundaryStencil2))  #two stencils are calculated

   end

   return BoundaryStencil

end

function apply_end_boundary_conditions(A, EndBoundaryConditions, NthDerivative, dz)

   #left end
   h = dz[1]
   BCFlag = EndBoundaryConditions[1]
   BoundaryStencil = calculateBoundaryStencils(BCFlag,h,NthDerivative)

   A[1,:] .= 0.0
   A[2,:] .= 0.0

   if (BCFlag == 1) | (BCFlag == 2)   #make this cleaner, combine
      A[2,1:length(BoundaryStencil[1])] = BoundaryStencil[1]
   else
      A[1,1:length(BoundaryStencil[1])] = BoundaryStencil[1]
      A[2,1:length(BoundaryStencil[2])] = BoundaryStencil[2]
   end

   #right end
   h = dz[end]
   BCFlag = EndBoundaryConditions[2]
   BoundaryStencil = calculateBoundaryStencils(BCFlag,h,NthDerivative)

   A[end,:] .= 0.0
   A[end-1,:] .= 0.0

   if (BCFlag == 1) | (BCFlag == 2)
      A[end-1,(end-(length(BoundaryStencil[1])-1)):end] = reverse(BoundaryStencil[1])
   else
      A[end,end-(length(BoundaryStencil[1])-1):end] = reverse(BoundaryStencil[1])
      A[end-1,end-(length(BoundaryStencil[2])-1):end] = reverse(BoundaryStencil[2])
   end

   return A

end

mutable struct Properties

      A::Array{Float64,1}
      Ix::Array{Float64,1}
      Iy::Array{Float64,1}
      J::Array{Float64,1}
      Cw::Array{Float64,1}
      xc::Array{Float64,1}
      yc::Array{Float64,1}
      xs::Array{Float64,1}
      ys::Array{Float64,1}
      xo::Array{Float64,1}
      yo::Array{Float64,1}
      Io::Array{Float64,1}
      E::Array{Float64,1}
      ν::Array{Float64,1}
      G::Array{Float64,1}
      kx::Float64
      ky::Float64
      kϕ::Float64
      hx::Float64
      hy::Float64
      uo::Array{Float64,1}
      vo::Array{Float64,1}
      ϕo::Array{Float64,1}
      P::Array{Float64,1}
      dz::Array{Float64,1}
      z::Array{Float64,1}
      dm::Array{Float64,1}
      Azz::Array{Float64,2}
      Azzzz::Array{Float64,2}
      supports::Array{Float64,1}
      xcoords::Array{Float64,1}
      ycoords::Array{Float64,1}
      Ai::Array{Float64,1}
      cxi::Array{Float64,1}
      cyi::Array{Float64,1}
      
end

   

function initialize(member_definitions, section_properties, material_properties, loads, springs, end_boundary_conditions, supports, imperfections)

   dz, z, dm = Mesh.define_line_element(member_definitions)

   NumberOfNodes=length(dz)+1


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

   #define imperfections 
   uo = imperfections[1]
   vo = imperfections[2]
   ϕo = imperfections[3] 

   #define load
   P =  loads

   #calculate derivative operators
   Azzzz,Azz = calculate_derivative_operators(dz) 

   #apply left and right end boundary condition stencils to derivative operators
   NthDerivative = 4
   Azzzz = apply_end_boundary_conditions(Azzzz, end_boundary_conditions, NthDerivative, dz)
    
   NthDerivative = 2
   Azz = apply_end_boundary_conditions(Azz, end_boundary_conditions, NthDerivative, dz)

   #initialize cross-section geometry and triangulation, add them later if needed
   xcoords = []
   ycoords = []
   Ai = []
   cxi = []
   cyi = []

   properties = Properties(A, Ix, Iy, J, Cw, xc, yc, xs, ys, xo, yo, Io, E, ν, G, kx, ky, kϕ, hx, hy, uo, vo, ϕo, P, dz, z, dm, Azz, Azzzz, supports, 
   xcoords, ycoords, Ai, cxi, cyi)

   return properties

end



function equations(properties)


   NumberOfNodes=length(properties.dz)+1

   #build identity matrix for ODE operations
   AI = Matrix(1.0I,NumberOfNodes,NumberOfNodes)

   #build operator matrix that doesn't update with load

   #start building operator matrix, LHS first
   A11 = zeros(NumberOfNodes,NumberOfNodes)
   A12 = zeros(NumberOfNodes,NumberOfNodes)
   A13 = zeros(NumberOfNodes,NumberOfNodes)
   A21 = zeros(NumberOfNodes,NumberOfNodes)
   A22 = zeros(NumberOfNodes,NumberOfNodes)
   A23 = zeros(NumberOfNodes,NumberOfNodes)
   A31 = zeros(NumberOfNodes,NumberOfNodes)
   A32 = zeros(NumberOfNodes,NumberOfNodes)
   A33 = zeros(NumberOfNodes,NumberOfNodes)

   #calculate operator quantities on LHS  AU=B
   for i = 1:NumberOfNodes
      A11[i,:] = properties.E .* properties.Iy .* properties.Azzzz[i,:] .+ properties.P .* properties.Azz[i,:] .+ properties.kx .* AI[i,:]
      A13[i,:] = properties.kx .*(properties.yo .- properties.hy) .* AI[i,:] .+ properties.P .* properties.yo .* properties.Azz[i,:]
      A22[i,:] = properties.E .* properties.Ix .* properties.Azzzz[i,:] .+ properties.P .* properties.Azz[i,:] .+ properties.ky .* AI[i,:]
      A23[i,:] = -properties.ky .* (properties.xo .-properties.hx) .* AI[i,:] - properties.P .* properties.xo .* properties.Azz[i,:]
      A31[i,:] = properties.kx.*(properties.yo .- properties.hy) .* AI[i,:] .+ properties.P .* properties.yo .* properties.Azz[i,:]
      A32[i,:] = -properties.ky.*(properties.xo .- properties.hx) .* AI[i,:] .- properties.P .* properties.xo .* properties.Azz[i,:]
      A33[i,:] = properties.E .* properties.Cw .* properties.Azzzz[i,:] .-(properties.G .* properties.J .- (properties.P .* properties.Io ./ properties.A)) .* properties.Azz[i,:] 
      .+ properties.kx .* (properties.yo .- properties.hy).^2 .*AI[i,:] .+ properties.ky.*(properties.xo .- properties.hx).^2 .*AI[i,:] .+ properties.kϕ .* AI[i,:]
   end

   #calculate RHS of AU=B
   B1 = -properties.P .* properties.Azz * properties.uo - properties.P .* properties.yo .* properties.Azz * properties.ϕo
   B2 = -properties.P .* properties.Azz * properties.vo + properties.P .* properties.xo .* properties.Azz * properties.ϕo
   B3 = -properties.P .* properties.yo .* properties.Azz * properties.uo + properties.P .* properties.xo .* properties.Azz * properties.vo - properties.P .* (properties.Io ./ properties.A) .* properties.Azz * properties.ϕo

   #reduce problem to free dof
   FixedDOF = [findall(x->abs(x-properties.supports[i])<=10e-6, properties.z) for i=1:length(properties.supports)]
   FixedDOF = VectorOfArray(FixedDOF)
   FixedDOF  = convert(Array,FixedDOF)
   FreeDOF = setdiff(1:NumberOfNodes,FixedDOF)

   Am = [A11[FreeDOF,FreeDOF] A12[FreeDOF,FreeDOF] A13[FreeDOF,FreeDOF];
      A21[FreeDOF,FreeDOF] A22[FreeDOF,FreeDOF] A23[FreeDOF,FreeDOF];
      A31[FreeDOF,FreeDOF] A32[FreeDOF,FreeDOF] A33[FreeDOF,FreeDOF]]

   Bm = [B1[FreeDOF]; B2[FreeDOF]; B3[FreeDOF]]

   return Am, Bm, FreeDOF

end

function residual!(R, U, K, F, properties, free_dof, inelasticity_flag)

   if inelasticity_flag == 1

      K, F = inelasticity_update(U, properties, free_dof)

   end

   for i=1:length(F)

      R[i] = transpose(K[i,:]) * (U) - F[i]
   
   end


   return R

end


# function inelasticity_update(U, properties, free_dof)


#    num_sections = length(properties.z)

#    for i=1:num_sections

#       Myy = InternalForces.moment(properties.z, properties.dm, -U, properties.E, properties.Iy)
#       ϵ_axial, ϵ_flexure, ϵ_total = calculate_strain(P, Myy[i], A[i], Iyy[i], cx, cxi, Es)

#       if abs(ϵ_total) < ϵ_limit


#    end

   #calculate moment from deformation
   #get cross-section cell discretization
   #calculate strains from P+M along the column
   #assign strain to each cell
   #reduce area of cell based on tangent elastic modulus,  Acell * Etangent/Ei  
   #calculate updated I




function solve(properties, inelasticity_flag)


   K, F, free_dof = equations(properties)

   num_nodes=length(properties.z)

   u=zeros(num_nodes)
   v=zeros(num_nodes)
   ϕ=zeros(num_nodes)

   deformation_guess = K \ F    #consider revising this for large systems, it might be slow...

   solution = nlsolve((R, U) ->residual!(R, U, K, F, properties, free_dof, inelasticity_flag), deformation_guess)

   u[free_dof] = solution.zero[1:length(free_dof)]
   v[free_dof] = solution.zero[length(free_dof)+1:2*length(free_dof)]
   ϕ[free_dof] = solution.zero[2*length(free_dof)+1:3*length(free_dof)]

   return u, v, ϕ

end



function deformed_shape(xcoords, ycoords, zcoords, u, v, ϕ, w, scale_u, scale_v, scale_ϕ, scale_w)

   x_deformed = []
   y_deformed = []
   z_deformed = []
   num_sections = length(zcoords)

   u = scale_u .* u
   v = scale_v .* v
   ϕ = scale_ϕ .* ϕ
   w = scale_w .* w

   for i = 1:num_sections

       if i == 1
           x_deformed = xcoords .+ (u[i] .+xcoords.*cos.(ϕ[i]).-ycoords.*sin.(ϕ[i]))
           y_deformed = ycoords .+ (v[i] .+xcoords.*sin.(ϕ[i]).+ycoords.*cos.(ϕ[i]))
           z_deformed = zcoords[i] .+ w[i, :]
       else
           x_deformed = [x_deformed; xcoords .+ (u[i] .+xcoords.*cos.(ϕ[i]).-ycoords.*sin.(ϕ[i]))]
           y_deformed = [y_deformed; ycoords .+ (v[i] .+xcoords.*sin.(ϕ[i]).+ycoords.*cos.(ϕ[i]))]
           z_deformed = [z_deformed; zcoords[i] .+ w[i, :]]
       end

   end

   return x_deformed, y_deformed, z_deformed

end


function axial_stress(P, A)

   σ_axial = P / A

end

function flexural_stress(M, y, I)

   σ_flexural = M * y / I

end


function normal_stresses(xcoords, ycoords, zcoords, u, v, properties)

   num_sections = length(zcoords)
   num_cross_section_nodes = length(xcoords)

   σ_axial = zeros(Float64, (num_sections, num_cross_section_nodes))
   σ_flexural_xx = zeros(Float64, (num_sections, num_cross_section_nodes))
   σ_flexural_yy = zeros(Float64, (num_sections, num_cross_section_nodes))
   σ_total_normal = zeros(Float64, (num_sections, num_cross_section_nodes))

   Mxx = InternalForces.moment(properties.z, properties.dm, -v, properties.E, properties.Ix)
   Myy = InternalForces.moment(properties.z, properties.dm, -u, properties.E, properties.Iy)

   for i=1:num_sections

       A=properties.A[i]
       P=properties.P[i]
       Ix = properties.Ix[i]
       Iy = properties.Iy[i]

       #need a better way to include axial stresses in the deformed configuration  
       # σ_axial[i, :] = axial_stress(-properties.P[i], properties.A[i]) * ones(Float64, num_cross_section_nodes)
       σ_flexural_xx[i, :] = flexural_stress.(Mxx[i], -(ycoords .- properties.yc[i]), properties.Ix[i]) 
       σ_flexural_yy[i, :] = flexural_stress.(Myy[i], -(xcoords .- properties.xc[i]), properties.Iy[i]) 
       #add warping stresses here one day...
       σ_total_normal[i,:] = σ_axial[i, :] + σ_flexural_xx[i, :] + σ_flexural_yy[i, :]

   end

   return σ_total_normal

end


function normal_stains(properties, σ_total_normal)
   
   ϵ_total_normal = σ_total_normal ./ properties.E

end

using NumericalIntegration

#integrate strains to calculate warping displacements

function warping_displacements(zcoords, ϵ_total_normal)

   num_sections = length(zcoords)
   num_cross_section_nodes = length(ϵ_total_normal[1,:])

   w = zeros(Float64, (num_sections, num_cross_section_nodes))
   w_total = zeros(Float64, (num_sections, num_cross_section_nodes))

   for i=1:num_sections

       for j = 1:num_cross_section_nodes

           w_total[i,j] = integrate(zcoords[1:end], ϵ_total_normal[1:end,j])
           
           w[i, j] = w_total[i,j]/2 - integrate(zcoords[1:i], ϵ_total_normal[1:i,j])  #not sure if this is general or not

       end

   end

   return w

end





end #module
