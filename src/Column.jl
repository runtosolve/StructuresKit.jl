module Column

using DiffEqOperators: CenteredDifference
using LinearAlgebra
using RecursiveArrayTools: VectorOfArray
using NLsolve

using ..Mesh
using ..InternalForces


export solve, normal_stresses, normal_strains, warping_displacements, deformed_shape



function calculateDerivativeOperators(dz)

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

function applyEndBoundaryConditions(A, EndBoundaryConditions, NthDerivative, dz)

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


function define(MemberDefinitions, SectionProperties, MaterialProperties, Loads, Springs, EndBoundaryConditions, Supports, imperfections)


   dz, z, dm = Mesh.define_line_element(MemberDefinitions)

   NumberOfNodes=length(dz)+1


   #define property vectors
   A =  Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, SectionProperties, 3, 1)
   Ix = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, SectionProperties, 3, 2)
   Iy = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, SectionProperties, 3, 3)
   J = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, SectionProperties, 3, 4)
   Cw = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, SectionProperties, 3, 5)
   xc = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, SectionProperties, 3, 6)
   yc = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, SectionProperties, 3, 7)
   xs = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, SectionProperties, 3, 8)
   ys = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, SectionProperties, 3, 9)

   xo = -(xc .- xs)
   yo = yc .- ys

   Io = Ix .+ Iy .+ A .* (xo.^2 + yo.^2)

   E = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, MaterialProperties, 4, 1)
   ν = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, MaterialProperties, 4, 2)
   G = E./(2 .*(1 .+ ν))

   kx = Springs[1]
   ky = Springs[2]
   kϕ = Springs[3]

   hx = Springs[4]
   hy = Springs[5]

   #define imperfections 
   uo = imperfections[1]
   vo = imperfections[2]
   ϕo = imperfections[3] 


   Azzzz,Azz = calculateDerivativeOperators(dz) #calculate derivative operators

   #apply left and right end boundary condition stencils to derivative operators
   NthDerivative = 4
   Azzzz = applyEndBoundaryConditions(Azzzz,EndBoundaryConditions,NthDerivative,dz)

   NthDerivative = 2
   Azz = applyEndBoundaryConditions(Azz,EndBoundaryConditions,NthDerivative,dz)

   #define load
   P =  Loads


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
      A11[i,:] = E.*Iy.*Azzzz[i,:] .+ P .* Azz[i,:] .+kx.*AI[i,:]
      A13[i,:] = kx.*(yo .-hy).*AI[i,:] .+ P.*yo .* Azz[i,:]
      A22[i,:] = E.*Ix.*Azzzz[i,:] .+ P.*Azz[i,:] .+ ky.*AI[i,:]
      A23[i,:] = -ky.*(xo .-hx).*AI[i,:] - P .* xo .* Azz[i,:]
      A31[i,:] = kx.*(yo .-hy).*AI[i,:] .+ P .* yo .* Azz[i,:]
      A32[i,:] = -ky.*(xo .-hx).*AI[i,:] .- P .* xo .* Azz[i,:]
      A33[i,:] = E.*Cw.*Azzzz[i,:] .-(G.*J .- (P.*Io./A)).*Azz[i,:] .+kx.*(yo .-hy).^2 .*AI[i,:] .+ky.*(xo .- hx).^2 .*AI[i,:] .+kϕ.*AI[i,:]
   end

   #calculate RHS of AU=B
   B1 = -P .* Azz*uo - P .* yo .* Azz*ϕo
   B2 = -P .* Azz*vo + P .* xo .* Azz*ϕo
   B3 = -P .* yo .* Azz*uo + P .* xo .* Azz*vo - P .* (Io ./ A) .* Azz*ϕo

   #reduce problem to free dof
   FixedDOF = [findall(x->abs(x-Supports[i])<=10e-6,z) for i=1:length(Supports)]
   FixedDOF = VectorOfArray(FixedDOF)
   FixedDOF  = convert(Array,FixedDOF)
   FreeDOF = setdiff(1:NumberOfNodes,FixedDOF)

   Am = [A11[FreeDOF,FreeDOF] A12[FreeDOF,FreeDOF] A13[FreeDOF,FreeDOF];
      A21[FreeDOF,FreeDOF] A22[FreeDOF,FreeDOF] A23[FreeDOF,FreeDOF];
      A31[FreeDOF,FreeDOF] A32[FreeDOF,FreeDOF] A33[FreeDOF,FreeDOF]]

   Bm = [B1[FreeDOF]; B2[FreeDOF]; B3[FreeDOF]]

   properties = NamedTuple{(:dm, :z, :P, :A, :Ix, :Iy, :J, :Cw, :xc, :yc, :xs, :ys, :xo, :yo, :Io, :E, :ν, :G, :kx, :ky, :kϕ, :hx, :hy, :uo, :vo, :ϕo)}((dm, z, P, A, Ix, Iy, J, Cw, xc, yc, xs, ys, xo, yo, Io, E, ν, G, kx, ky, kϕ, hx, hy, uo, vo, ϕo))


   return Am, Bm, FreeDOF, properties

end

function residual!(R, U, K, F)

   for i=1:length(F)
    R[i] = transpose(K[i,:]) * (U) - F[i]
   end

   return R

end


function solve(member_definitions, section_properties, material_properties, loads, springs, end_boundary_conditions, supports, imperfections)

   K, F, free_dof, properties = define(member_definitions, section_properties, material_properties, loads, springs, end_boundary_conditions, supports, imperfections)

   num_nodes=length(properties.z)

   u=zeros(num_nodes)
   v=zeros(num_nodes)
   ϕ=zeros(num_nodes)

   deformation_guess = K \ F    #consider revising this for large systems, it might be slow...

   solution = nlsolve((R, U) ->residual!(R, U, K, F), deformation_guess)

   u[free_dof] = solution.zero[1:length(free_dof)]
   v[free_dof] = solution.zero[length(free_dof)+1:2*length(free_dof)]
   ϕ[free_dof] = solution.zero[2*length(free_dof)+1:3*length(free_dof)]

   return u, v, ϕ, properties

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
