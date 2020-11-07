module Beam

using DiffEqOperators: CenteredDifference
using LinearAlgebra
using RecursiveArrayTools: VectorOfArray
using NLsolve

using ..Mesh


export solve



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


function define(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, SpringStiffness, EndBoundaryConditions, Supports, UniformLoad)


   dz, z, dm = Mesh.define_line_element(MemberDefinitions)

   NumberOfNodes=length(dz)+1


   #define property vectors
   Ix = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, SectionProperties, 3, 1)
   Iy = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, SectionProperties, 3, 2)
   Ixy = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, SectionProperties, 3, 3)
   J = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, SectionProperties, 3, 4)
   Cw = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, SectionProperties, 3, 5)

   E = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, MaterialProperties, 4, 1)
   ν = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, MaterialProperties, 4, 2)
   G = E./(2 .*(1 .+ ν))

   ax = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, LoadLocation, 5, 1)
   ay = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, LoadLocation, 5, 2)

   kx = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, SpringStiffness, 6, 1)
   kϕ = Mesh.create_line_element_property_array(MemberDefinitions, dm, dz, SpringStiffness, 6, 2)


   Azzzz,Azz = calculateDerivativeOperators(dz) #calculate derivative operators

   #apply left and right end boundary condition stencils to derivative operators
   NthDerivative = 4
   Azzzz = applyEndBoundaryConditions(Azzzz,EndBoundaryConditions,NthDerivative,dz)

   NthDerivative = 2
   Azz = applyEndBoundaryConditions(Azz,EndBoundaryConditions,NthDerivative,dz)

   #define load
   qx = UniformLoad[1].*ones(NumberOfNodes)
   qy = UniformLoad[2].*ones(NumberOfNodes)

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
      A11[i,:] = E.*Iy.*Azzzz[i,:] .+kx.*AI[i,:]
      A12[i,:] = E.*Ixy.*Azzzz[i,:]
      A13[i,:] = kx.*ay.*AI[i,:]
      A21[i,:] = E.*Ixy.*Azzzz[i,:]
      A22[i,:] = E.*Ix.*Azzzz[i,:]
      A31[i,:] = kx.*ay.*AI[i,:]
      A33[i,:] = E.*Cw.*Azzzz[i,:] .-G.*J.*Azz[i,:] .+kx.*ay.*ay.*AI[i,:] .+kϕ.*AI[i,:] .+qx.*ax.*AI[i,:] .-qy.*ay.*AI[i,:]
   end

   #calculate RHS of AU=B
   B1 = qx
   B2 = qy
   B3 = qx.*ay .+qy.*ax

   #reduce problem to free dof
   FixedDOF = [findall(x->abs(x-Supports[i])<=10e-6,z) for i=1:length(Supports)]
   FixedDOF = VectorOfArray(FixedDOF)
   FixedDOF  = convert(Array,FixedDOF)
   FreeDOF = setdiff(1:NumberOfNodes,FixedDOF)

   A = [A11[FreeDOF,FreeDOF] A12[FreeDOF,FreeDOF] A13[FreeDOF,FreeDOF];
      A21[FreeDOF,FreeDOF] A22[FreeDOF,FreeDOF] A23[FreeDOF,FreeDOF];
      A31[FreeDOF,FreeDOF] A32[FreeDOF,FreeDOF] A33[FreeDOF,FreeDOF]]

   B = [B1[FreeDOF]; B2[FreeDOF]; B3[FreeDOF]]

   return A, B, FreeDOF, Ix, Iy, Ixy, J, Cw, E, ν, G, ax, ay, kx, kϕ, z, dm

end

function residual!(R,U,K,F)

   for i=1:length(F)
    R[i] = transpose(K[i,:])*U-F[i]
   end

   return R

end


function solve(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, SpringStiffness, EndBoundaryConditions, Supports, UniformLoad)

   K, F, FreeDOF, Ix, Iy, Ixy, J, Cw, E, ν, G, ax, ay, kx, kϕ, z, dm = define(MemberDefinitions, SectionProperties, MaterialProperties, LoadLocation, SpringStiffness, EndBoundaryConditions, Supports, UniformLoad)


   NumberOfNodes=length(z)

   u=zeros(NumberOfNodes)
   v=zeros(NumberOfNodes)
   ϕ=zeros(NumberOfNodes)

   Uguess = K \ F

   solution = nlsolve((R,U) ->residual!(R,U,K,F),Uguess)

   u[FreeDOF] = solution.zero[1:length(FreeDOF)]
   v[FreeDOF] = solution.zero[length(FreeDOF)+1:2*length(FreeDOF)]
   ϕ[FreeDOF] = solution.zero[2*length(FreeDOF)+1:3*length(FreeDOF)]

   properties = NamedTuple{(:dm, :Ix, :Iy, :Ixy, :J, :Cw, :E, :ν, :G, :ax, :ay, :kx, :kϕ)}((dm, Ix, Iy, Ixy, J, Cw, E, ν, G, ax, ay, kx, kϕ))

   return z, u, v, ϕ, properties

end

end #module
