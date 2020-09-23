module BeamColumn

using DiffEqOperators: CenteredDifference
using LinearAlgebra
using RecursiveArrayTools: VectorOfArray
using NLsolve

using ..BeamMesh


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


function define(MemberDefinitions, SectionProperties, MaterialProperties, Loads, Springs, EndBoundaryConditions, Supports)


   dz, z, dm = BeamMesh.define(MemberDefinitions)

   NumberOfNodes=length(dz)+1


   #define property vectors
   A =  BeamMesh.propvector(MemberDefinitions, dm, dz, SectionProperties, 3, 1)
   Ix = BeamMesh.propvector(MemberDefinitions, dm, dz, SectionProperties, 3, 2)
   Iy = BeamMesh.propvector(MemberDefinitions, dm, dz, SectionProperties, 3, 3)
   J = BeamMesh.propvector(MemberDefinitions, dm, dz, SectionProperties, 3, 4)
   Cw = BeamMesh.propvector(MemberDefinitions, dm, dz, SectionProperties, 3, 5)

   xo = BeamMesh.propvector(MemberDefinitions, dm, dz, SectionProperties, 3, 6)
   yo = BeamMesh.propvector(MemberDefinitions, dm, dz, SectionProperties, 3, 7)

   Io = Ix .+ Iy .+ A .* (xo.^2 + yo.^2)

   E = BeamMesh.propvector(MemberDefinitions, dm, dz, MaterialProperties, 4, 1)
   ν = BeamMesh.propvector(MemberDefinitions, dm, dz, MaterialProperties, 4, 2)
   G = E./(2 .*(1 .+ ν))

   kx = BeamMesh.propvector(MemberDefinitions, dm, dz, Springs, 5, 1)
   ky = BeamMesh.propvector(MemberDefinitions, dm, dz, Springs, 5, 2)
   kϕ = BeamMesh.propvector(MemberDefinitions, dm, dz, Springs, 5, 3)

   hx = BeamMesh.propvector(MemberDefinitions, dm, dz, Springs, 5, 4)
   hy = BeamMesh.propvector(MemberDefinitions, dm, dz, Springs, 5, 5)

   Azzzz,Azz = calculateDerivativeOperators(dz) #calculate derivative operators

   #apply left and right end boundary condition stencils to derivative operators
   NthDerivative = 4
   Azzzz = applyEndBoundaryConditions(Azzzz,EndBoundaryConditions,NthDerivative,dz)

   NthDerivative = 2
   Azz = applyEndBoundaryConditions(Azz,EndBoundaryConditions,NthDerivative,dz)

   #define load and load locations
   P =  Loads[1]
   qx = Loads[2]
   qy = Loads[3]
   ax = Loads[4]
   ay = Loads[5]

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
   FixedDOF = [findall(x->abs(x-Supports[i])<=10e-6,z) for i=1:length(Supports)]
   FixedDOF = VectorOfArray(FixedDOF)
   FixedDOF  = convert(Array,FixedDOF)
   FreeDOF = setdiff(1:NumberOfNodes,FixedDOF)

   A = [A11[FreeDOF,FreeDOF] A12[FreeDOF,FreeDOF] A13[FreeDOF,FreeDOF];
      A21[FreeDOF,FreeDOF] A22[FreeDOF,FreeDOF] A23[FreeDOF,FreeDOF];
      A31[FreeDOF,FreeDOF] A32[FreeDOF,FreeDOF] A33[FreeDOF,FreeDOF]]

   B = [B1[FreeDOF]; B2[FreeDOF]; B3[FreeDOF]]


   properties = NamedTuple{(:dm, :z, :P, :Ix, :Iy, :J, :Cw, :xo, :yo, :Io, :E, :ν, :G, :ax, :ay, :kx, :ky, :kϕ, :hx, :hy, :qx, :qy)}((dm, z, P, Ix, Iy, J, Cw, xo, yo, Io, E, ν, G, ax, ay, kx, ky, kϕ, hx, hy, qx, qy))


   return A, B, FreeDOF, properties

end

function residual!(R,U,K,F)

   for i=1:length(F)
    R[i] = transpose(K[i,:])*U-F[i]
   end

   return R

end


function solve(MemberDefinitions, SectionProperties, MaterialProperties, Loads, Springs, EndBoundaryConditions, Supports)

   K, F, FreeDOF, properties = define(MemberDefinitions, SectionProperties, MaterialProperties, Loads, Springs, EndBoundaryConditions, Supports)


   NumberOfNodes=length(properties.z)

   u=zeros(NumberOfNodes)
   v=zeros(NumberOfNodes)
   ϕ=zeros(NumberOfNodes)

   Uguess = K \ F

   solution = nlsolve((R,U) ->residual!(R,U,K,F),Uguess)

   u[FreeDOF] = solution.zero[1:length(FreeDOF)]
   v[FreeDOF] = solution.zero[length(FreeDOF)+1:2*length(FreeDOF)]
   ϕ[FreeDOF] = solution.zero[2*length(FreeDOF)+1:3*length(FreeDOF)]

   return u, v, ϕ, properties

end

end #module
