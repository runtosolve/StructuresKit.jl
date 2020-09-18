module BeamMesh


export define, propvector


function define(memberDefinitions)

   #mesh along the length
   for i in eachindex(memberDefinitions)

      L = memberDefinitions[i][1]
      dL = memberDefinitions[i][2]
      numSegments = round(Int64, L/dL)

      if i == 1
         dz = ones(numSegments)*dL  #member discretization
      else
         dz = [dz; ones(numSegments)*dL]
      end

   end
   dz = dz
   z = [0; cumsum(dz)]

   #define what beam properties to apply at each node
   nodeprops = beamnodeprops(memberDefinitions)

   return dz, z, nodeprops

end

function beamnodeprops(memberDefinitions)

   #define member properties to use at each node

   if length(memberDefinitions) == 1

      L = memberDefinitions[1][1]
      dL = memberDefinitions[1][2]
      numSegments = round(Int64, L/dL)
      numNodes = numSegments+1   #number of nodes in a segment

      dm = ones(Int8, numNodes)

   elseif iseven(length(memberDefinitions))

      for i in eachindex(memberDefinitions)

         L = memberDefinitions[i][1]
         dL = memberDefinitions[i][2]
         numSegments = floor(Int64, L/dL)
         numNodes = numSegments+1   #number of nodes in a segment

         if i == 1
            dm = ones(Int8, numNodes-1)*i
         elseif i == length(memberDefinitions)
            dm = [dm; ones(Int8, numNodes)*i]
         else
            dm = [dm; ones(Int8, numNodes-1)*i]
         end

      end

   elseif isodd(length(memberDefinitions))

      middleSegment = (length(memberDefinitions) - 1)/2 + 1

      for i in eachindex(memberDefinitions)

         L = memberDefinitions[i][1]
         dL = memberDefinitions[i][2]
         numSegments = floor(Int64, L/dL)
         numNodes = numSegments+1   #number of nodes in a segment

         if i == 1
            dm = ones(Int8, numNodes-1)*i
         elseif i < middleSegment
            dm = [dm; ones(Int8, numNodes-1)*i]
         elseif i == middleSegment
            dm = [dm; ones(Int8, numNodes)*i]
         elseif i > middleSegment
            dm = [dm; ones(Int8, numNodes-1)*i]
         end

      end

   end

   dm = dm

   return dm

end


function propvector(MemberDefinitions, dm, dz, Property, PropertyOrder, PropertyType)


   z = [0; cumsum(dz)]

   NumberOfNodes = length(dm)

   A = zeros(NumberOfNodes)

   for i=1:NumberOfNodes

      PropertyIndex = MemberDefinitions[dm[i]][PropertyOrder]  #maps properties to each member along beam
      A[i] = Property[PropertyIndex][PropertyType]
   end

   return A

end

end #module
