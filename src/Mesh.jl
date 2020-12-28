module Mesh

using Triangle

export define_line_element, create_line_element_property_vector, extrude_nodes, extrude_elements, open_cross_section_tessellation, surface,
   calculate_weights


"""
    define_line_element(member_definitions)

Accepts `member_definitions` tuple and returns element discretization `dz`, line coordinates 'z', and nodal properties at each node in the line `node_props`.
"""

function define_line_element(member_definitions)

   #mesh along the length
   for i in eachindex(member_definitions)

      L = member_definitions[i][1]
      dL = member_definitions[i][2]
      num_segments = round(Int64, L/dL)

      if i == 1
         dz = ones(num_segments)*dL  #member discretization
      else
         dz = [dz; ones(num_segments)*dL]
      end

   end
   dz = dz
   z = [0; cumsum(dz)]

   #define what properties to apply at each node
   node_props = assign_line_element_nodal_properties(member_definitions)

   return dz, z, node_props

end

"""
   assign_line_element_nodal_properties(member_definitions)

Accepts `member_definitions` tuple and returns a properties assignment for each node line segment `dm`.   

There are some tricks to defining the nodal properties. If the number of `member_definitions` is even, the first and last segments are added as `num_nodes - 1`.  If the number of `member_definitions` is odd, middle segment is added as 'num_nodes - 1'.

"""


function assign_line_element_nodal_properties(member_definitions)

   #define member properties to use at each node

   if length(member_definitions) == 1

      L = member_definitions[1][1]
      dL = member_definitions[1][2]
      num_segments = round(Int64, L/dL)
      num_nodes = num_segments + 1   #number of nodes in a segment

      dm = ones(Int8, num_nodes)

   elseif iseven(length(member_definitions))

      for i in eachindex(member_definitions)

         L = member_definitions[i][1]
         dL = member_definitions[i][2]
         numSegments = floor(Int64, L/dL)
         num_nodes = num_segments+1   

         if i == 1
            dm = ones(Int8, num_nodes-1)*i
         elseif i == length(member_definitions)
            dm = [dm; ones(Int8, num_nodes)*i]
         else
            dm = [dm; ones(Int8, num_nodes-1)*i]
         end

      end

   elseif isodd(length(member_definitions))

      middle_segment = (length(member_definitions) - 1)/2 + 1

      for i in eachindex(member_definitions)

         L = member_definitions[i][1]
         dL = member_definitions[i][2]
         num_segments = floor(Int64, L/dL)
         num_nodes = num_segments + 1   

         if i == 1
            dm = ones(Int8, num_nodes-1)*i
         elseif i < middle_segment
            dm = [dm; ones(Int8, num_nodes-1)*i]
         elseif i == middle_segment
            dm = [dm; ones(Int8, num_nodes)*i]
         elseif i > middle_segment
            dm = [dm; ones(Int8, num_nodes-1)*i]
         end

      end

   end

   dm = dm

   return dm

end

"""
   create_line_element_property_vector(member_definitions, dm, dz, property, property_order, property_type)

Accepts the `member_definitions` tuple, property array 'dm' assigning a property to each node, discretization along the line segment `dz`, the `property` magnitude, `property_order` which defines how the properties are ordered in `member_definitions` , and `property_type` which is the flavor of the property (e.g., `Ix` or `Iy` or `J` or ...) .  The function returns `property_array` with the magnitude of a property assigned node by node.   

"""


function create_line_element_property_array(member_definitions, dm, dz, property, property_order, property_type)


   z = [0; cumsum(dz)]

   num_nodes = length(dm)

   property_array = zeros( num_nodes)

   for i=1:num_nodes
      property_index = member_definitions[dm[i]][property_order]  #maps properties to each segment along a line element
      property_array[i] = property[property_index][property_type]
   end

   return property_array

end

"""
   extrude_nodes(x, y, z)

Accepts the `x` and `y` cross-section coordinates and `z` coordinates along the extrusion. Returns the 3D `coordinates` in the form `[x y z]` from the extrusion.
"""

function extrude_nodes(x, y, z)

   num_nodes_section = length(x)
   coordinates = zeros(num_nodes_section, 3)   #initialize coordinates array
   
   for i in 1:length(z)
      z_section = z[i] .* ones(length(x))
      if i==1
         coordinates .= [x y z_section]   #start with first cross-section
      else
         coordinates = vcat(coordinates, [x y z_section])  #add following cross-sections along z
      end
   end

   return coordinates

end

"""
   extrude_elements(num_cross_section_nodes, num_cross_section_elements, num_extruded_layers)

Returns a triangular surface mesh `connectivity` matrix as `[node_i node_j node_k]`.   If `num_cross_section_nodes = num_cross_section_elements`, then the cross-section is closed or multi-branched.

"""

function extrude_elements(num_cross_section_nodes, num_cross_section_elements, num_extruded_layers)

   num_elements = num_cross_section_elements * (num_extruded_layers - 1)

   uppoly = zeros(Int, num_elements, 3)
   downpoly = zeros(Int, num_elements, 3)

   for j = 1:num_extruded_layers - 1

      node_offset =  num_cross_section_nodes * (j - 1)
      elem_offset =  (j - 1) * num_cross_section_elements

      for i = 1:num_cross_section_elements

         element_count = i + elem_offset
      
         #triangular element 
                     #c
                  #  #
               #     #
            #        #
         #           #
        #a###########b
      
         a_up = i
         b_up = i + 1
         c_up = i + 1 + num_cross_section_nodes  

         #address closed or multi-branch cross-sections
         if i == num_cross_section_nodes
            b_up = 1 + 1
            c_up = 1 + 1 + num_cross_section_nodes  
         end
            
         uppoly[element_count, :] = [a_up b_up c_up] .+ node_offset

         #triangular element 
         #b#########c
         #        #
         #      #
         #    #
         #  #
         #a

         a_down = i
         b_down = i + num_cross_section_nodes
         c_down = i + num_cross_section_nodes + 1

         #address closed or multi-branch cross-sections 
         if i == num_cross_section_nodes
            c_down = 1 + num_cross_section_nodes + 1
         end
      
         downpoly[element_count, :] = [a_down b_down c_down] .+ node_offset
      
      end
   
   end

   connectivity=vcat(uppoly,downpoly)

   return connectivity

end

"""
   open_cross_section_tessellation(x, y)

Accepts `x` and `y` coordinates for an open cross-section and determines a through-thickness triangular mesh `shape_mesh`.  The Array `cross_section_edges = [node(i) node(j)]` is also output, providing the nodal connectivity that describes the cross-section boundary.

"""

function open_cross_section_tessellation(x, y)

    #define cross-section edges
    cross_section_edges = zeros(Int, length(x), 2)
    for i = 1:length(x) - 1
        cross_section_edges[i, 1] = i
        cross_section_edges[i, 2] = i + 1 
    end
    cross_section_edges[end, :] = [length(x) 1]  #for open cross-section only
  
    points = [x y]
    points_map = [i for i=1:length(x)]
    shape_mesh = Triangle.constrained_triangulation(points, points_map, cross_section_edges)

    return shape_mesh, cross_section_edges

end


"""
   surface(xcoords, ycoords, zcoords)

Accepts 3D surface coordinates of a structural member and defines a triangular surface mesh with `coordinates` and `connectivity`.  

"""

function surface(xcoords, ycoords, zcoords)

   #generate triangular mesh at the ends
   section_mesh, cross_section_edges = open_cross_section_tessellation(xcoords, ycoords)

   #get member 3D coordinates via extrusion
   coordinates = extrude_nodes(xcoords, ycoords, zcoords)

   #define triangular mesh along extrusion
   num_cross_section_nodes = length(xcoords)
   num_extruded_layers = length(zcoords)
   num_cross_section_elements = size(cross_section_edges)[1]
   extrusion_connectivity = extrude_elements(num_cross_section_nodes, num_cross_section_elements, num_extruded_layers)

   # mesh beginning and end faces of extrusion
   #beginning face
   begin_connectivity = zeros(Int, length(section_mesh), 3)
   for i=1:length(section_mesh)
       begin_connectivity[i, :] = section_mesh[i]
   end

   #end face
   end_connectivity = begin_connectivity .+ num_cross_section_nodes * (num_extruded_layers - 1)

   #combine to define the full surface with end caps   
   connectivity = [extrusion_connectivity; begin_connectivity; end_connectivity]

   return coordinates, connectivity

end


#############################################################
# Fornberg algorithm

# This implements the Fornberg (1988) algorithm (https://doi.org/10.1090/S0025-5718-1988-0935077-0)
# to obtain Finite Difference weights over arbitrary points to arbitrary order.

function calculate_weights(order::Int, x0::T, x::AbstractVector) where T<:Real
   #=
       order: The derivative order for which we need the coefficients
       x0   : The point in the array 'x' for which we need the coefficients
       x    : A dummy array with relative coordinates, e.g., central differences
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


end #module
