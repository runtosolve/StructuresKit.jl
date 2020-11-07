module Mesh


export define_line_element, create_line_element_property_vector, extrude, extruded_surface


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
   extrude(x, y, z)

Accepts the `x` and `y` cross-section coordinates and `z` coordinates along the extrusion. Returns the 3D `coordinates` in the form `[x y z]` from the extrusion.
"""

function extrude(x, y, z)

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
   extruded_surface(num_cross_section_nodes, num_cross_section_elements, num_extruded_layers)

Returns a triangular surface mesh `connectivity` matrix as `[node_i node_j node_k]`.   If `num_cross_section_nodes = num_cross_section_elements`, then the cross-section is closed or multi-branched.

"""

function extruded_surface(num_cross_section_nodes, num_cross_section_elements, num_extruded_layers)

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

end #module
