module Truss

using LinearAlgebra

export define, solve


"""
        truss_object

This object holds all the definitions (material properties, section properties, geometry) and the solution (displacements, reactions, internal forces) for a system of truss elements.

num_dof_per_node::Int64
members::Array{NTuple{4,Int64},1}
node_geometry::Array{Float64,2}
section_properties::Array{Float64,1}
material_properties::Array{Float64,1}
supports::Array{Int64,1}
external_forces::Array{Float64,1}
L::Array{Float64,1}
E::Array{Float64,1}
A::Array{Float64,1}
θ::Array{Float64,1}
K::Array{Float64,2}
Kff::Array{Float64,2}
Ksf::Array{Float64,2}
Kfs::Array{Float64,2}
Kss::Array{Float64,2}
F::Array{Float64,1}
Fs::Array{Float64,1}
Ff::Array{Float64,1}
s::Array{Int64,1}
f::Array{Int64,1}
T::Array{Array{Float64,2},1}
k_element_local::Array{Array{Float64,2},1}
k_element_global::Array{Array{Float64,2},1}  
u::Array{Float64,1}
uf::Array{Float64,1}
f_element::Array{Array{Float64,1},2}  

"""


mutable struct truss_object

    num_dof_per_node::Int64
    members::Array{NTuple{4,Int64},1}
    node_geometry::Array{Float64,2}
    section_properties::Array{Float64,1}
    material_properties::Array{Float64,1}
    supports::Array{Int64,1}
    external_forces::Array{Float64,1}
    L::Array{Float64,1}
    E::Array{Float64,1}
    A::Array{Float64,1}
    θ::Array{Float64,1}
    K::Array{Float64,2}
    Kff::Array{Float64,2}
    Ksf::Array{Float64,2}
    Kfs::Array{Float64,2}
    Kss::Array{Float64,2}
    F::Array{Float64,1}
    Fs::Array{Float64,1}
    Ff::Array{Float64,1}
    s::Array{Int64,1}
    f::Array{Int64,1}
    T::Array{Array{Float64,2},1}
    k_element_local::Array{Array{Float64,2},1}
    k_element_global::Array{Array{Float64,2},1}  
    u::Array{Float64,1}
    uf::Array{Float64,1}
    f_element::Array{Array{Float64,1},2}  

    #This allows for the construction of a new object.
    truss_object() = new()
    
end


"""
        calculate_truss_element_lengths(members, node_geometry)

Calculate the truss member lengths in a structural system.

`members` is a NTuple{4, Int64} containing node i, node j, A, E assignments
`node_geometry` is array with the global x-y coordinates

"""

function calculate_truss_element_lengths(members, node_geometry)

    #Initialize array to hold truss lengths.
    truss_lengths=zeros(Float64, length(members))

    #Loop over each truss element.
    for i=1:length(members)

        #Define truss node numbers. 
        node_i=members[i][1]
        node_j=members[i][2]

        #Define truss node coordinates.
        node_i_xy=node_geometry[node_i,:]
        node_j_xy=node_geometry[node_j,:]

        #Calculate the truss element length.
        truss_lengths[i] =norm(node_i_xy - node_j_xy)

    end

    return truss_lengths

end

"""
    define_local_truss_element_stiffness_matrix(E, A, L)

Define a truss element stiffness matrix in its local coordinate system.

`E` is elastic modulus.
`A` is the cross-sectional area.
`L` is truss length.

There are 4 degrees of freedom for the truss element.

      2                 4
    1>^--------------->3^

"""

function define_local_truss_element_stiffness_matrix(E, A, L)
	
	k_truss_element_local = E * A/ L * [1  0 -1 0
		                  0  0 0 0
		                  -1 0  1 0
	                      0  0  0 0]
	
	return k_truss_element_local
	
end

"""
        define_truss_element_orientations(members, node_geometry)

Define the truss element orientation angle.

`members` is a NTuple{4, Int64} containing node i, node j, A, E assignments
`node_geometry` is array with the global x-y coordinates

A positive angle is defined counterclockwise from the horizon, rotating about node i.  

"""

function define_truss_element_orientations(members, node_geometry)

    # Initialize the orientation array for all truss elements in a structural system.
    θ = zeros(Float64, length(members))

    # Loop over the truss elements.
    for i=1:length(members)

        point_i = node_geometry[members[i][1], :]
        point_j = node_geometry[members[i][2], :]

        element_vector = point_j - point_i

        #Calculate the orientation angle.
        #https://en.wikipedia.org/wiki/Atan2
        θ[i] = atan(element_vector[2], element_vector[1])
    
    end

    return θ

end

"""
        vector_rotation_operator(θ)

Define an operator matrix that converts truss element displacements to a new, rotated coordinate system.   

`θ` is the orientation angle of a truss element.

A positive angle is defined counterclockwise from the horizon, rotating about node i.  

"""

function vector_rotation_operator(θ)
	
	T = [cos(θ) sin(θ) 0 0
		 -sin(θ) cos(θ) 0 0
		 0 0 cos(θ) sin(θ)
		 0 0  -sin(θ) cos(θ)]
	
	return T
	
end
	
"""
    assign_member_properties(members, member_property, property_order, property_type)

Define an array containing the member property values for each element in a structural system.   For example, define all the cross-sectional areas for all the members in a truss. 

`members` is a NTuple{4, Int64} containing node i, node j, A, E assignments.
`member_property` is an array containing all the available properties that a member can take on. 
`property order` is the index in the `members` NTuple where the specific property is assigned.

"""


function assign_member_properties(members, member_property, property_order, property_type)

    #Initialize the property array
    property = zeros(Float64, length(members))

    #Loop over all the members.
    for i=1:length(members)

        #Assign a specific property (area, elastic modulus) to each element in a structural system.
        property[i] = member_property[members[i][property_order]][property_type]

    end

    return property

end


"""
    calculate_global_element_stiffness_matrix(k_element_local, θ)

Transform a local element stiffness matrix to global coordinates.

`k_element_local` is the element stiffness matrix in its local coordinate system.
`θ` is the orientation angle of the element.

"""

function calculate_global_element_stiffness_matrix(k_element_local, T)

    #Transform the local element stiffness matrix to global coordinates.
    k_element_global = T' * k_element_local * T

    return k_element_global

end


"""
    assemble_global_stiffness_matrix(node_geometry, members, k_element_global, num_dof_per_node)

Assemble all the element stiffness matrices into the global system stiffness matrix.

`node_geometry` is an array with the global node x-y coordinates.
`members` is a NTuple{4, Int64} containing node i, node j, and property assignments.
`k_element_global` are all the element stiffness matrices in a 3D array.
`num_dof_per_node` defines the number of degrees of freedom per element node.  For a truss, this is 2, for a beam it is 3, ...

"""

function assemble_global_stiffness_matrix(node_geometry, members, k_element_global, num_dof_per_node)

    num_nodes = size(node_geometry)[1]

    k_system_global = zeros(Float64, num_nodes * num_dof_per_node, num_nodes * num_dof_per_node)

    for i=1:length(members)

        node_i = members[i][1]
        node_j = members[i][2]

        node_i_dof = [1;2] .+ (node_i - 1) * num_dof_per_node
        node_j_dof = [1;2] .+ (node_j - 1) * num_dof_per_node

        global_dof = [node_i_dof; node_j_dof]

        k_system_global[global_dof, global_dof] = k_system_global[global_dof, global_dof] + k_element_global[i]

    end

    return k_system_global

end


"""
    define(members, section_properties, material_properties, node_geometry, supports, external_forces)

Define and partition the global system stiffness matrix and external force vector.

`members` is a NTuple{4, Int64} containing node i, node j, and property assignments.
'section_properties' defines a library of section properties.  
'material_properties' defines a library of material properties.  
`node_geometry` is an array with the global node x-y coordinates.
`supports` defines the degrees of freedom that are free (0) and fixed(1).
`external_forces` is an array of the external forces applied to the nodes, in global coordinates.

Partition the stiffness matrix into Kff (free free), Kfs (free supported), Ksf (supported free), and Kss (supported supported).  Partition the external forces at the free dof into Ff. Also return `f` array of free global degrees of freedom, `s` fixed global degrees of freedom, `T` all the rotation matrices for each truss element, and the `k_element_local` and `k_element_global` matrices.

"""

function define(members, section_properties, material_properties, node_geometry, supports, external_forces)

    truss = truss_object()

    truss.members = members
    truss.section_properties = section_properties
    truss.material_properties = material_properties
    truss.node_geometry = node_geometry
    truss.supports = supports

    num_nodes = size(node_geometry)[1]
    truss.external_forces = zeros(Float64, num_nodes * 2)
    truss.external_forces .= external_forces

    #Define the number of degrees of freedom per node.
    truss.num_dof_per_node = 2

    #Calculate truss member lengths.
    truss.L = calculate_truss_element_lengths(members, node_geometry)

    #Calculate truss member orientations.
    truss.θ = define_truss_element_orientations(members, node_geometry)

    #Define truss properties for stiffness matrix calculations.
    truss.A = assign_member_properties(members, section_properties, 3, 1)
    truss.E = assign_member_properties(members, material_properties, 4, 1)

    #Calculate the local stiffness matrix for each element.
    truss.k_element_local = define_local_truss_element_stiffness_matrix.(truss.A, truss.E, truss.L)

    #Define rotation matrix for each element.
    truss.T = vector_rotation_operator.(truss.θ)

    #Rotate local stiffness matrix into global coordinates.
    truss.k_element_global = calculate_global_element_stiffness_matrix.(truss.k_element_local, truss.T)

    #Assemble the global stiffness matrix.
    truss.K = assemble_global_stiffness_matrix(node_geometry, members, truss.k_element_global, truss.num_dof_per_node)

    #Define free degrees of freedom.  
    truss.f = findall(x-> x==0, supports)

    #Define number of free dof.
    num_free_dof = length(truss.f)

    #Define fixed degrees of freedom.
    truss.s = findall(x-> x==1, supports)

    #Define number of fixed dof.
    num_fixed_dof = length(truss.s)

    #Initialize partitioned stiffness matrices.
    truss.Kff = zeros(Float64, (num_free_dof, num_free_dof))
    truss.Ksf = zeros(Float64, (num_fixed_dof, num_free_dof))
    truss.Kfs = zeros(Float64, (num_free_dof, num_fixed_dof))
    truss.Kss = zeros(Float64, (num_fixed_dof, num_fixed_dof))

    #Partition the global stiffness matrix.
    truss.Kff .= truss.K[truss.f, truss.f]
    truss.Ksf .= truss.K[truss.s, truss.f]
    truss.Kfs .= truss.K[truss.f, truss.s]
    truss.Kss .= truss.K[truss.s, truss.s]

    #Partition the external force vectors.
    truss.F = external_forces
    truss.Ff = zeros(Float64, length(truss.f))
    truss.Ff .= truss.F[truss.f]

    return truss

end



"""
    solve(K, Kff, Ksf, F, Ff, s, f, T, k_element_local)

Solve for system displacements, reactions, and internal forces.

`K` is the global stiffness matrix.
'Kff' is the partitioned portion of the global stiffness matrix, free by free.
'Ksf' is the partitioned portion of the global stiffness matrix, support by free. 
`F` is an array of the external forces.
`s` is an array of the fixed degrees of freedom.
`f` is an array of the support degrees of freedom.
`T` are all the rotation matrices for each element.
`k_element_local` are all the element local stiffness matrices.

Calculate the global displacement vector `uf` for the free degrees of freedom, provide the full global displacement vector `u`, and calculate all the element internal forces as an array of arrays `f_element`. 

"""

function solve(truss)

    #Define the total number of degrees of freedom in the system.
    num_dof = size(truss.K)[1]

    #Initialize the global displacement vector.
    truss.u = zeros(Float64, num_dof)

    #Calculate the global displacements at the free degrees of freedom.
    truss.uf = truss.Kff^-1 * truss.Ff

    #Insert uf into u.
    truss.u[truss.f] .= truss.uf

    #Calculate the global system reactions.
    truss.Fs = truss.Ksf * truss.uf

    #Insert Fs into F.
    truss.F[truss.s] .= truss.Fs

    #Define the number of elements.
    num_elem = length(truss.A)

    #Initialize an array of element force vectors.
    truss.f_element = fill(Float64[], num_elem, 1) 

    #Calculate the internal forces in each element.
    for i = 1:num_elem

        #Define node i and j for an element.
        node_i = truss.members[i][1]
        node_j = truss.members[i][2]

        #Calculate the global dof numbers at each element nodes i and j.
        node_i_dof = [1;2] .+ (node_i - 1) * truss.num_dof_per_node
        node_j_dof = [1;2] .+ (node_j - 1) * truss.num_dof_per_node

        #Assemble all the global dof numbers for the element.
        global_dof = [node_i_dof; node_j_dof]

        #Find the element displacements in global coordinates.
        u_element_global = zeros(Float64, length(global_dof))
        u_element_global .= truss.u[global_dof]

        #Transform the element displacements into local coordinates.
        u_element_local = truss.T[i] * u_element_global

        #Calculate the element local forces.
        truss.f_element[i] = truss.k_element_local[i] * u_element_local

    end

    return truss

end

end #module





