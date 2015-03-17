## gmsh.jl
#
# Provides functionality for integrating gmsh meshes into Julia
# Currently this module does not differentiate between element 
# types.
#
# TODO: 
#    Update to handle PhysicalNames and ElementNodeData
#    Test how bad this is for really big files
#
# @author: arbennett
#
module gmsh


## Mesh
#
#  The mesh datatype holds the information needed to define the mesh.
#
type Mesh
    n_nodes::Int64
    n_elements::Int64
    nodes::Array{Float64, 2}
    boundary_nodes::Array{Int64}
    elements::Array{Float64, 2}

    Mesh(
         n_nodes::Int64,
         n_elements::Int64,
         nodes::Array{Float64,2},
         boundary_nodes::Array{Int64},
         elements::Array{Float64, 2}
        ) = new(n_nodes, n_elements, nodes, boundary_nodes, elements)
end


## elementTypes
#
#  An enum like structure for mapping the element types to their descriptions
#
baremodule elementTypes 
    const LINE = 1
    const TRIANGLE = 2
    const QUADRANGLE = 3
    const TETRAHEDRON = 4
    const HEXAHEDRON = 5
    const PRISM = 6
    const PYRAMID = 7
    const LINE_3_ND = 8
    const TRIANGLE_6_ND = 9
    const QUADRANGLE_9_ND = 10
    const TET_10_ND = 11
    const HEX_27_ND = 12
    const PRISM_18_ND = 13
    const PYRAMID_14_ND = 14
    const POINT = 15
    const QUADRANGLE_8_ND = 16
    const HEX_20_ND = 17
    const PRISM_15_ND = 18
    const PYRAMID_13_ND = 19
end


## read
#
#  Reads in the file from the given fileName and returns
#  the associated mesh representation.
#
function read(file_name::ASCIIString)
    # Make sure that the file exists
    if !isfile(file_name)
        println("Error in gmsh.read():")
        println("  Could not find " * file_name * "!")
        exit(1) 
    end

    # Open the file and begin grabbing data from it
    file_lines = readlines(open(file_name))
    idx::Int64 = 1
    n_lines::Int64 = length(file_lines)
    
    # Read until the mesh format is found
    while file_lines[idx] != "\$MeshFormat\n" && idx < n_lines
        idx += 1
    end
    
    # Make sure that we didn't reach the end of the file
    if idx == n_lines
        println("Error in gmsh.read():")
        println("  Could not determine mesh format!")
        exit(1)
    end
    
    # Grab the mesh format
    idx += 1
    mesh_format = file_lines[idx]
    idx += 1
    
    # Read until the node section is found
    while file_lines[idx] != "\$Nodes\n" && idx < n_lines
        idx += 1
    end

    # Make sure that we didn't reach the end of the file
    if idx == n_lines
        println("Error in gmsh.read():")
        println("  Could not read node data!")
        exit(1)
    end

    # Get the number of nodes
    idx += 1
    n_nodes = int(file_lines[idx])
    
    # Read in the node data
    nodes = Array(Float64, n_nodes, 3)
    for i = 1:n_nodes
       idx += 1
       nodes[i,:] = float(split(file_lines[idx])[2:end])
    end
    
    # Read until the element section is found
    while file_lines[idx] != "\$Elements\n" && idx < n_lines
        idx += 1
    end
    
    # Make sure that we didn't reach the end of the file
    if idx == n_lines
        println("Error in gmsh.read():")
        println("  Could not read element data!")
        exit(1)
    end

    # Get the number of elements
    idx += 1
    n_elem_entries::Int64 = int(file_lines[idx])
    
    # Determine how many columns of space needed to store element data
    # and the number of elements/boundaries on the mesh
    n_cols::Int64 = 0
    n_elems::Int64 = 0
    n_bc_nodes::Int64 = 0
    n_bc_edges::Int64 = 0
    for i = 1:n_elem_entries
        split_line = split(file_lines[idx + i])
        # Determine what kind of element we are on 
        if split_line[2] == "1"
            # Edge with boundary condition
            n_bc_edges += 1
        elseif split_line[2] == "15"
            # Node with boundary condition
            n_bc_nodes += 1
        else
            # Regular element
            n_elems += 1
            n_cols = max(n_cols, length(split_line)-1)
        end    
    end

    # Pull the element and boundary data
    elements::Array{Float64, 2} = zeros(n_elems, n_cols)
    boundary_nodes::Array{Int64} = zeros(Int64, n_bc_nodes)
    boundary_edges::Array{Int64, 2} = zeros(Int64, n_bc_edges, 2)
    elem_idx::Int64, bc_node_idx::Int64, bc_edge_idx::Int64 = 0, 0, 0 
    for i = 1:n_elem_entries
        idx += 1
        split_line = split(file_lines[idx])
        # Determine what kind of element we are on 
        if split_line[2] == "1"
            # Edge with boundary condition
            bc_edge_idx += 1
            boundary_edges[bc_edge_idx] = int(split_line)[end-1,end]
        elseif split_line[2] == "15"
            # Node with boundary condition
            bc_node_idx += 1
            boundary_nodes[bc_node_idx] = int(split_line)[end]
        else
            # Regular element
            elem_idx += 1           
            elements[elem_idx,1:length(split_line) - 1 ] = int(split_line)[2:end]
        end    
    end
    
    # Create the mesh and return it
    return Mesh(n_nodes, n_elems, nodes, boundary_nodes, elements)
end


end # End module gmsh
