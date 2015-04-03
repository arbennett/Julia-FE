## solver.jl
#
#  The solver for the finite element method.  Currently this
#  solver implements the METHOD TYPE HERE on triangular elements
#  with a single quadrature point.
#
#  @author Andrew Bennett
#

## Global Definitions
#----------------------

## Finite element space
#
#  Definition of the C0 finite element space
#
type Solution_Space
    # Build the initial solution space with boundary conditions imposed
    function Solution_Space(mesh::Mesh, boundaryCondition::Function)
        n_nodes::Int64 = mesh.n_nodes
        n_elements::Int64 = mesh.n_elements
        n_boundary::Int64 = length(mesh.boundary_nodes)
        bc_count::Int64 = 0
        
        boundary_nodes::Array{Int64} = mesh.boundary_nodes
        
        # Go through and set up the boundary conditions
        for i = 1:n_nodes
            x = mesh.nodes[i,1]
            y = mesh.nodes[i,2]
            if i in boundary_nodes
                bc_count += 1
                node_vals[bc_count] = boundaryCondition(x,y)
            end
        end
        
        new(n_nodes, n_elements, n_boundary, boundary_nodes)
    end
end


## Quadrature points and weight
const xi = 1/3
const eta = 1/3
const weight = 1/2

## Shape functions & derivatives along axes
psi1(xi,eta) = 1 - xi - eta
const dpsi1_dxi = -1.0
const dpsi1_deta = -1.0

psi2(xi,eta) = xi
const dpsi2_dxi = 1.0
const dpsi2_deta = 0.0

psi3(xi,eta) = eta
const dpsi3_dxi = 0.0
const dpsi3_deta = 1.0

## Solver functions
#--------------------

## The top level solve function
#
function solve(mesh::Mesh, stiffness::Function, bc::Function, ext_force::Function)
	
	u = Solution_Space(mesh, boundaryCondition)
	size = u.n_nodes - u.n_boundary
	
	nodes_i, nodes_j, vals, load = assemble_matrices(mesh, u, stiffness, ext_force)
	
	stiffness = sparse(nodes_i, nodes_j, vals, size, size)
	
	U = K\F
	
	return U		
end


## Build the stiffness matrix and load vector
#
assemble_matrices(mesh::Mesh, field::Solution_Space, stiffness::Function, externalForce::Function, F::Array{Float64,1})
	size = field.n_boundary
	b_nodes = field.boundary_nodes
	global_load = zeros(field.n_nodes - b_nodes)
	count = 0
	
	nodes_i = zeros(Int, size)
	nodes_j = zeros(Int, size)
	vals = zeros(size)
	
	psi = [psi1(xi,eta), psi2(xi,eta), psi3(xi,eta)]
	Dpsi = zeros(3,2)
	
	# Stiffness and load for an element
	stiff = zeros(3,3)
	load = zeros(3)
	
	for elemIdx = 1:field.n_elements
		# Get element nodes, and their locations
		n = mesh.elements[elemIdx,:]
		x1, y1 = mesh.nodes[n[1],:]
		x2, y2 = mesh.nodes[n[2],:]
		x3, y3 = mesh.nodes[n[3],:]
		
		# Find the barycentric coordinates
		x = x1*psi1(xi,eta) + x2*psi2(xi,eta) + x3*psi3(xi,eta)
		y = y1*psi1(xi,eta) + y2*psi2(xi,eta) + y3*psi3(xi,eta)
		
		# Area and some slopes
    	area = (x1-x3)*(y2-y1) - (x1-x2)*(y3-y1)
    	dxi_dx, dxi_dy = (y3-y1)/area, -(x3-x1)/area
   		deta_dx, deta_dy = -(y2-y1)/area, (x2-x1)/area

	    # Compute Jacobian matrix
	    j11 = x1*dpsi1_dxi + x2*dpsi2_dxi + x3*dpsi3_dxi
	    j12 = x1*dpsi1_deta + x2*dpsi2_deta + x3*dpsi3_deta
    	j21 = y1*dpsi1_dxi + y2*dspi2_dxi + y3*dpsi3_dxi
    	j22 = y1*dpsi1_deta + y2*dpsi2_deta + y3*dpsi3_deta
    	
    	# And the value of the Jacobian
    	jacobian = abs(j11*j22 - j12*j21)
    	
	    Dpsi[1,1] = dpsi1_dxi*dxi_dx + dpsi1_deta*deta_dx
	    Dpsi[1,2] = dpsi1_dxi*dxi_dy + dpsi1_deta*deta_dy

	    Dpsi[2,1] = dpsi2_dxi*dxi_dx + dpsi2_deta*deta_dx
	    Dpsi[2,2] = dpsi2_dxi*dxi_dy + dpsi2_deta*deta_dy

	    Dpsi[3,1] = dpsi3_dxi*dxi_dx + dpsi3_deta*deta_dx
	    Dpsi[3,2] = dpsi3_dxi*dxi_dy + dpsi3_deta*deta_dy
	    
	    # Assemble the stiffness and load for this element
	    for i = 1:3
	    	load[i] = weight * jacobian * externalForce(x,y) * psi[i]
	    	for j = 1:3
	    		stiff[i,j] = weight * jacobian * stiffness(dpsi,i,j)
	    	end
	    end
	    
	    # Use local stiffness and load to build the global versions
	    for i = 1:3
	    	if in(n[i], b_nodes)
	    		global_load[n[i]] += load[i]
	    		for j = 1:3
	    			if in(n[j], b_nodes)
	    				count += 1
	    				nodes_i[count] = n[i]
	    				nodes_j[count] = n[j]
	    				vals[count] = stiff[i,j]
	    			else
	    				global_load[n[i]] -= stiff[i,j] * field.node_vals[n[j]]
	    			end
	    		end
	    	end
	    end
	end
	
	return nodes_i, nodes_j, vals, global_load
end


