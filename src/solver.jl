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
type solution_space
    # Build the initial solution space with boundary conditions imposed
    function solution_space(mesh::Mesh, boundaryCondition::Function)
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

#  TODO

