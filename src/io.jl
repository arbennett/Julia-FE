## io.jl
#
#  Provides io routines for writing solutions to a file
#

## Write the solution to a files
function writeSolution(mesh::Mesh, u::solution_space)
	f = open("u.dat","w")
	for i = 1:mesh.n_nodes:
		println(f, mesh.nodes[i,1], ",", mesh.nodes[i,2], "," u.node_vals[i])
	end
	close(f)
end

