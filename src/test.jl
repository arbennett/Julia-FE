# Making this work:
#  Sparse matrix construction needs to have every combo of internal nodes
#  Contributions of each combo of internal nodes need to be calculated in
#

require("gmsh.jl")
require("solver.jl")
require("poisson.jl")

println("Performing regular test...")
poisson("../meshes/test_mesh.msh")
println()

println("Performing internal BC test...")
poisson("../meshes/internal_bc.msh")
println()

println("Performing small test...")
poisson("../meshes/small_test.msh")
println()

