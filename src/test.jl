# Making this work:
#  Sparse matrix construction needs to have every combo of internal nodes
#  Contributions of each combo of internal nodes need to be calculated in
#

require("gmsh.jl")
require("solver.jl")
require("poisson.jl")

#poisson("../meshes/small_test.msh")
#poisson("../meshes/test_mesh.msh")
poisson("../meshes/internal_bc.msh")
