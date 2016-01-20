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

