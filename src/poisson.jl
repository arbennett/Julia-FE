## poisson.jl
#
#  IN PROGRESS!
#
#  An example script for solving a pde via the finite element method
#
#  @author Andrew Bennett
#
cwd = dirname(Base.source_path())
cd(cwd)

require("gmsh.jl")

## The stiffness for Poisson's equation
function stiffness(dphi, i, j)
    return dphi[i,1]*dphi[j,1] + dphi[i,2]*dphi[j,2]
end

## Definition of the boundary conditions
function boundaryCondition(x,y)
    return sin(x/pi) + sin(y/pi)
end

## An external force applied to the domain
function externalForce(x,y)
    return 10.0
end

## Solve Poisson's equation
#
function poisson(file_path, size=20)
  t0 = time()
  mesh = gmsh.read(file_path)
  #u = solver.solve(mesh, stiffness, boundaryCondition, externalForce)
  t1 = time()
  
  println("Time elapsed: ", t1-t0, "s")
  # WRITE THE SOLUTION TO A FILE FOR VIS
end

