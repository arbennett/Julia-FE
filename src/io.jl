## io.jl
#
#  Provides io routines for writing solutions to a file
#

require("gmsh.jl")

## Write the solution to a files
function write_solution(coords::Array{Float64,2}, u::Array{Float64,1})
    println("Writing solution to file....")
    f = open("solution.dat", "w")
    for i=1:length(u)
        println(f, "$(coords[i,1]), $(coords[i,2]), $(u[i])")
    end
end

## Draw a quick plot
function plot_solution(coords::Array{Float64,2}, u::Array{Float64,1})
    x = coords[:,1]
    y = coords[:,2]
    plot_trisurf(x,y,u) 
end

