module GalaxySim

export main, yr, pc, Msun


using Printf
using Glob
using JLD2


include("constants.jl")
include("particles.jl")
include("init.jl")
include("density.jl")
include("physics.jl")
include("tree.jl")
include("evolve.jl")

using .Constants
using .Particles
using .Init
using .Density
using .Physics
using .Evolve
using .Tree

# Lengths are pc, times are years
# Masses in solar mass



function main(config_file="static_eq.toml")
    if isfile("result.jld")
        rm("result.jld")
    end
    params = get_params(config_file)

    soln = evolve(params)
    println("saving")
    jldsave("result.jld"; u=soln.u, t=soln.t)
    println("saved")
    return soln
end


end
