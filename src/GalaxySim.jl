module GalaxySim

export main, yr, pc, Msun, G
export Particle


using Printf
using Glob
using JLD2


include("constants.jl")
include("particles.jl")
include("gal_files.jl")
include("init.jl")
include("density.jl")
include("physics.jl")
include("tree.jl")
include("evolve.jl")

using .Constants
using .Particles
using .GalFiles
using .Init
using .Density
using .Physics
using .Evolve
using .Tree

# Lengths are pc, times are years
# Masses in solar mass



function main(config_file="static_eq.toml")
    dir = @__DIR__
    file = joinpath(dir, config_file)
    params = get_params(file)
    soln = evolve(params)
end


end
