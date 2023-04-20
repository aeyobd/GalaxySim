module GalaxySim

export main, yr, pc, Msun, G
export Particle


using Printf
using Glob
using JLD2


include("constants.jl")
include("particles.jl")
include("gal_files.jl")
include("density.jl")
include("tree.jl")
include("init.jl")
include("physics.jl")
include("evolve.jl")

using .Constants
using .Particles
using .GalFiles
using .Density
using .Tree
using .Init
using .Physics
using .Evolve

# Lengths are pc, times are years
# Masses in solar mass



function main(config_file="static_eq.toml")
    dir = @__DIR__
    file = joinpath(dir, config_file)
    params = get_params(file)
    soln = evolve(params)
end


end
