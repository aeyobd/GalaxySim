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
include("evolve.jl")

using .Constants
using .Particles
using .Init
using .Density
using .Physics
using .Evolve

# Lengths are pc, times are years
# Masses in solar mass



function main(N=100, t_end=1e6yr)
    if isfile("result.jld")
        rm("result.jld")
    end

    soln = evolve(N, t_end)
    println("saving")
    jldsave("result.jld"; u=soln.u, t=soln.t)
    println("saved")
    return soln
end


end
