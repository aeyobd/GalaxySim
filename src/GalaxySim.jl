module GalaxySim

export main


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
    rm("result.jld")

    soln = evolve(N, t_end)
    println("saving")
    jldsave("result.jld"; u=soln.u, t=soln.t)
    println("saved")
end

function print_time(t, t_end)
    if t < 1e6*yr
        s = @sprintf("%4.0f yr", t/yr)
    elseif t < 1e9*yr
        s = @sprintf("%4.0f Myr", t/1e6yr)
    else
        s = @sprintf("%4.0f Gyr", t/1e9yr)
    end 

    p = t/t_end*100
    sp = @sprintf("%2.2f %% complete", p)

    print("t =\t" * s * ", " * sp * "\r")
end


end
