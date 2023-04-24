# sedov.jl
#
# This runs a sedov-taylor-von-newman blast wave
#
# Created 18-04-2023
# Author Daniel Boyea (boyea.2@osu.edu)
#
# Unfortunantly, this doesn't work well with the simulation as is
# so there is some incorrect physics in the viscosity prescriptions
# or I am not running this at high enough resolution
#
# The particles start off in a mostly regular grid, and the center
# particle is set at a high temperature.

using Pkg
Pkg.activate(".")

using GalaxySim
using LinearAlgebra

include("init/init.jl")
using .Init

function setup()
    params = Params("init/sedov.toml")

    T = 1e-4
    Th = 10^4

    R_max = 2*pc
    M_tot = 1Msun
    m = M_tot/(params.N)

    ρ_mean = M_tot/R_max^3 / m_p

    ps = Particle[]

    # create a grid of particles this time
    # but add a little scatter as well
    L = 3 
    for i in -L:L, j in -L:L, k in -L:L
        if i==j==k==0
            push!(ps, Particle(x=zeros(3), v=zeros(3), m=M_tot, T=Th, id=0))
        else
            x = R_max/2L * ( [i,j,k] .+ 0.3*randn(3))
            v = zeros(3)
            push!(ps, Particle(x=x, v=v, m=m, T=T))
        end
    end

    return ps, params
end


function run()
    ps, params = setup()
    evolve!(ps, params)
end


run()



