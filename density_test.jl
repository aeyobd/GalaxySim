# density_test.jl 
#
# initial_version
# author: Daniel Boyea (boyea.2@osu.edu)
#
# This file tests if the density estimator is good
# The script creates particals randomly spaced
# in a box, then calculates the density of a 
# particle placed at the center. 
#
# The script writes the data to `density.dat` 
# and `density_f.dat`, each of which are csv files

using Pkg
Pkg.activate(".")

using GalaxySim
using LinearAlgebra
using QuadGK
using Random

import GalaxySim.Constants: m_p # proton mass
include("init/init.jl")

const k_B = 1.3807e-16 # kelven-boltzman constant (cgs)


function setup()
    params = Params("init/static_eq.toml")

    Random.seed!(127)

    R_max = 100*pc
    M_tot = 1000Msun
    ρ_mean = M_tot/R_max^3 / m_p / 2^3
    println("ρ_mean = ", ρ_mean)

    m = M_tot/params.N # per-particle mass

    ps = Particle[]
    for i in 1:params.N
        if i == 1 # let the first particle be at the origin
            x = zeros(3) 
            v = zeros(3)
        else # make a random particle in the box (-R_max, R_max)^3
            x = R_max*(2*randn(3) .- 1)
            v = 2σ * Init.rand_tangent(x)
        end
        p = Particle(x=x, v=v, m=m, T=T, id=i)
        push!(ps, p)
    end

    return ps, params
end


function run()
    ps, params = setup()

    # need to find nearest neighbors to calculate density
    GalaxySim.Density.find_neighbors!(ps, params)
    for p in ps # other particles help calculate the current particles density
        GalaxySim.Density.solve_ρ!(p, params)
    end
    GalaxySim.Density.solve_ρ!(ps[1], params, save=true)

    print("ρ = ")
    println(ps[1].ρ/m_p)
    print("h = ")
    println(ps[1].h/pc)
    println("data saved at density.dat and density_f.dat")
end


run()



