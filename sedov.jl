# This scripts lets us run a 
#

using Pkg
Pkg.activate(".")

using GalaxySim
using LinearAlgebra

include("init/init.jl")
using .Init

function setup()
    params = Params("init/sedov.toml")

    T = 0
    Th = 10^4

    R_max = 5*pc
    M_tot = 10Msun
    m = M_tot/(params.N)

    ρ_mean = M_tot/R_max^3 / m_p

    ps = Particle[]

    L = 5 # needs to be consistant with N
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



