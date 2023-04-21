using Pkg
Pkg.activate(".")

using GalaxySim
using LinearAlgebra
using QuadGK
import GalaxySim.Constants: R_ig, m_p, G, k_B

function setup()
    params = Params("src/sedov.toml")


    T = 10
    Th = 1e4

    R_max = 10*pc
    M_tot = 10Msun
    m = M_tot/(params.N)

    ρ_mean = M_tot/R_max^3 / m_p
    println(ρ_mean)

    ps = Particle[]
    push!(ps, Particle(x=zeros(3), v=zeros(3), m=m, T=Th, id=0))

    for i in 1:params.N
        x = R_max * sqrt(rand()) * GalaxySim.Init.rand_unit_vector()
        v = zeros(3)
        push!(ps, Particle(x=x, v=v, m=m, T=T, id=i))
    end

    return ps, params
end


function run()
    ps, params = setup()
    evolve!(ps, params)
    println("finished")
end


run()



