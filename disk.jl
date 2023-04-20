using Pkg
Pkg.activate(".")

using GalaxySim
using LinearAlgebra
using QuadGK
import GalaxySim.Constants: R_ig, m_p, G, k_B

function setup()
    params = GalaxySim.Constants.get_params("src/disk.toml")

    T = params.T0

    R_max = params.R_bary*10
    m = params.M_bary/(params.N)

    ps = Particle[]
    v_vir = sqrt(G * params.M_bary/R_max)/4

    for i in 1:params.N
        vec = randn(2)
        push!(vec, 0.1*randn())
        x = vec * R_max 
        v = [0,0,1.0] × vec * v_vir * norm(vec)
        push!(ps, Particle(x=x, v=v, m=m, T=T, id=i))
    end

    return ps, params
end


function run()
    ps, params = setup()
    GalaxySim.Evolve.evolve!(ps, params)
    println("finished")
end


run()



