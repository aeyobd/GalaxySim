using Pkg
Pkg.activate(".")

using GalaxySim
using LinearAlgebra
using QuadGK

function setup()
    params = Params("init/disk.toml")

    T = 10

    R_bary = 1000pc
    R_max = R_bary*10

    M_bary = 1000Msun
    m = M_bary/params.N

    ps = Particle[]
    v_vir = sqrt(G * M_bary/R_max)/4

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
    evolve!(ps, params)
    println("finished")
end


run()



