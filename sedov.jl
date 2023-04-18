using Pkg
Pkg.activate(".")

using GalaxySim
using LinearAlgebra
using QuadGK
import GalaxySim.Constants: R_ig, m_p, G, k_B

function setup()
    params = GalaxySim.Constants.get_params("src/sedov.toml")

    T = params.T0
    Th = 10^5

    R_max = 5*pc
    M_tot = 1000Msun
    m = M_tot/(params.N)

    ρ_mean = M_tot/R_max^3 / m_p
    println(ρ_mean)

    ps = Particle[]
    push!(ps, Particle(x=zeros(3), v=zeros(3), m=m, T=Th, id=0))

    Ni = 3
    N = (2Ni)^3
    jitter = 0.1*R_max/Ni

    for i in -Ni:Ni, j in -Ni:Ni, k in -Ni:Ni
        if i==j==k==0
            continue
        end
        x = [i, j, k] * R_max/Ni + jitter*randn(3)
        v = zeros(3)
        push!(ps, Particle(x=x, v=v, m=m, T=T, id=2))
    end

    return ps, params
end


function run()
    ps, params = setup()
    GalaxySim.Evolve.evolve!(ps, params)
    println("finished")
end


run()



