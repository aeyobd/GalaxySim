using Pkg
Pkg.activate(".")

using GalaxySim
using LinearAlgebra
using QuadGK
import GalaxySim.Constants: R_ig, m_p, G

# cgs
const k_B = 1.3807e-16

function setup()
    params = GalaxySim.Constants.get_params("src/static_eq.toml")

    T = params.T0
    σ = sqrt(k_B*T/m_p)
    C = σ^2/(2π * G)

    ρ(r) = C / r^2
    M(r) = 4*π * C * r

    R_max = 100*pc
    M_tot = M(R_max)
    println("total mass")
    println(M_tot/Msun)
    println(σ/1e5)

    m = M_tot/params.N

    ps = Particle[]

    for i in 1:params.N
        r = rand()*R_max
        x = r * GalaxySim.Init.rand_unit_vector()
        v = 2σ * GalaxySim.Init.rand_tangent(x)
        push!(ps, Particle(x=x, v=v, m=m, T=T))
    end

    return ps, params
end


function run()
    ps, params = setup()
    GalaxySim.Evolve.evolve!(ps, params)
    println("finished")
end


run()



