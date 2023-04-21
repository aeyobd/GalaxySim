using Pkg
Pkg.activate(".")

using GalaxySim
using LinearAlgebra
using QuadGK
import GalaxySim.Constants: R_ig, m_p, G
include("init/init.jl")

# cgs
const k_B = 1.3807e-16

function setup()
    params = Params("init/static_eq.toml")

    T = 10
    σ = sqrt(k_B*T/m_p)
    C = σ^2/(2π * G)

    ρ(r) = C / r^2
    M(r) = 4*π * C * r

    R_max = 200pc
    M_tot = M(R_max)
    print("total mass = ")
    println(M_tot/Msun)
    print("mean density = ")
    ρ_mean = M_tot/R_max^3 / m_p
    println(ρ_mean)
    print("v disp = ")
    println(σ/1e5)

    m = M_tot/params.N

    ps = Particle[]

    for i in 1:params.N
        r = rand()*R_max
        x = r * Init.rand_unit_vector()
        v = 2σ * Init.rand_tangent(x)
        p = Particle(x=x, v=v, m=m, T=T, id=i)
        push!(ps, p)
    end

    return ps, params
end


function run()
    ps, params = setup()
    evolve!(ps, params)
end


run()



