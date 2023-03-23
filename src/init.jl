module init

import LinearAlgebra: norm, normalize
import SpecialFunctions: gamma
import Roots: find_zero
using Distributions
import Base: rand

include("particle.jl")

const G = 4.4985e-15


M_tot = 1e8
M_bary = 1e6 #M_sun
c = 10
R_virial = 1e3
Rs = R_virial/c

A_NFW = (log(1+c) - c/(1+c))
ρc = M_tot / ( 4π*R_virial^3 * A_NFW)

function ρ_DM(r)
    x = r / Rs
    return ρc / (x * (1+x^2) )
end

v0_virial = √(G*M_tot/R_virial)
function v_virial(r)
    x = r/Rs
    return v0_virial * √( 1/x * (log(1+c*x) - (c*x)/(1+c*x)) /A_NFW )
end

function a_DM(r::Real)
    G * (M_tot/A_NFW) * (r/(r+Rs) - log(1+r/Rs))/r^2
end

function a_DM(x::Vector)
    return a_DM(norm(x)) * normalize(x)
end


const Rp = Rs
const γ = 1.5
ρ_bary(r) = (3-γ)/4π * (M_bary*Rp)/(r^γ * (r+Rp)^(4-γ))

# radius sampler helper
∫ρ_bary(r) = r^(3-γ) * (r+Rp)^(γ-3)


function rand_r()
    p = 0.866rand() + 0.001
    find_zero(x->p-∫ρ_bary(x), (0, 2*R_virial))
end

const N_particles = 1000

function rand_m()
    return rand(Normal(M_bary/N_particles, M_bary/N_particles * 0.05))
end

function rand_unit_vector()
    θ = rand() * 2π
    ϕ = acos(2*rand() - 1)

    x = sin(ϕ) * cos(θ)
    y = cos(ϕ) * cos(θ)
    z = sin(θ)
    return [x,y,z]
end

function rand_x(R)
    return R .* rand_unit_vector()
end

function rand_speed()
    rand(Normal(1, 0.3))
end

function rand_v(R)
    return rand_speed()*v_virial(R)*rand_unit_vector() 
end

function rand_particle()
    R = rand_r()
    m = rand_m()

    x = rand_x(R)
    v = rand_v(R)

    return particle.Particle(x=x, v=v, m=m)
end

function rand_particles()
    return [rand_particle() for _ in 1:N_particles]
end


end
