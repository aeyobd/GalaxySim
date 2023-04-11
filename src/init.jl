module Init

export rand_particles, a_DM

import LinearAlgebra: norm, normalize
import SpecialFunctions: gamma
import Roots: find_zero
import Base: rand

using ..Constants
using ..Particles


const γ = 1.5 # dimensionless

const M_tot = 1e8 * Msun
const M_bary = 1e6 * Msun

const R_virial = 1e3 * pc
const c = 10 # dimensionless
const Rs = R_virial/c
const Rp = Rs

A_NFW = (log(1+c) - c/(1+c))
ρc = M_tot / ( 4π*R_virial^3 * A_NFW)


function ρ_DM(r)
    if r == 0
        return 0
    end

    x = r / Rs
    return ρc / (x * (1+x^2) )
end

v0_virial = √(G*M_tot/R_virial)
function v_virial(r)
    x = r/Rs
    return v0_virial * √( 1/x * (log(1+c*x) - (c*x)/(1+c*x)) /A_NFW )
end

function a_DM(r::Real)
    if r == 0
        return 0
    end
    G * (M_tot/A_NFW) * (r/(r+Rs) - log(1+r/Rs))/r^2
end

function a_DM(x::Vector)
    if norm(x) == 0
        return zeros(3)
    end
    return a_DM(norm(x)) * normalize(x)
end


ρ_bary(r) = (3-γ)/4π * (M_bary*Rp)/(r^γ * (r+Rp)^(4-γ))

# radius sampler helper
∫ρ_bary(r) = r^(3-γ) * (r+Rp)^(γ-3)


function rand_r()
    p = 0.866rand() + 0.001
    find_zero(x->p-∫ρ_bary(x), (0, 2*R_virial))
end


function rand_m(N)
    return M_bary/N*(1 + randn()*0.05)
end

function rand_unit_vector()
    r = randn(3)
    return r / norm(r)
end

function rand_x(R)
    return R .* rand_unit_vector()
end

function rand_speed()
    1 + 0.3*randn()
end

function rand_v(R)
    return rand_speed()*v_virial(R)*rand_unit_vector() 
end

function rand_particle(i=0, N=1)
    R = rand_r()
    m = rand_m(N)

    x = rand_x(R)
    v = rand_v(R)

    return Particle(x=x, v=v, m=m, id=i)
end

function rand_particles(N)
    return [rand_particle(i, N) for i in 1:N]
end


end
