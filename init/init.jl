module Init

export rand_particles, rand_unit_vector

using LinearAlgebra
import Roots: find_zero
import Base: rand
using Random

using StaticArrays

using GalaxySim.Constants
using GalaxySim.Particles





"""
    rand_particles(params)

Generates a list of random particles using the configuration stored in params

See CONFIG.md for more information about params

Returns
    particles: Vector{Particle}

"""
function rand_particles(params)
    set_seed!(params)
    return [rand_particle(i, params) for i in 1:params.N]
end



""" Sets the RNG seed if seed is set in params.  """
function set_seed!(params)
    if params.seed > 0
        Random.seed!(params.seed)
    end
end



""" Creates a random particle with identifier i"""
function rand_particle(i, params)
    r = rand_r(params)
    m = rand_m(params)
    x = rand_x(r, params)
    v = rand_v(r, params)

    return Particle(x=x, v=v, m=m, id=i)
end



""" The initial baryonic density function at r """
function ρ_bary(r::F, R_bary, M_bary, γ) 
    C = (3-γ)/4π * (M_bary*R_bary)

    return C/(r^γ * (r + R_bary)^(4 - γ))
end



""" The integral of `ρ_bary(r, params)` to calculate distribution in r"""
function ∫ρ_bary(x, R_bary, γ)
    x^(3 - γ) * (x + R_bary)^(γ-3)
end



""" A random radius of a particle, sampled from the distribution `ρ_bary(r, params)`"""
function rand_r(R_bary, γ=1.5)
    p = 0.98*rand() + 0.001
    # use the integrated distribution to get parameter
    find_zero(x->p - ∫ρ_bary(x, R_bary, γ), [1e-3R_bary, 1e3R_bary])
end



""" A random particle mass """
function rand_m(params)
    return params.M_bary/params.N * (1 + randn()*params.sigma_M)
end




""" A random particle position vector """
function rand_x(r::F, params)
    return r .* rand_unit_vector()
end




""" A random 3-dimensional unit vector """
function rand_unit_vector()
    r = randn(3)
    return r / norm(r)
end




""" A random 3-dimensional vector tangent to x """
function rand_tangent(x)
    y = randn(3)
    tang = x × y
    return normalize(tang)
end




""" A random 3-dimensional velocity at `r` """
function rand_v(r::F, params)
    return params.v_0*(1 + params.sigma_v*randn()) * v_virial(r, params) * rand_unit_vector() 
end




""" The virial velocity at r """
function v_virial(r::F, config)
    v0_virial = √(G*config.M_tot/config.R_virial)
    x = r/config.R_s
    return v0_virial * √( 1/x * (log(1 + config.c*x) - (config.c*x)/(1 + config.c*x)) 
                         /config.A_NFW )
end


# The dark matter distribution
# function ρ_DM(r::F, params)
#     ρ_c = params["M_tot"] / ( 4π*params["R_virial"]^3 * params["A_NFW"])
# 
#     if r == 0
#         return 0
#     end
# 
#     x = r / params.Rs
#     return ρc / (x * (1+x^2) )
# end


end
