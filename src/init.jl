module Init

export rand_particles

import LinearAlgebra: norm, normalize
import Roots: find_zero
import Base: rand

using StaticArrays

using ..Constants
using ..Particles



"""
Generates a list of random particles using
the configuration stored in params

Parameters
----------
params: the parameters dictionary
    in the dictionary:

    N: the number of particles (int)
    R_virial: the virial radius
    R_bary: the radius of the baryonic profile
    c: NFW profile compactness parameter
    gamma: baryonic compactness parameter

    M_tot: the total (DM) mass 
    M_bary: the total baryonic (gas) mass 

Returns
    particles: Vector{Particle}

"""
function rand_particles(params)
    return [rand_particle(i, params) for i in 1:params.N]
end



function rand_particle(i, params)
    r = rand_r(params)
    m = rand_m(params)
    x = rand_x(r, params)
    v = rand_v(r, params)

    return Particle(x=x, v=v, m=m, id=i)
end


# radius sampler helper
∫ρ_bary(r::F, params) = r^(3 - params.gamma) * (r + params.R_bary)^(params.gamma-3)
# ρ_bary(r) = (3-γ)/4π * (M_bary*Rp)/(r^γ * (r+Rp)^(4-γ))
#
function rand_r(params)
    p = 0.98*rand() + 0.001
    # use the integrated distribution to get parameter
    find_zero(x->p - ∫ρ_bary(x, params), [1e-3*params.R_bary, 1e3*params.R_bary])
end


function rand_m(params)
    return params.M_bary/params.N * (1 + randn()*params.sigma_M)
end


function rand_x(r::F, params)
    return r .* rand_unit_vector()
end


function rand_unit_vector()
    r = randn(3)
    return r / norm(r)
end


function rand_v(r::F, params)
    return params.v_0*(1 + params.sigma_v*randn()) * v_virial(r, params) * rand_unit_vector() 
end


function v_virial(r::F, params)
    v0_virial = √(G*params.M_tot/params.R_virial)
    x = r/params.Rs
    return v0_virial * √( 1/x * (log(1 + params.c*x) - (params.c*x)/(1 + params.c*x)) 
                         /params.A_NFW )
end

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
