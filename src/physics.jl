module Physics

export dm_star!
export du_cond!, du_P!
export dv_DM!, dv_P!
export cs!


using LinearAlgebra
using Debugger

using ..Particles
using ..Init
using ..Constants
using ..Density

# Lengths are pc, times are years
# Masses in solar mass




function dv_DM!(p, params)
    if !params.phys_DM
        return zeros(3)
    end
    p.dv_DM .= a_DM(norm(p.x), params) * normalize(p.x)
    return p.dv_DM
end

function a_DM(r::F, params)
    if r == 0
        return 0
    end
    G*params.M_tot/params.A_NFW * 1/r^2 * (r/(r + params.Rs) - log(1 + r/params.Rs))
end

function dv_P!(p::Particle, params)
    if !params.phys_pressure
        return 0.
    end

    p.dv_P .= zeros(3)
    for q in p.neighbors
        p.dv_P .+= -q.m*(p.P/p.ρ^2 + q.P/q.ρ^2 + Π(p, q, params)) * ∇W(p, q) 
    end
    return p.dv_P
end


function du_cond!(p::Particle, params)
    if !params.phys_conduction
        return 0.
    end

    p.du_cond = 0.

    # K = params.K_cond/(1 + params.rho_cond*m_p/p.ρ)
    kp = params.K_cond/p.ρ
    for q in p.neighbors
        kq = params.K_cond/q.ρ
        ρ_pq = (p.ρ + q.ρ)/2
        p.du_cond += - q.m * (kp+kq) * (p.u-q.u) * (p.x-q.x) ⋅ ∇W(p, q) / (
                                ρ_pq * dist(p, q)^2 + params.eta_visc^2*p.h^2)
    end

    return p.du_cond
end


function du_P!(p, params)
    if !params.phys_pressure
        return 0.
    end

    p.du_P = 0

    for q in p.neighbors
        p.du_P += 1/2*q.m*(p.P/p.ρ^2 + q.P/q.ρ^2 + Π(p, q, params)) * (p.v-q.v) ⋅ ∇W(p, q)
    end

    return p.du_P
end


function Π(p, q, params)
    vr = (p.x - q.x) ⋅ (p.v - q.v)
    if !params.phys_visc || vr > 0 
        return 0.
    end

    α = params.alpha_visc
    β = params.beta_visc
    ρ_pq = (p.ρ + q.ρ)/2
    c_pq = (p.c + q.c)/2

    μ_pq = p.h*vr/(dist(p, q)^2 + params.eta_visc^2*p.h^2)

    return ρ_pq \ (-α*c_pq*μ_pq + β*μ_pq^2)
end


function cs!(p)
    if p.T >= 0
        p.c = sqrt(5/3 * R_ig * p.T/p.μ)
    else
        throw(DomainError(p.T, "argument must be positive!"))
    end
    return p.c
end




"""
energy per unit pass
e = u + v^2/2
P=(γ-1)ρ u
γ=5/3
"""
function dm_star!(p::Particle, params)
    if !params.phys_star_formation
        return 0.
    end
    t_ff = sqrt(3π/(32 * G * p.ρ_gas))
    p.dm_star = params.eta_eff * p.m_gas/t_ff
    return p.dm_star
end



end
