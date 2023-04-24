# physics
#
# contains all other hydrodynamic/star formation physics
#
# Created 03-22-2023
# Author Daniel Boyea (boyea.2@osu.edu)



module Physics

export dm_star!
export du_cond!, du_P!
export dv_DM!, dv_P!
export du_visc!, dv_visc!
export c_sound, pressure, temp


using LinearAlgebra
using Logging

using ..Particles
using ..Constants
using ..Density


# pressure and temperature of the particles
pressure(p) = 2/3 * p.u * p.ρ
temp(p) = 2p.μ/(3R_ig) * p.u


"""
Acceleration due to pressure 
"""
function dv_P!(p::Particle, params)
    p.dv_P .= zeros(3)
    for q in p.neighbors
        p.dv_P .+= -q.m .* (
            -p.P/p.ρ^2/p.Ω .* ∇W(p, q)
            .+ q.P/q.ρ^2/q.Ω .* ∇W(q, p)
           ) 
    end
    return p.dv_P
end



"""
Change in energy due to pressure
"""
function du_P!(p, params)
    s = 0

    for q in p.neighbors
        s += q.m * (q.v .- p.v) ⋅ ∇W(p, q)
    end

    p.du_P = p.P/p.Ω/p.ρ^2 * s
    return p.du_P
end



"""
Sound speed in the gas of p 
"""
function c_sound(p)
    if p.T >= 0
        return sqrt(5/3 * R_ig * p.T/p.μ)
    else
        @debug p
        throw(DomainError(p.T, "T must be positive!"))
    end
end




"""
Signal speed between two particles
"""
function v_sig(p::Particle, q::Particle, params)
    v_r = (q.v .- p.v) ⋅ normalize(q.x .- p.x)

    if v_r <= 0
        return 1/2 * (p.c + q.c - params.beta*v_r)
    else # the particles are moving away from eachother
        return 0.
    end
end





"""
change in energy due to viscosity
"""
function du_visc!(p::Particle, params)
    p.du_visc = 0.

    for q in p.neighbors
        ρ_pq = (p.ρ + q.ρ)/2
        p.du_visc += -q.m/2ρ_pq*params.alpha*v_sig(p, q, params)^2*(
            dW(p, q) + dW(q, p))/2
    end

    return p.du_visc
end



"""
Acceleration due to viscosity
"""
function dv_visc!(p::Particle, params)
    p.dv_visc .= zeros(3)
    for q in p.neighbors
        ρ_pq = (p.ρ + q.ρ)/2
        v_r = (p.v - q.v) ⋅ normalize(p.x - q.x)
        p.dv_visc .+= -q.m/ρ_pq * v_sig(p, q, params) .* v_r .* ( ∇W(p, q)  .- ∇W(q, p) )
    end

    return p.dv_visc
end


function K(p::Particle, params)
    return params.K0 / (1 + params.n_K/(p.ρ/(m_p*p.μ)) )
end


""" Change in thermal energy due to heat conduction """
function du_cond!(p::Particle, params)
    p.du_cond = 0.

    kp = K(p, params) * p.ρ/(p.μ*m_p)
    for q in p.neighbors
        kq = K(q, params)* p.ρ/(p.μ*m_p)
        k = 4*kp*kq/(kp + kq)
        p.du_cond += q.m/(p.ρ + q.ρ) * k * (p.T - q.T) * dW(p, q)
    end
    # params.alpha*v_sig_u(p, q, params) * (p.u - q.u)
    return p.du_cond
end



"""
Energy signal speed between two particles
"""
function v_sig_u(p::Particle, q::Particle, params)
    return sqrt(2*abs(p.P - q.P)/(p.ρ + q.ρ))
end



"""
The change in star formation mass (in-place) for p
"""
function dm_star!(p::Particle, params)
    t_ff = sqrt(3π/(32 * G * p.ρ_gas))

    p.dm_star = params.eta_eff * p.m_gas/t_ff
end



end
