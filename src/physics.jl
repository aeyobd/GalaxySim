module Physics

export dm_star!
export du_cond!, du_P!
export dv_DM!, dv_P!
export c_sound, pressure, temp


using LinearAlgebra
using Debugger

using ..Particles
using ..Constants
using ..Density

# Lengths are pc, times are years
# Masses in solar mass



"""
dv_DM!(p::Particle, params::Params) -> Vector{3, F}

Acceleration due to dark matter
Modifies the dv_DM attribute of the
particle in-place.

I assume a NFW profile.
"""
function dv_DM!(p::Particle, params)
    if !params.phys_DM
        return zeros(3)
    end
    p.dv_DM .= a_DM(norm(p.x), params) * normalize(p.x)
    return p.dv_DM
end



"""
Helper function for dv_DM!, returns scalar part of acceleration due
to DM
"""
function a_DM(r::F, params)
    if r == 0
        return 0
    end
    G*params.M_tot/params.A_NFW * 1/r^2 * (r/(r + params.Rs) - log(1 + r/params.Rs))
end



"""
Acceleration due to pressure (and viscosity)
"""
function dv_P!(p::Particle, params)
    if !params.phys_pressure
        return 0.
    end

    p.dv_P .= zeros(3)
    for q in p.neighbors
        p.dv_P .+= -q.m .* (
            -p.P/p.ρ^2/p.Ω .* ∇W(p, q)
            .+ q.P/q.ρ^2/q.Ω .* ∇W(q, p)
            .+ Π(p, q, params)
           ) 
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
        p.du_cond += - q.m * (kp+kq) * (p.u-q.u) * (q.x-p.x) ⋅ ∇W(p, q) / (
                                ρ_pq * dist(p, q)^2 + params.eta_visc^2*p.h^2)
    end

    if !isfinite(p.du_cond)
        println("nan du")
        println(p)
        exit()
    end
    return p.du_cond
end


function du_P!(p, params)
    if !params.phys_pressure
        return 0.
    end

    s = 0

    for q in p.neighbors
        if q != p
            s += q.m * (q.v .- p.v) ⋅ ∇W(p, q)
        end
    end

    p.du_P = p.P/p.Ω/p.ρ^2 * s
    return p.du_P
end



"""
Viscosity function between two particles
"""
function Π(p, q, params)
    # calculate the dot product of v and r
    vr = (q.x - p.x) ⋅ (q.v - p.v)
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



"""
Sound speed in the gas of p 
"""
function c_sound(p)
    if p.T >= 0
        return sqrt(5/3 * R_ig * p.T/p.μ)
    else
        println(p)
        throw(DomainError(p.T, "T must be positive!"))
    end
end


"""
Pressure of particle p
(from current energy and μ)
"""
pressure(p) = 2/3 * p.u * p.ρ
temp(p) = 2p.μ/(3R_ig) * p.u



"""
The change in star formation mass (in-place) for p
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
