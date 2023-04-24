# Gravity
#
# methods to calculate gravitational force
#
# Created 11-04-2023
# Updated to use SPH consistant gravity 20-04-2023
#
# Author Daniel Boyea (boyea.2@osu.edu)


module Gravity

export dv_G!, L_grav, dv_DM!


using ..Constants
using ..Particles
using ..Density




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
Acceleration due to gravity
This is calculated following 
Price and Monaghan 2007.
"""
function dv_G!(p::Particle, params)
    p.dv_G .= zeros(3)

    for (q, d) in zip(p.neighbors, p.distances)
        if q == p # skip over self
            continue
        end
        p.dv_G .+= -G*q.m*(dϕ(d, q.h) + dϕ(d, p.h))/2 * (p.x-q.x)/d
        p.dv_G .+= -G*q.m/2*(ζ(p)/p.Ω*∇W(p, q) - ζ(q)/q.Ω * ∇W(q, p))
    end

end



"""
Helper function to calculate gravity
"""
function ζ(p::Particle)
    s = 0
    for (q, d) in zip(p.neighbors, p.distances)
        s += q.m * dϕ_dh(d, p.h)
    end
    s *= dh_dρ(p)
    return s
end




"""
Local softened gravitational potential (negative)
"""
function ϕ(r::F, h)
    q = r/h
    if 0 ≤ q < 1
        1/h*(2/3*q^2 - 3/10*q^4 + 1/10*q^5 - 7/5)
    elseif 1 ≤ q < 2
        1/h*(4/3*q^2 - q^3 + 3/10*q^4 - 1/30*q^5 - 8/5 + 1/(15q))
    elseif 2 ≤ q
        return -1/r
    else
        throw(DomainError(q, "argument must be real"))
    end
end



"""
Equal to ϕ'(r, h), positive, gravitational force kernel
"""
function dϕ(r::F, h)
    q = r/h
    if 0 ≤ q < 1
        1/h^2*(4/3*q - 6/5*q^3 + 1/2*q^4)
    elseif 1 ≤ q < 2
        1/h^2*(8/3*q - 3*q^2 + 6/5*q^3 - 1/6*q^4 - 1/(15*q^2))
    elseif 2 ≤ q
        return 1/r^2
    else
        throw(DomainError(q, "argument must be real"))
    end
end



"""
derivative of ϕ with respect to h
"""
function dϕ_dh(r::F, h)
    q = r/h
    if 0 ≤ q < 1
        1/h^2*(-2q^2 + 3/2*q^4 - 3/5*q^5 + 7/5)
    elseif 1 ≤ q < 2
        1/h^2*(-4q^2 + 4q^3 - 3/2*q^4 + 1/5*q^5 + 8/5)
    elseif 2 ≤ q
        return 0
    else
        throw(DomainError(q, "argument must be real"))
    end
end


ϕ(p::Particle, q::Particle) = ϕ(dist(p, q), p.h)

"""
Pairwise gravitational potential (halved)
"""
function L_grav(p::Particle, q::Particle)
    return G/2 * p.m * q.m * (ϕ(p, q) + ϕ(q, p))/2
end

"""
Gravitational potential of particles ps
"""
function L_grav(ps, params)
    if !params.phys_gravity
        return 0.
    end
    s = 0
    for p in ps
        for q in ps
            s += L_grav(p, q)
        end
    end
    return s
end

end
