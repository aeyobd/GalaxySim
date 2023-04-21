module Gravity

export dv_G!, L_grav


using ..Constants
using ..Particles
using ..Density


function dv_G!(p::Particle, params)
    p.dv_G .= zeros(3)
    if !params.phys_gravity
        return p.dv_G
    end

    for (q, d) in zip(p.neighbors, p.distances)
        if q == p
            continue
        end
        p.dv_G .+= -G*q.m*(dϕ(d, q.h) + dϕ(d, p.h))/2 * (p.x-q.x)/d
        p.dv_G .+= -G*q.m/2*(ζ(p)/p.Ω*∇W(p, q) - ζ(q)/q.Ω * ∇W(q, p))
    end

end



function ζ(p::Particle)
    s = 0
    for (q, d) in zip(p.neighbors, p.distances)
        s += q.m * dϕ_dh(d, p.h)
    end
    s *= dh_dρ(p)
    return s
end




"""
Local potential, negative
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
Equal to ϕ'(r, h), >0
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



# >0
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

function L_grav(p::Particle, q::Particle)
    return G/2 * p.m * q.m * (ϕ(p, q) + ϕ(q, p))/2
end

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
