module Gravity

export dv_G!


using ..Constants
using ..Particles
using ..Density


function dv_G!(p::particle, params)
    p.dv_G .= zeros(3)

    for (q, d) in zip(p.neighbors, p.distances)
        p.dv_G .+= -q.m*(dϕ(d, q.h) + dϕ(d, p.h))/2 * (p.x-q.x)/dist
        p.dv_G .+= -q.m/2*(ζ(p)/p.Ω*∇W(p, q) - ζ(q)/q.Ω * ∇W(q, p))
    end

end



function ζ(p::Particle)
    s = 0
    for (q, d) in (p.neighbors, p.distances)
        s += q.m * dϕ_dh(d, p.h)
    end
    s *= dh_dρ(p)
    return s
end




function ϕ(r::F, h)
    q = r/h
    if 0 ≤ q < 1
        1/h*(2/3*q^2 - 3/10*q^4 + 1/10*q^5 - 7/5)
    elseif
        1/h*(4/3*q^2 - q^3 + 3/10*q^4 - 1/30*q^5 - 8/5 + 1/(15q))
    elseif 2 ≤ q
        return -1/r
    else
        throw(DomainError(r, "argument must be real"))
    end
end



function dϕ(r::F, h)
    q = r/h
    if 0 ≤ q < 1
        1/h^2*(4/3*q - 6/5*q^3 + 1/2*q^4)
    elseif
        1/h^2*(8/3*q - 3*q^2 + 6/5*q^3 - 1/6*q^4 - 1/(15*q^2))
    elseif 2 ≤ q
        return -1/r^2
    else
        throw(DomainError(r, "argument must be real"))
    end
end




function dϕ_dh(r::F, h)
    q = r/h
    if 0 ≤ q < 1
        1/h^2*(-2q^2 + 3/2*q^4 - 3/5*q^5 + 7/5)
    elseif
        1/h^2*(-4q^2 + 4q^3 - 3/2*q^4 + 1/5*q^5 + 8/5)
    elseif 2 ≤ q
        return 0
    else
        throw(DomainError(r, "argument must be real"))
    end
end


end
