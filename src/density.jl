module Density

export dρ!, dh!, W, F_ab, ∇W, ρ, h, dist


using LinearAlgebra

using ..Constants
using ..Particles



function dρ!(p::Particle, params)
    p.dρ = 0
    for q in p.neighbors
        p.dρ += q.m*(q.v .- p.v) ⋅ ∇W(p, q)
    end
    return p.dρ
end


function dh!(p::Particle, params)
    p.dh = 0
    for q in p.neighbors
        p.dh += -p.h/(3p.ρ) * q.m*(q.v .- p.v) ⋅ ∇W(p, q)
    end
    return p.dh
end


function W(p::Particle, q::Particle)
    return (W(dist(p, q), p.h) + W(dist(p, q), q.h))/2
end

"""
The kernel function
"""
function W(r, h)
    if isnan(r) || isnan(h)
        throw(DomainError("r, h should be real, got $r, $h"))
    end

    return w(abs(r/h))/h^3
end


function ∇W(a::Particle, b::Particle)
    r_hat = normalize(a.x .- b.x)
    return F_ab(a, b) * r_hat
end


function F_ab(a::Particle, b::Particle)
    r = dist(a, b)
    q1 = abs(r/a.h)
    q2 = abs(r/b.h)
    F1 = dw(q1)/a.h^4
    F2 = dw(q2)/b.h^4
    return (F1 + F2)/2
end


"""
Wendland kernel function
see https://academic.oup.com/mnras/article/425/2/1068/1187211
"""
function w(q)
    σ = 495/32π
    c1 = max(0., 1-q)^6
    return σ * c1 * (1 + 6q + 35/3*q^2)
end


function dw(q)
    σ = -495/32π
    c2 = max(0, 1-q)^5
    return σ * 56/3 * c2 * (q + 5*q^2)
end









"""
Calculates the distance
between two positions
"""
function dist(p::Particle, q::Particle)
    return norm(p.x .- q.x)
end


# absolute calculations
#
"""
used for initial calculation of ρ
"""
function ρ(p::Particle, params)
    if length(p.neighbors) < 1
        println("warning, no neighbors")
        return params.rho_min
    end

    soln =  itersolve(p.ρ, [params.rho_min, params.rho_max], params.rho_maxiter, params.tol
                     ) do x
        h1 = h(x, p.m, params.eta)
        return ρ(p, h1, p.neighbors)
    end

    return soln
end



"""
x, m should be arraysh x is 3xN and m is N
ρ = ∑_b m_b W(r_a - r_b; h); (h smoothing length)
"""
function ρ(p0, h::Real, particles)
    s = 0
    for p in particles
        s += p.m * W(dist(p, p0), h)
    end
    s += p0.m * W(0, h)
    return s 
end




h(p::Particle, params) = h(p.ρ, p.m, params.eta)
h(ρ1, m, η) = η*(m/abs(ρ1))^(1/3)


function itersolve(f, x0, range, maxiter=100, tol=1e-3)
    x = x0
    for i in 1:maxiter
        dx = f(x) - x
        x += dx
        if abs(dx/x) < tol
            return x
        end

        # apply range limits
        x = max(x, range[1])
        x = min(x, range[2])
    end
    return x
end






end
