module Density

export ρ, h, W, ∇W, dist


import LinearAlgebra: norm, normalize
using ..Particles
using ..Constants


"""
Calculates the distance
between the positions of two particles
"""
function dist(p::Particle, q::Particle)
    return norm(p.x .- q.x)
end


function ρ(p::Particle)
    return solve(p, p.neighbors, p.distances)
end

h(p::Particle) = h(p.ρ, p.m)

h(ρ1, m) = η * (m/abs(ρ1))^(1/3)

"""
x, m should be arraysh x is 3xN and m is N
ρ = ∑_b m_b W(r_a - r_b; h); (h smoothing length)
"""
function ρ(p0, h::Real, particles, distances)
    s = 0
    for (p, d) in zip(particles, distances)
        s += p.m * W(d, h)
    end
    return s
end


function solve(p::Particle, particles, distances)
    soln =  itersolve(p.ρ, [ρ_min, ρ_max]) do x
        h1 = h(x, p.m)
        return ρ(p, h1, particles, distances)
    end

    return soln
end

function itersolve(f, x0, range, maxiter=100, tol=1e-4)
    x = x0
    for i in 1:maxiter
        dx = f(x) - x
        x += dx
        if abs(dx/x) < tol
            return x
        end
        x = max(x, range[1])
        x = min(x, range[2])
    end
    return x
end



function W(r, h)
    if isnan(r) || isnan(h)
        throw(DomainError("r, h should be real, got $r, $h"))
    end

    return w(abs(r/h))/h^3
end


"""
Welemnd kernel function
see https://academic.oup.com/mnras/article/425/2/1068/1187211
"""
function w(q)
    σ = 495/32π
    c1 = max(0., 1-q)^6
    return σ * c1 * (1 + 6q + 35/3 * q^2)
end

function dw(q)
    σ = 495/32π
    c2 = max(0, 1-q)^5
    return σ * 56/3 * c2 * (q + 5q^2)
end


function ∇W(a::Particle, b::Particle)
    r = dist(a, b)
    if r < 1e-6
        return zeros(3)
    end

    r_hat = normalize(a.x .- b.x)

    h = a.h
    q = abs(r/h)
    return dw(q) * r_hat / h^4
end



end
