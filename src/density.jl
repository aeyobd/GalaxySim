module Density

export solve_ρ!, Ω, W, ∇W, F_ab, dW_dh


using LinearAlgebra
using Logging

using ..Constants
using ..Particles


function Ω(p::Particle)
    s = 0.
    for q in p.neighbors
        s += q.m * dW_dh(p, q)
    end
    
    return 1 - s * dh_dρ(p)
end



"""
The SPH kernel function between two particles
"""
function W(p::Particle, q::Particle)
    return W(dist(p, q), p.h)
end



"""
The kernel gradient between a and b

Returns a vector pointing towards p

"""
function ∇W(p::Particle, q::Particle)
    r_hat = normalize(q.x .- p.x)
    return F_ab(p, q) * r_hat
end



"""
The scalar part of the (symmetriezed) scalar gradient
(is negative)
"""
function F_ab(a::Particle, b::Particle)
    r = dist(a, b)
    return dW(r, a.h)
end




"""
The kernel function for a distance r and smoothing length h
"""
W(r::F, h) = 1/h^3 * w(r/h)

dW_dh(p::Particle, q::Particle) = dW_dh(dist(p, q), p.h)

dW_dh(r::F, h) = -3/h^4*w(r/h) + 1/h^4*dw(r/h)



"""
cubic spline kernel
see documentation
"""
function w(q)
    σ = 1/π
    if 0 ≤ q < 1
        return σ*(1-3/2*q^2 + 3/4*q^3)
    elseif 1 ≤ q < 2
        return σ/4 * (2-q)^3
    else
        return 0.
    end
end



dW(r::F, h) =  1/h^4*dw(r/h)

function dw(q)
    σ = 1/π
    if 0 ≤ q < 1
        return σ*(-3*q + 9/4*q^2)
    elseif 1 ≤ q < 2
        return -3σ/4 * (2-q)^2
    else
        return 0.
    end
end




"""
Calculates the distance
between two particles
"""
dist(p::Particle, q::Particle) = norm(p.x .- q.x)



"""
solve_ρ!(p, params)

solves for the density of ρ using Newton-Raphson's method
sets ρ, h, and Ω for the particle
"""
function solve_ρ!(p, params)
    h0 = p.h + dh(p, params)

    for i in 1:params.rho_maxiter
        p.ρ = ρ(p, params)
        dh0 = -f(p, params)/df(p, params)
        p.h += dh0

        if abs(dh0/h0) < params.tol
            p.ρ = ρ(p, p.h, p.neighbors)
            p.Ω = Ω(p)
            return p
        end

        p.h = max(p.h, params.h_min)
        p.h = min(p.h, params.h_max)
    end

    @debug "failed to converge"
    return p
end



# h is the kernel density smoothing length
# h(p::Particle, params) = h(p.ρ, p.m, params.eta)
# h(ρ1, m, η) = η*(m/abs(ρ1))^(1/3)

ρ_h(p, params) = p.m/(p.h/params.eta)^(3)


f(p::Particle, params) = ρ_h(p, params) - p.ρ
df(p::Particle, params) = -3p.ρ/p.h * Ω(p)

function ρ(p::Particle, params)
    s = 0.
    for q in p.neighbors
        s += q.m * W(p, q)
    end

    return s 
end




"""
Calculates the change in h (density smoothing length) of particle p (in-place)
"""
dh(p::Particle, params) = dh_dρ(p) * dρ(p, params)
dh_dρ(p) = -p.h/3p.ρ


"""
Calculates the change in density of particle p (in-place)
"""
function dρ(p::Particle, params)
    s = 0
    for q in p.neighbors
        s += q.m*(q.v .- p.v) ⋅ ∇W(p, q)
    end

    return s/p.Ω
end





end
