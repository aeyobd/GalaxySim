module density
export ∇, w

import LinearAlgebra: norm
import NLsolve: nlsolve, newton
using IntervalRootFinding
using StaticArrays


const ρ_max = 1e-19
const ρ_min = 1e-26
const η = 0.1



function itersolve(f, x0, range, maxiter=100, tol=1e-8)
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


function solve(p, particles, distances)
    function fi(x)
        hi = h(x, p.m)
        return ρ(p, hi, particles, distances)
    end

    soln = itersolve(fi, p.ρ, [ρ_min, ρ_max])
    
    return soln
end


function ρ(p, particles, distances)
    return solve(p, particles, distances)
end

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

function ρ(p0, h::Real, particles)
    s = 0
    for p in particles
        d = dist(p.x, p0.x)
        s += p.m * W(d, h)
    end
    return s
end

function dist(ra, rb)
    return norm(ra .- rb)
end



function W(r, h)
    if isnan(r) || isnan(h)
        throw(DomainError("r, h should be real, got $r, $h"))
    end
    return w(abs(r/h))/h^3
end

"""
the weight function for 
calculation ρ(q) where q=r/h


normalized to 1 for a sphere of radius q=1
"""
function w(q)
    σ=1/120π #3-d normalization
    if q>=3
        return 0
    elseif 2<=q< 3
        return σ * (3-q)^5
    elseif 1<=q<2
        return σ *( (3-q)^5 - 6*(2-q)^5 )
    elseif 0 <=q<1
        return σ * ( (3-q)^5 - 6*(2-q)^5 + 15*(1-q)^5 )
    else
        throw(DomainError(w, "argument must be >= 0, got $q"))
    end
end


function ∂W_∂h(r, h)
    q = r/h
    σ = 1/24π * r/h^5

    if q>=3
        w1 = 0
    elseif 2<=q< 3
        w1 = σ * (3-q)^4
    elseif 1<=q<2
        w1 = σ *( (3-q)^4 - 6*(2-q)^4 )
    elseif 0 <=q<1
        w1 = σ * ( (3-q)^4 - 6*(2-q)^4 + 15*(1-q)^4 )
    else
        throw(DomainError(w, "argument must be >= 0, got r=$r, h=$h"))
    end

    w2 = -3/h * W(r, h)
    return w1 + w2
end

function ∇W(a, b)
    r_vec = b.x - a.x
    r = norm(r_vec)
    if r < 1e-6
        return [0.,0.,0.]
    end

    h = a.h

    q = abs(r/h)
    σ = -1/24π * 1/h^4

    if q < 1
        w1 = σ * ( (3-q)^4 - 6*(2-q)^4 + 15*(1-q)^4 )
    elseif 1<=q<2
        w1 = σ *( (3-q)^4 - 6*(2-q)^4 )
    elseif 2<=q< 3
        w1 = σ * (3-q)^4
    elseif 3 <= q
        w1 = 0
    else
        throw(DomainError(w, "argument must be >= 0, got r=$r, h=$h, $r_vec\n $a, $b"))
    end

    return w1 * r_vec /r
end
# functions for NLsolve
#

function f!(F, x, p, particles, distances)
    F[1] = ρ(p, x[2], particles, distances) - x[1]
    F[2] = h(x[1], p.m) - x[2]
end


function j!(J, x, p, particles, distances)
    # implement 
    #
    J[1,1] = -1
    J[1,2] = sum([p1.m * ∂W_∂h(d, x[2]) for (p1, d) in zip(particles, distances)])
    J[2,1] = -η*(p.m)^(1/3) * 1/3 * (abs(x[1]))^(-4/3)
    J[2,2] = -1
end

function h(ρ1, m)
    return η * (m/abs(ρ1))^(1/3)
end



end
