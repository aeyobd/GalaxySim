module calc
export ∇, w

import ForwardDiff: gradient
import LinearAlgebra: norm
import Roots: find_zero

mutable struct Particle
    m::AbstractFloat #static

    x::Vector
    v::Vector

    ρ::Real = nothing
    h::Real = nothing
    W::Vector = nothing
    
    u = nothing
    stars::Vector = nothing
end

function add_density!(p, particles)
    p.ρ = ρ(p, particles)
    p.h = h(p.ρ, p.m)
end

"""
Vectorized form can be used
as divergance operator
"""
function ∇(f, x)
    return gradient(f, x)
end

const ρ_init = 1

function ρ(p, particles)
    return find_zero(ρ_init) do ρ_a
        ρ_a - ρ(p.x, h(ρ_a, p.m), particles)
    end
end

"""
x, m should be arrays, x is 3xN and m is N
"""
function ρ(r, h, particles)
    s = 0
    for p in particles
        s += p.m * W(dist(r,p.x), h)
    end
    return s
end


function dist(ra, rb)
    return norm(ra .- rb)
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
        throw(DomainError(w, "argument must be >= 0"))
    end
end


function ∂W_∂h(r, h)
    q = r/h
    σ=1/24π * r/h^2
    if q>=3
        return 0
    elseif 2<=q< 3
        return σ * (3-q)^4
    elseif 1<=q<2
        return σ *( (3-q)^4 - 6*(2-q)^4 )
    elseif 0 <=q<1
        return σ * ( (3-q)^4 - 6*(2-q)^4 + 15*(1-q)^4 )
    else
        throw(DomainError(w, "argument must be >= 0"))
    end
end


function W(r, h)
    return w(r/h)
end

const η = 1
function h(ρ1, m)
    return η * ∛(m/ρ1)
end

end
