module density
export ∇, w

import LinearAlgebra: norm
import NLsolve: nlsolve, newton

const ρ_max = 100

function ρ(p, particles, ρ0=100, h0=0.1)

    f1!(F, x) = f!(F, x, p, particles)
    j1!(J, x) = j!(J, x, p, particles)

    solution = nlsolve(f1!, j1!, [ρ0, h0], xtol=2e-4, iterations=200)

    # println("iterations\t", solution.iterations)
    # println("f calls   \t", solution.f_calls)
    # println("J calls   \t", solution.g_calls)
    # println("converged \t", solution.x_converged)

    return solution.zero
end

"""
x, m should be arrays, x is 3xN and m is N
ρ = ∑_b m_b W(r_a - r_b; h); (h smoothing length)
"""
function ρ(p0, h::Real, particles)
    s = 0
    for p in (particles[findall(x->x!==p0, particles)])
        s += p.m * W(dist(p0.x,p.x), h)
    end
    return s
end

function dist(ra, rb)
    return norm(ra .- rb)
end


function W(r, h)
    if isnan(r) || isnan(h)
        return NaN
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
        throw(DomainError(w, "argument must be >= 0"))
    end
end


function ∂W_∂h(r, h)
    if isnan(r) || isnan(h)
        return NaN
    end
    q = abs(r/h)
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
        throw(DomainError(w, "argument must be >= 0"))
    end

    w2 = -3/h * W(r, h)
    return w1 + w2
end

# functions for NLsolve
#

function f!(F, x, p, particles)
    F[1] = ρ(p, x[2], particles) - x[1]
    F[2] = h(x[1], p.m) - x[2]
end

function j!(J, x, p, particles)
    # implement 
    #
    J[1,1] = -1
    J[1,2] = sum([p1.m * ∂W_∂h(dist(p1.x, p.x), x[2]) for p1 in particles])
    J[2,1] = -η*(p.m)^(1/3) * 1/3 * (abs(x[1]))^(-4/3)
    J[2,2] = -1
end

const η = 1.
function h(ρ1, m)
    return η * (m/abs(ρ1))^(1/3)
end



end
