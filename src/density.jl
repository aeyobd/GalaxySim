module Density

export solve_ρ!, find_neighbors!
export Ω, dh_dρ, dist
export W, ∇W, dW_dh, dW


using LinearAlgebra
using Logging
using NearestNeighbors
using Printf

using ..Constants
using ..Particles





"""
get the nearist parlticles 
"""
function find_neighbors!(ps, params)
    x = [p.x for p in ps]
    tree = BallTree(x)
    idxs, dists = knn(tree, x, params.NN, true)

    for i in 1:params.N
        idx = idxs[i]
        ds = dists[i]
        if any(idx .> params.N) || any(idx .<= 0)
            println(ps)
            print(ps[i])
            throw(error("neighbor search failed"))
        end
        filt = ds .> 0

        ps[i].neighbors = ps[idx[filt]]
        ps[i].distances = ds[filt]
    end

    return ps
end



function Ω(p::Particle)
    s = 0.
    for q in p.neighbors
        s += q.m * dW_dh(p, q)
    end
    
    return 1 - s * dh_dρ(p)
end




"""
Calculates the distance
between two particles
"""
dist(p::Particle, q::Particle) = norm(p.x .- q.x)



W(p::Particle, q::Particle) = W(dist(p, q), p.h)
W(r, h) = 1/h^3 * w(r/h)


∇W(p::Particle, q::Particle) = dist(p, q)==0 ? zeros(3) : normalize(q.x .- p.x) * dW(p, q)


# note all derivatives are negative
dW_dh(p::Particle, q::Particle) = dW_dh(dist(p, q), p.h)
dW_dh(r, h) = -3/h^4*w(r/h) - r/h^5*dw(r/h)


dW(p::Particle, q::Particle) = dW(dist(p, q), p.h)
dW(r, h) = 1/h^4*dw(r/h)


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
solve_ρ!(p, params)

solves for the density of ρ using Newton-Raphson's method
sets ρ, h, and Ω for the particle
"""
function solve_ρ!(p, params; save=false)
    if save
        file = open("density.dat", "w")
        println(file, "h,ρ,Ω")
        println(file, "$(p.h), $(p.ρ), $(p.Ω)")
    end

    h0 = p.h + dh(p, params)*p.dt
    h0 = max(h0, params.h_min)
    h0 = min(h0, params.h_max)
    if length(p.neighbors) < 1
        return 
    end

    if save
        file2 = open("density_f.dat", "w")
        println(file2, "h,f,df")
        for i in LinRange(-2, 4, 30)
            h = 10^i * pc
            p.h = h
            x = f(p, params)
            dx = df(p, params)
            println(file2, "$(h/pc), $(x/m_p), $(dx/m_p*pc)")
        end
        close(file2)
    end

    h1 = h0
    p.h = h0


    for i in 1:params.h_maxiter
        p.ρ = ρ(p, params)
        p.Ω = Ω(p)
        p.h -= f(p, params)/df(p, params)/2
        p.h = max(p.h, params.h_min)
        p.h = min(p.h, params.h_max)
        if save
            println(file, "$(p.h/pc), $(p.ρ*m_p), $(p.Ω)")
        end

        if abs((h1-p.h)/h0) < params.tol
            p.ρ = ρ(p, params)
            p.Ω = Ω(p)
            return p
        end
        h1 = p.h
    end

    @debug "failed to converge"

    if save
        close(file)
    end
    return p
end



# h is the kernel density smoothing length
# h(p::Particle, params) = h(p.ρ, p.m, params.eta)
# h(ρ1, m, η) = η*(m/abs(ρ1))^(1/3)



f(p::Particle, params) = ρ_h(p, params) - p.ρ
df(p::Particle, params) = p.Ω/dh_dρ(p)


ρ_h(p, params) = p.m/(p.h/params.eta)^(3)

function ρ(p::Particle, params)
    s = 0.
    for q in p.neighbors
        s += q.m * W(p, q)
    end
    s += p.m * W(0.0, p.h)
    return s 
end



"""
Calculates the change in h (density smoothing length) of particle p (in-place)
"""
dh(p::Particle, params) = dh_dρ(p) * dρ(p, params) * p.dt
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
