# density.jl
#
# Created 03-19-2023
# Revised 04-16-2023 
# Author: Daniel Boyea (boyea.2@osu.edu)
#
# This file contains the routines for
# density estimation.
# See the background.pdf for more detail, 
# but essential this solves a system
# for density ρ and the smoothing length h
# using Newton-Raphson's method.



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
Given a list of particles, finds the 
nearest neighbors to the particle and sets the
Particle.neighbors attribute equal to this list
"""
function find_neighbors!(ps, params)
    # This code uses the Nearest Neighbors package
    # to compute the nearest neighbors
    x = [p.x for p in ps]
    tree = BallTree(x) 
    idxs, dists = knn(tree, x, params.NN, true)

    # idx and dists are lists for each point, 
    # so lets iterate through them
    for i in 1:params.N
        idx = idxs[i]
        ds = dists[i]

        # If the index is nonsense, throw an error
        if any(idx .> params.N) || any(idx .<= 0)
            println(ps)
            print(ps[i])
            throw(error("neighbor search failed"))
        end

        # zeros cause singularities, but most importaintly
        # we don't want to count a particle as its own neighbor
        filt = ds .> 0

        ps[i].neighbors = ps[idx[filt]]
        ps[i].distances = ds[filt]
    end

    return ps
end


"""
The density correction term for particle ρ due to variable
smoothing length
"""
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


"""
W is the kernel function, determining the weights
to use in the calculation of density (and other parameters), i.e.
ρ_p = ∑ m_q W_pq
"""
W(p::Particle, q::Particle) = W(dist(p, q), p.h)
W(r, h) = 1/h^3 * w(r/h) # r is distance, h is smoothing length


"""
The gradient of the density kernel with respect to p
This should point towards p from q.
"""
∇W(p::Particle, q::Particle) = dist(p, q)==0 ? zeros(3) : normalize(q.x .- p.x) * dW(p, q)


"""
Derivative of the density kernel W with respect to h, 
the smoothing length
"""
dW_dh(p::Particle, q::Particle) = dW_dh(dist(p, q), p.h)
dW_dh(r, h) = -3/h^4*w(r/h) - r/h^5*dw(r/h)


"""
Derivative of the density kernel W with respect to r, 
the distance between particles
"""
dW(p::Particle, q::Particle) = dW(dist(p, q), p.h)
dW(r, h) = 1/h^4*dw(r/h)


"""
cubic spline kernel formula
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


"""
derivative of w(q)
"""
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
    if length(p.neighbors) < 1 
        return  # this shouldn't happen
    end

    if save # open up a file to write to if saving for the demonstration
        file = open("density.dat", "w")
        println(file, "h,ρ,Ω,f")
        println(file, "$(p.h/pc),$(p.ρ/m_p),$(p.Ω),$(f(p,params)/m_p)")
    end

    # initial guess of h
    h0 = p.h + dh(p, params)*p.dt

    if save
        # iterate over values of h
        # so we can see what f() looks like 
        # during the density test
        file2 = open("density_f.dat", "w")
        println(file2, "h,f,df")
        for i in LinRange(0.5, 2, 30)
            h = 10^i * pc
            p.h = h
            x = f(p, p.h, params)
            dx = df(p, params)
            println(file2, "$(p.h/pc),$(x/m_p),$(dx/m_p*pc)")
        end
        close(file2)
    end

    h1 = h0
    p.h = h0

    for i in 1:params.h_maxiter
        p.ρ = ρ(p, params)
        p.Ω = Ω(p)

        # newton-raphson's step, we want f() -> 0
        p.h -= f(p, params)/df(p, params)

        if !(params.h_min<=p.h<=params.h_max)
            solve_ρ_bisection!(p, params) # resort to simpler method 
            return 
        end

        if save # record the file
            x = f(p, p.h, params)
            println(file, "$(p.h/pc),$(p.ρ/m_p),$(p.Ω),$(x/m_p)")
        end

        if abs((h1-p.h)/h0) < params.tol # stop when the value converges (relative)
            p.ρ = ρ(p, params)
            p.Ω = Ω(p)
            return p
        end
        h1 = p.h
    end

    if save
        close(file)
    end

    # should have exited by now, resorting to next method
    return solve_ρ_bisection!(p, params)
end


"""
Instead solves for ρ-h by bisection method
"""
function solve_ρ_bisection!(p, params)
    # @info "newton's method failed, trying bisection"

    h_l = params.h_min # lower bound
    h_h = params.h_max # upper bound
    f_l = f(p, h_l, params)
    f_h = f(p, h_h, params)
    if sign(f_l) == sign(f_h) # no zero??
        if abs(f_l) < abs(f_h)
            p.h = h_l
        else 
            p.h = h_h
        end

        p.ρ = ρ(p, p.h, params)
        p.ρ = ρ(p, p.h, params)
        p.Ω = Ω(p)
        return p
    end

    for i in 1:params.h_maxiter
        h_mid = (h_l + h_h)/2
        p.h = h_mid
        p.ρ = ρ(p, p.h, params)
        f_mid = f(p, h_mid, params)

        if sign(f_mid) == 0 || abs(h_l-h_h)/h_mid < params.tol # converged
            p.ρ = ρ(p, params)
            p.Ω = Ω(p)
            return p
        elseif sign(f_mid) == sign(f_l)
            h_l = h_mid
            f_l = f_mid
        else 
            h_h = h_mid
            f_h = f_mid
        end
    end

    p.Ω = Ω(p)
    return p
end



# h is the kernel density smoothing length
# h(p::Particle, params) = h(p.ρ, p.m, params.eta)
# h(ρ1, m, η) = η*(m/abs(ρ1))^(1/3)


"""
f is the function who's root is found above,
it is the difference between the value of ρ
based on h and based on the neighbors of p.
"""
f(p::Particle, params) = ρ_h(p, params) - p.ρ
f(p::Particle, h, params) = ρ_h(p, h, params) - ρ(p, h, params)

df(p::Particle, params) = p.Ω/dh_dρ(p)


"""
the value of ρ given only h and the current mass
"""
ρ_h(p, params) = p.m/(p.h/params.eta)^(3)
ρ_h(p, h, params) = p.m/(h/params.eta)^(3)


"""
The density of the particle p 
using smoothing length h
"""
function ρ(p::Particle, h, params)
    s = 0.
    for (q, d) in zip(p.neighbors, p.distances)
        s += q.m * W(d, h)
    end
    s += p.m * W(0.0, h)
    return s 
end


"""
The density of the particle p 
"""
function ρ(p::Particle, params)
    s = 0.
    for (q, d) in zip(p.neighbors, p.distances)
        s += q.m * W(d, p.h)
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
