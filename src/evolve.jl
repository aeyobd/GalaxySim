module Evolve
export update!, evolve

using ..Physics
using ..Particles
using ..Constants
using ..Density
using ..Init

using NearestNeighbors
using Printf
using JLD2

import DifferentialEquations: ODEProblem, solve, init, TRBDF2
#
#
#
nvar = 3 + 3 + 1 + 1 


function particle_system!(dU, U, p, t)
    masses = p
    ps = matrix_to_particles(U, masses)
    update!(ps)
    dU .= particles_to_matrix(ps)
    dU .-= U
end


function evolve(N=100, t_end=1e9*yr)
    ps0 = rand_particles(N)

    U0 = particles_to_matrix(ps0)
    masses = [p.m for p in ps0]
    @save "mass.jld" masses

    # Create an ODEProblem using the particle_system! function
    prob = ODEProblem(particle_system!, U0, (0., t_end), masses)
    #
    # # Solve the problem
    integrator = init(prob, reltol=1e-3, saveat=LinRange(0., t_end, Nt))

    for i in integrator
        dt = i.t
        t_now = integrator.t[end]
        print_time(t_now, t_end)
    end

    return integrator.sol
end

function matrix_to_particles(U, masses)
    N = length(masses)

    ps =  [Particle(x=U[i, 1:3],
                    v=U[i, 4:6],
                    u=U[i, 7],
                    mstar=U[i, 8],
                    m=masses[i])
            for i in 1:N]

    calc!(ps)
    return ps
end

function particles_to_matrix(ps)
    N = length(ps)
    U = zeros(N, nvar)

    for i in 1:N
        p = ps[i]
        U[i, 1:3] .= p.x
        U[i, 4:6] .= p.v
        U[i, 7] = p.u
        U[i, 8] = p.mstar
    end

    return U
end

function calc!(particles::Vector{Particle})
    idxs, dists = neighbors(particles)

    for i in 1:length(particles)
        p = particles[i]

        p.neighbors = particles[idxs[i][2:end]]
        p.distances = dists[i][2:end]

        p.ρ = ρ(p)
        p.h = h(p)

        if p.u < 0
            p.u = 0
            @debug "negative T warning"
        end
        if p.mstar < 0
            p.mstar = 0
        end

        p.mgas = p.m - p.mstar
        p.ρgas = p.ρ * p.mgas/p.m
        if p.ρgas < 0
            p.ρgas = 0
            p.mgas = 0
        end

        p.T = p.u / (3/2 * R/μ)
        p.P = R/μ * p.T * p.ρ 
    end
    return particles
end

function update!(particles::Vector{Particle})

    for p in particles
        # dynamics
        p.x .+= p.v
        p.v .+= dv_G(p)
        p.v .+= dv_DM(p)

        # Hydrodynamics
        p.v .+= dv_P(p)
        p.u += du_P(p)

        # stellar/astro
        p.mstar += dm_star(p)
        p.u += du_cool(p)
        p.u += du_cond(p)
    end

    return particles
end


function neighbors(particles::Vector)
    x =  hcat(map(p->p.x, particles)...)
    nearest_neighbors(x, 20)
end




"""
get the nearist parlticles 
x is array ndxnp
"""
function nearest_neighbors(x, k=10)
    tree = KDTree(x)
    idxs, dists = knn(tree, x, k, true)
    return idxs, dists
end


function print_time(t, t_end)
    if t < 1e3*yr
        s = @sprintf("%4.0f kyr", t/1e3yr)
    elseif t < 1e6yr
        s = @sprintf("%4.0f yr", t/yr)
    elseif t < 1e9*yr
        s = @sprintf("%4.0f Myr", t/1e6yr)
    else
        s = @sprintf("%4.0f Gyr", t/1e9yr)
    end 

    p = t/t_end*100
    sp = @sprintf("%2.2f %% complete", p)

    print("t =\t" * s * ", " * sp * "\r")
end


end
