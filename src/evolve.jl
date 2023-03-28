module Evolve
export update!, evolve

using ..Physics
using ..Particles
using ..Constants
using ..Density
using ..Init

using NearestNeighbors
using Printf
using Glob

import DifferentialEquations: ODEProblem, solve, init
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


function evolve(N=100, t_end=1e6*yr)
    ps0 = rand_particles(N)

    U0 = particles_to_matrix(ps0)
    masses = [p.m for p in ps0]

    # Create an ODEProblem using the particle_system! function
    prob = ODEProblem(particle_system!, U0, (0., t_end), masses)
    #
    # # Solve the problem
    integrator = init(prob)

    for i in integrator
        dt = i.t
        t_now = integrator.t[end] / yr
        print("t = $t_now \r")
    end

    return integrator.sol
end


function open_files(N)
    # delete files
    for file in glob("data/mass*.dat")
        rm(file)
    end

    files = Vector()
    for i in 1:N
        fname = "data/mass$i.dat"
        touch("data/mass$i.dat")
        f = open("data/mass$i.dat", "w")
        println(f, "x1,x2,x3,v1,v2,v3,ρ,h,T,ms,f")
        push!(files, f)
    end
    println("Opened files")

    return files
end

function record_masses(files, masses)
    for (file, mass) in zip(files, masses)
        x1 = mass.x[1]/pc
        x2 = mass.x[2]/pc
        x3 = mass.x[3]/pc
        v1 = mass.v[1] /100_000 # km/s
        v2 = mass.v[2] /100_000 
        v3 = mass.v[3] /100_000 
        ρ = mass.ρ / Msun * pc^3
        h = mass.h / pc
        T = mass.T
        mstar = mass.mstar/Msun
        fstar = mass.mstar/mass.m
        write(file, "$x1,$x2,$x3,$v1,$v2,$v3,$ρ,$h,$T,$mstar,$fstar\n")
        flush(file)
    end
end

function close_files(files)
    close.(files)
end


function matrix_to_particles(U, masses)
    N = length(masses)

    ps =  [Particle(x=U[i, 1:3],
                    v=U[i, 4:6],
                    u=U[i, 7],
                    mstar=U[i, 8],
                    m=masses[i])
            for i in 1:N]
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

function next_ps(particles::Vector{Particle})
    p2 = deepcopy(particles)
    update!(p2)
    return p2
end

function update!(particles::Vector{Particle})
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

    for p in particles
        p.x .+= p.v
        p.v .+= dv_G(p)
        p.v .+= dv_P(p)

        p.mstar += dm_star(p)
        p.u += du_cool(p)
        p.u += du_cond(p)
        p.u += du_P(p)
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



end
