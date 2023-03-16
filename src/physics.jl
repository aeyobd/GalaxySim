module Physics
export main

using LinearAlgebra
import ForwardDiff: derivative, gradient
using StaticArrays
using NearestNeighbors

# Lengths are pc, times are years
# Masses in solar mass

const G = 4.4985e-15

mutable struct Mass
    x::SVector{3, AbstractFloat}
    v::SVector{3, AbstractFloat}
    m::AbstractFloat
end

function main(dt=10, t_end=500_000)
    m1 = Mass([0.01, 0, 0], [0, 6.7038e-8, 0], 1)
    m2 = Mass([-0.01, 0, 0], [0, -6.7038e-8, 0], 1)
    masses = [m1, m2]

    files = open_files(length(masses))
    # open_files
    for t in 0:dt:t_end
        record_masses(files, masses)
        update_particles!(masses, dt)
    end

    close_files(files)
end


function Φ(x, masses)
    s = 0
    for mass in masses
        s += mass.m / norm(mass.x - x)
    end

    s *= -G
    return s
end

function open_files(N)
    files = Vector()
    for i in 1:N
        f = open("data/mass$i.dat", "w")
        println(f, "x1,x2,x3,v1,v2,v3")
        push!(files, f)
    end

    return files
end

function record_masses(files, masses)
    for (file, mass) in zip(files, masses)
        x1 = mass.x[1]
        x2 = mass.x[2]
        x3 = mass.x[3]
        v1 = mass.v[1]
        v2 = mass.v[2]
        v3 = mass.v[3]
        write(file, "$x1,$x2,$x3,$v1,$v2,$v3\n")
    end
end

function close_files(files)
    close.(files)
end

function dΦ(x, others)
    gradient(x) do x 
        Φ(x, others)
    end
end


"""
get the nearist parlticles 
x is array ndxnp
"""
function nearist_neighbors(x, k=10)
    tree = KDTree(x)
    idxs, dists = knn(tree, data, k, true)
    return idxs, dists
end


"""
get density at r given particles and smoothing length h

ρ = ∑_b m_b W(r_a - r_b; h); (h smoothing length)
dv/dt = 
"""
function ρ(r, particles, h)
    # weighted sum over nearby particles
    # # e.g. quintic spline
end


"""
The ideal smoothing length at r
"""
function h(r, particles)

end


"""
Dρ/Dt = -ρ∇⋅v
Dv/Dt = -1/ρ ∇P
De/Dt = -1/ρ∇⋅Pv

energy per unit pass
e = u + v^2/2
P=(γ-1)ρ u
γ=5/3
"""
function update_particle!(particle, dt)

    particle.x .+= particle.v * dt 
    particle.v .-= dΦ(particle.x)*dt
end

function update_particles!(particles, dt)
    for i in 1:length(particles)
        particles[i].x .+= particles[i].v * dt
        particles[i].v .-= dΦ(particles[i].x, particles[1:end .!=i])
    end

end


end
