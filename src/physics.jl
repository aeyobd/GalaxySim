module Physics
export main


using LinearAlgebra
import ForwardDiff: derivative, gradient
using StaticArrays
using NearestNeighbors

include("init.jl")
include("particle.jl")

# Lengths are pc, times are years
# Masses in solar mass


const G = init.G




function main(dt=100e3, t_end=100e6)
    masses = init.rand_particles()

    files = open_files(length(masses))

    for t in 0:dt:t_end
        record_masses(files, masses)
        update_particles!(masses, dt)
    end

    close_files(files)
end


function open_files(N)
    files = Vector()
    for i in 1:N
        fname = "data/mass$i.dat"
        println(fname)
        touch("data/mass$i.dat")
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


function Φ(x, masses)
    s = 0
    for mass in masses
        s += mass.m / norm(mass.x - x)
    end

    s *= -G
    return s
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
function nearest_neighbors(x, k=10)
    tree = KDTree(x)
    idxs, dists = knn(tree, data, k, true)
    return idxs, dists
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
function update_particles!(particles, dt)
    for i in 1:length(particles)
        ps = particles[1:end .!=i]
        p = particles[i]

        p.x .+= p.v * dt
        p.v .+= init.a_DM(p.x) * dt

        # p.ρ, p.h = p # add this in

        # for q in ps
            # p.v -= q.m*(p.P/p.ρ^2 + q.P/q.ρ^2) * ∇W(p, q) * dt
            # p.u -= p.P/p.ρ^2 * q.m * (v.*∇W.(p, q)) * dt
        # end
        # p.P = K*p.ρ^γ
    end

end


end
