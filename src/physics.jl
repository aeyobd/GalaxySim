module Physics
export main


using LinearAlgebra
import ForwardDiff: derivative, gradient
using StaticArrays
using NearestNeighbors
using Printf
using Debugger
using Glob

include("init.jl")
include("particle.jl")
include("density.jl")

# Lengths are pc, times are years
# Masses in solar mass


const G = init.G
const R = init.R
const yr = init.yr
const pc = init.pc
μ = 1



function main(dt=100e3*yr, t_end=100e6*yr)
    masses = init.rand_particles()

    files = open_files(length(masses))

    for t in 0:dt:t_end

        print_time(t, t_end)

        record_masses(files, masses)
        update_particles!(masses, dt)
    end

    close_files(files)

    println("closed files, completed!")
    return 
end

function print_time(t, t_end)
    if t < 1e6*yr
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
        println(f, "x1,x2,x3,v1,v2,v3,ρ,h,T")
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
        v1 = mass.v[1]
        v2 = mass.v[2]
        v3 = mass.v[3]
        ρ = mass.ρ
        h = mass.h
        T = mass.T
        write(file, "$x1,$x2,$x3,$v1,$v2,$v3,$ρ,$h,$T\n")
        flush(file)
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
    idxs, dists = neighbors(particles)

    for i in 1:length(particles)
        p = particles[i]

        p.x .+= p.v * dt
        p.v .+= init.a_DM(p.x) * dt
        p.neighbors = idxs[i][2:end]
        p.distances = dists[i][2:end]

        p.ρ = density.ρ(p, particles[p.neighbors], p.distances)
        p.h = density.h(p.ρ, p.m)


        for (q, dist) in zip(particles[p.neighbors], p.distances)
            # pressure
            p.v .-= q.m*(p.P/p.ρ^2 + q.P/q.ρ^2) * density.∇W(p, q) * dt
            #gravity
            if dist != 0
                p.v .-= G*q.m/(dist)^3 * (q.x - p.x)
            end

            p.u += p.P/p.ρ^2 * q.m * sum((p.v-q.v) .* density.∇W(p, q)) * dt

            p.T = p.u / (3/2 * R/μ)

            p.P = R/μ * p.T * p.ρ 
        end

    end

end


end
