module Physics
export main

using LinearAlgebra
import ForwardDiff: derivative


const G = 4.30091727e−3 # pc⋅M⊙^-1  ⋅(km/s)2 

mutable struct Mass
    x::Vector{AbstractFloat}
    v::Vector{AbstractFloat}
    m::AbstractFloat
end

function main(dt=0.1, t_end=10)
    m1 = Mass([1], [0], 1)

    open("gal.dat", "w") do f
        write(f, "time\tx\tv\n")

        for t in 0:dt:t_end
            update_particle!(m1, dt)
            x1 = m1.x[begin]
            v1 = m1.v[begin]

            write(f, "$t\t$x1\t$v1\n")

        end
    end
end


function Φ(x, masses)
    s = 0
    for mass in masses
        s += mass.m / norm(masses.x - x)
    end

    s *= -G
    return s
end


const k = 1
V(x) = k*sum(x.^2)

function dΦ(x)
    return derivative.(V, x)
end

function update_particle!(particle, dt)
    particle.x .+= particle.v*dt
    particle.v .+= -dΦ(particle.x)/particle.m
end



end
