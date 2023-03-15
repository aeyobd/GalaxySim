module Physics
export main

using LinearAlgebra
import ForwardDiff: derivative, gradient

# Lengths are pc, times are years
# Masses in solar mass

const G = 4.4985e-15

mutable struct Mass
    x::Vector{AbstractFloat}
    v::Vector{AbstractFloat}
    m::AbstractFloat
end

function main(dt=100, t_end=100_000)
    m1 = Mass([0.01, 0, 0], [0, 6.7038e-7, 0], 0.0001)

    open("gal.dat", "w") do f
        write(f, "time\tx1\tx2\tx3\tv1\tv2\ta\n")

        for t in 0:dt:t_end
            x1 = m1.x[1]
            x2 = m1.x[2]
            x3 = m1.x[3]

            v1 = m1.v[1]
            v2 = m1.v[2]
            v3 = m1.v[3]

            d = dΦ(m1.x)
            a1 = d[1]
            write(f, "$t\t$x1\t$x2\t$x3\t$v1\t$v2\t$a1\n")

            update_particle!(m1, dt)
        end
    end
end


function Φ(x, masses)
    s = 0
    for mass in masses
        s += mass.m / norm(mass.x - x)
    end

    s *= -G
    return s
end


function dΦ(x)
    m0 = Mass([0,0,0], [0,0,0], 1)
    return gradient(x->Φ(x, [m0]), x)
end

function update_particle!(particle, dt)
    particle.x .+= particle.v * dt 
    particle.v .-= dΦ(particle.x)*dt
end



end
