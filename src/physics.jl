module Physics
export dm_star, du_cool, du_cond, dv_G, dv_P, du_P, dv_DM


using LinearAlgebra
using Debugger

using ..Particles
using ..Init
using ..Constants
using ..Density

# Lengths are pc, times are years
# Masses in solar mass

"""
Dρ/Dt = -ρ∇⋅v
Dv/Dt = -1/ρ ∇P
De/Dt = -1/ρ∇⋅Pv

energy per unit pass
e = u + v^2/2
P=(γ-1)ρ u
γ=5/3
"""
function dm_star(p::Particle, config)
    t_ff = sqrt(3π/(32 * G * p.ρgas))
    ρmin = (p.T/6000)^3 * (p.mgas/1.3e6Msun)^(-2)
    dms = p.ρgas > ρmin ? ϵ_eff * p.mgas/t_ff : 0
    return dms
end

f_rad = 1
Jfuv = 1
Γ0 = 2e-26
function du_cool(p::Particle)
    dms = dm_star(p)
    ΔA = (p.m/p.ρ)^(2/3)
    dt = 1
    Σsfr = dms/dt / ΔA
    Γ = Γ0 * (f_rad*Σsfr*(2.5e-3*Msun/1e6pc^2/yr) + Jfuv/0.0024)
    Λ = 2e-19 * exp(-1.184e-5/(p.T + 1000)) + 2.8e-28*sqrt(p.T)*exp(-92/p.T)

    n = p.ρgas/p.mgas
    du_cool = -n*(n*Λ - Γ)
    return du_cool
end

function dv_P(p::Particle)
    dv = zeros(3)
    for q in p.neighbors
        dv .+= q.m*(p.P/p.ρ^2 + q.P/q.ρ^2) * ∇W(p, q) 
    end
    return dv
end

function dv_G(p::Particle)
    dv = zeros(3)
    for (q, dist) in zip(p.neighbors, p.distances)
        if dist != 0
            dv .-= G*q.m/(dist)^3 * (q.x .- p.x)
        else
            @debug "particles atop eachother"
        end
    end

    return dv
end

dv_DM(p::Particle) = a_DM(p.x)

function du_cond(p::Particle)

    ΔT = 0.
    for (q, dist) in zip(p.neighbors, p.distances)
        ΔT += 2 * q.m/q.ρ * (p.T - q.T) * norm(∇W(p, q))/dist
    end
    K = K0/(1 + ρ0/p.ρ)
    return K * ΔT
end


function du_P(p)
    du = 0

    for q in p.neighbors
        du -= p.P/p.ρ/p.ρgas * q.m * (p.v-q.v) ⋅ ∇W(p, q)
    end
    return du
end



end
