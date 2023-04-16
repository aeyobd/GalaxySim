module Physics
export dm_star, du_cool, du_cond, dv_P, du_P
export a_DM


using LinearAlgebra
using Debugger

using ..Particles
using ..Init
using ..Constants
using ..Density

# Lengths are pc, times are years
# Masses in solar mass


function a_DM(r::F, params)
    if r == 0
        return 0
    end
    G*params.M_tot/params.A_NFW * 1/r^2 * (r/(r + params.Rs) - log(1 + r/params.Rs))
end


function a_DM(x, params)
    if !params.phys_dm
        return zeros(3)
    end
    if norm(x) == 0
        return zeros(3)
    end
    return a_DM(norm(x), params) * normalize(x)
end


function du_cool(p::Particle, params)
    if !params.phys_cooling
        return 0.
    end
    dms = dm_star(p, params)
    ΔA = (p.m/p.ρ)^(2/3)
    dt = 1
    Σsfr = dms/dt / ΔA

    Γ = params.Gamma_0 * (params.f_rad*Σsfr*(2.5e-3*Msun/1e6pc^2/yr) + params.Jfuv/0.0024)
    Λ = 2e-19 * exp(-1.184e-5/(p.T + 1000)) + 2.8e-28*sqrt(p.T)*exp(-92/p.T)

    n = p.ρgas/p.mgas
    du_cool = -n*(n*Λ - Γ)
    return du_cool
end

function dv_P(p::Particle, params)
    if !params.phys_pressure
        return 0.
    end
    dv = zeros(3)
    for q in p.neighbors
        dv .+= q.m*(p.P/p.ρ^2 + q.P/q.ρ^2) * ∇W(p, q) 
    end
    return dv
end


function du_cond(p::Particle, params)
    if !params.phys_conduction
        return 0.
    end
    ΔT = 0.
    for (q, dist) in zip(p.neighbors, p.distances)
        ΔT += 2 * q.m/q.ρ * (p.T - q.T) * norm(∇W(p, q))/dist
    end
    K = params.K_cond/(1 + params.rho_cond*m_p/p.ρ)
    return K * ΔT
end


function du_P(p, params)
    if !params.phys_pressure
        return 0.
    end
    du = 0

    for q in p.neighbors
        du -= p.P/p.ρ/p.ρgas * q.m * (p.v-q.v) ⋅ ∇W(p, q)
    end
    return du
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
function dm_star(p::Particle, params)
    if !params.phys_star_formation
        return 0.
    end
    t_ff = sqrt(3π/(32 * G * p.ρgas))
    # ρmin = (p.T/6000)^3 * (p.mgas/1.3e6Msun)^(-2)
    dms = params.eta_eff * p.mgas/t_ff 
    return dms
end



end
