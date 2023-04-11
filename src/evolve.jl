module Evolve
export update!, evolve

using ..Physics
using ..Constants
using ..Density
using ..Init
using ..Particles

using NearestNeighbors
using Printf
using JLD2

function RK4()
end

function evolve(params)
    ps = rand_particles(params["Np"])
    t = 0

    while t < params["Tmax"]
        update_particles!(ps, t)
        dt_min = get_dt(ps)

        t += dt_min
        print_time(t_now, t_end)
    end

    return integrator.sol
end

function update_particles!(ps, t)
    balltree = BallTree([p.xs for p in ps])

    for p in ps
        if p.t + p.dt < t
            continue
        end

        idxs = inrange(balltree, p.xs, p.h, true)
        p.neighbors = ps[idxs]

        p.ρ = ρ(p)
        p.h = h(p)

        # if p.u < 0
        #     p.u = 0
        # end
        # if p.mstar < 0
        #     p.mstar = 0
        # end

        # p.mgas = p.m - p.mstar
        # p.ρgas = p.ρ * p.mgas/p.m
        # if p.ρgas < 0
        #     p.ρgas = 0
        #     p.mgas = 0
        # end

    end

end

function constraint!(U)
    for i in 1:N(U)
        if u(U) < 0
            u!(U)[i] = 0
            @debug "negative T warning"
        end

    end
end


function dp(p, dt)
    ρ1 = calc_dρ(x(particle), U, masses, params)
end


function particle_system!(dU, U, p, t)
    params, masses, dt = p

    # needs something 

    constraint!(U)

    for i in 1:N(U)

    end

    xs!(dU, p.v)
    vs!(dU, dv_G() + a_DM(x) + dv_P(x))
    u!(dU, du_P() + du_cool() + du_cond())
    ms!(dU, md_star)

    return dU
end

function neighbors(particles::Vector)
    x =  hcat(map(p->p.x, particles)...)
    nearest_neighbors(x, Nn)
end


"""
get the nearist parlticles 
x is array ndxnp
"""

function print_time(t, t_end)
    if t < 1e3*yr
        s = @sprintf("%4.0f yr", t/yr)
    elseif t < 1e6yr
        s = @sprintf("%4.0f kyr", t/1e3yr)
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
