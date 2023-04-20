module Evolve

export evolve, evolve!, setup!

using ..Init
using ..Physics
using ..Density
using ..Particles
using ..Tree
using ..Constants
using ..GalFiles

using LinearAlgebra


function evolve(params)
    ps = rand_particles(params)

    evolve!(ps, params)
end


function evolve!(ps, params)
    files = open_files(params)
    set_densities!(ps, params)
    t = 0
    i = 0
    while t < params.t_end
        update_particles!(ps, t, params)

        # only save once every so many frames
        print_time(t, params.t_end)
        if i % params.save_skip == 0
            record_particles(files, ps, params)
        end

        t += get_dt(ps, t, params)
        i += 1
    end

    close_files(files)
end



function get_dt(ps, t, params)
    if !params.adaptive
        return params.dt_min
    end

    tm = minimum([(t - p.t) for p in ps])
    tm = max(params.dt_min, tm)
    return tm
end



function update_particles!(ps, t, params)
    tree = make_tree(ps, params)

    update_ps = which_update(ps, t, params)

    d1 = d3 = 1.35120719
    d2 = -1.70241438
    c1 = c4 = d1/2
    c2 = c3 = (d1 + d2)/2
    # where to calculate accel.
    e1 = c1+c4
    e2 = c2
    e3 = c3

    # use 4th order yoshida integrator
    for p in update_ps
        setup!(p, tree, params)

        p.x .+= c1*p.v*p.dt

        update!(p, tree, params)
        advance!(p, e1*p.dt, params)
    end

    for p in update_ps
        p.x .+= c2*p.v*p.dt
        # advances dv approprietly
        update!(p, tree, params)
        advance!(p, e2*p.dt, params)
    end

    for p in update_ps
        p.x .+= c3*p.v*p.dt

        update!(p, tree, params)
        advance!(p, e3*p.dt, params)
    end

    for p in update_ps
        p.x .+= c4*p.v*p.dt
        # last step doesn't increase v
        update_m_star!(p, params)
        update_dt!(p, t, params)
    end

end


function which_update(ps, t, params)
    if !params.adaptive
        for p in ps
            p.dt = params.dt_min
        end
        return ps
    end

    ps_new = Particle[]

    for p in ps
        if p.t + p.dt < t
            push!(ps_new, p)
        end
    end

    return ps_new
end


function setup!(p::Particle, tree, params)
    p.neighbors = find_within_r(p, tree, p.h)
    p.t += p.dt
end



function update!(p::Particle, tree, params)
    dh!(p, params)
    dρ!(p, params)

    du_P!(p, params)
    du_cond!(p, params)
    p.du = p.du_P + p.du_cond

    dv_DM!(p, params)
    dv_P!(p, params)
    if params.phys_gravity
        dv_G!(p, tree, params)
    end
    p.dv = @. p.dv_G + p.dv_P 
end


function advance!(p::Particle, dt, params)
    p.h = max(p.h +p.dh * dt, 1)
    p.ρ = max(p.ρ + p.dρ * dt, 0)
    p.v .+= p.dv * dt
    p.u += p.u * dt
    p.ρ_gas = p.ρ * p.m_gas/p.m
    cs!(p)

    p.T = 2p.μ/(3R_ig) * p.u
    p.P = R_ig/p.μ * p.ρ * p.T
end


function update_m_star!(p::Particle, params)
    if p.t > p.dt
        dm_star!(p, params)
        if p.m_star > p.m
            p.m_star = p.m
            p.dm_star = p.m_gas
        end
        p.m_star += p.dm_star * p.dt

        p.m_gas = p.m - p.m_star
    end
end



function update_dt!(p::Particle, t, params)
    dt = params.tol * 3/sqrt(8π * G * p.ρ)
    r = (p.m/p.ρ)^(1/3)
    dt_P = r / p.c
    dt = min(dt, dt_P)

    dtn = minimum(q.dt for q in p.neighbors)
    dt = min(dt, params.dt_rel_max * dtn)
    dt = min(dt, params.dt_max)
    dt = max(dt, params.dt_min)

    p.dt = dt
end




end
