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

    update_ps = which_update(ps, t, params)

    # leapfrog integration
    for p in update_ps
        p.v .+= p.dv * p.dt/2
        p.x .+= p.v * p.dt
    end

    tree = make_tree(ps, params)
    for p in update_ps
        # recalculate dv
        setup!(p, tree, params)
    end

    for p in update_ps
        update!(p, tree, params)
        p.v .+= p.dv * p.dt/2
    end

    for p in update_ps
        update_m_star!(p, params)
        update_dt!(p, t, params)
        p.t += p.dt
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
    dh!(p, params)
    dρ!(p, params)
    p.h += p.dh * p.dt
    p.ρ += p.dρ * p.dt

    p.P = R_ig/p.μ * p.ρ * p.T
    p.T = 2p.μ/(3R_ig) * p.u

    p.ρ_gas = p.ρ * p.m_gas/p.m

    cs!(p)
end



function update!(p::Particle, tree, params)
    du_P!(p, params)
    du_cond!(p, params)

    dv_DM!(p, params)
    dv_P!(p, params)
    if params.phys_gravity
        p.dv_G = zeros(3)
        dv_G!(p, tree, params)
    end

    p.du = p.du_P + p.du_cond
    p.u += p.du * p.dt
    p.dv .= p.dv_G .+ p.dv_P 
    p.t += p.dt

end


function update_m_star!(p::Particle, params)
    dm_star!(p, params)

    # check for positive...
    if p.m_star > p.m
        p.m_star = p.m
        p.dm_star = p.m_gas
    end
    p.m_star += p.dm_star * p.dt

    p.m_gas = p.m - p.m_star
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

    if dt ===NaN
        println("dt")
        exit()
    end
    p.dt = dt
end




end
