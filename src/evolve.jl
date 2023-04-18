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

    for p in update_ps
        update_h_position!(p)
    end
    
    for p in update_ps
        setup!(p, tree, params)
        update_energy!(p, params)
    end

    for p in update_ps
        update_velocity!(p, tree, params)
        update_h_position!(p)
        update_mstar!(p, params)
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
    p.distances = [dist(p, q) for q in p.neighbors]
    p.t += p.dt

    p.ρ = ρ(p, params)
    p.h = h(p, params)
    p.c = cs(p)
end



function update_h_position!(p::Particle)
    p.x .+= p.v * p.dt/2
end



function update_energy!(p::Particle, params)
    p.du_P = du_P(p, params) 
    p.du_cond = du_cond(p, params) 
    p.du_visc = du_visc(p, params) 
    ducool = du_cool(p, params)

    p.u += (p.du_P + p.du_cond + p.du_visc + ducool) * p.dt

    if p.u < 0
        p.u = 0
        # eventually steal energy from neighbors
    end

    p.T = 2/(3R_ig) * p.u
    p.P = R_ig/params.mu_0 * p.ρ * p.T
end



function update_velocity!(p::Particle, tree, params)
    dv_G = a_DM(p.x, params) 
    p.dv_P .= dv_P(p, params) 
    p.dv_visc .= dv_visc(p, params) 

    if params.phys_gravity
        a_G!(dv_G, p, tree, params)
    end

    dv = @. dv_G + p.dv_P + p.dv_visc

    p.v .+= dv * p.dt
end



function update_mstar!(p::Particle, params)
    p.mstar += dm_star(p, params)
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
