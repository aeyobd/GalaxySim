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
using Printf
using Logging


function evolve(params)
    ps = rand_particles(params)

    evolve!(ps, params)
end


function evolve!(ps, params)
    files = open_files(params)
    e_file = open("energy.dat", "w")
    log_file = open("log.txt", "w")
    global_logger(SimpleLogger(log_file))

    println(e_file,     "thermal,kinetic,grav,tot")
    set_densities!(ps, params)
    t = 0
    i = 0
    while t < params.t_end
        update_particles!(ps, t, params)

        # only save once every so many frames
        print_time(t, params.t_end)
        if i % params.save_skip == 0
            record_particles(files, ps, params)
            total_energy(ps, e_file, params)
        end

        t += get_dt(ps, t, params)
        i += 1
    end

    close_files(files)
    close(e_file)
    close(log_file)
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
        update_dt!(p, t, params)
        p.v .+= p.dv * p.dt/2
        p.x .+= p.v * p.dt

        p.t += p.dt
    end

    tree = make_tree(ps, params)

    for p in update_ps
        update!(p, tree, params)
    end

    for p in update_ps
        p.v .+= p.dv * p.dt/2
        update_m_star!(p, params)
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



function update!(p::Particle, tree, params)
    find_neighbors!(p, tree, params)

    dh!(p, params)
    dρ!(p, params)

    du_P!(p, params)
    du_cond!(p, params)
    dv_DM!(p, params)
    dv_P!(p, params)
    if params.phys_gravity
        p.dv_G = zeros(3)
        dv_G!(p, tree, params)
    end


    p.du = p.du_P + p.du_cond

    if p.du*p.dt + p.u < 0
        p.du = p.u/p.dt/2
        @debug "energy is near 0"
    end

    p.u += p.du * p.dt
    p.dv .= p.dv_G .+ p.dv_P 
    p.t += p.dt


    p.P = R_ig/p.μ * p.ρ * p.T
    p.T = 2p.μ/(3R_ig) * p.u

    # constraints on density evolution
    if p.ρ + p.dρ * p.dt < params.rho_min
        p.dρ = (params.rho_min-p.ρ)/p.dt
        @debug "density hit limit"
    end
    if p.h + p.dh * p.dt < params.h_min
        p.dh = (params.h_min - p.h)/p.dt
        @debug "h hit limit"
    end

    p.h += p.dh * p.dt
    p.ρ += p.dρ * p.dt
    p.ρ_gas = p.ρ * p.m_gas/p.m

    cs!(p)
end


function total_energy(ps, e_file, params)
    tree = make_tree(ps, params)

    thermal = 0.
    kinetic = 0.
    grav = 0.
    

    for p in ps
        thermal += p.u*p.m
        kinetic += 1/2*norm(p.v)^2*p.m
        # diviede by two since we double count pairs
        grav += U_G(p, tree, params)/2
    end
    tot = thermal + kinetic + grav

    @printf e_file "%8.4e, %8.4e, %8.4e, %8.4e\n" thermal kinetic grav tot

    if ps[1].t == 0
        println("initial energy: $tot")
    elseif ps[1].t + ps[1].dt >= params.t_end
        println("final energy: $tot")
    end

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
    if !params.adaptive
        return p.dt = params.dt_min
    end
    dt_f = p.h/norm(p.dv)
    dt_c = p.h/p.c/(1+0.6*params.alpha_visc)
    dt_g = 1/sqrt(8π*G*p.ρ)

    dt = 0.25*min(dt_f, dt_c, dt_g)

    dtn = minimum(q.dt for q in p.neighbors)
    dt = min(dt, params.dt_rel_max * dtn)
    dt = min(dt, params.dt_max)
    dt = max(dt, params.dt_min)

    p.dt = dt
end




end
