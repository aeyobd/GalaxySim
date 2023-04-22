"""
module Evolve


Contains the core of the simulation.

exports evolve!(ps::Vector{Particle}, params::Params)

"""
module Evolve

export evolve!

using ..Physics
using ..Density
using ..Particles
using ..Constants
using ..GalFiles
using ..Gravity

using LinearAlgebra
using Printf
using Logging



function evolve!(ps::Vector{Particle}, params)
    # open up files to write to 
    log_file = open("$(params.name)/log.txt", "w")
    global_logger(SimpleLogger(log_file, Logging.Debug))

    files = open_files(params)
    e_file = open("$(params.name)/energy.dat", "w")
    println(e_file,     "thermal,kinetic,grav,tot")

    find_neighbors!(ps, params)
    for p in ps
        solve_ρ!(p, params)
        p.T = temp(p)
        p.P = pressure(p)
        p.c = c_sound(p)
    end

    _, _, _, tot = energy(ps, params)
    println("initial energy: $tot")


    t = 0
    i = 0 # keep track of the number of frames
    while t < params.t_end
        update_particles!(ps, t, params)

        # only save once every so many frames
        if i % params.save_skip == 0
            record_particles(files, ps, params)
            save_energy(ps, e_file, params)
        end

        print_time(t, params.t_end)
        t += get_dt(ps, t, params)
        i += 1
    end

    _, _, _, tot = energy(ps, params)
    println("final energy: $tot")

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



"""
Updates the quantities for all particles

"""
function update_particles!(ps, t, params)
    # we only update particles as needed
    update_ps = which_update(ps, t, params)

    for p in update_ps
        p.t += p.dt
        p.x .+= p.v * p.dt/2
    end


    for p in update_ps
        p.v .+= p.dv * p.dt/2
    end


    find_neighbors!(ps, params)

    for p in update_ps
        solve_ρ!(p, params)
        update!(p, params)
    end

    for p in update_ps
        p.v .+= p.dv * p.dt/2
        p.x .+= p.v * p.dt/2
    end

    for p in update_ps
        # finally update other parameters
        update_dt!(p, t, params)
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
        if p.t + p.dt/2 < t # update particle when t = t_i+1/2
            push!(ps_new, p)
        end
    end

    return ps_new
end



function update!(p::Particle, params)
    if params.phys_pressure
        du_P!(p, params)
        dv_P!(p, params)
    end

    if params.phys_DM
        dv_DM!(p, params)
    end

    if params.phys_conduction
        du_cond!(p, params)
    end

    if params.phys_gravity
        dv_G!(p, params)
    end

    if params.phys_visc
        dv_visc!(p, params)
        du_visc!(p, params)
    end

    p.du = p.du_P + p.du_cond + p.du_visc
    @. p.dv = p.dv_P + p.dv_G + p.dv_DM + p.dv_visc

    # prevent energy from going below zero
    # This does violate energy conservation, 
    # but I haven't implemented a better scheme yet
    if p.du*p.dt + p.u < 0
        p.du = p.u/p.dt/2
        @debug "energy is near 0"
    end

    p.u += p.du * p.dt
    p.dv .= p.dv_G .+ p.dv_P 

    p.t += p.dt

    p.T = temp(p)
    p.P = pressure(p)
    p.c = c_sound(p)
end


"""
calculates the current system energy
"""
function energy(ps, params)
    grav = L_grav(ps, params)

    thermal = 0.
    kinetic = 0.
    for p in ps
        thermal += p.u*p.m
        kinetic += 1/2*norm(p.v)^2*p.m
    end

    tot = thermal + kinetic + grav

    return thermal, kinetic, grav, tot
end



function save_energy(ps, e_file, params)
    thermal, kinetic, grav, tot = energy(ps, params)
    @printf e_file "%8.4e, %8.4e, %8.4e, %8.4e\n" thermal kinetic grav tot
end




function update_m_star!(p::Particle, params)
    dm_star!(p, params)

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
    dt_c = p.h/p.c

    dt_g = 1/sqrt(8π*G*p.ρ)

    dt = 0.25*min(dt_f, dt_c, dt_g)

    factor = 2^ceil(log2(params.dt_max/dt))
    dt = params.dt_max/factor

    dtn = minimum(q.dt for q in p.neighbors)
    dt = min(dt, params.dt_rel_max * dtn)
    dt = min(dt, params.dt_max)
    dt = max(dt, params.dt_min)

    p.dt = dt
end




end
