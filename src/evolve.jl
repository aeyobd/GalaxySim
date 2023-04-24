# modeul Evolve
#
# Contains the core time-stepping loop
#
# Created 10-04-2021
# Author Daniel Boyea (boyea.2@osu.edu)
#
# This contains the main integration loop of the 
# simulator and controls timestepping
#
# I chose to use a leapfrog integration scheme,
# which should keep energy/momentum relatively
# conserved without being overly complex.
#
# I also allow the timestep of each particle to vary, 
# only updating when the simulation time passes that of the 
# particle. However, I do restrict the timesteps to powers
# of 2 (which is said to improve numerical stability).
#
# The timesteps are not perfectly syncronized however,
# so this is something a future version could improve
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



"""
Evolves the given vector of particles
with the given params until `params.t_end`
"""
function evolve!(ps::Vector{Particle}, params)
    files = open_files(params)
    setup!(ps, params)
    println("initial energy: ", energy(ps, params)[end])
    println()

    try
        main_loop!(ps, files, params)
    finally
        close_files(files)
        println()
        println("final energy: ", energy(ps, params)[end])
    end
end


""" calculate initial density temp and pressure for each particle"""
function setup!(ps, params)
    find_neighbors!(ps, params)
    for p in ps
        solve_ρ!(p, params)
        p.T = temp(p)
        p.P = pressure(p)
        p.c = c_sound(p)
    end
end


"""Main loop which updates particles"""
function main_loop!(ps, files, params)
    t = 0
    i = 0 # keep track of the number of loops
    while t < params.t_end
        # timestep the particles
        update_particles!(ps, t, params)
        
        if i % params.save_skip == 0 # only save once every so many frames
            record_particles(files, ps, t, params)
        end

        print_time(t, params.t_end)
        t += get_dt(ps, t, params)
        i += 1
    end
end



""" Determines the next simulation timestep """
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

Uses leapfrog integration.

"""
function update_particles!(ps, t, params)
    update_ps = which_update(ps, t, params) # we only update particles as needed

    for p in update_ps
        p.x .+= p.v * p.dt/2
    end


    for p in update_ps
        p.v .+= p.dv * p.dt/2
    end

    # recalculate ρ to find acceleration
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
        p.t += p.dt
        update_dt!(p, t, params)
        update_m_star!(p, params)
    end
end



"""
Determins which of the list of particles, ps, 
need updated if the current simulation time is t
"""
function which_update(ps, t, params)
    if !params.adaptive # use parameters to switch off as needed
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


"""
Updates the acceleration values,
and energy, temperature, and pressure of p
"""
function update!(p::Particle, params)
    # switches for each process
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

    # add together all energy and acceleration parameters
    p.du = p.du_P + p.du_cond + p.du_visc
    @. p.dv = p.dv_P + p.dv_G + p.dv_DM + p.dv_visc


    # prevent energy from going below zero
    # This does in fact violate energy conservation, 
    # but I haven't implemented a better scheme yet
    if p.du*p.dt + p.u < 0
        p.du = p.u/p.dt/2
        @debug "energy is near 0"
    end

    p.u += p.du * p.dt
    p.dv .= p.dv_G .+ p.dv_P 

    p.T = temp(p)
    p.P = pressure(p)
    p.c = c_sound(p)
end






""" Updates the star formation mass of particle p """
function update_m_star!(p::Particle, params)
    if !params.phys_star_formation
        return 0.
    end
    dm_star!(p, params)

    if p.m_star > p.m
        p.m_star = p.m
        p.dm_star = p.m_gas
    end

    p.m_star += p.dm_star * p.dt

    p.m_gas = p.m - p.m_star
    p.ρ_gas = p.m_gas/p.m * p.ρ
end




"""
Updates the particles timestep dt
Timesteps are rounded down to powers of 2 
of the maximum timestep
"""
function update_dt!(p::Particle, t, params)
    if !params.adaptive
        return p.dt = params.dt_min
    end
    # timestep based on acceleration
    dt_f = p.h/norm(p.dv)
    # sound timestep
    dt_c = p.h/p.c

    # gravitational timestep
    dt_g = 1/sqrt(8π*G*p.ρ)

    dt = 0.25*min(dt_f, dt_c, dt_g)

    # round down
    factor = 2^ceil(log2(params.dt_max/dt))
    dt = params.dt_max/factor

    # make sure timestp isn't too high or low
    dtn = minimum(q.dt for q in p.neighbors)
    dt = min(dt, params.dt_rel_max * dtn)
    dt = min(dt, params.dt_max)
    dt = max(dt, params.dt_min)

    p.dt = dt
end




end
