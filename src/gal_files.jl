# gal_files.jl
#
# Methods to output results
#
# Created 17-04-2023
# Author Daniel Boyea

module GalFiles

export print_time
export open_files, record_particles, close_files
export energy

using Printf
using Glob
using LinearAlgebra
using Logging

using ..Constants
using ..Gravity


# Function to nicely print the time 
# in human readable format
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

    print("t =\t$s, $sp\r")
end



# this should have been a dictionary :(
# File names to store data in
FILE_NAMES = [
    ("t", p->p.t/yr,"# the particle time in years"),
    ("x1", p->p.x[1]/pc, "# The particle position x1 in pc"),
    ("x2",  p->p.x[2]/pc, "# The particle position x2 in pc"),
    ("x3", p->p.x[3]/pc,"# The particle position x3 in pc"),
    ("v1",  p->p.v[1]/1e5, "# The particle velocity v1 in km/s"),
    ("v2",  p->p.v[2]/1e5, "# The particle velocity v2 in km/s"),
    ("v3", p->p.v[3]/1e5,    "# The particle velocity v3 in km/s"),
    ("rho", p->p.ρ/m_p,"# density in cm^-3"),
    ("h", p->p.h/pc,"# smoothing length in pc"),
    ("T", p->p.T,"# temperatures in K"),
    ("mstar", p->p.m_star/Msun,"# mass of star formation"),
    ("N_neighbors", p->length(p.neighbors),"# number of neighbors"),
    ("dt", p->p.dt/yr,"# particle timesteps in years"),
    ("P", p->p.P,"# particle pressure in erg/cm"),
    ("du_P", p->p.du_P, "# change in energy due to pressure in ergs"),
    ("du_C", p->p.du_cond, "# change in energy due to conduction"),
    ("du_visc", p->p.du_visc,"# change in energy due to dissipation"),
    ("dv_P", p->norm(p.dv_P),"# acceleration due to pressure"),
    ("dv_G", p->norm(p.dv_G),"# acceleration due to gravity"),
    ("dv_visc", p->norm(p.dv_visc),"# acceleration due to viscosity "),
    ("dv_DM", p->norm(p.dv_DM),"# acceleration due to dark matter (cm/s^2)"),
   ]


"""
Creates and opens the files to write the simulation output to
"""
function open_files(params)
    mkpath(params.name)

    # delete old
    for file in glob("$(params.name)/*.dat")
        rm(file)
    end

    log_file = open("$(params.name)/log.txt", "w")
    global_logger(SimpleLogger(log_file, Logging.Debug))

    e_file = open("$(params.name)/energy.dat", "w")
    println(e_file,     "thermal,kinetic,grav,tot")

    var_files = Vector()

    for (basename, _, header) in FILE_NAMES
        fname = "$(params.name)/$(basename).dat"
        file = open(fname, "w")
        println(file, header)
        println(file, "# each row is a timestep and each column is a particle")
        push!(var_files, file)
    end

    println("Opened files")

    return var_files, e_file, log_file
end


"""
Records the current values to files
"""
function record_particles(files, particles, params)
    var_files, e_file, log_file = files

    for (file, names) in zip(var_files, FILE_NAMES)
        var = names[2]

        for p in particles
            val = var(p)
            @printf file "%12.8e    " val
        end

        println(file)
    end

    save_energy(particles, e_file, params)
end


"""
Closes the simulation files
"""
function close_files(files)
    var_files, e_file, log_file = files
    close.(var_files)
    close(e_file)
    close(log_file)
end



"""
calculates the current system thermal, kinetic, gravitational, and total energy
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



"""
Saves the energy value to a file
"""
function save_energy(ps, e_file, params)
    thermal, kinetic, grav, tot = energy(ps, params)
    @printf e_file "%8.4e, %8.4e, %8.4e, %8.4e\n" thermal kinetic grav tot
end



end
