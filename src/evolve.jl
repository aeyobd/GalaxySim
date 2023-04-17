module Evolve
export evolve

using ..Init
using ..Physics
using ..Density
using ..Particles
using ..Tree
using ..Constants

using Printf
using JLD2
using LinearAlgebra
using Glob


function evolve(params)
    ps = rand_particles(params)
    files = open_files(params)

    t = 0

    while t < params.t_end
        update_particles!(ps, int_t, dt, params)
        print_time(t, params.t_end)
        record_particles(files, ps, params)

        t += get_dt(ps, params)
    end

    close_files(files)

    return
end


function evolve!(ps, params)
    files = open_files(params)

    t = 0
    while t < params.t_end
        update_particles!(ps, t, params)
        print_time(t, params.t_end)

        record_particles(files, ps, params)

        t += get_dt(ps, params)
    end

    close_files(files)
end



function get_dt(ps, params)
    if !params.adaptive
        return params.dt_min
    end

    tm = minimum([p.dt for p in ps])
    return params.dt_min
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
end



function update_h_position!(p::Particle)
    p.x .+= p.v * p.dt/2
end



function update_energy!(p::Particle, params)
    p.u += du_P(p, params) * p.dt
    p.u += du_cond(p, params) * p.dt
    p.u += du_cool(p, params) * p.dt

    p.T = 2/(3R_ig) * p.u
    p.P = R_ig/params.mu_0 * p.ρ * p.T
end



function update_velocity!(p::Particle, tree, params)
    a = a_DM(p.x, params)
    a .+= dv_P(p, params)
    a_G!(a, p, tree, params)

    p.v .+= a * p.dt
end



function update_mstar!(p::Particle, params)
    p.mstar += dm_star(p, params)
end



function update_dt!(p::Particle, t, params)
    if length(p.distances) < 1
        p.dt = params.dt_min
        p.t += p.dt
        return
    end

    dt = max(params.tol * 3/sqrt(8π * G * p.ρ), params.dt_min)

    dtn = minimum(q.dt for q in p.neighbors)
    dt = min(dt, params.dt_rel_max * dtn)
    dt = min(dt, params.dt_max)
    p.t += p.dt
end




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



FILE_NAMES = (
    "t",
    "x1", "x2", "x3",
    "v1", "v2", "v3",
    "rho",
    "h",
    "T",
    "mstar"
   )

VAR_NAMES = (p->p.t/yr,
             p->p.x[1]/pc, p->p.x[2]/pc, p->p.x[3]/pc,
             p->p.v[1], p->p.v[2], p->p.v[3],    
             p->p.ρ,
             p->p.h/pc,
             p->p.T,
             p->p.mstar,
            )

function open_files(params)
    mkpath(params.name)

    # delete files
    for file in glob("$(params.name)/*.dat")
        rm(file)
    end

    files = Vector()

    for basename in FILE_NAMES
        fname = "$(params.name)/$(basename).dat"
        file = open(fname, "w")
        push!(files, file)
    end

    println("Opened files")

    return files
end




function record_particles(files, particles, params)
    for (file, var) in zip(files, VAR_NAMES)
        for p in particles
            val = var(p)
            print(file, "$val,")

        end

        println(file)
    end

end



function close_files(files)
    close.(files)
end




end
