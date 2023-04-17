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
    files = open_files(params.N)

    int_t = 0
    dt = params.dt_min
    t = 0

    while t < params.t_end
        update_particles!(ps, int_t, dt, params)

        print_time(t, params.t_end)

        if int_t % 64 == 0
            record_particles(files, ps, params)
        end

        int_t += get_dt(ps, params)

        t = int_t * params.dt_min
    end

    close_files(files)

    return
end


function evolve!(ps, params)
    files = open_files(params.N)

    int_t = 0
    t = 0
    dt = params.dt_min
    while t < params.t_end
        update_particles!(ps, int_t, dt, params)
        print_time(t, params.t_end)

        if int_t % 64 == 0
            record_particles(files, ps, params)
        end

        t = int_t * params.dt_min
        int_dt = get_dt(ps, params)
        int_t += int_dt
        dt = int_dt * params.dt_min
    end

    close_files(files)
end



function get_dt(ps, params)
    tm = minimum([p.dt for p in ps])/params.dt_min
    tm = min(tm, 2^32)

    if !params.adaptive
        return 1
    end
    return 1
end



function update_particles!(ps, t_int, dt, params)
    tree = make_tree(ps, params)

    update_ps = which_update!(ps, t_int, dt, params)
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
        update_dt!(p, t_int + 1, params)
    end

end


function which_update(ps, t_int, dt, params)
    if !params.adaptive
        for p in ps
            p.dt = params.dt_min
        end
        return ps
    end

    ps_new = []

    for p in ps
        if p.dt < 2*params.dt_min
            p.dt = params.dt_min
            push!(ps_new, p)
        else
            if factor > 2^params.dt_pow_max
                factor = 2^params.dt_pow_max
            end

            if t_int % factor == 0
                p.dt = factor*params.dt_min
                push!(ps_new, p)
            end
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


function update_dt!(p::Particle, t_int, params)
    if length(p.distances) < 1
        p.dt = params.dt_min
        return
    end

    dt = max(params.tol * 3/sqrt(8π * G * p.ρ), params.dt_min)

    dtn = minimum(q.dt for q in p.neighbors)
    dt = min(dt, 2^params.dt_rel_pow_max * dtn)

    if dt < 2*params.dt_min
        dt = params.dt_min
        return 
    end

    factor = 2^floor(Int, log2(dt/params.dt_min))
    if factor > 2^params.dt_pow_max
        factor = 2^params.dt_pow_max
    end

    if t_int % factor == 0
        p.dt = factor * params.dt_min
    else
        factor = 2^Base.trailing_zeros(t_int % factor)
        p.dt = factor * params.dt_min
    end

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


function open_files(N)
    # delete files
    for f in glob("data/*.dat")
        rm(f)
    end
    files = Vector()
    for i in 1:N
        fname = "data/mass$i.dat"
        f = open(fname, "w")
        println(f, "t,x1,x2,x3,v1,v2,v3,ρ,h,T,ms,f")
        push!(files, f)
    end
    println("Opened files")

    return files
end


function record_particles(files, masses, params)
    for (file, mass) in zip(files, masses)
        x1 = mass.x[1]/pc
        x2 = mass.x[2]/pc
        x3 = mass.x[3]/pc
        v1 = mass.v[1] /100_000 # km/s
        v2 = mass.v[2] /100_000 
        v3 = mass.v[3] /100_000 
        ρ = mass.ρ / Msun * pc^3
        h = mass.h / pc
        T = mass.T
        t = mass.t * params.dt_min / yr
        mstar = mass.mstar/Msun
        fstar = mass.mstar/mass.m
        write(file, "$t,$x1,$x2,$x3,$v1,$v2,$v3,$ρ,$h,$T,$mstar,$fstar\n")
        flush(file)
    end
end

function close_files(files)
    close.(files)
end




end
