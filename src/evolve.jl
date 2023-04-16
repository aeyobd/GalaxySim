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

    t = 0
    dt = params.dt_min

    while t < params.t_end/params.dt_min

        update_particles!(ps, t, dt, params)

        print_time(t*params.dt_min, params.t_end)

        if t % 64 == 0
            record_particles(files, ps)
        end

        t += get_dt(ps, params)
    end

    close_files(files)

    return
end



function get_dt(ps, params)
    tm = minimum([p.dt for p in ps])/params.dt_min
    tm = min(tm, 2^32)

    return max(1, tm)
end



function update_particles!(ps, t_int, dt, params)
    tree = make_tree(ps, params)

    update_ps = which_update(ps, t_int, dt, params)
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
        update_other!(p, params)
    end

end


function which_update(ps, t_int, dt, params)
    ps_new = []
    for p in ps
        if p.dt < 2*params.dt_min
            p.dt = params.dt_min
            push!(ps_new, p)
        else
            factor = 2^floor(Int, log2(p.dt/params.dt_min))
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
    p.t += p.t + p.dt

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


function update_other!(p::Particle, params)
    p.mstar += dm_star(p, params)
    if length(p.distances) > 0
        p.dt = max(params.tol * 3/sqrt(8π * G * p.ρ), params.dt_min)
    else
        p.dt = params.dt_min
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


function record_particles(files, masses)
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
        t = mass.t / yr
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
