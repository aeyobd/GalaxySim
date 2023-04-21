module GalFiles

export print_time
export open_files, record_particles, close_files

using Printf
using Glob
using LinearAlgebra

using ..Constants



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
    "mstar",
    "N_neighbors",
    "dt",
    "P",
    "du_P",
    "du_C",
    "dv_P",
    "dv_G",
   )

VAR_NAMES = (p->p.t/yr,
             p->p.x[1]/pc, p->p.x[2]/pc, p->p.x[3]/pc,
             p->p.v[1]/1e5, p->p.v[2]/1e5, p->p.v[3]/1e5,    
             p->p.ρ/m_p,
             p->p.h/pc,
             p->p.T,
             p->p.m_star/Msun,
             p->length(p.neighbors),
             p->p.dt/yr,
             p->p.P,
             p->p.du_P,
             p->p.du_cond,
             p->norm(p.dv_P),
             p->norm(p.dv_G),
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
