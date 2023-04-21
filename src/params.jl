# params.jl
#
# Revision history
#       2023-04-21: Created 
#                   Daniel Boyea (boyea.2@osu.edu)
#
# This file contains the definition of the Params struct
# This struct holds all simulation parameters and options.
# Initialization is treated seperately (see `../init/`)
#
# The params struct includes a option to parse
# options from a TOML file, but all defaults are declared
# here. 

using ..Constants
using TOML # to read in options from file


@Base.kwdef struct Params
    name::String = "data"           # the name of the directory to place data files in

    N::Int = 100                    # The number of simulation particles
    NN::Int = 20                    # The number of neighbors to use for each particle

    t_end::F = 1e8                  # How long to integrate for
    save_skip::Int = 25             # Number of timesteps between writing results to file
    adaptive::Bool = true           # Whether or not to use adaptive timestepping

    dt_min::F = 1e3yr               # minimum value for dt, or dt when not adaptive
    dt_max::F = 1024e3yr            # maximum allowed timestep
    dt_rel_max::F = 16              # maximum dt relative to neighbors
    
    eta::F = 0.74                    # density scale parameters
    h_max::F = 1000pc                # smoothing length 
    h_min::F = 0.01pc
    h_maxiter::Int = 20             # maximum iterations when solving for h
    tol::F = 1e-3                   # required accuracy for h



    ######### Physics ###########
    
    phys_gravity::Bool = false
    
    # Dark matter profile
    phys_DM::Bool = false
    c::F = 10
    M_tot::F = 1e8Msun
    R_virial::F = 1e3pc
    A_NFW::F = log(1+c) - c/(1+c)   # helper constant
    Rs::F = R_virial/c              # helper constant


    # pressure and viscosity
    phys_pressure::Bool = false
    phys_visc::Bool = false
    alpha_visc::F = 1
    beta_visc::F = 2
    eta_visc::F = 0.1


    # thermal conduction
    phys_conduction::Bool = false
    K_cond::F = 4e7
    n_K::F = 0.05


    # star formation
    phys_star_formation::Bool = false
    sfe::F = 0.01
    rho_0::F = 1m_p


    # feedback
    phys_feedback::Bool = false
    Jfuv::F = 0.0014
    Gamma_0::F = 0
    f_rad::F = 2

end


"""
Takes a filename and returns params
"""
function Params(filename) 
    # read the file in through TOML
    df = TOML.parsefile(filename)

    # convert natural units in file to 
    # cgs
    for (name, _) in df
        if name in ["M_bary", "M_tot"]
            df[name] *= Msun
        elseif name in ["h_max", "h_min",
                       "R_virial", "R_bary"]
            df[name] *= pc
        elseif name in ["rho_max", "rho_min", "rho_0"]
            df[name] *= m_p
        elseif name in ["dt_min", "dt_max", "t_end"]
            df[name] *= yr
        end
    end

    # parse the options into keyword arguments to pass
    # to params
    df_symb = Dict(Symbol(k)=>v for (k, v) in df)
    return Params(;df_symb...)
end


