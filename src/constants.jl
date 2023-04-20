module Constants
using TOML

export G, R_ig, yr, pc, Msun, m_p, k_B
export get_params
export F


# which type to use for floats?
F = Float64


# We use CGS for everything
const R_ig = 8.314e7   #erg/K/mol
const G = 6.67e-8 
const Msun = 1.989e33  # grams
const pc = 3.086e18    # cm
const yr = 3.15e7      # seconds
const m_p = 1.6726e-24 # grams
const k_B = 1.3806e-16 # erg/K




@Base.kwdef struct Params
    name::String = "sim"

    # simulation/timestepping 
    N::Int = 100

    seed::Int = -1 # to turn off the sed
    save_skip::Int = 25

    t_end::F = 1e8

    adaptive::Bool = true

    dt_min::F = 1e3yr
    dt_max::F = 1e7yr
    dt_rel_max::F = 20
    
    r_plummer::F = 0.2 # grav softening
    theta::F = 0.05 # BH theta

    # density estimation...
    # probably don't need all of these
    eta::F = 2
    h_max::F = 100pc
    h_min::F = 0.1pc
    rho_min::F = 1e-20m_p
    rho_max::F = 1e20m_p
    tol::F = 1e-3
    rho_maxiter::Int = 40



    # initialization parameters
    # Will likely move out of here
    # and into script which calls evolve!
    T0::F = 10
    mu_0::F = 1.4
    M_bary::F = 1e6Msun
    R_bary::F = 100pc
    sigma_M::F = 0
    sigma_v::F = 0
    v_0::F = 0
    gamma::F = 1.5


    ######### Physics ###########
    
    phys_gravity::Bool = false
    
    # Dark matter profile
    phys_DM::Bool = false
    c::F = 10
    M_tot::F = 1e8Msun
    R_virial::F = 1e3pc
    A_NFW::F = log(1+c) - c/(1+c)
    Rs::F = R_virial/c


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


function Params(filename)
    df =  TOML.parsefile(filename)

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
        elseif name == "c"
        end
    end


    df_symb = Dict(Symbol(k)=>v for (k, v) in df)
    return Params(;df_symb...)
end




end
