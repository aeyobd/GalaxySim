module Constants
using TOML

export G, R_ig, yr, pc, Msun, m_p, k_B
export get_params
export F


# We use CGS for everything
const R_ig = 8.314e7   #erg/K/mol

const G = 6.67e-8 
const Msun = 1.989e33  # grams
const pc = 3.086e18    # cm
const yr = 3.15e7      # seconds
const m_p = 1.6726e-24 # grams
const k_B = 1.3806e-16 # erg/K

F = Float64

# except these are in natural units
function get_params(filename)
   params =  TOML.parsefile(filename)
   params["M_bary"] *= Msun
   params["M_tot"] *= Msun
   params["R_virial"] *= pc
   params["R_bary"] *= pc
   params["A_NFW"] = (log(1+params["c"]) - params["c"]/(1+params["c"]))
   params["Rs"] = params["R_virial"]/params["c"]
   params["t_end"] *= yr
   params["dt_min"] *= yr
   params["dt_max"] *= yr
   params["h_max"] *= pc
   params["rho_min"] *= m_p
   params["rho_max"] *= m_p
   params["rho_0"] *= m_p

   q = NamedTuple{ Tuple(Symbol.(keys(params))) }(values(params))
   return q
end



end
