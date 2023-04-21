module Constants

export G, R_ig, yr, pc, Msun, m_p, k_B
export F


# which type to use for floats?
F = Float64


# All constants in cgs
const R_ig = 8.314e7   #erg/K/mol
const G = 6.67e-8 
const Msun = 1.989e33  # grams
const pc = 3.086e18    # cm
const yr = 3.15e7      # seconds
const m_p = 1.6726e-24 # grams
const k_B = 1.3806e-16 # erg/K


end
