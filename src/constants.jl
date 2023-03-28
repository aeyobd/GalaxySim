module Constants

export G, R, yr, pc, Msun
export ρ_max, ρ_min, η
export μ, ρ0, K0, ϵ_eff


const R = 8.314e7 #erg/K/mol

const G = 6.67e-8 # We use CGS for everything
const Msun = 1.989e33
const pc = 3.086e18
const yr = 3.15e7


# limits on density to help solve
const ρ_max = 1e-19
const ρ_min = 1e-26
const η = 0.1

dt::Real = 1e3yr


# mean mass in mp per particle of gas
μ = 1

# constants for state equtions
ρ0 = 1.169e-25
K0 = 4e7
ϵ_eff = 0.01

end
