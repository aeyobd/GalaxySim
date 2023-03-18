module init

import LinearAlgebra: norm, normalize
import SpecialFunctions: gamma

const G = 4.4985e-15


M_tot = 1e8
M_bary = 1e6 #M_sun
c = 10
R_virial = 1e3
Rs = R_virial/c

A_NFW = (log(1+c) - c/(1+c))
ρc = M_tot / ( 4π*R_virial^3 * A_NFW)

function ρ_DM(r)
    x = r / Rs
    return ρc / (x * (1+x^2) )
end

v0_virial = √(G*M_tot/R_virial)
function v_virial(r)
    x = r/Rs
    return v0_virial * √( 1/x * (log(1+c*x) - (c*x)/(1+c*x)) /A_NFW )
end

function a_DM(r::Real)
    G * (M_tot/A_NFW) * (r/(r+Rs) - log(1+r/Rs))/r^2
end

function a_DM(x::Vector)
    return a_DM(norm(x)) * normalize(x)
end


const β=5
const k_ρ = gamma(β/2)/(gamma((β-3)/2) * π^(3/2))
const Rp = Rs
function ρ_bary(r)
    3*M_bary/(4π * Rp^3) / (1 + (r/Rp)^2)^(β/2)
end


end
