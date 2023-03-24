
module particle
using Printf

const G = 6.67e-8 # We use CGS for everything
const Msun = 1.989e33
const pc = 3.086e18
const yr = 3.15e7
const R = 8.314e7 #erg/K/mol
const μ = 1 #g/mol

Base.@kwdef mutable struct Particle
    x::Vector
    v::Vector
    m::Real

    ρ::Real = 0.1 * Msun/pc^3
    h::Real = 100 * pc
    W::Vector = []
    # 
    # u::Vector = []
    # stars::Vector = []
    neighbors::Vector = []
    distances::Vector = []
    T::Real = 1000
    u::Real = 3/2*R*T
    P::Real = R/μ * ρ * T
    id::Real = 0
end


function Base.show(io::IO, p::Particle)
    println(io)

    @printf io "x\t%0.2f\t%0.2e\t%0.2e\n"     (p.x/pc)...
    @printf io "v\t%0.2e\t%0.2e\t%0.2e\n"     p.v...
    @printf io "ρ\t%0.2e\n"     p.ρ
    @printf io "T\t%0.3e\n"     p.T
    @printf io "u\t%0.2e\n"     p.u
    @printf io "P\t%0.2e\n"     p.P
    @printf io "ID\t%d\n"       p.id
end

end
