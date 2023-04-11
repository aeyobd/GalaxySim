module Particles
export Particle, DParticle

using Printf

using ..Constants

const μ = 1.4

Base.@kwdef mutable struct Particle
    x::Vector
    v::Vector
    m::Real

    ρ::Real = 0.1 * Msun/pc^3
    h::Real = 100 * pc
    W::Vector = []
    # 
    # stars::Vector = []
    #
    neighbors::Vector = []
    distances::Vector = []
    T::Real = 10000
    u::Real = 3/2*R_ig*T
    P::Real = R_ig/μ * ρ * T
    id::Real = 0

    ρgas::Real = ρ
    mgas::Real = m
    mstar::Real = 0
end


function Base.show(io::IO, p::Particle)
    println(io)

    @printf io "x\t%0.2f\t%0.2e\t%0.2e\n"     (p.x/pc)...
    @printf io "v\t%0.2e\t%0.2e\t%0.2e\n"     p.v...
    @printf io "ρ\t%0.2e\n"                     p.ρ
    @printf io "T\t%0.3e\n"                     p.T
    @printf io "u\t%0.2e\n"                     p.u
    @printf io "P\t%0.2e\n"                     p.P
    @printf io "M_stars\t%0.2e\n"               p.mstar
    @printf io "ID\t%d\n"                       p.id
end

function Base.copy(p::Particle)
    p1 = particle(x=copy(p.x),
                  v=copy(p.v),
                  m=copy(p.m),
                  ρ=copy(p.ρ))
    return p1
end

function Base.show(io::IO, ::MIME"text/plain", p::Particle)
    print(io, "particle at ")
    @printf io "(%8.1f,%8.1f,%8.1f)    "     (p.x/pc)...
    @printf io "T=%8.1e\t"                     p.T
    @printf io "ID=%d"                       p.id
    return io
end


end
