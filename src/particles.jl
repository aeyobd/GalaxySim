module Particles
export Particle
using StaticArrays

using Printf

using ..Constants

const μ = 1.4


Base.@kwdef mutable struct Particle
    x::MVector{3,F}
    v::MVector{3,F}
    m::F


    ρ::F = 1 * m_p
    h::F = 5 * pc
    c::F = 0
    μ::F = 1.4

    # stars::Vector = []
    neighbors::Vector{Particle} = []

    T::F = 0.1
    u::F = 3/2*R_ig*T
    P::F = R_ig/μ * ρ * T

    ρ_gas::F = ρ
    m_gas::F = m
    m_star::F = 0

    t::F = 0
    dt::F = 0

    id::Int = 0

    # derivatives :)
    dv::MVector{3, F} = zeros(3)
    dm_star::F = 0
    dρ::F = 0
    du::F = 0
    dh::F = 0

    # for tracking each process
    du_P::F = 0
    du_cond::F = 0

    dv_P::MVector{3, F} = zeros(3)
    dv_G::MVector{3, F} = zeros(3)
    dv_DM::MVector{3, F} = zeros(3)
end


Base.@kwdef mutable struct NParticle
    x::MVector{3,F}
    v::MVector{3,F}
    m::F
    ρ::F = 1 * m_p
    h::F = 5 * pc
    c::F = 0
    μ::F = 1.4  # mean molecular mass in m_p
    T::F = 0.1
    u::F = 3/2*R_ig*T
    P::F = R_ig/μ * ρ * T
end

function interpolate(p::Particle, symbol, t)
    f = getproperty(p, :symbol)
    df_dt = getproperty(p, :"d$(symbol)_dt")
    dt = t - p.t
    return f + df_dt * dt
end


function Base.show(io::IO, p::Particle)
    println(io)

    @printf io "x (pc)\t%8.2f\t%8.2f\t%8.2f\n"     (p.x/pc)...
    @printf io "v (km/s)\t%8.2f\t%8.2f\t%8.2f\n"     (p.v/1e5)...
    @printf io "ρ (m_p/cc)\t%8.2e\n"                     p.ρ
    @printf io "h (pc)\t%8.2e\n"                     (p.h/pc)
    @printf io "T (K) \t%8.3e\n"                     p.T
    @printf io "u\t%8.2e\n"                     p.u
    @printf io "P\t%8.2e\n"                     p.P
    @printf io "M_stars\t%0.2f\n"               (p.m_star/Msun)
    @printf io "ID\t%d\n"                       p.id
    println(io)
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
