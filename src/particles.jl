module Particles
export Particle, NParticle, AParticle, interpolate
using StaticArrays

using Printf

using ..Constants


abstract type AParticle end


Base.@kwdef mutable struct Particle <: AParticle
    x::MVector{3,F}
    v::MVector{3,F}
    m::F


    ρ::F = 1 * m_p
    h::F = 5 * pc
    c::F = 0
    μ::F = 1.4

    # stars::Vector = []
    neighbors::Vector{AParticle} = []

    T::F = 0.1
    u::F = 3/2*R_ig/μ*T
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


"""
Used as neighbors for each particle
"""
Base.@kwdef mutable struct NParticle <: AParticle
    x::MVector{3,F} = zeros(3)
    v::MVector{3,F} = zeros(3)

    m::F = 0
    ρ::F = 1 * m_p
    h::F = 0

    c::F = 0
    μ::F = 1.4  # mean molecular mass in m_p

    T::F = 0.1
    u::F = 3/2*R_ig*T
    P::F = R_ig/μ * ρ * T
    dt::F = 0

    # weights
    w::F = 0
    dw::MVector{3, F} = zeros(3)
end


function interpolate(q::Particle, t)
    q1 = NParticle()
    dt = t - q.t

    @. q1.x = q.x + q.v*dt
    @. q1.v = q.v + q.dv*dt
    q1.m = q.m
    q1.ρ = q.ρ + q.dρ*dt
    q1.h = q.h + q.dh*dt

    q1.dt = q.dt
    q1.μ = q.μ

    q1.u = q.u + q.du*dt
    q1.T = 2*q1.μ/(3R_ig) * q1.u
    q1.P = R_ig/q1.μ * q1.ρ * q1.T

    return q1
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



function Base.show(io::IO, ::MIME"text/plain", p::Particle)
    print(io, "particle at ")
    @printf io "(%8.1f,%8.1f,%8.1f)    "     (p.x/pc)...
    @printf io "T=%8.1e\t"                     p.T
    @printf io "ID=%d"                       p.id
    return io
end


end
