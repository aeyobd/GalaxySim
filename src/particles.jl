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
    h::F = 3 * pc
    c::F = 0
    μ::F = 1.4
    Ω::F = 1

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


function interpolate(q::Particle, t, params)
    q1 = NParticle()
    dt = t - q.t

    @. q1.x = q.x + q.v*dt
    @. q1.v = q.v + q.dv*dt
    q1.m = q.m
    q1.ρ = max(q.ρ + q.dρ*dt, params.rho_min)
    q1.h = max(q.h + q.dh*dt, params.h_min)

    q1.dt = q.dt
    q1.μ = q.μ

    q1.u = max(q.u + q.du*dt, 0)
    q1.T = 2*q1.μ/(3R_ig) * q1.u
    q1.P = R_ig/q1.μ * q1.ρ * q1.T

    # TODO: q1.c = ...

    return q1
end


function Base.show(io::IO, p::Particle)
    println(io)

    @printf io "x (pc)          %8.2f\t%8.2f\t%8.2f\n"     (p.x/pc)...
    @printf io "v (km s⁻¹)      %8.2f\t%8.2f\t%8.2f\n"      (p.v/1e5)...
    @printf io "ρ (mₚ/cc)       %8.2e\t%8.2e\n"             (p.ρ/m_p) (p.dρ/m_p)
    @printf io "h (pc)          %8.2e\t%8.2e\n"             (p.h/pc) (p.dh/pc)
    @printf io "u (erg/g)       %8.2e\t%8.2e\n"              p.u p.du
    @printf io "T (K)           %8.3e\n"                     p.T
    @printf io "P (erg cm⁻³)    %8.2e\n"                     p.P
    @printf io "M_⋆ (M_⊙)       %0.2f\n"               (p.m_star/Msun)
    @printf io "t (yr)          %0.2f\n"               (p.t/yr)
    @printf io "dt (yr)         %0.2f\n"               (p.dt/yr)
    @printf io "ID              %d\n"                       p.id
    println(io)
end


function Base.show(io::IO, p::NParticle)
    println(io)

    @printf io "x (pc)          %8.2f\t%8.2f\t%8.2f\n"     (p.x/pc)...
    @printf io "v (km s⁻¹)      %8.2f\t%8.2f\t%8.2f\n"      (p.v/1e5)...
    @printf io "ρ (mₚ/cc)       %8.2e\n"             (p.ρ/m_p) 
    @printf io "h (pc)          %8.2e\n"             (p.h/pc)
    @printf io "u (erg/g)       %8.2e\n"              p.u 
    @printf io "T (K)           %8.3e\n"                     p.T
    @printf io "P (erg cm⁻³)    %8.2e\n"                     p.P
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
