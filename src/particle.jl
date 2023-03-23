
module particle


Base.@kwdef mutable struct Particle
    x::Vector
    v::Vector
    m::Real

    ρ::Real = 0.1
    h::Real = 100
    W::Vector = []
    # 
    # u::Vector = []
    # stars::Vector = []
    neighbors::Vector = []
    distances::Vector = []
    P::Real = 0
    T::Real = 10
    u::Real = 0
end


end
