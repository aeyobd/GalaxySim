
module particle


Base.@kwdef mutable struct Particle
    x::Vector
    v::Vector
    m::Real

    ρ::Real = NaN
    h::Real = NaN
    W::Vector = []
    # 
    # u::Vector = []
    # stars::Vector = []
    neighbors::Vector = []
    distances::Vector = []
end


end
