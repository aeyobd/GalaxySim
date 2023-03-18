
module particle


Base.@kwdef mutable struct Particle
    x
    v
    m

    # ρ::Real = NaN
    # h::Real = NaN
    # W::Vector = NaN
    # 
    # u::Vector = []
    # stars::Vector = []
    # neighbors::Vector = []
    # distances::Vector = []
end

end
