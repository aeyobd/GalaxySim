# GalaxySim.jl
#
# This file simply exports all the other files
#
# Creatd 15-03-2023
# Revised to add other files
# Author Daniel Boyea (boyea.2@osu.edu)


module GalaxySim

export yr, pc, Msun, G, m_p
export Particle, Params, evolve!

include("constants.jl")
include("params.jl")
include("particles.jl")
include("density.jl")
include("gravity.jl")
include("gal_files.jl")
include("physics.jl")
include("evolve.jl")


using .Constants
using .Particles
using .Density
using .Gravity
using .GalFiles
using .Physics
using .Evolve


end
