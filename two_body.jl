# Created 
#
#
# Author: Daniel Boyea (boyea.2@osu.edu)
#
# This contains a simple test of the local gravity 
# Two particles are created which co-orbit the origin
# The energy is conserved to <1% over 10Gyr, and 
# the path only deviates slightly.
# The smoothing length is intentionally set smaller than the 
# distance to test just the n-body portion. 
#
# I have also tested elliptical orbits, which also conserve energy, 
# which is also an easy modification of the code (just change v)
#
# Essentially, gravity performs with minimal errors.


using Pkg

Pkg.activate(".")

using GalaxySim
using LinearAlgebra


function setup()
    params = Params("init/two_body.toml")

    # initial radii, masses, and velocities
    r_0 = 10pc
    M_0 = 10Msun
    v_0 = sqrt(G*M_0/r_0)/2 # divide by two because of reduced mass

    v_vec_1 = [0,1,0] * v_0
    v_vec_2 = [0,-1,0] * v_0

    r_vec_1 = [-1,0,0] * r_0
    r_vec_2 = [1,0,0] * r_0

    p1 = Particle(
        x = r_vec_1,
        v = v_vec_1,
        m = M_0,
        dt = params.dt_min
    )

    p2 = Particle(
        x = r_vec_2,
        v = v_vec_2,
        m = M_0,
        dt = params.dt_min
    )

    # return a list of the particles
    ps = [p1, p2]
    return params, ps
end


function run()
    params, ps = setup()
    evolve!(ps, params)
end


run()



