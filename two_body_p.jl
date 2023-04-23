# two_body_p.jl
#
# A test of the pressure forces between two bodies
#
# Created 19-04-2023
# Author Daniel Boyea (boyea.2@osu.edu)
#
# Honeslty, this was a test just to make sure I didn't
# add a minus sign in the pressure equations. 
# See the associated plots, the pressure and density smoothly increase,
# until the particles reach their nearest point, then the particles 
# move away from eachother.
# Energy is also successfuly conserved.


using Pkg

Pkg.activate(".")

using GalaxySim
using LinearAlgebra


function setup()
    params = Params("init/two_body_p.toml")

    # create two particles 
    # with a head on collision
    r_0 = 10pc
    M_0 = 10Msun
    v_0 = 3*sqrt(G*M_0/r_0)

    v_vec_1 = [1,0,0] * v_0
    v_vec_2 = [-1,0,0] * v_0

    r_vec_1 = [-1,0,0] * r_0
    r_vec_2 = [1,0,0] * r_0

    p1 = Particle(
        x = r_vec_1,
        v = v_vec_1,
        m = M_0,
        T = 10
    )

    p2 = Particle(
        x = r_vec_2,
        v = v_vec_2,
        m = M_0,
        T=10
    )

    ps = [p1, p2]
    return params, ps
end


function run()
    params, ps = setup()
    GalaxySim.Evolve.evolve!(ps, params)
end


run()



