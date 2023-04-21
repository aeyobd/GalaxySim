using Pkg

Pkg.activate(".")

using GalaxySim
using LinearAlgebra


function setup()
    params = Params("init/two_body.toml")

    r_0 = 10pc
    M_0 = 10Msun
    v_0 = sqrt(G*M_0/r_0)/4

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

    ps = [p1, p2]
    return params, ps
end


function run()
    params, ps = setup()
    evolve!(ps, params)
    println("finished")
end


run()



