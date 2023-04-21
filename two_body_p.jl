using Pkg

Pkg.activate(".")

using GalaxySim
using LinearAlgebra


function setup()
    params = Params("init/two_body_p.toml")

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
    println(ps)

    GalaxySim.Evolve.evolve!(ps, params)

    println(ps)
    println("finished")
end


run()



