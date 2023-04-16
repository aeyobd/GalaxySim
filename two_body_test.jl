using Pkg

Pkg.activate(".")

using GalaxySim
using LinearAlgebra


function setup()
    params = GalaxySim.Constants.get_params("src/two_body.toml")

    r_0 = 10pc
    M_0 = 10Msun
    v_0 = sqrt(G*M_0/r_0)

    v_vec_1 = [0,0,0] * v_0
    v_vec_2 = [0,-1,0] * v_0

    r_vec_1 = [0,0,0] * r_0
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
        m = 0.001M_0,
        dt = params.dt_min
    )

    ps = [p1, p2]
    return params, ps
end


function run()
    params, ps = setup()
    println(ps)
    println("initial energy")
    e = sum(0.5*p.m*norm(p.v)^2 for p in ps)
    println(e)

    t_int = 0
    file = open("two_body.dat", "w")
    println(file, "t,x1,y1,x2,y2")
    dt = params.dt_min

    while t_int < params.t_end/params.dt_min
        t = t_int * params.dt_min

        GalaxySim.Evolve.update_particles!(ps, t, dt, params)
        x1 = ps[1].x[1] / pc
        y1 = ps[1].x[2] / pc
        x2 = ps[2].x[1] / pc
        y2 = ps[2].x[2] / pc

        println(file, "$(t/yr),$x1,$y1,$x2,$y2")
        print("t = $(t/yr) years\r")

        t_int += GalaxySim.Evolve.get_dt(ps, params)
    end


    close(file)
    println()
    println(ps)
    println("final energy")
    e = sum(0.5*p.m*norm(p.v)^2 for p in ps)
    println(e)

    println("finished")
end


run()



