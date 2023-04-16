using GalaxySim
using Test
import ForwardDiff: derivative
using QuadGK
using LinearAlgebra

using GalaxySim


@testset "density.jl" begin
    # test that W is normalized
    @test quadgk(x->GalaxySim.Density.W(x,1)*4π*x^2, 0, 5)[1] ≈ 1 atol = 0.0001

    # Test that the analytic derivative works
    s = 0
    for i in 0:0.5:10
        x = i*pc
        p1 = Particle(x=zeros(3), v=zeros(3), m=1Msun, h=1pc)
        p2 = Particle(x=[x, 0, 0], v=zeros(3), m=1Msun)

        expected = derivative(x->GalaxySim.Density.W(x, 1pc), x) 
        predicted = norm(GalaxySim.Density.∇W(p1, p2))
        s += abs(expected - predicted)
    end
    @test s ≈ 0 atol = 1e-5


    # test density calculation for particles in a box
    params = GalaxySim.Constants.get_params("../src/static_eq.toml")
    N = 1000
    m = 1
    ps = Particle[]
    for i in 1:1000
        x = rand(3) * pc
        push!(ps, Particle(x=x, v=zeros(3), m=1))
    end

    tree = GalaxySim.Tree.make_tree(ps, params)
    ρ_predicted = 0
    for p in ps
        GalaxySim.Evolve.setup!(p, tree, params)
        ρ_predicted += p.ρ
    end

    ρ_predicted /= N
    ρ_expected = N*m/pc^3

    @test ρ_predicted ≈ ρ_expected atol = 1e-2


end


@testset "init.jl" begin
    params = GalaxySim.Constants.get_params("../src/static_eq.toml")

    # test the density integral
    @test GalaxySim.Init.∫ρ_bary(1000*params.R_virial, params) ≈ 1 atol = 1e-3
    @test GalaxySim.Init.∫ρ_bary(0., params) ≈ 0 atol = 1e-3


    # test density calculation for initial distribution
    ps = GalaxySim.Init.rand_particles(params)
    tree = GalaxySim.Tree.make_tree(ps, params)

    s = 0.
    for p in ps
        GalaxySim.Evolve.setup!(p, tree, params)
        predicted = p.ρ
        r = norm(p.x)
        expected = GalaxySim.Init.ρ_bary(r, params)
        s += expected - predicted
    end
    s /= params.N

    @test s ≈ 0 atol=1e-2

    @test sum(p.m for p in ps) ≈ params.M_bary rtol = 1e-4
    # @test quadgk(x->GalaxySim.Init.ρ_bary(x)*4π*x^2, 0, GalaxySim.Init.R_virial)[1] ≈ GalaxySim.Init.M_bary atol=1
end
