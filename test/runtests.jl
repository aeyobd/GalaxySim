# runtests.jl
#
# This file contains some basic unit tests
#
# Created 16-04-2023
# Author Daniel Boyea (boyea.2@osu.edu)
#
# To run this file, 
# activate this package from the julia command line
# by entering package mode (by typing `]`)
# and typing `activate .` from the main project directory,
# then type `test`
#
# I'm sure there's another way but I don't know it yet :/


using GalaxySim
using Test
import ForwardDiff: derivative
using QuadGK
using LinearAlgebra

using GalaxySim


@testset "density.jl" begin
    # test that W adds up to 1
    @test quadgk(x->GalaxySim.Density.W(x,2)*4π*x^2, 0, 5)[1] ≈ 1 atol = 0.0001

    # Test that the analytic derivatives are as expected
    s = 0 # residual sums
    s1 = 0
    for i in 0.5:0.5:10 
        x = i*pc # just test over a reasonable range of xs
        p1 = Particle(x=zeros(3), v=zeros(3), m=1Msun, h=1pc)
        p2 = Particle(x=[x, 0, 0], v=zeros(3), m=1Msun)

        expected = derivative(x->GalaxySim.Density.W(x, 1pc), x) 
        predicted = norm(GalaxySim.Density.∇W(p1, p2))
        s += abs(expected - predicted)

        expected = derivative(x->GalaxySim.Density.W(GalaxySim.Density.dist(p1, p2), x), x) 
        predicted = GalaxySim.Density.dW_dh(p1, p2)
        s1 += abs(expected - predicted)
    end

    @test s ≈ 0 atol = 1e-5
    @test s1 ≈ 0 atol = 1e-5


end


