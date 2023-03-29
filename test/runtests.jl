using GalaxySim
using Test
import ForwardDiff: derivative
using QuadGK
using LinearAlgebra

using GalaxySim

@testset "GalaxySim.jl" begin
    # Write your tests here.
end

@testset "density.jl" begin
    @test quadgk(x->GalaxySim.Density.W(x,1)*4π*x^2, 0, 5)[1] ≈ 1 atol = 0.0001

    # Test that the analytic derivative works
    s = 0
    for x in 0:0.05:1
        s += abs(derivative(x->GalaxySim.Density.W(x,1), x) - norm(GalaxySim.Density.∇W(x, 1)))
    end
    @test s ≈ 0 atol = 1e-5
end


@testset "init.jl" begin
    @test quadgk(x->GalaxySim.Init.ρ_bary(x)*4π*x^2, 0, GalaxySim.Init.R_virial)[1] ≈ GalaxySim.Init.M_bary atol=1


end
