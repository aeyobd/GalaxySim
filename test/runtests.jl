using GalaxySim
using Test
import ForwardDiff: derivative
using QuadGK

@testset "GalaxySim.jl" begin
    # Write your tests here.
end

@testset "calc.jl" begin
    include("../src/calc.jl")

    @test quadgk(x->calc.w(x)*4π*x^2, 0, 5)[1] ≈ 1 atol = 0.0001

    # Test that the analytic derivative works
    s = 0
    for x in 0.1:0.1:3
        s += abs(derivative(x->calc.W(1,x), x) - calc.∂W_∂h(1, x))
    end
    @test s ≈ 0 atol = 1e-5
end
