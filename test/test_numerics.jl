using Test
using GridGeneration
using Statistics

@testset "Numerical Methods" begin
    @testset "Linear Interpolation" begin
        @testset "Simple Linear Function" begin
            x = [0.0, 1.0, 2.0, 3.0]
            y = [0.0, 2.0, 4.0, 6.0]  # y = 2x
            
            interp = GridGeneration.LinearInterpolate(x, y)
            
            @test interp(0.0) ≈ 0.0
            @test interp(1.0) ≈ 2.0
            @test interp(0.5) ≈ 1.0
            @test interp(1.5) ≈ 3.0
            @test interp(2.5) ≈ 5.0
        end

        @testset "Constant Function" begin
            x = [0.0, 1.0, 2.0]
            y = [5.0, 5.0, 5.0]
            
            interp = GridGeneration.LinearInterpolate(x, y)
            
            @test interp(0.0) ≈ 5.0
            @test interp(0.5) ≈ 5.0
            @test interp(1.5) ≈ 5.0
        end
    end

    @testset "ODE Solver - Analytic Method" begin
        # Test with constant metric (should give uniform spacing)
        xs = range(0, 1, length=20)
        M(x) = 1.0  # Constant metric
        
        sol = GridGeneration.SolveODE(M, collect(xs); solver=:analytic)
        
        @test length(sol) == length(xs)
        @test sol[1] ≈ 0.0 atol=1e-10
        # With constant metric, spacing should be approximately uniform
        spacings = diff(sol)
        @test std(spacings) < 0.1  # Low variation in spacing
    end

    @testset "ODE Solver - Fixed N" begin
        xs = range(0, 1, length=15)
        M(x) = 1.0 + x  # Linearly varying metric
        
        N_target = 25
        sol = GridGeneration.SolveODEFixedN(M, collect(xs), N_target; solver=:analytic)
        
        @test length(sol) == N_target
        @test sol[1] ≈ 0.0 atol=1e-10
        @test sol[end] ≈ 1.0 atol=1e-3
    end

    @testset "Central Difference" begin
        # Test on f(x) = x^2, f'(x) = 2x
        f(x) = x^2
        xs = range(0, 1, length=21)
        
        df = GridGeneration.CentralDiff(f, xs)
        
        # Check at interior points (central diff most accurate there)
        mid_idx = 11
        @test df[mid_idx] ≈ 2 * xs[mid_idx] atol=0.1
    end
end
