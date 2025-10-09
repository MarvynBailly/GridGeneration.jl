using Test
using GridGeneration

@testset "Projection Functions" begin
    @testset "2D to 1D Projection - Straight Line" begin
        # Horizontal line from (0,0) to (1,0)
        curve = [range(0, 1, length=10) zeros(10)]'
        xs = GridGeneration.ProjectBoundary2Dto1D(curve)
        
        @test length(xs) == 10
        @test xs[1] ≈ 0.0 atol=1e-10
        @test xs[end] ≈ cumsum(curve, dims=1)[end] atol=1e-10 # End should be the length of the line
    end

    @testset "2D to 1D Projection - Vertical Line" begin
        # Vertical line from (0,0) to (0,1)
        curve = [zeros(10) range(0, 1, length=10)]'
        xs = GridGeneration.ProjectBoundary2Dto1D(curve)
        
        @test length(xs) == 10
        @test xs[1] ≈ 0.0 atol=1e-10
        @test xs[end] ≈ cumsum(curve, dims=1)[end] atol=1e-10 # End should be the length of the line    end
    end

    @testset "1D to 2D Projection Round-trip" begin
        # Create a curve and project it
        N = 8
        original = [range(0, 1, length=N) range(0, 0.5, length=N)]'
        xs = GridGeneration.ProjectBoundary2Dto1D(original)
        
        # Project back with same number of points
        reconstructed = GridGeneration.ProjectBoundary1Dto2D(original, xs)
        
        @test size(reconstructed) == size(original)
        # Should be very close to original since we used same parametrization
        @test reconstructed ≈ original atol=1e-10
    end
end
