using Test
using GridGeneration

@testset "Integration Tests" begin
    @testset "Complete Grid Generation - No Splitting" begin
        # Create initial grid
        N = 8
        top = [range(0, 1, length=N) ones(N)]
        right = [ones(N) range(1, 0, length=N)]
        bottom = [range(1, 0, length=N) zeros(N)]
        left = [zeros(N) range(0, 1, length=N)]
        
        initialGrid = TFI([top, right, bottom, left])
        
        # Constant metric
        M(x, y) = [1.0, 1.0]
        
        # Simple parameters - no splitting, no smoothing
        params = SimParams(
            useSplitting = false,
            useEdgeSolver = false,
            useSmoothing = false
        )
        
        # Should run without error
        result = @test_nowarn GenerateGrid(
            initialGrid, [], [], M; params=params
        )
        
        @test length(result) == 6  # Returns 6 values
    end

    @testset "Complete Grid Generation - With Splitting" begin
        # Create initial grid
        N = 12
        top = [range(0, 2, length=N) ones(N)]
        right = [fill(2, N) range(1, 0, length=N)]
        bottom = [range(2, 0, length=N) zeros(N)]
        left = [zeros(N) range(0, 1, length=N)]
        
        initialGrid = TFI([top, right, bottom, left])
        
        # Constant metric
        M(x, y) = [1.0, 1.0]
        
        # Parameters with splitting
        params = SimParams(
            useSplitting = true,
            splitLocations = [[6], [6]],
            useEdgeSolver = false,
            useSmoothing = false
        )
        
        # Should run without error
        result = @test_nowarn GenerateGrid(
            initialGrid, [], [], M; params=params
        )
        
        smoothBlocks, blocks, bndInfo, interInfo, finalErrors, finalIterations = result
        
        # Should have created 4 blocks
        @test length(blocks) == 4
    end

    @testset "Grid Generation - With Edge Solving" begin
        # Create smaller grid for faster testing
        N = 8
        top = [range(0, 1, length=N) ones(N)]
        right = [ones(N) range(1, 0, length=N)]
        bottom = [range(1, 0, length=N) zeros(N)]
        left = [zeros(N) range(0, 1, length=N)]
        
        initialGrid = TFI([top, right, bottom, left])
        
        # Simple metric
        M(x, y) = [1.0, 1.0]
        
        # Parameters with edge solving
        params = SimParams(
            useSplitting = false,
            useEdgeSolver = true,
            boundarySolver = :analytic,
            useSmoothing = false
        )
        
        # Should run without error
        result = @test_nowarn GenerateGrid(
            initialGrid, [], [], M; params=params
        )
        
        @test length(result) == 6
    end
end
