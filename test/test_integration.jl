using Test
using GridGeneration

function unit_square_grid(N; width=1.0)
    # TFI convention: every edge runs left→right / bottom→top
    top = [range(0, width, length=N) ones(N)]
    right = [fill(width, N) range(0, 1, length=N)]
    bottom = [range(0, width, length=N) zeros(N)]
    left = [zeros(N) range(0, 1, length=N)]
    return TFI([top, right, bottom, left])
end

@testset "Integration Tests" begin
    M(x, y) = [100.0, 100.0]

    @testset "GenerateGrid - All Stages Disabled" begin
        initialGrid = unit_square_grid(8)
        params = SimParams(useSplitting = false, useEdgeSolver = false, useSmoothing = false)

        result = GenerateGrid(initialGrid, [], [], M; params=params)

        @test length(result) == 6
        @test length(result.blocks) == 1
        @test result.blocks[1] == initialGrid
        @test result.smoothBlocks === result.blocks
        @test isempty(result.finalErrors)
    end

    @testset "GenerateGrid - Splitting Only" begin
        initialGrid = unit_square_grid(12; width=2.0)
        params = SimParams(useSplitting = true, splitLocations = [[6], [6]],
                           useEdgeSolver = false, useSmoothing = false)

        smoothBlocks, blocks, bndInfo, interInfo, finalErrors, finalIterations =
            GenerateGrid(initialGrid, [], [], M; params=params)

        @test length(blocks) == 4
        @test smoothBlocks === blocks
    end

    @testset "GenerateGrid - Edge Solving" begin
        initialGrid = unit_square_grid(8)
        for solver in (:analytic, :numeric)
            params = SimParams(useSplitting = false, useEdgeSolver = true,
                               boundarySolver = solver, useSmoothing = false)

            result = GenerateGrid(initialGrid, [], [], M; params=params)

            @test length(result.blocks) == 1
            blk = result.blocks[1]
            @test size(blk, 1) == 2
            # Constant metric of 100 on a unit square: about 1/sqrt(1/100) = 10 points per side
            @test 8 <= size(blk, 2) <= 12
            @test 8 <= size(blk, 3) <= 12
            # Corners are preserved
            @test blk[:, 1, 1] ≈ initialGrid[:, 1, 1]
            @test blk[:, end, end] ≈ initialGrid[:, end, end]
        end
    end

    @testset "GenerateGrid - Full Pipeline With Smoothing" begin
        initialGrid = unit_square_grid(10)
        params = SimParams(useSplitting = false, useEdgeSolver = true, boundarySolver = :analytic,
                           useSmoothing = true, smoothMethod = :ellipticSS,
                           elliptic = EllipticParams(max_iter = 2000, useBottomWall = false))

        result = GenerateGrid(initialGrid, [], [], M; params=params)

        @test length(result.smoothBlocks) == 1
        @test size(result.smoothBlocks[1]) == size(result.blocks[1])
        @test length(result.finalErrors) == 1
        @test result.finalIterations[1] >= 1
        # A uniform grid is already a solution of the unforced equations
        @test result.smoothBlocks[1] ≈ result.blocks[1] atol=1e-8
    end
end
