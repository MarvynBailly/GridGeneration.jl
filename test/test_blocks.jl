using Test
using GridGeneration

@testset "Block Operations" begin
    @testset "Block Splitting" begin
        # Create a simple 10×10 grid
        N = 10
        top = [range(0, 1, length=N) ones(N)]
        right = [ones(N) range(1, 0, length=N)]
        bottom = [range(1, 0, length=N) zeros(N)]
        left = [zeros(N) range(0, 1, length=N)]
        
        initialGrid = TFI([top, right, bottom, left])
        
        # Split at middle indices
        splitLocations = [[5], [5]]
        bndInfo = []
        interInfo = []
        
        blocks, newBndInfo, newInterInfo = GridGeneration.SplitBlock(
            initialGrid, splitLocations, bndInfo, interInfo
        )
        
        # Should create 4 blocks (2×2 split)
        @test length(blocks) == 4
        
        # Each block should be a 3D array
        for block in blocks
            @test ndims(block) == 3
            @test size(block, 1) == 2  # x and y coordinates
        end
        
        # Interface information should be created
        @test !isempty(newInterInfo)
    end

    @testset "GetNeighbors" begin
        # Create simple interface structure
        interInfo = [
            Dict("blockA" => 1, "blockB" => 2,
                 "start_blkA" => [5, 1, 1], "end_blkA" => [5, 5, 1]),  # Horizontal
            Dict("blockA" => 1, "blockB" => 3,
                 "start_blkA" => [1, 5, 1], "end_blkA" => [5, 5, 1])   # Vertical
        ]
        
        # Block 1 should have neighbors in both directions
        neighbors_h = GridGeneration.GetNeighbors(1, interInfo, 1; include_start=true)
        neighbors_v = GridGeneration.GetNeighbors(1, interInfo, 2; include_start=true)
        
        @test 1 in neighbors_h
        @test 1 in neighbors_v
        @test 2 in neighbors_h  # Horizontal neighbor
        @test 3 in neighbors_v  # Vertical neighbor
    end
end
