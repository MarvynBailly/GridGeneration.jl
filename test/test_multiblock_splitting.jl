using Test
using GridGeneration

# Comprehensive tests for multi-block splitting with automatic propagation
# Tests verify:
# 1. Correct number of resulting blocks after propagation
# 2. Exact block sizes match expected dimensions
# 3. Boundary information correctly updated with all affected blocks
# 4. Interface information includes all expected connections
# 5. No spurious interfaces created

@testset "Multi-Block Splitting" begin
    @testset "2×2 Grid with Propagating Horizontal Split" begin
        # Create a 2×2 multi-block grid matching the task.md scenario
        # Initial 2×2 layout: 
        #   ┌─────┬─────┐
        #   │  3  │  4  │  <- Top row
        #   ├─────┼─────┤
        #   │  1  │  2  │  <- Bottom row
        #   └─────┴─────┘
        #
        # Block 1 (lower-left) and Block 2 (lower-right) share a VERTICAL interface (side-by-side)
        # Block 1 (lower-left) and Block 3 (upper-left) share a HORIZONTAL interface (stacked vertically)
        #
        # When we split block 1 horizontally (j-direction split: [[], [split]]):
        # - This creates a horizontal line within block 1
        # - Block 2 is to the right of block 1, sharing a vertical interface
        # - For one-to-one grid correspondence across the vertical interface, 
        #   block 2 MUST also be split at the same j-location
        # - Blocks 3 and 4 remain unchanged (not affected by horizontal split)
        #
        # Final 6-block layout after split and propagation:
        #   ┌─────┬─────┐
        #   │  5  │  6  │  <- Top row (block 3 → 5; block 4 → 6)
        #   ├─────┼─────┤
        #   │  2  │  4  │  <- Bottom row (block 1 split → 1,2; block 2 split → 3,4)
        #   ├─────┼─────┤
        #   │  1  │  3  │  
        #   └─────┴─────┘
        #
        # Block 1 (old, lower-left) splits → blocks 1 (bottom) + 2 (top)
        # Block 2 (old, lower-right) splits → blocks 3 (bottom) + 4 (top) [propagated horizontally]
        # Block 3 (old, upper-left) → block 5 (unchanged, renumbered)
        # Block 4 (old, upper-right) → block 6 (unchanged, renumbered)
        
        # Start with a single 10×10 grid and split it into 2×2
        N = 10
        top = [collect(range(0.0, 1.0, length=N)) ones(N)]
        right = [ones(N) collect(range(1.0, 0.0, length=N))]
        bottom = [collect(range(1.0, 0.0, length=N)) zeros(N)]
        left = [zeros(N) collect(range(0.0, 1.0, length=N))]
        
        initialGrid = GridGeneration.TFI([top, right, bottom, left])
        
        # Split into 2×2 blocks
        splitLocs = [[5], [5]]
        bndInfo = [
            Dict("name" => "top", "faces" => [Dict("start" => [1, N, 1], "end" => [N, N, 1])]),
            Dict("name" => "right", "faces" => [Dict("start" => [N, 1, 1], "end" => [N, N, 1])]),
            Dict("name" => "bottom", "faces" => [Dict("start" => [1, 1, 1], "end" => [N, 1, 1])]),
            Dict("name" => "left", "faces" => [Dict("start" => [1, 1, 1], "end" => [1, N, 1])])
        ]
        
        blocks, bndInfo, interInfo = GridGeneration.SplitBlock(initialGrid, splitLocs, bndInfo, [])
        
        @test length(blocks) == 4
        @test length(interInfo) == 4  # 2 vertical + 2 horizontal
        
        # Now split block 1 horizontally (in j-direction) at midpoint
        # This should propagate to block 3 (connected via horizontal interface)
        block1_nj = size(blocks[1], 3)
        mid_j = div(block1_nj, 2)
        
        splitRequests = Dict(
            1 => [Int[], [mid_j]]  # Horizontal split in block 1
        )
        
        newBlocks, newBndInfo, newInterInfo = GridGeneration.SplitMultiBlock(
            blocks, splitRequests, bndInfo, interInfo
        )
        
        # Expected result: 6 blocks
        # Original blocks 1 and 2 each split into 2 -> blocks 3 and 4 remain unchanged
        # New numbering (row-major from SplitBlock):
        # Block 1 becomes: 1 (bottom part), 2 (top part)
        # Block 2 becomes: 3 (bottom part), 4 (top part) [propagated]
        # Block 3 remains:  5 (unchanged, renumbered)
        # Block 4 remains:  6 (unchanged, renumbered)
        @test length(newBlocks) == 6
        
        # Verify blocks maintain grid structure
        for block in newBlocks
            @test ndims(block) == 3
            @test size(block, 1) == 2
        end
        
        # Expected block sizes
        # Block 1 was 5x5, split at j=2 gives: 5x2 and 5x4
        # Block 2 was 6x5, split at j=2 gives: 6x2 and 6x4 (propagated)
        # Block 3 was 5x6, unchanged
        # Block 4 was 6x6, unchanged
        expected_sizes = [
            (2, 5, 2),  # Block 1 (new ID 1) - bottom part of old block 1
            (2, 5, 4),  # Block 2 (new ID 2) - top part of old block 1
            (2, 6, 2),  # Block 3 (new ID 3) - bottom part of old block 2 (propagated split)
            (2, 6, 4),  # Block 4 (new ID 4) - top part of old block 2 (propagated split)
            (2, 5, 6),  # Block 5 (new ID 5) - old block 3 unchanged
            (2, 6, 6)   # Block 6 (new ID 6) - old block 4 unchanged
        ]
        
        for (idx, expected_size) in enumerate(expected_sizes)
            @test size(newBlocks[idx]) == expected_size
        end
        
        # Check boundary information
        # Should have 4 named boundaries: top, right, bottom, left
        @test length(newBndInfo) == 4
        
        bnd_names = Set([b["name"] for b in newBndInfo])
        @test bnd_names == Set(["top", "right", "bottom", "left"])
        
        # Bottom boundary: After j-split, only jSeg=1 touches the original bottom edge
        # Block 1: bottom part of old block 1 (jSeg=1)
        # Block 3: bottom part of old block 2 (jSeg=1, propagated split)
        bottom_bnd = filter(b -> b["name"] == "bottom", newBndInfo)[1]
        bottom_blocks = sort([f["block"] for f in bottom_bnd["faces"]])
        @test bottom_blocks == [1, 3]
        
        # Top boundary: Unchanged blocks 5 and 6 (old blocks 3 and 4)
        # Block 5: old block 3 (top-left, unchanged)
        # Block 6: old block 4 (top-right, unchanged)
        # Top boundary: Unchanged blocks 5 and 6 (old blocks 3 and 4)
        # Block 5: old block 3 (top-left, unchanged)
        # Block 6: old block 4 (top-right, unchanged)
        top_bnd = filter(b -> b["name"] == "top", newBndInfo)[1]
        top_blocks = sort([f["block"] for f in top_bnd["faces"]])
        @test top_blocks == [5, 6]
        
        # Left boundary: All blocks on the left edge (i=1)
        # After j-split, both sub-blocks of old block 1 have i=1 edge
        # Block 1: bottom part of old block 1 (jSeg=1, left edge at i=1)
        # Block 2: top part of old block 1 (jSeg=2, left edge at i=1)
        # Block 5: old block 3 (left edge at i=1)
        left_bnd = filter(b -> b["name"] == "left", newBndInfo)[1]
        left_blocks = sort([f["block"] for f in left_bnd["faces"]])
        @test left_blocks == [1, 2, 5]
        
        # Right boundary: All blocks on the right edge (i=max)
        # After j-split (propagated to block 2), both sub-blocks have i=max edge
        # Block 3: bottom part of old block 2 (jSeg=1, right edge at i=max)
        # Block 4: top part of old block 2 (jSeg=2, right edge at i=max)
        # Block 6: old block 4 (right edge at i=max)
        right_bnd = filter(b -> b["name"] == "right", newBndInfo)[1]
        right_blocks = sort([f["block"] for f in right_bnd["faces"]])
        @test right_blocks == [3, 4, 6]
        
        # Check interface information
        # Expected 7 interfaces total:
        # - 2 internal interfaces from the horizontal splits (1-2, 3-4)
        # - 5 cross-block interfaces (2-5, 4-6 horizontal; 1-3, 2-4, 5-6 vertical)
        @test length(newInterInfo) == 7
        
        # Helper to find interface
        function find_interface(interfaces, blockA, blockB)
            findfirst(i -> (i["blockA"] == blockA && i["blockB"] == blockB) || 
                          (i["blockA"] == blockB && i["blockB"] == blockA), interfaces)
        end
        
        # Verify all 7 expected interfaces exist:
        
        # 1. Internal horizontal interfaces from splits
        @test !isnothing(find_interface(newInterInfo, 1, 2))  # Internal in old block 1 (bottom-left quadrant)
        @test !isnothing(find_interface(newInterInfo, 3, 4))  # Internal in old block 1 (top-left quadrant, propagated)
        
        # 2. Vertical interfaces within old block 1 (connecting bottom and top quadrants)
        @test !isnothing(find_interface(newInterInfo, 1, 3))  # Left edge of old block 1
        @test !isnothing(find_interface(newInterInfo, 2, 4))  # Right edge of old block 1
        
        # 3. Horizontal interfaces between old blocks 1 and 2
        @test !isnothing(find_interface(newInterInfo, 2, 5))  # Bottom-right of old block 1 to bottom of old block 2
        @test !isnothing(find_interface(newInterInfo, 4, 6))  # Top-right of old block 1 to top of old block 2
        
        # 4. Vertical interface within old block 2
        @test !isnothing(find_interface(newInterInfo, 5, 6))  # Old block 2 bottom-top connection
    end
    
    @testset "No Propagation When No Interface" begin
        # Create 2 independent blocks (no shared interface)
        N = 10
        top = [collect(range(0.0, 1.0, length=N)) ones(N)]
        right = [ones(N) collect(range(1.0, 0.0, length=N))]
        bottom = [collect(range(1.0, 0.0, length=N)) zeros(N)]
        left = [zeros(N) collect(range(0.0, 1.0, length=N))]
        
        block1 = GridGeneration.TFI([top, right, bottom, left])
        block2 = GridGeneration.TFI([top, right, bottom, left])
        
        blocks = [block1, block2]
        bndInfo = []
        interInfo = []  # No interfaces between blocks
        
        splitRequests = Dict(
            1 => [Int[], [5]]  # Vertical split in block 1 only
        )
        
        newBlocks, _, _ = GridGeneration.SplitMultiBlock(blocks, splitRequests, bndInfo, interInfo)
        
        # Block 1 should split into 2, block 2 unchanged -> 3 total
        @test length(newBlocks) == 3
    end
    
    @testset "Mixed Splits in 2×2 Grid" begin
        # Create 2×2 grid and apply splits in multiple directions
        # Initial 2×2 layout: Lower-left=1, lower-right=2, upper-left=3, upper-right=4
        #
        # Split requests:
        # - Block 1 (lower-left): horizontal split (j-direction)
        # - Block 2 (lower-right): vertical split (i-direction)
        #
        # Propagation effects:
        # - Block 1 horizontal split → propagates to block 3 (vertically aligned)
        # - Block 2 vertical split → propagates to block 4 (horizontally aligned)
        # - Cross-propagation creates additional splits
        #
        # Final 9-block layout:
        #   ┌──────┬──────┬──────┐
        #   │   6  │      │      │  <- Top row
        #   ├──────┤   8  │  9   │
        #   │   5  │      │      │
        #   ├──────┼──────┼──────┤
        #   │   4  │      │      │  <- Bottom row  
        #   │      ├──────┤      │
        #   ├──────┤   2  │  7   │
        #   │   3  ├──────┤      │
        #   │      │   1  │      │
        #   └──────┴──────┴──────┘
        
        N = 10
        top = [collect(range(0.0, 1.0, length=N)) ones(N)]
        right = [ones(N) collect(range(1.0, 0.0, length=N))]
        bottom = [collect(range(1.0, 0.0, length=N)) zeros(N)]
        left = [zeros(N) collect(range(0.0, 1.0, length=N))]
        
        initialGrid = GridGeneration.TFI([top, right, bottom, left])
        
        # Split into 2×2
        splitLocs = [[5], [5]]
        bndInfo = []
        blocks, bndInfo, interInfo = GridGeneration.SplitBlock(initialGrid, splitLocs, bndInfo, [])
        
        @test length(blocks) == 4
        
        # Apply both horizontal and vertical splits
        # Block 1: horizontal split
        # Block 2: vertical split
        block1_nj = size(blocks[1], 3)
        block2_ni = size(blocks[2], 2)
        
        splitRequests = Dict(
            1 => [Int[], [div(block1_nj, 2)]],  # Horizontal in block 1
            2 => [[div(block2_ni, 2)], Int[]]   # Vertical in block 2
        )
        
        newBlocks, newBndInfo, newInterInfo = GridGeneration.SplitMultiBlock(
            blocks, splitRequests, bndInfo, interInfo
        )
        
        # With corrected propagation logic:
        # - Block 1 horizontal split [[], [j]] → propagates horizontally to block 2 (side-by-side)
        # - Block 2 vertical split [[i], []] → propagates vertically to block 4 (stacked)
        # This creates cascading effects resulting in 9 blocks
        @test length(newBlocks) == 9
        
        # Verify all blocks are valid
        for block in newBlocks
            @test ndims(block) == 3
            @test size(block, 1) == 2
        end
        
        # Verify interfaces exist (at least the internal ones from splits)
        @test length(newInterInfo) >= 4  # At minimum, internal interfaces from the splits
        
        # Verify no boundary info (empty bndInfo passed in)
        @test isempty(newBndInfo)
    end
    
    @testset "Complete 2×2→9 Block Split - Ground Truth Test" begin
        # This test represents the complete ground truth scenario:
        # 1. Start with a unit square grid
        # 2. Split into 2×2 grid (4 blocks)
        # 3. Split block 1 with both i and j splits
        # 4. Propagation creates 9 blocks total
        # 
        # Block layout after all splits:
        #   ┌─────┬─────┬─────┐
        #   │  7  │  8  │  9  │  <- Top row
        #   ├─────┼─────┼─────┤
        #   │  3  │  4  │  6  │  <- Middle row
        #   ├─────┼─────┼─────┤
        #   │  1  │  2  │  5  │  <- Bottom row
        #   └─────┴─────┴─────┘
        
        N = 50
        top = [collect(range(0.0, 1.0, length=N)) ones(N)]
        right = [ones(N) collect(range(0.0, 1.0, length=N))]
        bottom = [collect(range(0.0, 1.0, length=N)) zeros(N)]
        left = [zeros(N) collect(range(0.0, 1.0, length=N))]

        initialGrid = GridGeneration.TFI([top, right, bottom, left])

        # Split into 2×2 blocks
        splitLocsInitial = [[N÷2], [N÷2]]
        bndInfoInitial = [
            Dict("name" => "top", "faces" => [Dict("start" => [1, N, 1], "end" => [N, N, 1])]),
            Dict("name" => "right", "faces" => [Dict("start" => [N, 1, 1], "end" => [N, N, 1])]),
            Dict("name" => "bottom", "faces" => [Dict("start" => [1, 1, 1], "end" => [N, 1, 1])]),
            Dict("name" => "left", "faces" => [Dict("start" => [1, 1, 1], "end" => [1, N, 1])])
        ]

        blocks, bndInfo, interInfo = GridGeneration.SplitBlock(initialGrid, splitLocsInitial, bndInfoInitial, [])

        # Now split block 1 with both i and j splits
        block1_nj = size(blocks[1], 3)
        mid_j = div(block1_nj, 2)
        mid_i = div(size(blocks[1], 2), 2)

        splitRequests = Dict(
            1 => [[mid_i], [mid_j]]
        )

        newBlocks, newBndInfo, newInterInfo = GridGeneration.SplitMultiBlock(
            deepcopy(blocks), deepcopy(splitRequests), deepcopy(bndInfo), deepcopy(interInfo)
        )
        
        # Verify 9 blocks created
        @test length(newBlocks) == 9
        
        # Verify block sizes
        @test size(newBlocks[1]) == (2, 12, 12)
        @test size(newBlocks[2]) == (2, 14, 12)
        @test size(newBlocks[3]) == (2, 12, 14)
        @test size(newBlocks[4]) == (2, 14, 14)
        @test size(newBlocks[5]) == (2, 26, 12)
        @test size(newBlocks[6]) == (2, 26, 14)
        @test size(newBlocks[7]) == (2, 12, 26)
        @test size(newBlocks[8]) == (2, 14, 26)
        @test size(newBlocks[9]) == (2, 26, 26)
        
        # Ground Truth: Boundary Information
        # Left boundary should have blocks 1, 3, 7
        left_bnd = findfirst(b -> b["name"] == "left", newBndInfo)
        @test !isnothing(left_bnd)
        left_faces = newBndInfo[left_bnd]["faces"]
        @test length(left_faces) == 3
        left_blocks = sort([f["block"] for f in left_faces])
        @test left_blocks == [1, 3, 7]
        # Verify exact coordinates for left boundary
        for face in left_faces
            @test face["start"][1] == 1  # i = 1 (left edge)
            @test face["end"][1] == 1
            block_id = face["block"]
            block_nj = size(newBlocks[block_id], 3)
            @test face["end"][2] == block_nj  # j spans full height
        end
        
        # Bottom boundary should have blocks 1, 2, 5
        bottom_bnd = findfirst(b -> b["name"] == "bottom", newBndInfo)
        @test !isnothing(bottom_bnd)
        bottom_faces = newBndInfo[bottom_bnd]["faces"]
        @test length(bottom_faces) == 3
        bottom_blocks = sort([f["block"] for f in bottom_faces])
        @test bottom_blocks == [1, 2, 5]
        # Verify exact coordinates
        for face in bottom_faces
            @test face["start"][2] == 1  # j = 1 (bottom edge)
            @test face["end"][2] == 1
            block_id = face["block"]
            block_ni = size(newBlocks[block_id], 2)
            @test face["end"][1] == block_ni  # i spans full width
        end
        
        # Right boundary should have blocks 5, 6, 9
        right_bnd = findfirst(b -> b["name"] == "right", newBndInfo)
        @test !isnothing(right_bnd)
        right_faces = newBndInfo[right_bnd]["faces"]
        @test length(right_faces) == 3
        right_blocks = sort([f["block"] for f in right_faces])
        @test right_blocks == [5, 6, 9]
        # Verify exact coordinates
        for face in right_faces
            block_id = face["block"]
            block_ni = size(newBlocks[block_id], 2)
            @test face["start"][1] == block_ni  # i = max (right edge)
            @test face["end"][1] == block_ni
            block_nj = size(newBlocks[block_id], 3)
            @test face["end"][2] == block_nj  # j spans full height
        end
        
        # Top boundary should have blocks 7, 8, 9
        top_bnd = findfirst(b -> b["name"] == "top", newBndInfo)
        @test !isnothing(top_bnd)
        top_faces = newBndInfo[top_bnd]["faces"]
        @test length(top_faces) == 3
        top_blocks = sort([f["block"] for f in top_faces])
        @test top_blocks == [7, 8, 9]
        # Verify exact coordinates
        for face in top_faces
            block_id = face["block"]
            block_nj = size(newBlocks[block_id], 3)
            @test face["start"][2] == block_nj  # j = max (top edge)
            @test face["end"][2] == block_nj
            block_ni = size(newBlocks[block_id], 2)
            @test face["end"][1] == block_ni  # i spans full width
        end
        
        # Ground Truth: Interface Information
        # Expected 12 interfaces total
        @test length(newInterInfo) == 12
        
        # Helper to find interface
        function find_interface(interfaces, blockA, blockB)
            findfirst(i -> (i["blockA"] == blockA && i["blockB"] == blockB) || 
                          (i["blockA"] == blockB && i["blockB"] == blockA), interfaces)
        end
        
        # Verify all expected interfaces exist
        expected_interfaces = [
            (1, 2), (1, 3), (2, 4), (3, 4), (5, 6),
            (7, 8), (2, 5), (4, 6), (3, 7), (4, 8),
            (6, 9), (8, 9)
        ]
        
        for (blockA, blockB) in expected_interfaces
            idx = find_interface(newInterInfo, blockA, blockB)
            @test !isnothing(idx)  # "Interface $blockA-$blockB should exist"
        end
        
        # Verify specific interface coordinates (spot check a few key ones)
        # Interface 1-2 (horizontal, within bottom-left quadrant)
        idx_1_2 = find_interface(newInterInfo, 1, 2)
        inter_1_2 = newInterInfo[idx_1_2]
        if inter_1_2["blockA"] == 1
            @test inter_1_2["start_blkA"] == [12, 1, 1]
            @test inter_1_2["end_blkA"] == [12, 12, 1]
            @test inter_1_2["start_blkB"] == [1, 1, 1]
            @test inter_1_2["end_blkB"] == [1, 12, 1]
        else
            @test inter_1_2["start_blkB"] == [12, 1, 1]
            @test inter_1_2["end_blkB"] == [12, 12, 1]
            @test inter_1_2["start_blkA"] == [1, 1, 1]
            @test inter_1_2["end_blkA"] == [1, 12, 1]
        end
        
        # Interface 2-5 (horizontal, bottom-left to bottom-right)
        idx_2_5 = find_interface(newInterInfo, 2, 5)
        inter_2_5 = newInterInfo[idx_2_5]
        if inter_2_5["blockA"] == 2
            @test inter_2_5["start_blkA"] == [14, 1, 1]
            @test inter_2_5["end_blkA"] == [14, 12, 1]
            @test inter_2_5["start_blkB"] == [1, 1, 1]
            @test inter_2_5["end_blkB"] == [1, 12, 1]
        else
            @test inter_2_5["start_blkB"] == [14, 1, 1]
            @test inter_2_5["end_blkB"] == [14, 12, 1]
            @test inter_2_5["start_blkA"] == [1, 1, 1]
            @test inter_2_5["end_blkA"] == [1, 12, 1]
        end
    end
end
