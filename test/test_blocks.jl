using Test
using GridGeneration

@testset "Block Operations" begin
    @testset "Block Splitting" begin
        # Build a deterministic square boundary and TFI grid
        N = 10
        top    = [collect(range(0.0, 1.0, length=N))  .+ 0.0  ones(N)]  # (N×2) but keep consistent layout
        top    = [collect(range(0.0, 1.0, length=N)) ones(N)]
        right  = [ones(N) collect(range(1.0, 0.0, length=N))]
        bottom = [collect(range(1.0, 0.0, length=N)) zeros(N)]
        left   = [zeros(N) collect(range(0.0, 1.0, length=N))]

        initialGrid = TFI([top, right, bottom, left])

        # Split in the middle -> 2×2 blocks
        splitLocations = [[5], [5]]

        # Define parent boundary info (global indices: 1..N)
        bndInfo = [
            Dict("name" => "top",    "faces" => [Dict("start" => [1, N, 1], "end" => [N, N, 1])]),
            Dict("name" => "right",  "faces" => [Dict("start" => [N, 1, 1], "end" => [N, N, 1])]),
            Dict("name" => "bottom", "faces" => [Dict("start" => [1, 1, 1], "end" => [N, 1, 1])]),
            Dict("name" => "left",   "faces" => [Dict("start" => [1, 1, 1], "end" => [1, N, 1])])
        ]

        interInfo = []

        blocks, newBndInfo, newInterInfo = GridGeneration.SplitBlock(
            initialGrid, splitLocations, bndInfo, interInfo
        )

        # Expect 4 blocks
        @test length(blocks) == 4

        # Expected local sizes due to inclusive indexing: segments are [1:5] (len 5) and [5:10] (len 6)
        expected_sizes = [(5,5), (6,5), (5,6), (6,6)]  # (ni,nj) for blocks 1..4 in enumeration order

        for (idx, (ni, nj)) in enumerate(expected_sizes)
            blk = blocks[idx]
            @test ndims(blk) == 3
            @test size(blk, 1) == 2
            @test size(blk, 2) == ni
            @test size(blk, 3) == nj
        end

        # Expect 4 internal interfaces: 2 vertical (1<->2, 3<->4) and 2 horizontal (1<->3, 2<->4)
        @test length(newInterInfo) == 4

        # Helper: find an interface that matches fields
        find_interface(interfaces, kv) = findfirst(i -> all(k -> i[k] == kv[k], keys(kv)), interfaces)

        # Vertical bottom row: block 1 <-> 2
        expect_v1 = Dict(
            "blockA" => 1, "blockB" => 2,
            "start_blkA" => [5,1,1], "end_blkA" => [5,5,1],
            "start_blkB" => [1,1,1], "end_blkB" => [1,5,1]
        )
        @test !isnothing(find_interface(newInterInfo, expect_v1))

        # Vertical top row: block 3 <-> 4
        expect_v2 = Dict(
            "blockA" => 3, "blockB" => 4,
            "start_blkA" => [5,1,1], "end_blkA" => [5,6,1],
            "start_blkB" => [1,1,1], "end_blkB" => [1,6,1]
        )
        @test !isnothing(find_interface(newInterInfo, expect_v2))

        # Horizontal left column: block 1 <-> 3
        expect_h1 = Dict(
            "blockA" => 1, "blockB" => 3,
            "start_blkA" => [1,5,1], "end_blkA" => [5,5,1],
            "start_blkB" => [1,1,1], "end_blkB" => [5,1,1]
        )
        @test !isnothing(find_interface(newInterInfo, expect_h1))

        # Horizontal right column: block 2 <-> 4
        expect_h2 = Dict(
            "blockA" => 2, "blockB" => 4,
            "start_blkA" => [1,5,1], "end_blkA" => [6,5,1],
            "start_blkB" => [1,1,1], "end_blkB" => [6,1,1]
        )
        @test !isnothing(find_interface(newInterInfo, expect_h2))

        # Check grouped boundary info: each outer boundary should be present and split into two faces
        names = [d["name"] for d in newBndInfo]
        @test all(x -> x in names, ["top", "right", "bottom", "left"]) 

        # bottom should have two faces (blocks 1 and 2)
        bottom_group = filter(d->d["name"]=="bottom", newBndInfo)[1]
        bottom_blocks = sort([f["block"] for f in bottom_group["faces"]])
        @test bottom_blocks == [1,2]
    end
end
