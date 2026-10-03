using Test
using GridGeneration

# Interfaces whose side A is the upper/right block, or whose two sides run in opposite
# directions, must still be remapped onto the touching edges after splitting.
@testset "Multi-block splitting: interface orientation" begin
    rect(x0, x1, y0, y1, ni, nj) = TFI([hcat(range(x0, x1, length=ni), fill(y1, ni)), hcat(fill(x1, nj), range(y0, y1, length=nj)),
                                        hcat(range(x0, x1, length=ni), fill(y0, ni)), hcat(fill(x0, nj), range(y0, y1, length=nj))])

    function interfaces_match(blocks, interInfo)
        all(interInfo) do it
            sa, ea, sb, eb = it["start_blkA"], it["end_blkA"], it["start_blkB"], it["end_blkB"]
            ra = (sa[1] <= ea[1] ? (sa[1]:ea[1]) : (sa[1]:-1:ea[1]), sa[2] <= ea[2] ? (sa[2]:ea[2]) : (sa[2]:-1:ea[2]))
            rb = (sb[1] <= eb[1] ? (sb[1]:eb[1]) : (sb[1]:-1:eb[1]), sb[2] <= eb[2] ? (sb[2]:eb[2]) : (sb[2]:-1:eb[2]))
            pa = blocks[it["blockA"]][:, ra...]
            pb = blocks[it["blockB"]][:, rb...]
            size(pa) == size(pb) && pa ≈ pb
        end
    end

    upper = rect(0, 2, 1, 2, 9, 5)   # block 1
    lower = rect(0, 2, 0, 1, 9, 5)   # block 2

    @testset "side A on the bottom edge of the upper block" begin
        itf = Any[Dict{String,Any}("blockA" => 1, "blockB" => 2,
            "start_blkA" => [1, 1, 1], "end_blkA" => [9, 1, 1],      # bottom of upper
            "start_blkB" => [1, 5, 1], "end_blkB" => [9, 5, 1],      # top of lower
            "offset" => [0.0, 0.0, 0.0], "angle" => 0.0)]
        @test interfaces_match([upper, lower], itf)

        blocks, _, interInfo = SplitMultiBlock([upper, lower], [(1, [[5], Int[]])], Any[], itf)
        @test length(blocks) == 4                      # split propagated into the lower block
        @test count(it -> it["blockA"] in (1, 2) && it["blockB"] in (3, 4), interInfo) == 2
        @test interfaces_match(blocks, interInfo)
    end

    @testset "sides running in opposite directions" begin
        itf = Any[Dict{String,Any}("blockA" => 1, "blockB" => 2,
            "start_blkA" => [9, 1, 1], "end_blkA" => [1, 1, 1],      # bottom of upper, right to left
            "start_blkB" => [9, 5, 1], "end_blkB" => [1, 5, 1],      # top of lower, right to left
            "offset" => [0.0, 0.0, 0.0], "angle" => 0.0)]
        @test interfaces_match([upper, lower], itf)

        blocks, _, interInfo = SplitMultiBlock([upper, lower], [(2, [[3], Int[]])], Any[], itf)
        @test length(blocks) == 4
        @test interfaces_match(blocks, interInfo)
    end

    @testset "side-by-side with side A on the right" begin
        left = rect(0, 1, 0, 2, 5, 9)
        right = rect(1, 2, 0, 2, 5, 9)
        itf = Any[Dict{String,Any}("blockA" => 2, "blockB" => 1,
            "start_blkA" => [1, 1, 1], "end_blkA" => [1, 9, 1],      # left edge of right block
            "start_blkB" => [5, 1, 1], "end_blkB" => [5, 9, 1],      # right edge of left block
            "offset" => [0.0, 0.0, 0.0], "angle" => 0.0)]
        blocks, _, interInfo = SplitMultiBlock([left, right], [(1, [Int[], [4]])], Any[], itf)
        @test length(blocks) == 4
        @test interfaces_match(blocks, interInfo)
    end
end
