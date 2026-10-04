using Test
using GridGeneration

@testset "Turtle domain loading" begin
    G = GridGeneration
    rect(x0, x1, y0, y1, ni, nj) = TFI([hcat(range(x0, x1, length=ni), fill(y1, ni)), hcat(fill(x1, nj), range(y0, y1, length=nj)),
                                        hcat(range(x0, x1, length=ni), fill(y0, ni)), hcat(fill(x0, nj), range(y0, y1, length=nj))])

    # Interface points on both sides must coincide exactly
    function interfaces_match(blocks, interInfo)
        all(interInfo) do it
            sa, ea, sb, eb = it["start_blkA"], it["end_blkA"], it["start_blkB"], it["end_blkB"]
            pa = blocks[it["blockA"]][:, sa[1]:ea[1], sa[2]:ea[2]]
            pb = blocks[it["blockB"]][:, sb[1]:eb[1], sb[2]:eb[2]]
            size(pa) == size(pb) && pa ≈ pb
        end
    end

    function roundtrip_load(blocks, bndInfo, interInfo; kwargs...)
        mesh3D, b3, i3 = G.convert_2D_to_3D(blocks, bndInfo, interInfo, 0.1, 2)
        file = tempname() * ".grid"
        try
            redirect_stdout(devnull) do
                G.write_turtle_grid(mesh3D, i3, b3, file)
                load_turtle_grid(file; kwargs...)
            end
        finally
            rm(file; force=true)
        end
    end

    # Block "big" (6×9) with block "small" (6×5) attached to the lower half of its right edge
    big = rect(0, 1, 0, 1, 6, 9)
    small = rect(1, 2, 0, 0.5, 6, 5)
    wall = Any[Dict{String,Any}("name" => "BCWall", "faces" => Any[
        Dict{String,Any}("block" => 1, "start" => [1, 1], "end" => [6, 1]),
        Dict{String,Any}("block" => 2, "start" => [1, 1], "end" => [6, 1])])]

    @testset "conversion to package conventions" begin
        even = rect(1, 2, 0, 1, 6, 9)
        itf = Any[Dict{String,Any}("blockA" => 1, "blockB" => 2, "start_blkA" => [6, 1], "end_blkA" => [6, 9],
                                   "start_blkB" => [1, 1], "end_blkB" => [1, 9], "offset" => [0.0, 0.0, 0.0], "angle" => 0.0)]
        blocks, bndInfo, interInfo, centers = roundtrip_load([big, even], wall, itf)

        @test length(blocks) == 2
        @test blocks[1] ≈ big
        @test length(centers) == 2
        # 1-based, faces under "faces", periodic k self-interfaces dropped
        @test haskey(bndInfo[1], "faces") && !haskey(bndInfo[1], "faceInfo")
        @test [f["block"] for f in bndInfo[1]["faces"]] == [1, 2]
        @test length(interInfo) == 1
        @test (interInfo[1]["blockA"], interInfo[1]["blockB"]) == (1, 2)
        @test interInfo[1]["start_blkA"][1:2] == [6, 1]
        @test interfaces_match(blocks, interInfo)
    end

    @testset "uneven interface on side $side" for side in (:A, :B)
        itf = side == :A ?
            Dict{String,Any}("blockA" => 1, "blockB" => 2, "start_blkA" => [6, 1], "end_blkA" => [6, 5],
                             "start_blkB" => [1, 1], "end_blkB" => [1, 5]) :
            Dict{String,Any}("blockA" => 2, "blockB" => 1, "start_blkA" => [1, 1], "end_blkA" => [1, 5],
                             "start_blkB" => [6, 1], "end_blkB" => [6, 5])
        merge!(itf, Dict("offset" => [0.0, 0.0, 0.0], "angle" => 0.0))

        blocks, bndInfo, interInfo, _ = roundtrip_load([big, small], wall, Any[itf])

        # big is split at j = 5 (the shared node line) into two 6×5 blocks
        @test size.(blocks) == [(2, 6, 5), (2, 6, 5), (2, 6, 5)]
        @test length(interInfo) == 2          # internal split interface + the original one
        @test interfaces_match(blocks, interInfo)
        # the bottom wall stays on the lower sub-block of big and on small
        @test sort([f["block"] for f in bndInfo[1]["faces"]]) == [1, 3]

        # without fixing, the uneven interface is left as is
        blocks2, _, interInfo2, _ = roundtrip_load([big, small], wall, Any[itf]; fix_uneven_interfaces=false)
        @test length(blocks2) == 2
        @test length(interInfo2) == 1
    end

    @testset "interface crossing a split is rejected" begin
        splitRequests = Dict(1 => [Int[], [5]])
        mapping = Dict(1 => [1, 2], 2 => [3])
        newBlocks = [big[:, :, 1:5], big[:, :, 5:9], small]
        itf = Dict{String,Any}("blockA" => 1, "blockB" => 2, "start_blkA" => [6, 1, 1], "end_blkA" => [6, 9, 1],
                               "start_blkB" => [1, 1, 1], "end_blkB" => [1, 5, 1])
        @test_throws ErrorException G.update_external_interface(itf, splitRequests, mapping, newBlocks)
    end
end
