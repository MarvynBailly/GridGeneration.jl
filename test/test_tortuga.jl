using Test
using GridGeneration

@testset "Tortuga grid I/O" begin
    G = GridGeneration
    N1, N2 = 6, 5
    unit_block(x0, x1) = TFI([hcat(range(x0, x1, length=N1), ones(N1)), hcat(fill(x1, N2), range(0, 1, length=N2)),
                              hcat(range(x0, x1, length=N1), zeros(N1)), hcat(fill(x0, N2), range(0, 1, length=N2))])
    blocks = [unit_block(0.0, 1.0), unit_block(1.0, 2.0)]

    # 1-based 2D boundary/interface info, as produced by splitting and edge solving
    bndInfo = Any[Dict{String,Any}("name" => "BCWall", "faces" => Any[
        Dict{String,Any}("block" => 1, "start" => [1, 1], "end" => [N1, 1]),
        Dict{String,Any}("block" => 2, "start" => [1, 1], "end" => [N1, 1])])]
    interInfo = Any[Dict{String,Any}("blockA" => 1, "blockB" => 2,
        "start_blkA" => [N1, 1], "end_blkA" => [N1, N2], "start_blkB" => [1, 1], "end_blkB" => [1, N2],
        "offset" => [0.0, 0.0, 0.0], "angle" => 0.0)]

    k_layers = 2
    mesh3D, bnd3D, int3D = G.convert_2D_to_3D(blocks, bndInfo, interInfo, 0.1, k_layers)
    file = tempname() * ".grid"
    try
        G.write_turtle_grid(mesh3D, int3D, bnd3D, file)
        Xno2D, _, _, interfaces, boundaries = G.ImportTurtleGrid(file)

        @test length(Xno2D) == 2
        @test all(size.(Xno2D) .== size.(blocks))
        @test all(Xno2D[b] ≈ blocks[b] for b in 1:2)

        # Boundary faces come back 0-based with the k-extent added
        @test length(boundaries) == 1
        @test strip(boundaries[1]["name"]) == "BCWall"
        faces = boundaries[1]["faceInfo"]
        @test [f["block"] for f in faces] == [0, 1]
        @test faces[1]["start"] == [0, 0, 0]
        @test faces[1]["end"] == [N1 - 1, 0, k_layers - 1]

        # The block-to-block interface survives (plus periodic k self-interfaces from the extrusion)
        itf = only(filter(i -> i["blockA"] != i["blockB"], interfaces))
        @test (itf["blockA"], itf["blockB"]) == (0, 1)
        @test itf["start_blkA"] == [N1 - 1, 0, 0]
        @test itf["end_blkB"] == [0, N2 - 1, k_layers - 1]
    finally
        rm(file; force=true)
    end
end
