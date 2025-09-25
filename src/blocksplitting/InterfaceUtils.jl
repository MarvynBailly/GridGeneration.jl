function UpdateInterInfo(interInfo, blocks; verbose=false)
    """
    Update interface endpoints to match current block sizes.

    For each interface side (A and B), any index != 1 is set to the block's
    maximum in that dimension: i->Nx, j->Ny. k is forced to 1.
    Assumes blocks are arrays with Nx = size(block, 2), Ny = size(block, 3).
    """
    # small helper: ensure we can mutate the vectors
    ensure_vec_int(x) = x isa Vector{Int} ? x : Vector{Int}(x)

    for inter in interInfo
        bA = inter["blockA"]; bB = inter["blockB"]
        blockA = blocks[bA]
        blockB = blocks[bB]

        NxA = size(blockA, 2); NyA = size(blockA, 3)
        NxB = size(blockB, 2); NyB = size(blockB, 3)

        sA = ensure_vec_int(inter["start_blkA"]); eA = ensure_vec_int(inter["end_blkA"])
        sB = ensure_vec_int(inter["start_blkB"]); eB = ensure_vec_int(inter["end_blkB"])

        # clamp to {1, Nx}/{1, Ny}; keep k == 1
        sA[1] = (sA[1] == 1) ? 1 : NxA;  sA[2] = (sA[2] == 1) ? 1 : NyA;  sA[3] = 1
        eA[1] = (eA[1] == 1) ? 1 : NxA;  eA[2] = (eA[2] == 1) ? 1 : NyA;  eA[3] = 1
        sB[1] = (sB[1] == 1) ? 1 : NxB;  sB[2] = (sB[2] == 1) ? 1 : NyB;  sB[3] = 1
        eB[1] = (eB[1] == 1) ? 1 : NxB;  eB[2] = (eB[2] == 1) ? 1 : NyB;  eB[3] = 1

        inter["start_blkA"] = sA; inter["end_blkA"] = eA
        inter["start_blkB"] = sB; inter["end_blkB"] = eB

        if verbose
            println("Updated inter A=$(bA) B=$(bB): ",
                    "A[start=", sA, ", end=", eA, "]  ",
                    "B[start=", sB, ", end=", eB, "]")
        end
    end
    return interInfo
end
