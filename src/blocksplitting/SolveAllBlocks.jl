"""
    SolveAllBlocks(metric, blocks, bndInfo, interInfo; solver=:analytic)
        -> (blocks, bndInfo, interInfo)

Redistribute the edge points of every block according to `metric(x, y) -> (M11, M22)` and
rebuild each block with [`TFI`](@ref). The optimal number of points in each direction is
computed per edge pair and shared across neighbouring blocks so interfaces stay conforming.
`solver` is `:analytic` or `:numeric`. Boundary and interface indices are updated to the
new block sizes.
"""
function SolveAllBlocks(metric, blocks, bndInfo, interInfo; solver =:analytic)
    # blockDirOptN = similar(blocks)
    blockDirOptN = Vector{Vector{Int}}(undef, length(blocks))
    
    for i in 1:length(blockDirOptN)
        blockDirOptN[i] = [-1,-1]
    end

    # Compute max number of points for each block
    for (blockId, block) in enumerate(blocks)
        for dir in 1:2  # 1=horizontal, 2=vertical
            blockNeighbors = GetNeighbors(blockId, interInfo, dir; include_start=true)

            if dir == 1  # horizontal
                left   = block[:, 1, :]
                right  = block[:, end, :]
                optN = GridGeneration.GetOptNEdgePair(left, right, metric; solver=solver)
            else  # vertical
                bottom   = block[:, :, 1]
                top  = block[:, :, end]
                optN = GridGeneration.GetOptNEdgePair(bottom, top, metric; solver=solver)
            end

            # Update blockDirOptN with max number of points
            for computeBlocks in blockNeighbors
                if blockDirOptN[computeBlocks][dir] < optN
                    blockDirOptN[computeBlocks][dir] = optN
                end
            end
        end
    end

    # Solve all blocks using optimal number
    computedBlocks = similar(blocks)
    for (blockId, block) in enumerate(blocks)
        optNs = blockDirOptN[blockId]
        computedBlock, bndInfo, interInfo = GridGeneration.SolveBlockFixedN(block, bndInfo, interInfo, metric, optNs; solver=solver)
        computedBlocks[blockId] = computedBlock
    end

    # Update boundary and interface information
    GridGeneration.UpdateBndInfo!(bndInfo, computedBlocks; verbose=false)
    updatedInterInfo = GridGeneration.UpdateInterInfo(interInfo, computedBlocks; verbose=false)

    return computedBlocks, bndInfo, updatedInterInfo
end
