"""
Solve all blocks with optimal point distribution based on metric field.
"""
function SolveAllBlocks(metric, blocks, bndInfo, interInfo; solver =:analytic)
    blockDirOptN = similar(blocks)
    
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
