"""
Split single block grid into multiple blocks based on specified split locations.

Inputs:
- block: 3D array representing grid (dimensions: (ni, nj, nk))
- splitLocations: Array of two arrays [x-splits, y-splits] as grid indices
- bndInfo: Boundary condition dictionary
- interInfo: Interface connectivity dictionary
"""
function SplitBlock(block, splitLocations, bndInfo, interInfo)
    blocks = []
    horzSplits = [1, splitLocations[1]..., size(block, 2)]
    vertSplits = [1, splitLocations[2]..., size(block, 3)]

    blockId = 1
    blockBoundaries = []
    interfaces = []
    for j in 1:length(vertSplits)-1
        for i in 1:length(horzSplits)-1
            ni = horzSplits[i+1] - horzSplits[i] + 1
            nj = vertSplits[j+1] - vertSplits[j] + 1
            nk = 1
            
            subblock = block[:, horzSplits[i]:horzSplits[i+1], vertSplits[j]:vertSplits[j+1]]
            push!(blocks, subblock) 

            # Get boundary info
            blockInfo = Dict(
                "block" => blockId,
                "start" => (horzSplits[i], vertSplits[j]),
                "end" => (horzSplits[i+1], vertSplits[j+1]),
            )

            boundaries = GridGeneration.GetTouchingBoundaries(blockInfo, bndInfo)
            append!(blockBoundaries, boundaries)
            
            # Get interface info - look forward to next block if not at end
            if i < length(horzSplits) - 1
                blockBId = blockId + 1
                push!(interfaces, Dict(
                    "blockA" => blockId, "start_blkA" => [ni,1,1], "end_blkA" => [ni,nj,nk],
                    "blockB" => blockBId, "start_blkB" => [1,1,1], "end_blkB" => [1,nj,nk],
                    "offset" => [0.0, 0.0, 0.0], "angle" => 0.0))
            end

            # if not at top, look up
            if j < length(vertSplits) - 1
                blockBId = blockId + length(horzSplits) - 1
                push!(interfaces, Dict(
                    "blockA" => blockId, "start_blkA" => [1,nj,1], "end_blkA" => [ni,nj,nk],
                    "blockB" => blockBId, "start_blkB" => [1,1,1], "end_blkB" => [ni,1,nk],
                    "offset" => [0.0, 0.0, 0.0], "angle" => 0.0))
            end

            blockId += 1
        end
    end
    updatedBndInfo = GridGeneration.GroupBoundariesByName(blockBoundaries)
    updatedInterInfo = interfaces

    return blocks, updatedBndInfo, updatedInterInfo
end