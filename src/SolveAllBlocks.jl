"""
Get all (transitively) neighboring block IDs reachable from `blockId`
through interfaces aligned with `dir` (1=horizontal, 2=vertical).
- `interInfo` entries must have keys: "blockA","blockB","start_blkA","end_blkA".
- Set `include_start=true` to include `blockId` in the result.
"""
function GetNeighbors(blockId::Int, interInfo, dir::Int; include_start::Bool=false)
    seen = Set{Int}()

    function visit(bid::Int)
        # already visited this block
        bid ∈ seen && return
        push!(seen, bid)

        for info in interInfo
            if info["blockA"] == bid || info["blockB"] == bid
                # Orientation is defined by the A-side span (shared for both sides)
                startA = info["start_blkA"]
                endA   = info["end_blkA"]
                is_vertical   = (startA[1] != endA[1])
                is_horizontal = (startA[2] != endA[2])

                if (dir == 2 && is_vertical) || (dir == 1 && is_horizontal)
                    other = (info["blockA"] == bid) ? info["blockB"] : info["blockA"]
                    visit(other)
                end
            end
        end
    end

    visit(blockId)

    return include_start ? collect(seen) : [b for b in seen if b != blockId]
end


function SolveAllBlocks(metric, blocks, bndInfo, interInfo; solver = "analytic")
    blockDirOptN = similar(blocks)

    for i in 1:length(blockDirOptN)
        blockDirOptN[i] = [-1,-1]
    end

    # compute max number of points for each block
    for (blockId, block) in enumerate(blocks)
        for dir in 1:2  # 1 for horizontal, 2 for vertical
            blockNeighbors = GetNeighbors(blockId, interInfo, dir; include_start=true)
            # println("blockId: ", blockId, " dir: ", dir, " neighbors: ", blockNeighbors)

            #############
            # method 1
            #############
            if dir == 1  # horizontal
                left   = block[:, 1, :]
                right  = block[:, end, :]
                optN = GridGeneration.GetOptNEdgePair(left, right, metric; solver=solver)
            else  # vertical
                bottom   = block[:, :, 1]
                top  = block[:, :, end]
                optN = GridGeneration.GetOptNEdgePair(bottom, top, metric; solver=solver)
            end

            # Update the blockDirOptN with the max number of points
            for computeBlocks in blockNeighbors
                if blockDirOptN[computeBlocks][dir] < optN
                    blockDirOptN[computeBlocks][dir] = optN
                end
            end
        end
    end

    # solve all blocks using the optimal number
    computedBlocks = similar(blocks)
    for (blockId, block) in enumerate(blocks)
        optNs = blockDirOptN[blockId]

        computedBlock, bndInfo, interInfo = GridGeneration.SolveBlockFixedN(block, bndInfo, interInfo, metric, optNs; solver=solver)

        # println("Solving block $blockId with ni=$(optNs[1]) and nj=$(optNs[2])")
        # println("Size of block $blockId: ", size(computedBlock))
        computedBlocks[blockId] = computedBlock
    end

    # update the boundary information and interface information 
    GridGeneration.UpdateBndInfo!(bndInfo, computedBlocks; verbose=false)
    updatedInterInfo = GridGeneration.UpdateInterInfo(interInfo, computedBlocks; verbose=false)

    return computedBlocks, bndInfo, updatedInterInfo
end
