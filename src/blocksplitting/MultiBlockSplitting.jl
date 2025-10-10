"""
Multi-block splitting with automatic propagation across interfaces.

When a block is split, the split must propagate to connected blocks:
- Horizontal splits (in j-direction) propagate through vertical interfaces
- Vertical splits (in i-direction) propagate through horizontal interfaces
"""

"""
Determine which direction an interface connects (1=horizontal edge, 2=vertical edge).
Returns the direction that varies along the interface.
"""
function GetInterfaceOrientation(interface::Dict)
    startA = interface["start_blkA"]
    endA = interface["end_blkA"]
    
    is_vertical = (startA[1] != endA[1])    # i varies -> vertical interface
    is_horizontal = (startA[2] != endA[2])  # j varies -> horizontal interface
    
    if is_vertical
        return 2  # vertical interface connects blocks in i-direction
    elseif is_horizontal
        return 1  # horizontal interface connects blocks in j-direction
    else
        error("Invalid interface: start and end are identical")
    end
end

"""
Map a split location from one block to a connected block through an interface.

Args:
- sourceSplit: split index in source block (in i or j depending on splitDir)
- sourceBlockId, targetBlockId: block IDs
- interface: interface dict connecting the blocks
- sourceBlock, targetBlock: the actual grid blocks
- splitDir: 1 for i-direction (vertical split line), 2 for j-direction (horizontal split line)

Returns the mapped split index in target block, or nothing if split doesn't propagate.

Propagation rules:
- A split in j-direction (horizontal line) propagates through horizontal interfaces (those with varying i)
- A split in i-direction (vertical line) propagates through vertical interfaces (those with varying j)
"""
function MapSplitAcrossInterface(sourceSplit::Int, sourceBlockId::Int, targetBlockId::Int, 
                                  interface::Dict, sourceBlock, targetBlock, splitDir::Int)
    # Determine which block is A and which is B
    if interface["blockA"] == sourceBlockId && interface["blockB"] == targetBlockId
        sourceStart = interface["start_blkA"]
        sourceEnd = interface["end_blkA"]
        targetStart = interface["start_blkB"]
        targetEnd = interface["end_blkB"]
    elseif interface["blockB"] == sourceBlockId && interface["blockA"] == targetBlockId
        sourceStart = interface["start_blkB"]
        sourceEnd = interface["end_blkB"]
        targetStart = interface["start_blkA"]
        targetEnd = interface["end_blkA"]
    else
        return nothing  # Interface doesn't connect these blocks
    end
    
    # Determine interface orientation
    is_vertical_interface = (sourceStart[1] != sourceEnd[1])    # i varies
    is_horizontal_interface = (sourceStart[2] != sourceEnd[2])  # j varies
    
    # CORRECTED LOGIC:
    # A j-split (horizontal line in block, splitDir=2) propagates through HORIZONTAL interfaces (blocks side-by-side)
    # An i-split (vertical line in block, splitDir=1) propagates through VERTICAL interfaces (blocks stacked vertically)
    #
    # This ensures one-to-one grid point correspondence across interfaces:
    # - Horizontal split [[], [j]] creates horizontal line → must align with side-by-side neighbor (horizontal interface)
    # - Vertical split [[i], []] creates vertical line → must align with vertically stacked neighbor (vertical interface)
    
    if splitDir == 2 && is_horizontal_interface
        # Horizontal split (j-direction) propagates through horizontal interface (side-by-side blocks)
        sourceNj = size(sourceBlock, 3)
        targetNj = size(targetBlock, 3)
        
        # Check if interface connects at left/right edges (constant i)
        if sourceStart[1] == sourceEnd[1]  # Interface at constant i
            sourceInterfaceI = sourceStart[1]
            targetInterfaceI = targetStart[1]
            
            if sourceInterfaceI == size(sourceBlock, 2) && targetInterfaceI == 1
                # Right of source connects to left of target
                return sourceSplit
            elseif sourceInterfaceI == 1 && targetInterfaceI == size(targetBlock, 2)
                # Left of source connects to right of target
                return sourceSplit
            end
        end
        
    elseif splitDir == 1 && is_vertical_interface
        # Vertical split (i-direction) propagates through vertical interface (stacked blocks)
        sourceNi = size(sourceBlock, 2)
        targetNi = size(targetBlock, 2)
        
        # Check if interface connects at top/bottom edges (constant j)
        if sourceStart[2] == sourceEnd[2]  # Interface at constant j
            sourceInterfaceJ = sourceStart[2]
            targetInterfaceJ = targetStart[2]
            
            if sourceInterfaceJ == size(sourceBlock, 3) && targetInterfaceJ == 1
                # Top of source connects to bottom of target
                return sourceSplit
            elseif sourceInterfaceJ == 1 && targetInterfaceJ == size(targetBlock, 3)
                # Bottom of source connects to top of target
                return sourceSplit
            end
        end
    end
    
    return nothing
end

"""
Propagate splits across all connected blocks recursively.

Args:
- blocks: Vector of grid blocks
- splitRequests: Dict mapping blockId => [[i_splits], [j_splits]]
- bndInfo: Boundary information
- interInfo: Interface information

Returns a Dict mapping blockId => [[i_splits], [j_splits]] with all propagated splits.
"""
function PropagateSplits(blocks::Vector, splitRequests::Dict{Int, Vector{Vector{Int}}}, 
                        bndInfo, interInfo)
    numBlocks = length(blocks)
    
    # Initialize split requirements for each block
    # splitReqs[blockId] = [[i_splits], [j_splits]]
    splitReqs = Dict{Int, Vector{Vector{Int}}}()
    for i in 1:numBlocks
        splitReqs[i] = [Int[], Int[]]  # [i-splits, j-splits]
    end
    
    # Add initial user-requested splits
    for (blockId, splits) in splitRequests
        if blockId > 0 && blockId <= numBlocks
            splitReqs[blockId][1] = sort(unique(vcat(splitReqs[blockId][1], splits[1])))
            splitReqs[blockId][2] = sort(unique(vcat(splitReqs[blockId][2], splits[2])))
        end
    end
    
    # Propagate splits iteratively until no new splits are added
    changed = true
    maxIter = 100  # prevent infinite loops
    iter = 0
    
    while changed && iter < maxIter
        changed = false
        iter += 1
        
        # For each block with splits, propagate to connected blocks
        for blockId in 1:numBlocks
            iSplits = splitReqs[blockId][1]
            jSplits = splitReqs[blockId][2]
            
            # Propagate vertical splits (in i-direction)
            for split in iSplits
                for interface in interInfo
                    if interface["blockA"] == blockId || interface["blockB"] == blockId
                        targetId = (interface["blockA"] == blockId) ? interface["blockB"] : interface["blockA"]
                        
                        mappedSplit = MapSplitAcrossInterface(split, blockId, targetId, 
                                                              interface, blocks[blockId], 
                                                              blocks[targetId], 1)
                        
                        if !isnothing(mappedSplit) && !(mappedSplit in splitReqs[targetId][1])
                            push!(splitReqs[targetId][1], mappedSplit)
                            changed = true
                        end
                    end
                end
            end
            
            # Propagate horizontal splits (in j-direction)
            for split in jSplits
                for interface in interInfo
                    if interface["blockA"] == blockId || interface["blockB"] == blockId
                        targetId = (interface["blockA"] == blockId) ? interface["blockB"] : interface["blockA"]
                        
                        mappedSplit = MapSplitAcrossInterface(split, blockId, targetId, 
                                                              interface, blocks[blockId], 
                                                              blocks[targetId], 2)
                        
                        if !isnothing(mappedSplit) && !(mappedSplit in splitReqs[targetId][2])
                            push!(splitReqs[targetId][2], mappedSplit)
                            changed = true
                        end
                    end
                end
            end
        end
        
        # Sort and deduplicate after each iteration
        for blockId in 1:numBlocks
            splitReqs[blockId][1] = sort(unique(splitReqs[blockId][1]))
            splitReqs[blockId][2] = sort(unique(splitReqs[blockId][2]))
        end
    end
    
    if iter >= maxIter
        @warn "Split propagation reached maximum iterations ($maxIter). Results may be incomplete."
    end
    
    return splitReqs
end

"""
Split multiple blocks with automatic propagation of split constraints across interfaces.

Inputs:
- blocks: Vector of 3D grid arrays
- splitRequests: Vector of tuples (blockId, [[i_splits], [j_splits]])
                or Dict mapping blockId => [[i_splits], [j_splits]]
- bndInfo: Boundary condition dictionary
- interInfo: Interface connectivity dictionary

Returns:
- newBlocks: Vector of split blocks (globally renumbered)
- newBndInfo: Updated boundary information
- newInterInfo: Updated interface information
"""
function SplitMultiBlock(blocks::Vector, splitRequests, bndInfo, interInfo)
    # Convert splitRequests to dict format if needed
    if splitRequests isa Vector
        splitDict = Dict{Int, Vector{Vector{Int}}}()
        for req in splitRequests
            if req isa Tuple && length(req) == 2
                splitDict[req[1]] = req[2]
            end
        end
        splitRequests = splitDict
    end
    
    # Propagate splits across interfaces
    allSplits = PropagateSplits(blocks, splitRequests, bndInfo, interInfo)
    
    # Split each block and collect results
    allNewBlocks = []
    blockMapping = Dict{Int, Vector{Int}}()  # oldBlockId => [newBlockIds...]
    nextBlockId = 1
    
    # We need to track how blocks map to track interfaces
    # Store for each original block: its sub-blocks and local interface info
    blockSplitData = []
    
    for (oldBlockId, block) in enumerate(blocks)
        splits = allSplits[oldBlockId]
        
        if isempty(splits[1]) && isempty(splits[2])
            # No splits for this block - keep it as is
            push!(allNewBlocks, block)
            blockMapping[oldBlockId] = [nextBlockId]
            
            push!(blockSplitData, Dict(
                "oldId" => oldBlockId,
                "newIds" => [nextBlockId],
                "splits" => splits,
                "subBlocks" => [block],
                "localInterfaces" => []
            ))
            
            nextBlockId += 1
        else
            # Split this block
            subBlocks, localBndInfo, localInterInfo = SplitBlock(block, splits, [], [])
            
            numSubs = length(subBlocks)
            newIds = collect(nextBlockId:(nextBlockId + numSubs - 1))
            blockMapping[oldBlockId] = newIds
            
            append!(allNewBlocks, subBlocks)
            
            push!(blockSplitData, Dict(
                "oldId" => oldBlockId,
                "newIds" => newIds,
                "splits" => splits,
                "subBlocks" => subBlocks,
                "localInterfaces" => localInterInfo
            ))
            
            nextBlockId += numSubs
        end
    end
    
    # Build new interface information
    newInterInfo = []
    
    # Add internal interfaces within each split block
    for splitData in blockSplitData
        for localInter in splitData["localInterfaces"]
            # Renumber the interface to global IDs
            localBlockA = localInter["blockA"]
            localBlockB = localInter["blockB"]
            baseId = splitData["newIds"][1]
            
            globalInter = copy(localInter)
            globalInter["blockA"] = baseId + localBlockA - 1
            globalInter["blockB"] = baseId + localBlockB - 1
            push!(newInterInfo, globalInter)
        end
    end
    
    # Map old interfaces to new interfaces between split blocks
    for oldInter in interInfo
        oldA = oldInter["blockA"]
        oldB = oldInter["blockB"]
        
        newIdsA = blockMapping[oldA]
        newIdsB = blockMapping[oldB]
        
        # If neither block was split, simple 1-to-1 mapping
        if length(newIdsA) == 1 && length(newIdsB) == 1
            newInter = copy(oldInter)
            newInter["blockA"] = newIdsA[1]
            newInter["blockB"] = newIdsB[1]
            
            # Update coordinates to match (potentially resized) blocks
            blockA = allNewBlocks[newIdsA[1]]
            blockB = allNewBlocks[newIdsB[1]]
            NxA = size(blockA, 2); NyA = size(blockA, 3)
            NxB = size(blockB, 2); NyB = size(blockB, 3)
            
            # Use same logic as UpdateInterInfo: clamp to edge coordinates
            sA = Vector{Int}(newInter["start_blkA"])
            eA = Vector{Int}(newInter["end_blkA"])
            sB = Vector{Int}(newInter["start_blkB"])
            eB = Vector{Int}(newInter["end_blkB"])
            
            sA[1] = (sA[1] == 1) ? 1 : NxA;  sA[2] = (sA[2] == 1) ? 1 : NyA;  sA[3] = 1
            eA[1] = (eA[1] == 1) ? 1 : NxA;  eA[2] = (eA[2] == 1) ? 1 : NyA;  eA[3] = 1
            sB[1] = (sB[1] == 1) ? 1 : NxB;  sB[2] = (sB[2] == 1) ? 1 : NyB;  sB[3] = 1
            eB[1] = (eB[1] == 1) ? 1 : NxB;  eB[2] = (eB[2] == 1) ? 1 : NyB;  eB[3] = 1
            
            newInter["start_blkA"] = sA; newInter["end_blkA"] = eA
            newInter["start_blkB"] = sB; newInter["end_blkB"] = eB
            
            push!(newInterInfo, newInter)
            continue
        end
        
        # At least one block was split - need to create multiple interface segments
        splitDataA = blockSplitData[oldA]
        splitDataB = blockSplitData[oldB]
        
        # Get split configuration for each block
        iSplitsA = splitDataA["splits"][1]  # i-direction splits
        jSplitsA = splitDataA["splits"][2]  # j-direction splits
        iSplitsB = splitDataB["splits"][1]
        jSplitsB = splitDataB["splits"][2]
        
        # Number of segments in each direction for each block
        numIsegsA = length(iSplitsA) + 1
        numJsegsA = length(jSplitsA) + 1
        numIsegsB = length(iSplitsB) + 1
        numJsegsB = length(jSplitsB) + 1
        
        # Determine interface orientation
        orientation = GetInterfaceOrientation(oldInter)
        
        if orientation == 1  # Horizontal interface (blocks side-by-side, j varies along interface)
            # The interface runs in the j-direction
            # Need to match j-segments between the two blocks
            # i-position is fixed at the interface
            
            # For horizontal interface, only the edge i-segments touch the interface:
            # - BlockA: rightmost i-segment (i = numIsegsA) touches the interface
            # - BlockB: leftmost i-segment (i = 1) touches the interface
            
            # Block A's j-segments should match Block B's j-segments
            if numJsegsA == numJsegsB
                # Iterate through all j-segments (these align along the interface)
                for jSeg in 1:numJsegsA
                    # Only the edge i-segments connect across the interface
                    iSegA = numIsegsA  # Rightmost segment of A touches the interface
                    iSegB = 1           # Leftmost segment of B touches the interface
                    
                    # Calculate sub-block index (row-major order: [i, j])
                    idxA = (jSeg - 1) * numIsegsA + iSegA
                    idxB = (jSeg - 1) * numIsegsB + iSegB
                    
                    if idxA <= length(newIdsA) && idxB <= length(newIdsB)
                        newInter = copy(oldInter)
                        newInter["blockA"] = newIdsA[idxA]
                        newInter["blockB"] = newIdsB[idxB]
                        
                        # Get dimensions of the new sub-blocks
                        blockA = allNewBlocks[newIdsA[idxA]]
                        blockB = allNewBlocks[newIdsB[idxB]]
                        NxA = size(blockA, 2); NyA = size(blockA, 3)
                        NxB = size(blockB, 2); NyB = size(blockB, 3)
                        
                        # For horizontal interface (j varies), i is fixed at edge
                        # BlockA is on the left (i=Nx), BlockB is on the right (i=1)
                        newInter["start_blkA"] = [NxA, 1, 1]
                        newInter["end_blkA"] = [NxA, NyA, 1]
                        newInter["start_blkB"] = [1, 1, 1]
                        newInter["end_blkB"] = [1, NyB, 1]
                        
                        push!(newInterInfo, newInter)
                    end
                end
            end
            
        else  # Vertical interface (blocks stacked vertically, i varies along interface)
            # The interface runs in the i-direction
            # Need to match i-segments between the two blocks
            # j-position is fixed at the interface
            
            # For vertical interface, only the edge j-segments touch the interface:
            # - BlockA: topmost j-segment (j = numJsegsA) touches the interface
            # - BlockB: bottommost j-segment (j = 1) touches the interface
            
            # Block A's i-segments should match Block B's i-segments
            if numIsegsA == numIsegsB
                # Iterate through all i-segments (these align along the interface)
                for iSeg in 1:numIsegsA
                    # Only the edge j-segments connect across the interface
                    jSegA = numJsegsA  # Topmost segment of A touches the interface
                    jSegB = 1           # Bottommost segment of B touches the interface
                    
                    # Calculate sub-block index (row-major order: [i, j])
                    idxA = (jSegA - 1) * numIsegsA + iSeg
                    idxB = (jSegB - 1) * numIsegsB + iSeg
                    
                    if idxA <= length(newIdsA) && idxB <= length(newIdsB)
                        newInter = copy(oldInter)
                        newInter["blockA"] = newIdsA[idxA]
                        newInter["blockB"] = newIdsB[idxB]
                        
                        # Get dimensions of the new sub-blocks
                        blockA = allNewBlocks[newIdsA[idxA]]
                        blockB = allNewBlocks[newIdsB[idxB]]
                        NxA = size(blockA, 2); NyA = size(blockA, 3)
                        NxB = size(blockB, 2); NyB = size(blockB, 3)
                        
                        # For vertical interface (i varies), j is fixed at edge
                        # BlockA is on the bottom (j=Ny), BlockB is on the top (j=1)
                        newInter["start_blkA"] = [1, NyA, 1]
                        newInter["end_blkA"] = [NxA, NyA, 1]
                        newInter["start_blkB"] = [1, 1, 1]
                        newInter["end_blkB"] = [NxB, 1, 1]
                        
                        push!(newInterInfo, newInter)
                    end
                end
            end
        end
    end
    
    # Build new boundary information
    # For each original boundary, distribute to new sub-blocks that touch it
    newBndInfo = []
    
    for bnd in bndInfo
        name = bnd["name"]
        newFaces = []
        
        for face in bnd["faces"]
            oldBlockId = face["block"]
            newIds = blockMapping[oldBlockId]
            
            if length(newIds) == 1
                # Block wasn't split - update block ID only
                newFace = Dict{String, Any}()
                newFace["block"] = newIds[1]
                newFace["start"] = copy(face["start"])
                newFace["end"] = copy(face["end"])
                # Copy any other fields from original face
                for (key, val) in face
                    if key ∉ ["block", "start", "end"]
                        newFace[key] = val
                    end
                end
                push!(newFaces, newFace)
            else
                # Block was split - determine which sub-blocks touch this boundary
                splitData = blockSplitData[oldBlockId]
                iSplits = splitData["splits"][1]
                jSplits = splitData["splits"][2]
                numIsegs = length(iSplits) + 1
                numJsegs = length(jSplits) + 1
                
                # Determine which edge this boundary is on by looking at start/end coordinates
                start = face["start"]
                endPt = face["end"]
                
                # Determine boundary edge type
                # Left edge: i=1, j varies
                # Right edge: i=Ni, j varies  
                # Bottom edge: j=1, i varies
                # Top edge: j=Nj, i varies
                
                isLeftEdge = (start[1] == 1 && endPt[1] == 1)
                isRightEdge = (start[1] != 1 && endPt[1] != 1 && start[2] != endPt[2])  # i is max, j varies
                isBottomEdge = (start[2] == 1 && endPt[2] == 1)
                isTopEdge = (start[2] != 1 && endPt[2] != 1 && start[1] != endPt[1])  # j is max, i varies
                
                # Determine which sub-blocks touch this boundary edge
                touchingSubBlocks = []
                
                if isLeftEdge
                    # Left edge: only leftmost i-segment (iSeg = 1) touches boundary
                    # All j-segments touch it
                    for jSeg in 1:numJsegs
                        idx = (jSeg - 1) * numIsegs + 1  # iSeg = 1
                        if idx <= length(newIds)
                            push!(touchingSubBlocks, (newIds[idx], idx))
                        end
                    end
                    
                elseif isRightEdge
                    # Right edge: only rightmost i-segment (iSeg = numIsegs) touches boundary
                    # All j-segments touch it
                    for jSeg in 1:numJsegs
                        idx = (jSeg - 1) * numIsegs + numIsegs  # iSeg = numIsegs
                        if idx <= length(newIds)
                            push!(touchingSubBlocks, (newIds[idx], idx))
                        end
                    end
                    
                elseif isBottomEdge
                    # Bottom edge: only bottommost j-segment (jSeg = 1) touches boundary
                    # All i-segments touch it
                    for iSeg in 1:numIsegs
                        idx = iSeg  # jSeg = 1, so index = (1-1)*numIsegs + iSeg = iSeg
                        if idx <= length(newIds)
                            push!(touchingSubBlocks, (newIds[idx], idx))
                        end
                    end
                    
                elseif isTopEdge
                    # Top edge: only topmost j-segment (jSeg = numJsegs) touches boundary
                    # All i-segments touch it
                    for iSeg in 1:numIsegs
                        idx = (numJsegs - 1) * numIsegs + iSeg  # jSeg = numJsegs
                        if idx <= length(newIds)
                            push!(touchingSubBlocks, (newIds[idx], idx))
                        end
                    end
                end
                
                # Create boundary faces for touching sub-blocks
                for (subBlockId, subBlockIdx) in touchingSubBlocks
                    newFace = Dict{String, Any}()
                    newFace["block"] = subBlockId
                    newFace["start"] = copy(face["start"])
                    newFace["end"] = copy(face["end"])
                    # Copy any other fields from original face
                    for (key, val) in face
                        if key ∉ ["block", "start", "end"]
                            newFace[key] = val
                        end
                    end
                    push!(newFaces, newFace)
                end
            end
        end
        
        if !isempty(newFaces)
            push!(newBndInfo, Dict("name" => name, "faces" => newFaces))
        end
    end
    
    # Update boundary info to match new block sizes
    UpdateBndInfo!(newBndInfo, allNewBlocks)
    
    return allNewBlocks, newBndInfo, newInterInfo
end
