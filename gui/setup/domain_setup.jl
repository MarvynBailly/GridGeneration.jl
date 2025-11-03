"""
Domain setup and initialization for the Grid Generation GUI.
"""

"""
    setup_rectangle_domain(width, height, num_points_width, num_points_height)

Create a simple rectangular domain with TFI initial grid.

# Arguments
- `width`: Width of the rectangle
- `height`: Height of the rectangle
- `num_points_width`: Number of grid points in the width direction
- `num_points_height`: Number of grid points in the height direction

# Returns
- `initialGrid`: Vector containing the initial block
- `initialBndInfo`: Boundary information dictionary
- `initialInterfaceInfo`: Interface information dictionary (empty for single block)
"""
function setup_rectangle_domain(width, height, num_points_width, num_points_height)
    # Define boundary curves
    top = hcat(range(0, stop=width, length=num_points_width), fill(height, num_points_width))
    bottom = hcat(range(0, stop=width, length=num_points_width), fill(0.0, num_points_width))
    right = hcat(fill(width, num_points_height), range(0, stop=height, length=num_points_height))
    left = hcat(fill(0.0, num_points_height), range(0, stop=height, length=num_points_height))
    
    # Create initial block using TFI
    initial_block = GridGeneration.TFI([top, right, bottom, left])
    
    # Setup boundary information
    initialBndInfo = []
    push!(initialBndInfo, Dict(
        "name" => "bottom", 
        "faces" => [Dict("block" => 1, "start" => [1,1,1], "end" => [num_points_width,1,1])]
    ))
    push!(initialBndInfo, Dict(
        "name" => "right", 
        "faces" => [Dict("block" => 1, "start" => [num_points_width,1,1], "end" => [num_points_width,num_points_height,1])]
    ))
    push!(initialBndInfo, Dict(
        "name" => "top", 
        "faces" => [Dict("block" => 1, "start" => [1,num_points_height,1], "end" => [num_points_width,num_points_height,1])]
    ))
    push!(initialBndInfo, Dict(
        "name" => "left", 
        "faces" => [Dict("block" => 1, "start" => [1,1,1], "end" => [1,num_points_height,1])]
    ))
    
    # No interfaces for single block
    initialInterfaceInfo = []
    
    return [initial_block], initialBndInfo, initialInterfaceInfo
end

"""
    detect_uneven_interfaces(blocks, interfaceInfo)

Detect interfaces that don't span the entire block edge and generate split requirements
to create aligned sub-blocks.

# Arguments
- `blocks`: Vector of grid blocks
- `interfaceInfo`: Interface information dictionary (must not contain self-referential interfaces)

# Returns
- `splitRequests`: Dict mapping blockId => [[i_splits], [j_splits]] for uneven interfaces

# Algorithm
For each interface, check if the start/end indices align with block boundaries:
- If an interface doesn't start at i=1 or j=1, or doesn't end at i=Ni or j=Nj,
  it's an uneven interface
- Generate split at the misaligned index to create separate sub-blocks
"""
function detect_uneven_interfaces(blocks, interfaceInfo)
    splitRequests = Dict{Int, Vector{Vector{Int}}}()
    
    for interface in interfaceInfo
        blockA_id = interface["blockA"]
        blockB_id = interface["blockB"]
        
        blockA = blocks[blockA_id]
        blockB = blocks[blockB_id]
        
        niA, njA = size(blockA, 2), size(blockA, 3)
        niB, njB = size(blockB, 2), size(blockB, 3)
        
        startA = interface["start_blkA"]
        endA = interface["end_blkA"]
        startB = interface["start_blkB"]
        endB = interface["end_blkB"]
        
        # Initialize split arrays for this block if not already present
        if !haskey(splitRequests, blockA_id)
            splitRequests[blockA_id] = [Int[], Int[]]
        end
        if !haskey(splitRequests, blockB_id)
            splitRequests[blockB_id] = [Int[], Int[]]
        end
        
        # Check Block A for uneven interface
        # Determine if this is a vertical (i varies) or horizontal (j varies) interface
        if startA[1] != endA[1]
            # Vertical interface on A (i varies, j is constant)
            # Check if the constant j is at a boundary
            if startA[2] != 1 && startA[2] != njA
                push!(splitRequests[blockA_id][2], startA[2])
            end
            # Check if i spans the full range
            if startA[1] != 1
                push!(splitRequests[blockA_id][1], startA[1])
            end
            if endA[1] != niA
                push!(splitRequests[blockA_id][1], endA[1] + 1)
            end
        else
            # Horizontal interface on A (j varies, i is constant)
            # Check if the constant i is at a boundary
            if startA[1] != 1 && startA[1] != niA
                push!(splitRequests[blockA_id][1], startA[1])
            end
            # Check if j spans the full range
            if startA[2] != 1
                push!(splitRequests[blockA_id][2], startA[2])
            end
            if endA[2] != njA
                push!(splitRequests[blockA_id][2], endA[2] + 1)
            end
        end
        
        # Check Block B for uneven interface
        if startB[1] != endB[1]
            # Vertical interface on B (i varies, j is constant)
            if startB[2] != 1 && startB[2] != njB
                push!(splitRequests[blockB_id][2], startB[2])
            end
            if startB[1] != 1
                push!(splitRequests[blockB_id][1], startB[1])
            end
            if endB[1] != niB
                push!(splitRequests[blockB_id][1], endB[1] + 1)
            end
        else
            # Horizontal interface on B (j varies, i is constant)
            if startB[1] != 1 && startB[1] != niB
                push!(splitRequests[blockB_id][1], startB[1])
            end
            if startB[2] != 1
                push!(splitRequests[blockB_id][2], startB[2])
            end
            if endB[2] != njB
                push!(splitRequests[blockB_id][2], endB[2] + 1)
            end
        end
    end
    
    # Deduplicate and sort split indices
    for blockId in keys(splitRequests)
        splitRequests[blockId][1] = unique(sort(splitRequests[blockId][1]))
        splitRequests[blockId][2] = unique(sort(splitRequests[blockId][2]))
    end
    
    return splitRequests
end

"""
    convert_turtle_to_1based_indexing!(bndInfo, interfaceInfo)

Convert Turtle grid boundary and interface data from 0-based to 1-based indexing.

# Arguments
- `bndInfo`: Boundary information dictionary (modified in-place)
- `interfaceInfo`: Interface information dictionary (modified in-place)
"""
function convert_turtle_to_1based_indexing!(bndInfo, interfaceInfo)
    # Rename and reformat boundary info for 1-based indexing
    for bnd in bndInfo
        # Rename faceInfo to faces
        if haskey(bnd, "faceInfo")
            bnd["faces"] = bnd["faceInfo"]
            delete!(bnd, "faceInfo")
        end
        
        # Remove redundant nbrFaces field
        if haskey(bnd, "nbrFaces")
            delete!(bnd, "nbrFaces")
        end
        
        # Convert 0-based to 1-based indexing
        for face in bnd["faces"]
            face["block"] = face["block"] + 1  # Convert block ID
            face["start"] = face["start"] .+ 1  # Convert start indices
            face["end"] = face["end"] .+ 1      # Convert end indices
        end
    end
    
    # Convert interface info from 0-based to 1-based
    for itf in interfaceInfo
        itf["blockA"] = itf["blockA"] + 1
        itf["blockB"] = itf["blockB"] + 1
        itf["start_blkA"] = itf["start_blkA"] .+ 1
        itf["end_blkA"] = itf["end_blkA"] .+ 1
        itf["start_blkB"] = itf["start_blkB"] .+ 1
        itf["end_blkB"] = itf["end_blkB"] .+ 1
    end
end

"""
    create_internal_interfaces_for_split_block(oldBlockId, splitRequests, blockMapping, newBlocks)

Create internal interfaces between sub-blocks within a split block.

# Arguments
- `oldBlockId`: Original block ID before splitting
- `splitRequests`: Dictionary of split requests
- `blockMapping`: Mapping from old block IDs to new block IDs
- `newBlocks`: Vector of new blocks after splitting

# Returns
- `internalInterfaces`: Vector of interface dictionaries for internal interfaces
"""
function create_internal_interfaces_for_split_block(oldBlockId, splitRequests, blockMapping, newBlocks)
    internalInterfaces = []
    
    if !haskey(splitRequests, oldBlockId) || length(blockMapping[oldBlockId]) <= 1
        return internalInterfaces
    end
    
    # This block was split - add internal interfaces
    iSplits = splitRequests[oldBlockId][1]
    jSplits = splitRequests[oldBlockId][2]
    numIsegs = length(iSplits) + 1
    numJsegs = length(jSplits) + 1
    newIds = blockMapping[oldBlockId]
    
    # Add vertical interfaces (between i-segments, constant i at split)
    for iSplitIdx in 1:length(iSplits)
        for jSeg in 1:numJsegs
            # Left sub-block index
            leftIdx = (jSeg - 1) * numIsegs + iSplitIdx
            # Right sub-block index
            rightIdx = leftIdx + 1
            
            if leftIdx <= length(newIds) && rightIdx <= length(newIds)
                leftBlock = newBlocks[newIds[leftIdx]]
                rightBlock = newBlocks[newIds[rightIdx]]
                niLeft = size(leftBlock, 2)
                njLeft = size(leftBlock, 3)
                njRight = size(rightBlock, 3)
                
                push!(internalInterfaces, Dict(
                    "blockA" => newIds[leftIdx],
                    "blockB" => newIds[rightIdx],
                    "start_blkA" => [niLeft, 1, 1],
                    "end_blkA" => [niLeft, njLeft, 1],
                    "start_blkB" => [1, 1, 1],
                    "end_blkB" => [1, njRight, 1],
                    "offset" => [0.0, 0.0, 0.0],
                    "angle" => 0.0
                ))
            end
        end
    end
    
    # Add horizontal interfaces (between j-segments, constant j at split)
    for jSplitIdx in 1:length(jSplits)
        for iSeg in 1:numIsegs
            # Bottom sub-block index
            bottomIdx = (jSplitIdx - 1) * numIsegs + iSeg
            # Top sub-block index
            topIdx = bottomIdx + numIsegs
            
            if bottomIdx <= length(newIds) && topIdx <= length(newIds)
                bottomBlock = newBlocks[newIds[bottomIdx]]
                topBlock = newBlocks[newIds[topIdx]]
                niBottom = size(bottomBlock, 2)
                njBottom = size(bottomBlock, 3)
                niTop = size(topBlock, 2)
                
                push!(internalInterfaces, Dict(
                    "blockA" => newIds[bottomIdx],
                    "blockB" => newIds[topIdx],
                    "start_blkA" => [1, njBottom, 1],
                    "end_blkA" => [niBottom, njBottom, 1],
                    "start_blkB" => [1, 1, 1],
                    "end_blkB" => [niTop, 1, 1],
                    "offset" => [0.0, 0.0, 0.0],
                    "angle" => 0.0
                ))
            end
        end
    end
    
    return internalInterfaces
end

"""
    update_external_interface(oldInter, splitRequests, blockMapping, newBlocks)

Update an existing interface to reference new block IDs and coordinates after splitting.

# Arguments
- `oldInter`: Original interface dictionary
- `splitRequests`: Dictionary of split requests
- `blockMapping`: Mapping from old block IDs to new block IDs
- `newBlocks`: Vector of new blocks after splitting

# Returns
- `newInter`: Updated interface dictionary
"""
function update_external_interface(oldInter, splitRequests, blockMapping, newBlocks)
    oldA = oldInter["blockA"]
    oldB = oldInter["blockB"]
    newIdsA = blockMapping[oldA]
    newIdsB = blockMapping[oldB]
    
    # Determine which sub-block of A contains the interface
    if length(newIdsA) == 1
        newBlockA = newIdsA[1]
    else
        # Block A was split - need to determine which sub-block has the interface
        iSplits = splitRequests[oldA][1]
        jSplits = splitRequests[oldA][2]
        numIsegs = length(iSplits) + 1
        numJsegs = length(jSplits) + 1
        
        # Find which segment contains the interface start point
        startA = oldInter["start_blkA"]
        
        # Determine i-segment
        iSeg = 1
        for (idx, split) in enumerate(iSplits)
            if startA[1] >= split
                iSeg = idx + 1
            end
        end
        
        # Determine j-segment  
        jSeg = 1
        for (idx, split) in enumerate(jSplits)
            if startA[2] >= split
                jSeg = idx + 1
            end
        end
        
        # Calculate sub-block index (row-major: [i, j])
        subIdx = (jSeg - 1) * numIsegs + iSeg
        newBlockA = newIdsA[subIdx]
    end
    
    # Block B should not be split for uneven interface fix
    newBlockB = newIdsB[1]
    
    # Create new interface with updated block IDs
    newInter = copy(oldInter)
    newInter["blockA"] = newBlockA
    newInter["blockB"] = newBlockB
    
    # Update coordinates to match new block sizes
    blockA = newBlocks[newBlockA]
    blockB = newBlocks[newBlockB]
    NxA = size(blockA, 2); NyA = size(blockA, 3)
    NxB = size(blockB, 2); NyB = size(blockB, 3)
    
    # For split blocks, the interface should now span the entire edge of the sub-block
    sA = Vector{Int}(oldInter["start_blkA"])
    eA = Vector{Int}(oldInter["end_blkA"])
    
    if length(newIdsA) > 1
        # Block A was split - interface is now on the edge of the sub-block
        if sA[1] == eA[1]
            # Vertical interface (constant i) - map to edge
            sA[1] = (sA[1] == oldInter["start_blkA"][1] && sA[1] > 1) ? NxA : sA[1]
            eA[1] = sA[1]
            # j spans the full range of the sub-block
            sA[2] = 1
            eA[2] = NyA
        else
            # Horizontal interface (constant j) - map to edge
            sA[2] = (sA[2] == oldInter["start_blkA"][2] && sA[2] > 1) ? NyA : sA[2]
            eA[2] = sA[2]
            # i spans the full range of the sub-block
            sA[1] = 1
            eA[1] = NxA
        end
    else
        # Clamp to valid range
        sA[1] = clamp(sA[1], 1, NxA)
        sA[2] = clamp(sA[2], 1, NyA)
        eA[1] = clamp(eA[1], 1, NxA)
        eA[2] = clamp(eA[2], 1, NyA)
    end
    sA[3] = 1; eA[3] = 1
    
    # Block B coordinates should remain valid (not split)
    sB = Vector{Int}(oldInter["start_blkB"])
    eB = Vector{Int}(oldInter["end_blkB"])
    sB[1] = clamp(sB[1], 1, NxB)
    sB[2] = clamp(sB[2], 1, NyB)
    eB[1] = clamp(eB[1], 1, NxB)
    eB[2] = clamp(eB[2], 1, NyB)
    sB[3] = 1; eB[3] = 1
    
    newInter["start_blkA"] = sA; newInter["end_blkA"] = eA
    newInter["start_blkB"] = sB; newInter["end_blkB"] = eB
    
    return newInter
end

"""
    update_boundaries_for_split_blocks(bndInfo, blocks, splitRequests, blockMapping, newBlocks, newInterfaceInfo)

Update boundary information after block splitting, distributing boundaries to correct sub-blocks.

# Arguments
- `bndInfo`: Original boundary information
- `blocks`: Original blocks before splitting
- `splitRequests`: Dictionary of split requests
- `blockMapping`: Mapping from old block IDs to new block IDs
- `newBlocks`: Vector of new blocks after splitting
- `newInterfaceInfo`: Updated interface information (used to detect which edges are interfaces)

# Returns
- `newBndInfo`: Updated boundary information
"""
function update_boundaries_for_split_blocks(bndInfo, blocks, splitRequests, blockMapping, newBlocks, newInterfaceInfo)
    # Helper function to check if an edge has an interface
    function has_interface_on_edge(blockId, edge_type)
        for itf in newInterfaceInfo
            if itf["blockA"] == blockId || itf["blockB"] == blockId
                # Get the relevant interface coordinates
                if itf["blockA"] == blockId
                    startItf = itf["start_blkA"]
                    endItf = itf["end_blkA"]
                else
                    startItf = itf["start_blkB"]
                    endItf = itf["end_blkB"]
                end
                
                subBlock = newBlocks[blockId]
                ni = size(subBlock, 2)
                nj = size(subBlock, 3)
                
                # Check if interface is on the specified edge
                if edge_type == :left && startItf[1] == 1 && endItf[1] == 1
                    return true
                elseif edge_type == :right && startItf[1] == ni && endItf[1] == ni
                    return true
                elseif edge_type == :bottom && startItf[2] == 1 && endItf[2] == 1
                    return true
                elseif edge_type == :top && startItf[2] == nj && endItf[2] == nj
                    return true
                end
            end
        end
        return false
    end
    
    newBndInfo = []
    for bnd in bndInfo
        name = bnd["name"]
        newFaces = []
        
        for face in bnd["faces"]
            oldBlockId = face["block"]
            newIds = blockMapping[oldBlockId]
            
            if length(newIds) == 1
                # Block wasn't split - keep as is
                newFace = copy(face)
                newFace["block"] = newIds[1]
                push!(newFaces, newFace)
            else
                # Block was split - determine which sub-blocks touch this boundary
                iSplits = splitRequests[oldBlockId][1]
                jSplits = splitRequests[oldBlockId][2]
                numIsegs = length(iSplits) + 1
                numJsegs = length(jSplits) + 1
                
                # Determine boundary edge type
                startFace = face["start"]
                endFace = face["end"]
                
                # Check which edge this boundary is on
                oldBlock = blocks[oldBlockId]
                niOld = size(oldBlock, 2)
                njOld = size(oldBlock, 3)
                
                isLeftEdge = (startFace[1] == 1 && endFace[1] == 1)
                isRightEdge = (startFace[1] == niOld && endFace[1] == niOld)
                isBottomEdge = (startFace[2] == 1 && endFace[2] == 1)
                isTopEdge = (startFace[2] == njOld && endFace[2] == njOld)
                
                # Distribute boundary to appropriate sub-blocks (excluding those with interfaces)
                if isLeftEdge
                    # Left edge: only leftmost i-segment (iSeg=1) touches boundary
                    for jSeg in 1:numJsegs
                        subIdx = (jSeg - 1) * numIsegs + 1
                        if subIdx <= length(newIds)
                            blockId = newIds[subIdx]
                            if !has_interface_on_edge(blockId, :left)
                                subBlock = newBlocks[blockId]
                                newFace = Dict(
                                    "block" => blockId,
                                    "start" => [1, 1, 1],
                                    "end" => [1, size(subBlock, 3), 1]
                                )
                                push!(newFaces, newFace)
                            end
                        end
                    end
                elseif isRightEdge
                    # Right edge: only rightmost i-segment (iSeg=numIsegs) touches boundary
                    for jSeg in 1:numJsegs
                        subIdx = (jSeg - 1) * numIsegs + numIsegs
                        if subIdx <= length(newIds)
                            blockId = newIds[subIdx]
                            if !has_interface_on_edge(blockId, :right)
                                subBlock = newBlocks[blockId]
                                ni = size(subBlock, 2)
                                nj = size(subBlock, 3)
                                newFace = Dict(
                                    "block" => blockId,
                                    "start" => [ni, 1, 1],
                                    "end" => [ni, nj, 1]
                                )
                                push!(newFaces, newFace)
                            end
                        end
                    end
                elseif isBottomEdge
                    # Bottom edge: only bottommost j-segment (jSeg=1) touches boundary
                    for iSeg in 1:numIsegs
                        subIdx = iSeg
                        if subIdx <= length(newIds)
                            blockId = newIds[subIdx]
                            if !has_interface_on_edge(blockId, :bottom)
                                subBlock = newBlocks[blockId]
                                ni = size(subBlock, 2)
                                newFace = Dict(
                                    "block" => blockId,
                                    "start" => [1, 1, 1],
                                    "end" => [ni, 1, 1]
                                )
                                push!(newFaces, newFace)
                            end
                        end
                    end
                elseif isTopEdge
                    # Top edge: only topmost j-segment (jSeg=numJsegs) touches boundary
                    for iSeg in 1:numIsegs
                        subIdx = (numJsegs - 1) * numIsegs + iSeg
                        if subIdx <= length(newIds)
                            blockId = newIds[subIdx]
                            if !has_interface_on_edge(blockId, :top)
                                subBlock = newBlocks[blockId]
                                ni = size(subBlock, 2)
                                nj = size(subBlock, 3)
                                newFace = Dict(
                                    "block" => blockId,
                                    "start" => [1, nj, 1],
                                    "end" => [ni, nj, 1]
                                )
                                push!(newFaces, newFace)
                            end
                        end
                    end
                end
            end
        end
        
        if !isempty(newFaces)
            push!(newBndInfo, Dict("name" => name, "faces" => newFaces))
        end
    end
    
    return newBndInfo
end

"""
    fix_uneven_interfaces!(blocks, interfaceInfo, bndInfo)

Fix uneven interfaces by splitting blocks to align interface boundaries.

# Arguments
- `blocks`: Vector of grid blocks (modified in-place)
- `interfaceInfo`: Interface information dictionary (modified in-place)
- `bndInfo`: Boundary information dictionary (modified in-place)

# Returns
- Nothing (modifies arguments in-place)
"""
function fix_uneven_interfaces!(blocks, interfaceInfo, bndInfo)
    splitRequests = detect_uneven_interfaces(blocks, interfaceInfo)
    
    # Remove empty split requests
    splitRequests = filter(p -> !isempty(p.second[1]) || !isempty(p.second[2]), splitRequests)
    
    if isempty(splitRequests)
        println("No uneven interfaces detected - all interfaces properly aligned")
        return
    end
    
    println("Detected uneven interface(s) - splitting blocks to align interfaces (NO propagation)...")
    for (blockId, splits) in splitRequests
        println("  Block $blockId: i-splits=$(splits[1]), j-splits=$(splits[2])")
    end
    
    # Split blocks WITHOUT propagation
    newBlocks = []
    blockMapping = Dict{Int, Vector{Int}}()  # oldBlockId => [newBlockIds...]
    nextBlockId = 1
    
    for oldBlockId in 1:length(blocks)
        if haskey(splitRequests, oldBlockId)
            # Split this block
            subBlocks, _, _ = GridGeneration.SplitBlock(blocks[oldBlockId], splitRequests[oldBlockId], [], [])
            numSubs = length(subBlocks)
            newIds = collect(nextBlockId:(nextBlockId + numSubs - 1))
            blockMapping[oldBlockId] = newIds
            append!(newBlocks, subBlocks)
            nextBlockId += numSubs
        else
            # Keep block as-is
            push!(newBlocks, blocks[oldBlockId])
            blockMapping[oldBlockId] = [nextBlockId]
            nextBlockId += 1
        end
    end
    
    # Update interfaces manually (no propagation logic)
    newInterfaceInfo = []
    
    # First, add internal interfaces within split blocks
    for oldBlockId in 1:length(blocks)
        internalInterfaces = create_internal_interfaces_for_split_block(
            oldBlockId, splitRequests, blockMapping, newBlocks
        )
        append!(newInterfaceInfo, internalInterfaces)
    end
    
    # Second, update existing interfaces between different blocks
    for oldInter in interfaceInfo
        newInter = update_external_interface(oldInter, splitRequests, blockMapping, newBlocks)
        push!(newInterfaceInfo, newInter)
    end
    
    # Update boundary info - distribute to correct sub-blocks
    newBndInfo = update_boundaries_for_split_blocks(
        bndInfo, blocks, splitRequests, blockMapping, newBlocks, newInterfaceInfo
    )
    
    # Update to final values (modify in-place via reassignment in caller)
    empty!(blocks)
    append!(blocks, newBlocks)
    
    empty!(interfaceInfo)
    append!(interfaceInfo, newInterfaceInfo)
    
    empty!(bndInfo)
    append!(bndInfo, newBndInfo)
    
    # Update boundary info to match new block sizes
    GridGeneration.UpdateBndInfo!(bndInfo, blocks)
    
    println("After fixing: $(length(blocks)) blocks, $(length(interfaceInfo)) interfaces")
end

"""
    setup_turtle_grid_domain(metricFieldFile, gridFolder; fix_uneven_interfaces=true)

Load domain from Turtle grid files and optionally fix uneven interfaces.

# Arguments
- `metricFieldFile`: Path to the metric field file
- `gridFolder`: Path to the folder containing grid files
- `fix_uneven_interfaces`: If true, automatically split blocks to align uneven interfaces (default: true)

# Returns
- `initialGrid`: Vector of initial blocks
- `initialBndInfo`: Boundary information dictionary
- `initialInterfaceInfo`: Interface information dictionary
- `M`: Metric function for the domain
"""
function setup_turtle_grid_domain(metricFieldFile, gridFile; fix_uneven_interfaces=true)
    # Load Turtle grid data
    metricData, datatype, gridfile = GridGeneration.readTurtleFields(metricFieldFile)
    # gridfile = joinpath(gridFolder, gridfile)
    blocks, centers, Xfa, interfaceInfo, bndInfo = GridGeneration.ImportTurtleGrid(gridFile)
    
    # Setup metric field
    tree, refs = GridGeneration.setup_metric_tree(centers)
    M = (x, y) -> GridGeneration.find_nearest_kd(metricData, tree, refs, x, y)

    # Convert from 0-based to 1-based indexing
    convert_turtle_to_1based_indexing!(bndInfo, interfaceInfo)
    
    # Filter out self-referential interfaces (periodic k-direction interfaces)
    # Keep only interfaces where blockA != blockB
    interfaceInfo = filter(itf -> itf["blockA"] != itf["blockB"], interfaceInfo)
    
    # Optionally detect and fix uneven interfaces
    if fix_uneven_interfaces
        fix_uneven_interfaces!(blocks, interfaceInfo, bndInfo)
    end
    
    return blocks, bndInfo, interfaceInfo, M
end

"""
    create_default_metric(s=1000)

Create a simple default metric function.

# Arguments
- `s`: Scaling factor for the metric

# Returns
- `M`: Metric function M(x,y) that returns a 2-element vector
"""
function create_default_metric(s=1000)
    return (x, y) -> [s * x^2, s]
end
