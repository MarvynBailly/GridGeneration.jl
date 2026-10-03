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
                push!(splitRequests[blockA_id][1], endA[1])  # split on the shared node line
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
                push!(splitRequests[blockA_id][2], endA[2])  # split on the shared node line
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
                push!(splitRequests[blockB_id][1], endB[1])  # split on the shared node line
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
                push!(splitRequests[blockB_id][2], endB[2])  # split on the shared node line
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
    newBlockA, sA, eA = remap_interface_side(oldInter["blockA"], oldInter["start_blkA"], oldInter["end_blkA"],
                                             splitRequests, blockMapping, newBlocks)
    newBlockB, sB, eB = remap_interface_side(oldInter["blockB"], oldInter["start_blkB"], oldInter["end_blkB"],
                                             splitRequests, blockMapping, newBlocks)

    newInter = copy(oldInter)
    newInter["blockA"] = newBlockA
    newInter["blockB"] = newBlockB
    newInter["start_blkA"] = sA; newInter["end_blkA"] = eA
    newInter["start_blkB"] = sB; newInter["end_blkB"] = eB
    return newInter
end

# Map one side of an interface (block id + start/end indices on the old block) onto the
# sub-block that contains it after splitting, with indices local to that sub-block.
function remap_interface_side(oldId, start, stop, splitRequests, blockMapping, newBlocks)
    newIds = blockMapping[oldId]
    s = Vector{Int}(start[1:2])
    e = Vector{Int}(stop[1:2])

    if length(newIds) == 1
        newId = newIds[1]
        blk = newBlocks[newId]
        s .= clamp.(s, 1, [size(blk, 2), size(blk, 3)])
        e .= clamp.(e, 1, [size(blk, 2), size(blk, 3)])
        return newId, [s; 1], [e; 1]
    end

    iSplits, jSplits = splitRequests[oldId]
    numIsegs = length(iSplits) + 1

    # Segment containing the interface: the last split at or below the start index.
    # Sub-block (iSeg, jSeg) starts at node iStart = iSplits[iSeg-1] (or 1) in the old block.
    iLo, iHi = minmax(s[1], e[1])
    jLo, jHi = minmax(s[2], e[2])
    if any(k -> iLo < k < iHi, iSplits) || any(k -> jLo < k < jHi, jSplits)
        error("Interface on block $oldId ($(start[1:2]) to $(stop[1:2])) crosses a split line " *
              "(i-splits $iSplits, j-splits $jSplits); splitting an interface into pieces is not supported")
    end
    iSeg = 1 + count(<=(iLo), iSplits)
    jSeg = 1 + count(<=(jLo), jSplits)
    newId = newIds[(jSeg - 1) * numIsegs + iSeg]

    iOffset = iSeg == 1 ? 0 : iSplits[iSeg - 1] - 1
    jOffset = jSeg == 1 ? 0 : jSplits[jSeg - 1] - 1
    s .-= [iOffset, jOffset]
    e .-= [iOffset, jOffset]

    blk = newBlocks[newId]
    s .= clamp.(s, 1, [size(blk, 2), size(blk, 3)])
    e .= clamp.(e, 1, [size(blk, 2), size(blk, 3)])
    return newId, [s; 1], [e; 1]
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
        @info("No uneven interfaces detected - all interfaces properly aligned")
        return
    end
    
    @info("Detected uneven interface(s) - splitting blocks to align interfaces (NO propagation)...")
    for (blockId, splits) in splitRequests
        @info("  Block $blockId: i-splits=$(splits[1]), j-splits=$(splits[2])")
    end
    
    # Split blocks WITHOUT propagation
    newBlocks = []
    blockMapping = Dict{Int, Vector{Int}}()  # oldBlockId => [newBlockIds...]
    nextBlockId = 1
    
    for oldBlockId in 1:length(blocks)
        if haskey(splitRequests, oldBlockId)
            # Split this block
            subBlocks, _, _ = SplitBlock(blocks[oldBlockId], splitRequests[oldBlockId], [], [])
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
    UpdateBndInfo!(bndInfo, blocks)
    
    @info("After fixing: $(length(blocks)) blocks, $(length(interfaceInfo)) interfaces")
end

"""
    load_turtle_grid(gridFile; fix_uneven_interfaces=true) -> (blocks, bndInfo, interInfo, centers)

Read a Tortuga grid file and convert it to the package's conventions: `[2, Ni, Nj]` blocks
(first k-layer), and 1-based boundary/interface information with boundary faces under
`"faces"`. Periodic k-direction self-interfaces are dropped. With `fix_uneven_interfaces=true`,
blocks whose interfaces do not span a whole edge are split so every interface is edge-to-edge.

Also returns the cell-centre coordinates `centers` (before any splitting), which index a metric
field read with [`readTurtleFields`](@ref).
"""
function load_turtle_grid(gridFile; fix_uneven_interfaces=true)
    blocks, centers, _, interfaceInfo, bndInfo = ImportTurtleGrid(gridFile)

    convert_turtle_to_1based_indexing!(bndInfo, interfaceInfo)

    # Drop self-referential interfaces (periodic k-direction interfaces)
    interfaceInfo = filter(itf -> itf["blockA"] != itf["blockB"], interfaceInfo)

    if fix_uneven_interfaces
        fix_uneven_interfaces!(blocks, interfaceInfo, bndInfo)
    end

    return blocks, bndInfo, interfaceInfo, centers
end

"""
    setup_turtle_grid_domain(metricFieldFile, gridFile; metricGridFile=nothing, fix_uneven_interfaces=true)
        -> (blocks, bndInfo, interInfo, M)

Load the Tortuga grid `gridFile` with [`load_turtle_grid`](@ref) and the metric field
`metricFieldFile` with [`readTurtleFields`](@ref). Returns the blocks, 1-based
boundary/interface information, and a metric function `M(x, y) -> [M11, M22]` that looks up
the metric at the nearest cell centre. The result can be passed straight to
[`GenerateGrid`](@ref) or [`SolveAllBlocks`](@ref).

The metric is stored per cell of the grid it was computed on, which need not be `gridFile`
(e.g. when regenerating a grid that was itself produced earlier). That grid is
`metricGridFile` if given; otherwise the grid named in the field file's header, if it exists
next to `metricFieldFile`; otherwise `gridFile`. Its cell layout must match the field.
"""
function setup_turtle_grid_domain(metricFieldFile, gridFile; metricGridFile=nothing, fix_uneven_interfaces=true)
    metricData, _, headerGrid = readTurtleFields(metricFieldFile)
    blocks, bndInfo, interfaceInfo, centers = load_turtle_grid(gridFile; fix_uneven_interfaces=fix_uneven_interfaces)

    if metricGridFile === nothing
        candidate = joinpath(dirname(metricFieldFile), strip(headerGrid))
        metricGridFile = (!isempty(strip(headerGrid)) && isfile(candidate)) ? candidate : gridFile
    end
    metricCenters = abspath(metricGridFile) == abspath(gridFile) ? centers : ImportTurtleGrid(metricGridFile)[2]
    check_metric_layout(metricData, metricCenters, metricFieldFile, metricGridFile, headerGrid)

    tree, refs = setup_metric_tree(metricCenters)
    M = (x, y) -> find_nearest_kd(metricData, tree, refs, x, y)

    return blocks, bndInfo, interfaceInfo, M
end

# The metric field must have one array per block of the metric grid, sized like its cells.
function check_metric_layout(metricData, centers, metricFieldFile, metricGridFile, headerGrid)
    M11 = metricData["M11"]
    cellSizes = [size(c)[2:3] for c in centers]
    fieldSizes = [size(m)[1:2] for m in M11]
    if cellSizes != fieldSizes
        error("Metric field $(basename(metricFieldFile)) (blocks of $(fieldSizes) cells, computed on " *
              "\"$(strip(headerGrid))\") does not match the cells of $(basename(metricGridFile)) " *
              "($(cellSizes)). Pass the grid the metric was computed on as `metricGridFile`.")
    end
end
