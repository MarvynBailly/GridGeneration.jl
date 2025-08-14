# Single Block with Splitting
Next, let's build `GridGeneration` functions that allow the user to add splits into the blocks. The flow I'm thinking will be
- Input single Tortuga block and desired split locations with the block
- break the block into sub-blocks at the split locations
  - Create new interfaces and add to `interInfo`
  - Update `bndInfo` for the blocks
  - Save block in blocks
- For (blockId, block) in enumerate(blocks)
  - for each dir in unvistedDirs[blockId]
    - get neighboring block ids
    - add this block and neighboring blocks to computeBlocks
    - for each edge in computeBlocks, solve for optimal $N$
    - take optimal $N$ for this direction to be the max of all opt $N$s
    - solve ODE with optimal $N$ along dir for each computeBlock
    - mark dir as visited for each computeBlock

Let's make this into two main `GridGeneration` functions: `GridGeneration.SplitBlock(block, splitLocations, bndInfo, interInfo)` and `GridGeneration.SolveAllBlocks(metric, blocks, bndInfo, interInfo).` And thus all we have to do is
- Input single Tortuga block and desired split locations with the block
- Split block into multi-block grid and update bndInfo and interInfo
- Solve for the optimal distribution in each block



## Split Block Algorithm
We want something along the lines:

```julia
splitLocations = [
    [ horizontal split indices... ],   # split along the x axis
    [ vertical split indices... ]      # split along the y axis
]

blocks, bndInfo, interInfo = GridGeneration.SplitBlock(block, splitLocations, bndInfo, interInfo)
```

where the `SplitBlock` function 

```julia
function SplitBlock(block, splitLocations, bndInfo, interInfo)
    blocks = []
    horzSplits = [1, splitLocations[1]..., size(block, 2)]
    vertSplits = [1, splitLocations[2]..., size(block, 3)]

    parentni = size(block, 2)
    parentnj = size(block, 3)
    parentnk = 1

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

            # get boundary info 
            blockInfo = Dict(
                "block" => blockId,
                "start" => (horzSplits[i], vertSplits[j]),
                "end" => (horzSplits[i+1], vertSplits[j+1]),
                "parentDims" => (parentni, parentnj, parentnk)
            )

            boundaries = GetTouchingBoundaries(blockInfo, bndInfo)
            append!(blockBoundaries, boundaries)
            
            # get interface info
            # if not at the end, look forward to the next block
            if i < length(horzSplits) - 1
                blockBId = blockId + 1
                # println("blockID: ", blockId, " next blockID: ", blockBId)
                push!(interInfo, Dict(
                    "blockA" => blockId, "start_blkA" => [ni,1,1], "end_blkA" => [ni,nj,nk],
                    "blockB" => blockBId, "start_blkB" => [1,1,1], "end_blkB" => [1,nj,nk],
                    "offset" => [0.0, 0.0, 0.0], "angle" => 0.0))
            end

            # if not at top, look up
            if j < length(vertSplits) - 1
                blockBId = blockId + length(horzSplits) - 1
                # println("blockID: ", blockId, " next blockID: ", blockBId)
                push!(interInfo, Dict(
                    "blockA" => blockId, "start_blkA" => [1,nj,1], "end_blkA" => [ni,nj,nk],
                    "blockB" => blockBId, "start_blkB" => [1,1,1], "end_blkB" => [ni,1,nk],
                    "offset" => [0.0, 0.0, 0.0], "angle" => 0.0))
            end

            blockId += 1
        end
    end
    updatedBndInfo = GroupBoundariesByName(blockBoundaries)
    updatedInterInfo = interfaces

    return blocks, updatedBndInfo, updatedInterInfo
end
```

with two helper functions `GetTouchingBoundaries` and `GroupBoundariesByName`. The former function will take in the block id and the starting and ending indices of the new subblock. 

Here we handle the interface info by looking to the right and above.

### Get Touching Boundaries
The function will check each of the edges of the child (newly created subblock) to see if they are on the boundary of their parent block (the original block being split). If the edge is on a boundary, inherit the boundary type from parent to child. We can check if the child is on the boundary of the parent by seeing if child start indices are equal to `1` or if the end indices are equal to size of the parent.

```julia
function GetTouchingBoundaries(blockInfo::Dict, bndInfo)
    # Parent extents
    N1, N2, N3 = blockInfo["parentDims"]

    # Child window in parent coordinates (inclusive)
    i0, j0 = blockInfo["start"]
    i1, j1 = blockInfo["end"]

    # Child-local sizes
    ni = i1 - i0 + 1
    nj = j1 - j0 + 1

    child_id = blockInfo["block"]
    out = Vector{Dict{String,Any}}()

    # LEFT side of child touches parent's LEFT boundary if i0 == 1
    if i0 == 1
        push!(out, Dict(
            "block" => child_id,
            "name"  => getBoundaryNameBySide(bndInfo; side=:left),
            "start" => [1, 1, 1],
            "end"   => [1, nj, 1],
        ))
    end

    # RIGHT side of child touches parent's RIGHT boundary if i1 == N1
    if i1 == N1
        push!(out, Dict(
            "block" => child_id,
            "name"  => getBoundaryNameBySide(bndInfo; side=:right),
            "start" => [ni, 1, 1],
            "end"   => [ni, nj, 1],
        ))
    end

    # BOTTOM side of child touches parent's BOTTOM boundary if j0 == 1
    if j0 == 1
        push!(out, Dict(
            "block" => child_id,
            "name"  => getBoundaryNameBySide(bndInfo; side=:bottom),
            "start" => [1, 1, 1],
            "end"   => [ni, 1, 1],
        ))
    end

    # TOP side of child touches parent's TOP boundary if j1 == N2
    if j1 == N2
        push!(out, Dict(
            "block" => child_id,
            "name"  => getBoundaryNameBySide(bndInfo; side=:top),
            "start" => [1, nj, 1],
            "end"   => [ni, nj, 1],
        ))
    end

    return out
end
```

with the helper function `GetBoundaryNameBySide(bndInfo; side=:side)` where we take advantage of the fact that bndInfo will only have information about the single inputted block. This function we can write as

```julia
function GetBoundaryNameBySide(bndInfo; side=:none)
    




end
```

### Group Boundaries by Name
Algorithm reformulates the list of boundaries into the Tortuga format, a list of dicts of dicts, where the first dict is the type of boundary and the sub dict are the blocks and edges.

```julia
function GroupBoundariesByName(faceList)
    boundaryGroups = Dict()

    for face in faceList
        name = face["name"]
        
        # Copy face and remove redundant name
        faceCopy = copy(face)
        delete!(faceCopy, "name")
        
        if !haskey(boundaryGroups, name)
            boundaryGroups[name] = []
        end
        push!(boundaryGroups[name], faceCopy)
    end

    groupedInfo = []
    for (name, faces) in boundaryGroups
        if !isempty(faces)
            push!(groupedInfo, Dict("name" => name, "faces" => faces))
        end
    end

    return groupedInfo
end
```


## Solve All Blocks Algorithm

