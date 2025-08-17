# Single Block with Splitting
 

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

### Get Touching Boundaries - Idea 1
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


### Get Touching Boundaries - Idea 2 (Current)
Another option is to use this old existing function of mine:
```julia
function GetTouchingBoundaries(block::Dict, bndInfo)
    # Child window in parent (global) coordinates, inclusive
    block_i1, block_j1 = block["start"]
    block_i2, block_j2 = block["end"]
    blockId = block["block"]

    # Global -> child-local index maps
    to_local_i(i) = i - block_i1 + 1
    to_local_j(j) = j - block_j1 + 1

    touchingFaces = Vector{Dict{String,Any}}()

    for bnd in bndInfo
        name = bnd["name"]
        for face in bnd["faces"]
            faceStart = face["start"];  faceEnd = face["end"]
            i1, j1 = faceStart[1], faceStart[2]
            i2, j2 = faceEnd[1],   faceEnd[2]

            # Vertical face on parent's left/right boundary?
            if i1 == i2 && (i1 == block_i1 || i1 == block_i2)
                jlo = max(min(j1, j2), block_j1)
                jhi = min(max(j1, j2), block_j2)
                if jhi > jlo
                    push!(touchingFaces, Dict(
                        "name"  => name,
                        "block" => blockId,
                        "start" => [to_local_i(i1), to_local_j(jlo), 1],
                        "end"   => [to_local_i(i2), to_local_j(jhi), 1],
                    ))
                end

            # Horizontal face on parent's bottom/top boundary?
            elseif j1 == j2 && (j1 == block_j1 || j1 == block_j2)
                ilo = max(min(i1, i2), block_i1)
                ihi = min(max(i1, i2), block_i2)
                if ihi > ilo
                    push!(touchingFaces, Dict(
                        "name"  => name,
                        "block" => blockId,
                        "start" => [to_local_i(ilo), to_local_j(j1), 1],
                        "end"   => [to_local_i(ihi), to_local_j(j2), 1],
                    ))
                end
            end
        end
    end

    return touchingFaces
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

### Split Examples

Inputting the following splits

```julia 
splitLocations = [
    [ 300, 500 ],   # split along the x axis
    [ 40 ]      # split along the y axis
]
```

yields:

![test](../../assets/images/SingleBlock/split_example.svg)

## Solve All Blocks Algorithm

With the single block splitting working, let's move to creating an algorithm to solve the ODE across all the edges of the blocks. The challenging part here comes from the required consistency across the edges of the blocks. Consider blocks $B_1$ with left/right edges $e$ and $f$ and neighboring block $B_2$ with left/right edges $f$ and $g$. Now if we solve the ODE along $e,f,g$ and find the optimal number of points (let's denote this via $N_e, N_f,$ and $N_g$) are such that $N_e > N_f, N_g > N_f$ with $N_e \neq N_g$, then we will have an issue solving along $e$ with the optimal number of points. To fix this issue, we can take the global optimal number $N_\text{opt}$ of the pairs, that is $N_\text{opt} = \max(N_e, N_f, N_g)$ and solve the $e,f,$ and $g$ edges with $N_\text{opt}$.

Let's have an array `blockInstructions` of size of blocks that will contain the max number of points in each direction. Now we loop through the blocks and do the following:
- for each dir
  - collect the neighboring block IDs using the interface info and save in array
  - Solve for the number of points for in the dir on both edges
  - Send the number of points to each neighbor with the following logic:
    - if incoming is greater than saved, overwrite with incoming, else keep the saved value. We can do this by looking at the correct spot in the `blockInstructions` using the neighboring block Id, and in the correct direction using the current dir. 

Once we loop through all blocks, each block should now contain the max number of points to use in each direction. We can finally loop through all the blocks and solve the ODE with the correct number of points.

### Algorithm

```julia
function SolveAllBlocks(metric, blocks, bndInfo, interInfo)
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
                optN = GridGeneration.GetOptNEdgePair(left, right, metric)
            else  # vertical
                bottom   = block[:, :, 1]
                top  = block[:, :, end]
                optN = GridGeneration.GetOptNEdgePair(bottom, top, metric)
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
    p1 = plot()
    for (blockId, block) in enumerate(blocks)
        optNs = blockDirOptN[blockId]

        computedBlock, bndInfo, interInfo = GridGeneration.SolveBlockFixedN(block, bndInfo, interInfo, metric, optNs)

        computedBlocks[blockId] = computedBlock
    end

    # update the boundary information and interface information 
    GridGeneration.UpdateBndInfo!(bndInfo, computedBlocks; verbose=false)
    updatedInterInfo = GridGeneration.UpdateInterInfo(interInfo, computedBlocks; verbose=false)

    return computedBlocks, bndInfo, updatedInterInfo
end

```


We can find the neighboring cells in a direction using this recursive function:

```julia
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
```


### Examples
Putting this all together and using our custom metric creator, we get the following results:

