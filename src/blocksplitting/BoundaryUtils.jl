"""
Update boundary information with new block sizes.
"""
function UpdateBndInfo!(bndInfo, blocks; verbose=false)
    for bc in bndInfo
        if verbose
            println("Updating boundary condition: ", bc["name"])
            println("Size of faces: ", length(bc["faces"]))
        end

        for face in bc["faces"]
            blockId = face["block"]
            block = blocks[blockId]
            nx = size(block, 2)
            ny = size(block, 3)
            nk = 1
            
            if verbose
                println("Block size: ", size(block))
                println("Block ID: ", blockId)
                println("Face start: ", face["start"], " end: ", face["end"])
            end

            # Update start and end positions
            start = face["start"]
            end_ = face["end"]

            if start[1] != 1 
                start[1] = nx
            end
            if start[2] != 1 
                start[2] = ny
            end
            if end_[1] != 1 
                end_[1] = nx
            end
            if end_[2] != 1 
                end_[2] = ny
            end
        end
    end
    if verbose
        println("Updated boundary conditions: ", bndInfo)
    end
end



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
        for face in bnd["faceInfo"]
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