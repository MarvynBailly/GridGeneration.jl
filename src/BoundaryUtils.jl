function UpdateBndInfo!(bndInfo, blocks; verbose=false)
    # Description: Update the boundary and interface information with the new sizes of the blocks
    # Input: bndInfo, interInfo, blocks
    # Output: updated bndInfo, interInfo
    
    
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

            # update the start and end positions
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