function convert_2D_to_3D(mesh2D, extrusion_length=0.1, n_layers=2)
    mesh3D = similar(mesh2D)

    for block2D in mesh2D
        # X2D, Y2D = block2D
        X2D, Y2D = block2D[1,:,:], block2D[2,:,:]
        ni, nj = size(X2D)
        nk = n_layers

        # Initialize 3D arrays
        X3D = zeros(ni, nj, nk)
        Y3D = zeros(ni, nj, nk)
        Z3D = zeros(ni, nj, nk)

        # loop over and extrude
        for k in 1:nk
            X3D[:,:,k] .= X2D
            Y3D[:,:,k] .= Y2D
            Z3D[:,:,k] .= (k-1)*extrusion_length
        end

        # Push the block into the 3D mesh
        # push!(mesh3D, [X3D, Y3D, Z3D])
        # save as [dim,x,y,z]
        push!(mesh3D, [X3D, Y3D, Z3D])  
    end
    return mesh3D
end

function generate_interfaces(mesh, nblocks_horizontal)
    interfaces = []

    num_blocks = length(mesh)

    # p = plot(legend=false)
    for blk_idx in 1:num_blocks
        X,Y,Z = mesh[blk_idx]
        ni, nj, nk = size(X).-1
        blk_num = blk_idx - 1 # zero-based indexing 

        # println("Block number: ", blk_num, " with size: ", ni, nj, nk)

        col = (blk_idx - 1) % nblocks_horizontal + 1  # column position
        row = (blk_idx - 1) ÷ nblocks_horizontal + 1  # row position

        # println("Row: ", row, " Column: ", col)

        # for k in 1:nk
            # look up 
            if col < nblocks_horizontal
                right_blk_num = blk_num + 1

                push!(interfaces, Dict(
                    "blockA"=>blk_num, "start_blkA"=>[0,nj,0], "end_blkA"=>[ni,nj,nk],
                    "blockB"=>right_blk_num, "start_blkB"=>[0,0,0], "end_blkB"=>[ni,0,nk],
                    "offset"=>[0.0,0.0,0.0], "angle"=>0.0))
            end
            # left
            if (blk_idx + nblocks_horizontal) <= num_blocks
                top_blk_num = blk_num + nblocks_horizontal
                push!(interfaces, Dict(
                    "blockA"=>blk_num, "start_blkA"=>[ni,0,0], "end_blkA"=>[ni,nj,nk],
                    "blockB"=>top_blk_num, "start_blkB"=>[0,0,0], "end_blkB"=>[0,nj,nk],
                    "offset"=>[0.0,0.0,0.0], "angle"=>0.0))
            end
        # end
    end

    return interfaces
end

function generate_boundaries(mesh, nblocks_horizontal)
    boundaries = Dict(
        "BCInflow" => [],
        "BCOutflow" => [],
        "BCWall" => [],
    )

    num_blocks = length(mesh)
    nblocks_vertical = num_blocks ÷ nblocks_horizontal

    # Add boundary condition across k=0
    for blk_idx in 1:4
        println("Block number: ", blk_idx)
        X, Y, Z = mesh[blk_idx]
        ni, nj, nk = size(X) .- 1
        blk_num = blk_idx - 1 
        push!(boundaries["BCWall"], Dict(
            "block" => blk_num,
            "start" => [0, 0, 0],
            "end"   => [ni, nj, 0]
        ))
    end

    # Add boundary condition across k=top
    for blk_idx in 1:4
        X, Y, Z = mesh[blk_idx]
        ni, nj, nk = size(X) .- 1
        blk_num = blk_idx - 1 
        push!(boundaries["BCWall"], Dict(
            "block" => blk_num,
            "start" => [0, 0, nk],
            "end"   => [ni, nj, nk]
        ))
    end


    # Explicitly identify blocks on right wall
    for blk_idx in (num_blocks - nblocks_horizontal + 1):num_blocks
        X, Y, Z = mesh[blk_idx]
        ni, nj, nk = size(X) .- 1
        blk_num = blk_idx - 1  # zero-based indexing

        # for k in 1:nk
            push!(boundaries["BCOutflow"], Dict(
                "block" => blk_num,
                "start" => [ni, 0, 0],
                "end"   => [ni, nj, nk]
            ))
        # end
    end

    # Explicitly identify blocks on left wall
    for blk_idx in 1:nblocks_horizontal
        X, Y, Z = mesh[blk_idx]
        ni, nj, nk = size(X) .- 1
        blk_num = blk_idx - 1  # zero-based indexing

        # for k in 1:nk
            push!(boundaries["BCInflow"], Dict(
                "block" => blk_num,
                "start" => [0, 0, 0],
                "end"   => [0, nj, nk]
            ))
        # end
    end

    # Explicitly identify blocks on bottom wall
    for blk_idx in 1:nblocks_horizontal:nblocks_horizontal*nblocks_vertical
        X, Y, Z = mesh[blk_idx]
        ni, nj, nk = size(X) .- 1
        blk_num = blk_idx - 1  # zero-based indexing

        # for k in 1:nk
            push!(boundaries["BCWall"], Dict(
                "block" => blk_num,
                "start" => [0, 0, 0],
                "end"   => [ni, 0, nk]
            ))
        # end
    end

    # Explicitly identify blocks on top wall
    for blk_idx in nblocks_horizontal:nblocks_horizontal:nblocks_horizontal*nblocks_vertical
        X, Y, Z = mesh[blk_idx]
        ni, nj, nk = size(X) .- 1
        blk_num = blk_idx - 1  # zero-based indexing

        # for k in 1:nk
            push!(boundaries["BCWall"], Dict(
                "block" => blk_num,
                "start" => [0, nj, 0],
                "end"   => [ni, nj, nk]
            ))
        # end
    end




    # --- Convert dictionary to list ---
    boundary_list = []
    for (name, faces) in boundaries
        if !isempty(faces)
            push!(boundary_list, Dict("name"=>name, "faces"=>faces))
        end
    end

    return boundary_list
end

function write_turtle_grid(mesh, interfaces, boundaries, filename)
    # Takes in mesh, a collection of blocks stored in (x,y,z) format, and writes the mesh to a .mesh file in the format of the turtle grid.
    
    open(filename, "w") do fid
        MAGICNUMBER = 751123
        VERSION = 3*100^2 + 1*100 + 47
        # VERSION =1

        IO_NAME_LENGTH = 200
        nbrBlocks = length(mesh)
        nbrInterfaces = length(interfaces)
        nbrBoundaries = length(boundaries)

        # write the header
        write(fid, Int32.([MAGICNUMBER, VERSION, nbrBlocks, nbrInterfaces, nbrBoundaries]))

        # size of each block
        for blk in mesh
            write(fid, collect(Int32, size(blk[1]) .- 1))
        end

        # Interfaces
        for iface in interfaces
            write(fid, Int32(iface["blockA"]), Int32.(iface["start_blkA"]), Int32.(iface["end_blkA"]))
            write(fid, Int32(iface["blockB"]), Int32.(iface["start_blkB"]), Int32.(iface["end_blkB"]))
            write(fid, Float64.(iface["offset"]), Float64(iface["angle"]))
        end

        # Write Boundaries
        for bnd in boundaries
            # Convert name to UTF-8 bytes
            name = bnd["name"]

            # PAD BEFORE CONVERTING
            if length(name) > IO_NAME_LENGTH
                name = name[1:IO_NAME_LENGTH]
            else
                name = rpad(name, IO_NAME_LENGTH)
            end
            name_bytes = collect(codeunits(name))

            # Ensure exact length of IO_NAME_LENGTH
            # name_bytes = collect(codeunits(name))
            # if length(name_bytes) > IO_NAME_LENGTH
            #     name_bytes = name_bytes[1:IO_NAME_LENGTH]
            # else
            #     name_bytes = vcat(name_bytes, zeros(UInt8, IO_NAME_LENGTH - length(name_bytes)))
            # end
            write(fid, name_bytes)
            write(fid, Int32(length(bnd["faces"])))
            for face in bnd["faces"]
                write(fid, Int32(face["block"]), Int32.(face["start"]), Int32.(face["end"]))
            end
        end
        # Node coordinates (k,j,i ordering, exactly like MATLAB)
        for d in 1:3
            for blk in mesh
                ni, nj, nk = size(blk[d])
                for i in 1:ni
                    for j in 1:nj
                        for k in 1:nk
                            write(fid, blk[d][i,j,k])
                        end
                    end
                end
            end
        end
        println("Grid written successfully to $filename")
    end
end

function create_turtle_grid(blocks, num_hor_blocks, filename; extrusion_length = 0.1, k_layers = 10, showPlots=false)

    mesh3d = convert_2D_to_3D(blocks, extrusion_length, k_layers)
    interfaces = generate_interfaces(mesh3d, num_hor_blocks)
    boundaries = generate_boundaries(mesh3d, num_hor_blocks)
    write_turtle_grid(mesh3d, interfaces, boundaries, filename)

    @info "Turtle grid written to $filename"

    # if showPlots
    #     visualizer_interface(blocks, interfaces)
    #     visualizer_boundaries(blocks, boundaries)
    #     p1 = visualize_grid(mesh3d, interfaces, boundaries)
    #     p2 = visualize_grid_3D(mesh3d, interfaces, boundaries)
    #     p = plot(p1, p2, layout = @layout [a b])
    #     display(p)
    # end
end