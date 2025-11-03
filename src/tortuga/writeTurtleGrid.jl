function convert_2D_to_3D(mesh2D, bndInfo2d, interfaceInfo2d, extrusion_length=0.1, n_layers=2)
    mesh3D = []

    # extrude each block
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

        # save as array [dim, X3D, Y3D, Z3D]
        block3d = zeros(3, ni, nj, nk)
        block3d[1, :, :, :] .= X3D
        block3d[2, :, :, :] .= Y3D
        block3d[3, :, :, :] .= Z3D
        push!(mesh3D, block3d)
    end

    # convert interface and bnd info to 0 based indexing for 3D
    # extrude bnd and interface info
    bndInfo3D = convert_bndInfo_to_0based_3D(bndInfo2d, n_layers)
    interfaceInfo3D = convert_interfaceInfo_to_0based_3D(interfaceInfo2d, n_layers)


    # add symmetric interfaces at top and bottom
    for (i, block) in enumerate(mesh3D)
        ni, nj, nk = size(block[1,:,:,:]) .- 1
        # Add interfaces at bottom (k=0) and top (k=nk)
        interface = Dict{String, Any}()
        interface["blockA"] = i - 1
        interface["blockB"] = i - 1
        interface["start_blkA"] = [0, 0, 0]
        interface["end_blkA"] = [ni, nj, 0]
        interface["start_blkB"] = [0, 0, nk]
        interface["end_blkB"] = [ni, nj, nk]
        # compute the offset from the top layer to the bottom layer
        offset = [0.0, 0.0, extrusion_length * (nk)]
        interface["offset"] = offset
        interface["angle"] = 0.0
        push!(interfaceInfo3D, interface)
    end

    return mesh3D, bndInfo3D, interfaceInfo3D
end

"""
    convert_bndInfo_to_0based_3D(bndInfo2d, n_layers)

Convert 2D boundary information from 1-based to 0-based indexing and add k-dimension.

# Arguments
- `bndInfo2d`: Boundary information with 1-based indexing (Julia format)
- `n_layers`: Number of layers in k-direction

# Returns
- `bndInfo3D`: Boundary information with 0-based indexing and k-dimension
"""
function convert_bndInfo_to_0based_3D(bndInfo2d, n_layers)
    bndInfo3D = deepcopy(bndInfo2d)
    
    for bnd in bndInfo3D
        for face in bnd["faces"]
            # Convert block ID from 1-based to 0-based
            face["block"] = face["block"] - 1
            
            # Convert start and end indices from 1-based to 0-based
            face["start"] = face["start"] .- 1
            face["end"] = face["end"] .- 1
            
            # Add k-dimension (from 0 to n_layers-1)
            face["start"] = [face["start"][1], face["start"][2], 0]
            face["end"] = [face["end"][1], face["end"][2], n_layers - 1]
        end
    end
    
    return bndInfo3D
end

"""
    convert_interfaceInfo_to_0based_3D(interfaceInfo2d, n_layers)

Convert 2D interface information from 1-based to 0-based indexing and add k-dimension.

# Arguments
- `interfaceInfo2d`: Interface information with 1-based indexing (Julia format)
- `n_layers`: Number of layers in k-direction

# Returns
- `interfaceInfo3D`: Interface information with 0-based indexing and k-dimension
"""
function convert_interfaceInfo_to_0based_3D(interfaceInfo2d, n_layers)
    interfaceInfo3D = deepcopy(interfaceInfo2d)
    
    for itf in interfaceInfo3D
        # Convert block IDs from 1-based to 0-based
        itf["blockA"] = itf["blockA"] - 1
        itf["blockB"] = itf["blockB"] - 1
        
        # Convert start and end indices from 1-based to 0-based
        itf["start_blkA"] = itf["start_blkA"] .- 1
        itf["end_blkA"] = itf["end_blkA"] .- 1
        itf["start_blkB"] = itf["start_blkB"] .- 1
        itf["end_blkB"] = itf["end_blkB"] .- 1
        
        # Add k-dimension (from 0 to n_layers-1)
        itf["start_blkA"] = [itf["start_blkA"][1], itf["start_blkA"][2], 0]
        itf["end_blkA"] = [itf["end_blkA"][1], itf["end_blkA"][2], n_layers - 1]
        itf["start_blkB"] = [itf["start_blkB"][1], itf["start_blkB"][2], 0]
        itf["end_blkB"] = [itf["end_blkB"][1], itf["end_blkB"][2], n_layers - 1]
    end
    
    return interfaceInfo3D
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
            write(fid, collect(Int32, size(blk)[2:4] .- 1))
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

            write(fid, name_bytes)
            write(fid, Int32(length(bnd["faces"])))
            for face in bnd["faces"]
                write(fid, Int32(face["block"]), Int32.(face["start"]), Int32.(face["end"]))
            end
        end
        # Node coordinates (k,j,i ordering for each dimension, matching Turtle format)
        for d in 1:3
            for blk in mesh
                ni, nj, nk = size(blk)[2:4]
                for i in 1:ni
                    for j in 1:nj
                        for k in 1:nk
                            write(fid, blk[d,i,j,k])
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