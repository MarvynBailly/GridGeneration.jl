function visualizer_interface(mesh, interfaces)
    # Visualize the mesh and interfaces
    for interface in interfaces
        blkA = interface["blockA"]
        blkB = interface["blockB"]
        start_blkA = interface["start_blkA"]
        end_blkA = interface["end_blkA"]
        start_blkB = interface["start_blkB"]
        end_blkB = interface["end_blkB"]

        # Extract the blocks from the mesh
        blockA = mesh[blkA + 1]
        blockB = mesh[blkB + 1]        

        println("Block A: $blkA, Block B: $blkB")
        println("Size of Block A: $(size(blockA[1]))")
        println("Size of Block B: $(size(blockB[1]))")
        println("Start Block A: $start_blkA, End Block A: $end_blkA")
        println("Start Block B: $start_blkB, End Block B: $end_blkB")

        X1, Y1, = blockA[1], blockA[2]
        X1_start = X1[start_blkA[1] + 1, start_blkA[2] + 1]
        Y1_start = Y1[start_blkA[1] + 1, start_blkA[2] + 1]
        X1_end = X1[end_blkA[1] + 1, end_blkA[2] + 1]
        Y1_end = Y1[end_blkA[1] + 1, end_blkA[2] + 1]

        X2, Y2, = blockB[1], blockB[2]
        X2_start = X2[start_blkB[1] + 1, start_blkB[2] + 1]
        Y2_start = Y2[start_blkB[1] + 1, start_blkB[2] + 1]
        X2_end = X2[end_blkB[1] + 1, end_blkB[2] + 1]
        Y2_end = Y2[end_blkB[1] + 1, end_blkB[2] + 1]

 
        l = @layout [a b]
        p = plot(layout = l, size=(1000, 500))#, legend=false)

        # Plot A: Full mesh with highlighted blocks
        for (i, blk) in enumerate(mesh)
            X, Y = blk[1], blk[2]
            clr = (i == blkA + 1 || i == blkB + 1) ? :red : :gray
            plotGrid!(p[1], X, Y, clr)
        end
        title!(p[1], "Full Mesh with Highlighted Blocks")
        # plot!(p[1], legend=false)

        # Plot B: Zoom-in on interface
        clr = :blue
        plotGrid!(p[2], X1, Y1, clr)
        plotGrid!(p[2], X2, Y2, clr)

        scatter!(p[2],[X1_start], [Y1_start], color=:red, label="Block A Interface start", markersize=8)
        scatter!(p[2],[X1_end], [Y1_end], color=:green, label="Block A Interface end", markersize=8)
        scatter!(p[2], [X2_start], [Y2_start], color=:red, label="Block B Interface start", markersize=5)
        scatter!(p[2], [X2_end], [Y2_end], color=:green, label="Block B Interface end", markersize=5)
        title!(p[2], "Interface Between Block A and B")

        display(p)
        readline()
    end    
end

function visualizer_boundaries(mesh, boundaries)
    for boundary in boundaries
        name = boundary["name"]
        faces = boundary["faces"]

        println("Boundary Name: $name")
        println("Number of Faces: $(length(faces))")

        for face in faces
            blk = face["block"]
            start = face["start"]
            end_ = face["end"]
            println("Face: Block $blk, Start $start, End $end_")

            block = mesh[blk + 1]
            X1, Y1 = block[1], block[2]
            X1_start = X1[start[1] + 1, start[2] + 1]
            Y1_start = Y1[start[1] + 1, start[2] + 1]
            X1_end = X1[end_[1] + 1, end_[2] + 1]
            Y1_end = Y1[end_[1] + 1, end_[2] + 1]

            # Layout with two plots
            l = @layout [a b]
            p = plot(layout = l, size=(1000, 500))

            # Plot A: Full mesh with highlighted block
            for (i, blkmesh) in enumerate(mesh)
                X, Y = blkmesh[1], blkmesh[2]
                color = (i == blk + 1) ? :red : :gray
                plotGrid!(p[1], X, Y, color)
            end
            title!(p[1], "Full Mesh with Highlighted Block")

            # Plot B: Zoomed-in block with face highlighted
            color = :blue
            plotGrid!(p[2], X1, Y1, color)
            scatter!(p[2], [X1_start], [Y1_start], color=:red, label="Face Start", markersize=8)
            scatter!(p[2], [X1_end], [Y1_end], color=:green, label="Face End", markersize=8)
            title!(p[2], "Boundary Face on Block $blk")

            display(p)
            readline()
        end
    end
end


function visualize_grid(mesh, interfaces, boundaries)
    p = plot(legend=false, size=(1000, 800), aspect_ratio=:equal)

    # Plot each block's grid in gray
    for (i, blk) in enumerate(mesh)
        X, Y, Z = blk[1], blk[2], blk[3]
        plotGrid!(p, X[:,:,1], Y[:,:,1], :gray)

        # Compute block center and plot block number
        center_x = mean(X[:,:,1])
        center_y = mean(Y[:,:,1])
        annotate!(p, center_x, center_y, text(string(i-1), :blue, 16)) # i-1 because you use 0-based indexing
    end

    # Plot interfaces (red lines)
    for iface in interfaces
        blkA = iface["blockA"]
        start_blkA = iface["start_blkA"]
        end_blkA = iface["end_blkA"]

        block = mesh[blkA+1]
        X, Y = block[1][:,:,1], block[2][:,:,1]

        x_start = X[start_blkA[1]+1, start_blkA[2]+1]
        y_start = Y[start_blkA[1]+1, start_blkA[2]+1]
        x_end   = X[end_blkA[1]+1,   end_blkA[2]+1]
        y_end   = Y[end_blkA[1]+1,   end_blkA[2]+1]

        plot!(p, [x_start, x_end], [y_start, y_end], color=:red, linewidth=2)
    end

    # Plot boundaries with labels
    for bnd in boundaries
        name = bnd["name"]
        faces = bnd["faces"]

        for face in faces
            blk = face["block"]
            start = face["start"]
            end_ = face["end"]

            block = mesh[blk+1]
            X, Y = block[1][:,:,1], block[2][:,:,1]

            x_start = X[start[1]+1, start[2]+1]
            y_start = Y[start[1]+1, start[2]+1]
            x_end   = X[end_[1]+1, end_[2]+1]
            y_end   = Y[end_[1]+1, end_[2]+1]

            # Choose color based on boundary name
            bcolor = :black
            if occursin("Inflow", name)
                bcolor = :green
            elseif occursin("Outflow", name)
                bcolor = :orange
            elseif occursin("Wall", name)
                bcolor = :blue
            end

            plot!(p, [x_start, x_end], [y_start, y_end], color=bcolor, linewidth=2, linestyle=:dash)

            # Put label near the middle of the boundary line
            mid_x = (x_start + x_end) / 2
            mid_y = (y_start + y_end) / 2
            annotate!(p, mid_x, mid_y, text(name, bcolor, 14))
        end
    end

    # display(p)
    return p
end

function visualize_grid_3D(mesh, interfaces, boundaries)
    p = plot(legend=true, size=(1000, 800), aspect_ratio=:equal, camera=(30, 10))

    # Plot each block's grid (gray surface)
    for (i, blk) in enumerate(mesh)
        X, Y, Z = blk[1], blk[2], blk[3]

        for k in 1:size(Z, 3)
            plotGrid3D!(p, X[:,:,k], Y[:,:,k], Z[:,:,k], :gray)

            # center_x = mean(X[:,:,k])
            # center_y = mean(Y[:,:,k])
            # z = Z[1, 1, k] 
            # annotate!(p, center_x, center_y, z, text(string(i-1), :blue, 16))
        end
    end

    # Plot interfaces (red lines)
    for iface in interfaces
        blkA = iface["blockA"]
        start_blkA = iface["start_blkA"]
        end_blkA = iface["end_blkA"]

        block = mesh[blkA+1]
        X, Y = block[1][:,:,1], block[2][:,:,1]

        x_start = X[start_blkA[1]+1, start_blkA[2]+1]
        y_start = Y[start_blkA[1]+1, start_blkA[2]+1]
        x_end   = X[end_blkA[1]+1,   end_blkA[2]+1]
        y_end   = Y[end_blkA[1]+1,   end_blkA[2]+1]

        for k in 1:size(block[3], 3)
            z = block[3][start_blkA[1]+1, start_blkA[2]+1, k]
            plot!(p, [x_start, x_end], [y_start, y_end], [z,z], color=:red, linewidth=2, label=false)
        end
    end


    # Plot interfaces (red lines)
    for (i, iface) in enumerate(interfaces)
        blkA = iface["blockA"]
        start_blkA = iface["start_blkA"]
        end_blkA = iface["end_blkA"]
        block = mesh[blkA+1]
        X, Y, Z = block[1], block[2], block[3]

        x_start = X[start_blkA[1]+1, start_blkA[2]+1, start_blkA[3]+1]
        y_start = Y[start_blkA[1]+1, start_blkA[2]+1, start_blkA[3]+1]
        z_start = Z[start_blkA[1]+1, start_blkA[2]+1, start_blkA[3]+1]

        x_end = X[end_blkA[1]+1, end_blkA[2]+1, end_blkA[3]+1]
        y_end = Y[end_blkA[1]+1, end_blkA[2]+1, end_blkA[3]+1]
        z_end = Z[end_blkA[1]+1, end_blkA[2]+1, end_blkA[3]+1]

        plot!(p, [x_start, x_end], [y_start, y_end], [z_start, z_end],
              color=:red, linewidth=2, label=(i == 1 ? "Interface" : false))
    end

    # Plot boundaries with labels
    for bnd in boundaries
        name = bnd["name"]
        faces = bnd["faces"]

        for face in faces
            blk = face["block"]
            start = face["start"]
            end_ = face["end"]

            block = mesh[blk+1]
            X, Y = block[1][:,:,1], block[2][:,:,1]

            x_start = X[start[1]+1, start[2]+1]
            y_start = Y[start[1]+1, start[2]+1]
            x_end   = X[end_[1]+1, end_[2]+1]
            y_end   = Y[end_[1]+1, end_[2]+1]

            # Choose color based on boundary name
            bcolor = :black
            if occursin("Inflow", name)
                bcolor = :green
            elseif occursin("Outflow", name)
                bcolor = :orange
            elseif occursin("Wall", name)
                bcolor = :blue
            end

            for k in 1:size(block[3], 3)
                z = block[3][start[1]+1, start[2]+1, k]
                plot!(p, [x_start, x_end], [y_start, y_end], [z,z], color=bcolor, linewidth=2, linestyle=:dash, label=false)
            end
        end
    end


    # Plot boundaries (colored dashed lines)
    has_label = Dict("Inflow" => false, "Outflow" => false, "Wall" => false)

    for bnd in boundaries
        name = bnd["name"]
        faces = bnd["faces"]

        for face in faces
            blk = face["block"]
            start = face["start"]
            end_ = face["end"]
            block = mesh[blk+1]
            X, Y, Z = block[1], block[2], block[3]

            x_start = X[start[1]+1, start[2]+1, start[3]+1]
            y_start = Y[start[1]+1, start[2]+1, start[3]+1]
            z_start = Z[start[1]+1, start[2]+1, start[3]+1]

            x_end = X[end_[1]+1, end_[2]+1, end_[3]+1]
            y_end = Y[end_[1]+1, end_[2]+1, end_[3]+1]
            z_end = Z[end_[1]+1, end_[2]+1, end_[3]+1]

            bcolor = :black
            bname = ""
            if occursin("Inflow", name)
                bcolor = :green
                bname = "Inflow"
            elseif occursin("Outflow", name)
                bcolor = :orange
                bname = "Outflow"
            elseif occursin("Wall", name)
                bcolor = :blue
                bname = "Wall"
            end

            plot!(p, [x_start, x_end], [y_start, y_end], [z_start, z_end],
                  color=bcolor, linewidth=2, linestyle=:dash,
                  label=(has_label[bname] ? false : bname))
            has_label[bname] = true
        end
    end

    # display(p)
    return p
end