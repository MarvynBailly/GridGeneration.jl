using Plots

function PlotGridAngleDeviation(blocks, angleDeviations)
    blockMax = [maximum(dev) for dev in angleDeviations]
    cmax = maximum(blockMax)
    
    
    # --- Create the Base Plot ---
    p = plot(
        aspect_ratio = :equal,
        title = "Angle Deviation from 90°",
        xlabel = "X",
        ylabel = "Y"
    )

    # --- Create the custom White-to-Red color map ---
    cmap = cgrad([:white, :red])
    
    # --- Optimization: Prepare arrays for shapes and colors ---
    shapes = Shape[]
    colors = RGBA[]
    
    # --- Loop through each block to GATHER data (no plotting yet) ---
    for block_index in 1:length(blocks)
        x = blocks[block_index][1, :, :]
        y = blocks[block_index][2, :, :]
        data = angleDeviations[block_index]
        nrows, ncols = size(x)
        
        for i in 1:(nrows-1)
            for j in 1:(ncols-1)
                # Define the cell's geometry
                cell_x = [x[i,j], x[i+1,j], x[i+1,j+1], x[i,j+1]]
                cell_y = [y[i,j], y[i+1,j], y[i+1,j+1], y[i,j+1]]
                push!(shapes, Shape(cell_x, cell_y))
                
                # Determine the cell's color
                cell_value = data[i, j]
                normalized_val = cmax > 0 ? clamp(cell_value / cmax, 0, 1) : 0
                push!(colors, get(cmap, normalized_val))
            end
        end
    end

    # --- Optimization: Plot all shapes in a SINGLE, fast command ---
    plot!(p, shapes,
        color = reshape(colors, 1, :), # Reshape color array for broadcasting
        linecolor = :black,
        lw = 0.05, 
        label = ""
    )

    # --- Add the colorbar once, after all shapes are plotted ---
    scatter!(p, [], [], zcolor=0:cmax, clims=(0,cmax),
             color=cmap, markeralpha=0, label="",
             colorbar_title="Deviation (°)")

    return p
end
