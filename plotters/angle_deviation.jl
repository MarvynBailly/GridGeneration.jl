using Plots


function PlotGridAngleDeviation(
    blocks,
    angleDeviations;
    cmax = nothing,
    title::String = "Angle Deviation PDF",
    colorbar_title::String = "Deviation (\$^{\\circ}\$)",
    color_scheme = :thermal, # Default to a built-in scheme for flexibility
    grid_color = :gray,
    grid_alpha::Real = 0.6,
    grid_linewidth::Real = 0.05,
    plt = nothing
)
    p = plt !== nothing ? p = plt : p = plot()

    if cmax === nothing
        cmax = maximum([maximum(dev) for dev in angleDeviations])
    end

    # --- Initialize Plot and Color Map ---
    plot!(p,
        aspect_ratio = :equal,
        title = title,
        xlabel = "X",
        ylabel = "Y",
        background_color = :white, # Set background for better grid line visibility,
        label = ""
    )

    # Use the user-defined color scheme
    cmap = cgrad(color_scheme)

    # --- Prepare arrays for shapes and colors ---
    shapes = Plots.Shape[]
    colors = [] 

    # --- Loop to GATHER cell data (Optimized for vectorization) ---
    for block_index in 1:length(blocks)
        # block[1, :, :] is X coordinates, block[2, :, :] is Y coordinates
        x = blocks[block_index][1, :, :]
        y = blocks[block_index][2, :, :]
        data = angleDeviations[block_index]
        nrows, ncols = size(x)

        # Use @inbounds for performance boost if array bounds are guaranteed to be safe
        @inbounds for j in 1:(ncols-1)
            @inbounds for i in 1:(nrows-1)
                # Create the quadrilateral shape coordinates
                cell_x = [x[i,j], x[i+1,j], x[i+1,j+1], x[i,j+1]]
                cell_y = [y[i,j], y[i+1,j], y[i+1,j+1], y[i,j+1]]
                push!(shapes, Plots.Shape(cell_x, cell_y))

                # Normalize and get color
                cell_value = data[i, j]
                normalized_val = clamp(cell_value / cmax, 0, 1)
                push!(colors, get(cmap, normalized_val))
            end
        end
    end

    # --- 2. Plot all shapes in a single, highly efficient call ---
    # The reshape(colors, 1, :) is key to mapping the color vector to the shape vector
    plot!(p, shapes,
        color = reshape(colors, 1, :),
        # linecolor = :transparent, # Removes borders between cells
        label = "",
        lw = grid_linewidth
    )

    # --- 3. Overlay the grid mesh as a separate, vectorized layer ---
    # This is done using the new customizable grid parameters
    for block in blocks
        x = block[1, :, :]
        y = block[2, :, :]
        
        # Plot lines along x-indices (vertical grid lines)
        plot!(p, x, y, legend=false, color=grid_color, lw=grid_linewidth, alpha=grid_alpha)
        
        # Plot lines along y-indices (horizontal grid lines)
        plot!(p, x', y', legend=false, color=grid_color, lw=grid_linewidth, alpha=grid_alpha)
    end

    # --- 4. Add the colorbar using a scatter trick ---
    # Use user-defined colorbar_title and the color_scheme
    scatter!(p, [Inf], [Inf], # Plots a single point off-screen to only show colorbar
             zcolor=[cmax], clims=(0,cmax),
             color=cmap, markeralpha=0, label="",
             colorbar_title=colorbar_title,
             colorbar_titlefontsize=10)

    return p
end
