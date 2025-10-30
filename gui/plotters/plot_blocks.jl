"""
Plotting utilities for the Grid Generation GUI.
"""

# ---------- Helper Functions ----------
# Ordered index range (inclusive)
@inline ordered_range(a::Int, b::Int) = a <= b ? (a:b) : (b:a)

"""
    extract_face_xy(block, start, endp)

Extract a face polyline (xs, ys) from a block given start/end index triplets [i,j,k].
The face must be axis-aligned in index space.
Returns empty arrays if indices are out of bounds.
"""
function extract_face_xy(block, start::AbstractVector{<:Integer}, endp::AbstractVector{<:Integer})
    i0, j0 = start[1], start[2]
    i1, j1 = endp[1],  endp[2]
    
    # Get block dimensions
    ni, nj = size(block, 2), size(block, 3)
    
    if i0 == i1
        # Constant i (vertical face)
        # Check bounds and clamp indices to valid range
        if i0 < 1 || i0 > ni
            @warn "Index i=$i0 out of bounds [1, $ni], clamping"
            i_clamped = clamp(i0, 1, ni)
        else
            i_clamped = i0
        end
        
        j0_clamped = clamp(j0, 1, nj)
        j1_clamped = clamp(j1, 1, nj)
        
        jrange = ordered_range(j0_clamped, j1_clamped)
        xs = collect(block[1, i_clamped, jrange])
        ys = collect(block[2, i_clamped, jrange])
        return xs, ys
    elseif j0 == j1
        # Constant j (horizontal face)
        # Check bounds and clamp indices to valid range
        if j0 < 1 || j0 > nj
            @warn "Index j=$j0 out of bounds [1, $nj], clamping"
            j_clamped = clamp(j0, 1, nj)
        else
            j_clamped = j0
        end
        
        i0_clamped = clamp(i0, 1, ni)
        i1_clamped = clamp(i1, 1, ni)
        
        irange = ordered_range(i0_clamped, i1_clamped)
        xs = collect(block[1, irange, j_clamped])
        ys = collect(block[2, irange, j_clamped])
        return xs, ys
    else
        # Not axis-aligned - return empty arrays
        @warn "Face is not axis-aligned in index space: start=$start end=$endp"
        return Float64[], Float64[]
    end
end

"""
    plot_blocks_with_highlighting!(ax, blocks_obs, highlight_obs, bndInfo_obs, interfaceInfo_obs; grid_stride=4)

Helper function to create a reactive plot of grid blocks with boundary and interface highlighting.

# Arguments
- `ax`: The axis to plot on
- `blocks_obs`: Observable containing vector of blocks
- `highlight_obs`: Observable{Bool} to toggle boundary highlighting
- `bndInfo_obs`: Observable containing boundary info
- `interfaceInfo_obs`: Observable containing interface info
- `grid_stride`: Plot every Nth grid line (default=4)

This function visualizes:
- Grid lines in light gray (strided for clarity)
- Interfaces as dashed red lines
- Boundaries color-coded by name when highlighting is enabled
"""
function plot_blocks_with_highlighting!(ax::Axis, blocks_obs, highlight_obs, 
                                       bndInfo_obs, interfaceInfo_obs; 
                                       grid_stride::Int=1)
    
    # Observables for plot data
    grid_segments_obs = Observable(Point2f[])
    boundary_segments_obs = Observable(Point2f[])
    boundary_colors_obs = Observable(RGBAf[])
    interface_segments_obs = Observable(Point2f[])
    
    # Color palette for boundary names
    name_colors = [
        RGBAf(1,0,0,1),      # red
        RGBAf(0,0,1,1),      # blue
        RGBAf(0,1,0,1),      # green
        RGBAf(1,0,1,1),      # magenta
        RGBAf(1,0.647,0,1),  # orange
        RGBAf(0,1,1,1),      # cyan
        RGBAf(0.502,0,0.502,1), # purple
    ]
    
    onany(blocks_obs, highlight_obs, bndInfo_obs, interfaceInfo_obs) do blocks, highlight_on, bndInfo, interfaceInfo
        # Reset all segments
        grid_points = Point2f[]
        boundary_points = Point2f[]
        boundary_colors = RGBAf[]
        interface_points = Point2f[]
        
        !isnothing(blocks) && isa(blocks, Vector) || return
        
        # 1) Plot block grids (light gray, strided)
        for block in blocks
            ni, nj = size(block, 2), size(block, 3)
            
            # ξ-lines (constant j)
            for j in 1:grid_stride:nj
                for i in 1:(ni-1)
                    push!(grid_points, Point2f(block[1, i, j], block[2, i, j]))
                    push!(grid_points, Point2f(block[1, i+1, j], block[2, i+1, j]))
                end
            end
            
            # η-lines (constant i)
            for i in 1:grid_stride:ni
                for j in 1:(nj-1)
                    push!(grid_points, Point2f(block[1, i, j], block[2, i, j]))
                    push!(grid_points, Point2f(block[1, i, j+1], block[2, i, j+1]))
                end
            end
        end
        
        # 2) Plot boundaries (colored by name if highlighting is on)
        if highlight_on && bndInfo !== nothing && !isempty(bndInfo)
            names = [bc["name"] for bc in bndInfo]
            uniq = unique(names)
            colmap = Dict(name => name_colors[1 + mod(i-1, length(name_colors))] 
                         for (i, name) in enumerate(uniq))
            
            for bc in bndInfo
                name = bc["name"]
                col = get(colmap, name, RGBAf(0,0,0,1))
                
                for face in bc["faces"]
                    blkid = face["block"]
                    if blkid <= length(blocks)
                        try
                            xs, ys = extract_face_xy(blocks[blkid], face["start"], face["end"])
                            # Skip if extraction failed (empty arrays)
                            if !isempty(xs) && !isempty(ys)
                                for i in 1:(length(xs)-1)
                                    push!(boundary_points, Point2f(xs[i], ys[i]))
                                    push!(boundary_points, Point2f(xs[i+1], ys[i+1]))
                                    push!(boundary_colors, col)
                                end
                            end
                        catch e
                            @warn "Failed to extract boundary face for block $blkid"#: $e"
                        end
                    end
                end
            end
        end
        
        # 3) Plot interfaces (dashed red lines)
        if interfaceInfo !== nothing && !isempty(interfaceInfo)
            for itf in interfaceInfo
                blkA = itf["blockA"]
                if blkA <= length(blocks)
                    try
                        xs, ys = extract_face_xy(blocks[blkA], itf["start_blkA"], itf["end_blkA"])
                        # Skip if extraction failed (empty arrays)
                        if !isempty(xs) && !isempty(ys)
                            for i in 1:(length(xs)-1)
                                push!(interface_points, Point2f(xs[i], ys[i]))
                                push!(interface_points, Point2f(xs[i+1], ys[i+1]))
                            end
                        end
                    catch e
                        @warn "Failed to extract interface for block $blkA"#: $e"
                    end
                end
            end
        end
        
        grid_segments_obs[] = grid_points
        boundary_segments_obs[] = boundary_points
        boundary_colors_obs[] = boundary_colors
        interface_segments_obs[] = interface_points
        autolimits!(ax)
    end
    
    # Plot layers (order matters for visibility)
    linesegments!(ax, grid_segments_obs, color=RGBAf(0.5, 0.5, 0.5, 0.5), linewidth=0.6)
    linesegments!(ax, interface_segments_obs, color=RGBAf(1,0,0,0.9), linewidth=2.0, linestyle=:dash)
    linesegments!(ax, boundary_segments_obs, color=boundary_colors_obs, linewidth=2.5)
    
    notify(blocks_obs) # Trigger initial plot
end
