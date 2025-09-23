using Plots

# ---------- helpers ----------
# Ordered index range (inclusive)
@inline ordered_range(a::Int, b::Int) = a <= b ? (a:b) : (b:a)

# Extract a face polyline (xs,ys) from a block given start/end index triplets [i,j,k]
function extract_face_xy(block, start::AbstractVector{<:Integer}, endp::AbstractVector{<:Integer})
    i0, j0 = start[1], start[2]
    i1, j1 = endp[1],  endp[2]
    if i0 == i1
        jrange = ordered_range(j0, j1)
        xs = collect(block[1, i0, jrange])
        ys = collect(block[2, i0, jrange])
        return xs, ys
    elseif j0 == j1
        irange = ordered_range(i0, i1)
        xs = collect(block[1, irange, j0])
        ys = collect(block[2, irange, j0])
        return xs, ys
    else
        error("Face is not axis-aligned in index space: start=$start end=$endp")
    end
end

# Midpoint of a face (for label placement)
function face_midpoint(block, start::AbstractVector{<:Integer}, endp::AbstractVector{<:Integer})
    xs, ys = extract_face_xy(block, start, endp)
    m = cld(length(xs), 2)
    return xs[m], ys[m]
end

# A small discrete palette to cycle through boundary names
name_colors = [:red, :blue, :green, :magenta, :orange, :cyan, :purple,
                              :brown, :olive, :pink, :teal, :navy, :darkorange, :darkgreen]

# ---------- main ----------
"""
    plot_blocks_interfaces_boundaries(blocks, interfaces, bndInfo;
        grid_stride=4, grid_lw=0.6, boundary_lw=2.5, interface_lw=2.0,
        show_block_ids=true, legend=false, aspect=:equal, boundary_colors=nothing)

Draws all blocks, boundaries (colored by name), and interfaces.

- `blocks`      :: Vector of (2, ni, nj) arrays
- `interfaces`  :: Vector of Dicts with keys like "blockA","start_blkA","end_blkA",...
- `bndInfo`     :: Vector of Dicts, each with "name"::String, "faces"::Vector{Dict}

Keyword args:
- `grid_stride`  : plot every `grid_stride`th line for clarity
- `boundary_colors` : optional Dict{String,Any} mapping boundary name -> Plots color
"""
function plot_blocks_interfaces_boundaries(blocks, interfaces, bndInfo;
        grid_stride::Int=4, grid_lw=0.6, boundary_lw=2.5, interface_lw=2.0,
        show_block_ids::Bool=true, legend::Bool=false, aspect=:equal,
        boundary_colors=nothing, titleName=nothing)

    plt = plot(; legend=legend, aspect_ratio=aspect,title=titleName)

    # 1) plot block grids (light)
    for (bid, blk) in enumerate(blocks)
        ni, nj = size(blk, 2), size(blk, 3)
        # ξ-lines
        for j in 1:grid_stride:nj
            plot!(plt, blk[1, :, j], blk[2, :, j], color=:gray, alpha=0.5, lw=grid_lw, label=false)
        end
        # η-lines
        for i in 1:grid_stride:ni
            plot!(plt, blk[1, i, :], blk[2, i, :], color=:gray, alpha=0.5, lw=grid_lw, label=false)
        end
        if show_block_ids
            ci = max(1, round(Int, ni/2)); cj = max(1, round(Int, nj/2))
            annotate!(plt, blk[1, ci, cj], blk[2, ci, cj], text("B$bid", 8, :black))
        end
    end

    # 2) boundaries colored by name
    if bndInfo !== nothing && !isempty(bndInfo)
        names = [bc["name"] for bc in bndInfo]
        uniq  = unique(names)
        colmap = boundary_colors === nothing ?
                 Dict(name => name_colors[1 + mod(i-1, length(name_colors))]
                      for (i,name) in enumerate(uniq)) :
                 boundary_colors

        for bc in bndInfo
            name = bc["name"]
            col  = get(colmap, name, :black)
            placed = false
            for face in bc["faces"]
                blkid = face["block"]
                xs, ys = extract_face_xy(blocks[blkid], face["start"], face["end"])
                plot!(plt, xs, ys; color=col, lw=boundary_lw, label=false)
                if !placed
                    mx, my = face_midpoint(blocks[blkid], face["start"], face["end"])
                    annotate!(plt, mx, my, text(name, 8, col))
                    placed = true
                end
            end
        end
    end

    # 3) interfaces (dashed)
    if interfaces !== nothing && !isempty(interfaces)
        for (idx, itf) in enumerate(interfaces)
            # Draw only once (blockA side). Adjust if you prefer both.
            blkA = itf["blockA"]
            xs, ys = extract_face_xy(blocks[blkA], itf["start_blkA"], itf["end_blkA"])
            plot!(plt, xs, ys; color=:red, lw=interface_lw, ls=:dash, alpha=0.9, label=false)
            scatter!(plt, [xs[1], xs[end]], [ys[1], ys[end]], color=:red, ms=3, label=false)
            mx, my = face_midpoint(blocks[blkA], itf["start_blkA"], itf["end_blkA"])
            annotate!(plt, mx, my, text("I$idx", 8, :black))
        end
    end

    xlabel!(plt, "x"); ylabel!(plt, "y")
    return plt
end
