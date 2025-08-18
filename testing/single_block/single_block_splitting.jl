include("../../src/GridGeneration.jl")

include("airfoil/GetAirfoilGrid.jl")
include("airfoil/GetBoundary.jl")
include("airfoil/metric/Metric.jl")
include("airfoil/metric/CustomMetric.jl")
include("../../plotter/metric_grid_plotter.jl")
include("../../plotter/blocks_interfaces_boundaries.jl")


using MAT


function SplitBlock(block, splitLocations, bndInfo, interInfo)
    blocks = []
    horzSplits = [1, splitLocations[1]..., size(block, 2)]
    vertSplits = [1, splitLocations[2]..., size(block, 3)]

    # parentni = size(block, 2)
    # parentnj = size(block, 3)
    # parentnk = 1

    blockId = 1
    blockBoundaries = []
    interfaces = []
    for j in 1:length(vertSplits)-1
        for i in 1:length(horzSplits)-1
            # println("Splitting block: ", blockId, " at i: ", i, " j: ", j)
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
                # "parentDims" => (parentni, parentnj, parentnk)
            )

            boundaries = GetTouchingBoundaries(blockInfo, bndInfo)
            append!(blockBoundaries, boundaries)
            
            # get interface info
            # if not at the end, look forward to the next block
            if i < length(horzSplits) - 1
                blockBId = blockId + 1
                # println("blockID: ", blockId, " next blockID: ", blockBId)
                push!(interfaces, Dict(
                    "blockA" => blockId, "start_blkA" => [ni,1,1], "end_blkA" => [ni,nj,nk],
                    "blockB" => blockBId, "start_blkB" => [1,1,1], "end_blkB" => [1,nj,nk],
                    "offset" => [0.0, 0.0, 0.0], "angle" => 0.0))
            end

            # if not at top, look up
            if j < length(vertSplits) - 1
                blockBId = blockId + length(horzSplits) - 1
                # println("blockID: ", blockId, " next blockID: ", blockBId)
                push!(interfaces, Dict(
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


#################
# Real Metric Data
#################

# load in the initial grid - no trailing edge 
initialGrid = GetAirfoilGrid(airfoilPath="examples/single_block_ns/airfoil/data/A-airfoil.txt", radius = 2)
# throw away trailing edge stuff
airfoilGrid = initialGrid[:, 101:end-100, :]
airfoil = airfoilGrid[:,:,1]


metricFunc = make_getMetric(airfoil;
    A_airfoil = 50.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
    A_origin  = 500.0,  ℓ_origin  = 0.2, p_origin  = 10,   
    floor     = 1e-4,  origin_center=(0.85, -0.1),
profile   = :rational)  # or :gauss


# define the boundary information
bndInfo = getBoundaryConditions(airfoilGrid)

# define interInfo
interInfo = Any[]


plt1 = plot_blocks_interfaces_boundaries([airfoilGrid], interInfo, bndInfo;
    grid_stride= 1,
    show_block_ids=true,
    legend=false,
    boundary_colors=Dict("BCWall"=>:purple, "BCInflow"=>:green, "BCOutflow"=>:blue) # optional
)



splitLocations = [
    [ 300, 500 ],   # split along the x axis
    [ 40 ]      # split along the y axis
]


blocks, bndInfo, interInfo = SplitBlock(airfoilGrid, splitLocations, bndInfo, interInfo)



plt2 = plot_blocks_interfaces_boundaries(blocks, interInfo, bndInfo;
    grid_stride= 1,
    show_block_ids=true,
    legend=false,
    boundary_colors=Dict("BCWall"=>:purple, "BCInflow"=>:green, "BCOutflow"=>:blue) # optional
)

p3 = plot(plt1, plt2, layout=@layout([a b]), size=(1000, 500),
    title="Single Block Splitting",
    xlabel="x", ylabel="y",
    aspect_ratio=:equal)


# #     boundary=nothing,
# #     # boundary=nothing,
# #     boundary_style=:fill,
# #     boundary_fillcolor=:black, boundary_fillalpha=1,
# #     mask_inside=true,
# #     title="Distance-based isotropic metric: large near airfoil and (0,0)",
# #     xlabel="x", ylabel="y",
# #     cb_label="M11(x,y)", clim=(0, 1.0),
# #     colormap=:viridis, equal_aspect=true)


# w = (x,y) -> first(metricFunc(x,y))   

# xl,xr = -2, 2 
# yl, yr = -2, 2
# xs = range(xl, xr, length=450)
# ys = range(yl, yr, length=300)

# p1, _ = plot_scalar_field(w, xs, ys; boundary=nothing,
#                          title="Distance-based isotropic metric",
#                          cb_label="M11(x,y)", equal_aspect=true, colormap =cgrad([RGB(1,1,1), RGB(1,0,0)]))#colormap=cgrd(:imola, rev=true))

# plot!(p1, xlims=(xl, xr), ylims=(yl, yr))

# X = block[1, :, :]
# Y = block[2, :, :]



# for j in 1:size(X, 1)
#     plot!(p1, X[j, :], Y[j, :], color=:black, lw=0.8, label=false, alpha=0.7)
# end

# for i in 1:size(X, 2)
#     plot!(p1, X[:, i], Y[:, i], color=:black, lw=0.8, label=false, alpha=0.7)
# end



# # add a zoomed in figure 
# p2 = deepcopy(p1)
# plot!(p2, xlims=(-0.15, 1.15), ylims=(-0.25, 0.25),
#         title="Zoomed in view",
#         xlabel="x", ylabel="y",
#         cb_label="M11(x,y)", equal_aspect=true, colormap=:imola)

# for j in 1:size(X, 1)
#     plot!(p2, X[j, :], Y[j, :], color=:black, lw=0.5, label=false, alpha=0.4)
# end

# for i in 1:size(X, 2)
#     plot!(p2, X[:, i], Y[:, i], color=:black, lw=0.5, label=false, alpha=0.4)
# end

    
# scatter!(p2, X[:, :], Y[:, :], color=:blue, label="", markersize=1.5, 
#          markerstrokewidth=0, shape=:diamond,)

# p3 = plot(p1, p2, layout=@layout([a b]), size=(1000, 500))