using Base.Threads

include("../StorStegSolver.jl")
include("../BlockSolver.jl")
include("../../examples/airfoil/data/GetAirfoilGrid.jl")
include("../../examples/airfoil/data/GetBoundary.jl")





##############################################
#################### Main ####################
##############################################
function AirfoilSolver_threaded(params)
    
    # preform block splitting, solve all blocks, and use TFI to generate the grids
    blocks, bndInfo, interInfo = BlockSplitting(
        initialGrid = params.initialGrid, 
        bndInfo = params.bndInfo, 
        interInfo = params.interInfo, 
        splitLocations = params.splitLocations, 
        M = params.M,
        showPlots = params.showPlots)


    # block = blocks[1]
    # smoothBlocks = similar(blocks)

    smoothBlocks = Vector{Array{Float64,3}}(undef, length(blocks))
    finalErrors = Vector{Float64}(undef, length(blocks))
    finalIterations = Vector{Int}(undef, length(blocks))

    println("Number of threads: ", nthreads())

    @threads for i in eachindex(blocks)
        println("Solving for block $i on thread $(threadid())...")
        xr, yr, finalError, finalIter = EllipticSolver(x = blocks[i][1, :, :], y = blocks[i][2, :, :],
                                max_iter = params.max_iter, tol = params.tol, ω = params.ω,
                                s_left = params.s_left, a_decay_left = params.a_decay_left, b_decay_left = params.b_decay_left,
                                s_right = params.s_right, a_decay_right = params.a_decay_right, b_decay_right = params.b_decay_right,
                                s_top = params.s_top, a_decay_top = params.a_decay_top, b_decay_top = params.b_decay_top,
                                s_bottom = params.s_bottom, a_decay_bottom = params.a_decay_bottom, b_decay_bottom = params.b_decay_bottom,
                                use_top_wall = params.use_top_wall, use_bottom_wall = params.use_bottom_wall, use_left_wall = params.use_left_wall, use_right_wall = params.use_right_wall,
                                verbose = params.EllipticSolver_verbose
                                    ) 

        smoothBlocks[i] = permutedims(cat(xr, yr, dims=3), (3,1,2))
        finalErrors[i] = finalError
        finalIterations[i] = finalIter
    end
    @info "Final errors for each block: $(finalErrors)"
    @info "Final iterations for each block: $(finalIterations)"
    @info "Finished solving all blocks."

    if params.showPlots
        p = plot()
        for block in smoothBlocks
            x, y = block[1,:,:], block[2,:,:]
            plot_grid(x,y,"Block Grid After Elliptic", plt=p, c = RGB(0.0, 0.0, 0.0))
        end
        display(p)
    end
end

function AirfoilSolver(params)
    println("using $(params.use_left_wall) $(params.use_right_wall) $(params.use_bottom_wall) $(params.use_top_wall)")
    # preform block splitting, solve all blocks, and use TFI to generate the grids
    blocks, bndInfo, interInfo = BlockSplitting(
        initialGrid = params.initialGrid, 
        bndInfo = params.bndInfo, 
        interInfo = params.interInfo, 
        splitLocations = params.splitLocations, 
        M = params.M,
        showPlots = params.showPlots)

    @info "Split block into $(length(blocks)) blocks."

    # block = blocks[1]
    # smoothBlocks = similar(blocks)

    smoothBlocks = Vector{Array{Float64,3}}(undef, length(blocks))
    finalErrors = Vector{Float64}(undef, length(blocks))
    finalIterations = Vector{Int}(undef, length(blocks))


    for i in eachindex(blocks)
        @info "Solving for block $i ..."
        xr, yr, finalError, finalIter = EllipticSolver(x = blocks[i][1, :, :], y = blocks[i][2, :, :],
                                max_iter = params.max_iter, tol = params.tol, ω = params.ω,
                                s_left = params.s_left, a_decay_left = params.a_decay_left, b_decay_left = params.b_decay_left,
                                s_right = params.s_right, a_decay_right = params.a_decay_right, b_decay_right = params.b_decay_right,
                                s_top = params.s_top, a_decay_top = params.a_decay_top, b_decay_top = params.b_decay_top,
                                s_bottom = params.s_bottom, a_decay_bottom = params.a_decay_bottom, b_decay_bottom = params.b_decay_bottom,
                                use_top_wall = params.use_top_wall, use_bottom_wall = params.use_bottom_wall, use_left_wall = params.use_left_wall, use_right_wall = params.use_right_wall,
                                verbose = params.EllipticSolver_verbose
                                    ) 

        smoothBlocks[i] = permutedims(cat(xr, yr, dims=3), (3,1,2))
        finalErrors[i] = finalError
        finalIterations[i] = finalIter
    end
    @info "Final errors for each block: $(finalErrors)"
    @info "Final iterations for each block: $(finalIterations)"
    @info "Finished solving all blocks."

    if params.showPlots
        p = plot(legend=false)
        for block in smoothBlocks
            x, y = block[1,:,:], block[2,:,:]
            plot_grid(x,y,"Block Grid After Elliptic", plt=p, c = RGB(0.0, 0.0, 0.0))
        end
        display(p)
    end
    return smoothBlocks, blocks#, bndInfo, interInfo
end





# @time main(params)
# main(params)





# plot_grid(xr, yr, "Block 2 Grid After Elliptic", c = RGB(0.0, 0.0, 0.0))
# plot_grid(block[1, :, :], block[2, :, :], "Block 2 Grid", c = RGB(0.0, 0.0, 0.0))

# blocks_new =[] 

# for (i,block) in enumerate(blocks)
#     println("Solving for block $block...")
#     xr, yr = main(x=blocks[i][1, :, :], y=blocks[i][2, :, :])

# end


# --- Plotting function ---


#     display(plot_grid(x,y,"o",plt = plot(
#     axis = nothing, legend = false, framestyle = :nothing, aspect_ratio = :equal,  border = :none, background_color = RGB(.27, .33, .27),
# ), c = RGB(0.9, 0.9, 0.9)
# ))