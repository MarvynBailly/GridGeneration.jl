using Plots
using LinearAlgebra
using MAT: matread


include("../src/GridGeneration.jl")
include("plotter/metric_grid_plotter.jl")
include("plotter/blocks_interfaces_boundaries.jl")




function BlockSplitting(;initialGrid=initialGrid, bndInfo=bndInfo, interInfo=interInfo, splitLocations=splitLocations, M=M, showPlots = false)




    plt1 = plot_blocks_interfaces_boundaries([initialGrid], interInfo, bndInfo;
        grid_stride= 1,
        show_block_ids=true,
        legend=false,
        titleName="Before Splitting",
        boundary_colors=Dict("BCWall"=>:purple, "BCInflow"=>:green, "BCOutflow"=>:blue)
    )

    # #################
    # ## Split the blocks
    # #################

    blocks, bndInfo, interInfo = GridGeneration.SplitBlock(initialGrid, splitLocations, bndInfo, interInfo)

    plt2 = plot_blocks_interfaces_boundaries(blocks, interInfo, bndInfo;
        grid_stride= 1,
        show_block_ids=true,
        legend=false,
        titleName="After Splitting",
        boundary_colors=Dict("BCWall"=>:purple, "BCInflow"=>:green, "BCOutflow"=>:blue) # optional
    )


    #################
    ##### Run the solver on them
    #################
    blocks, bndInfo, interInfo = GridGeneration.SolveAllBlocks(M, blocks, bndInfo, interInfo)


    plt3 = plot_blocks_interfaces_boundaries(blocks, interInfo, bndInfo;
        grid_stride= 1,
        show_block_ids=true,
        legend=false,
        titleName="After Solving",
        boundary_colors=Dict("BCWall"=>:purple, "BCInflow"=>:green, "BCOutflow"=>:blue) # optional
    )

    splitPlot = plot(plt1, plt2, plt3, layout=(1,3), size=(1800,600))
    if showPlots display(splitPlot) end
return blocks, bndInfo, interInfo
end

# include("../ellipitic/ss.jl")
# for block in blocks 
#     x, y = block[1,:,:], block[2,:,:]
#     display(plot_grid(x,y,"o",plt = plot(
#         axis = nothing, legend = false, framestyle = :nothing, aspect_ratio = :equal,  border = :none, background_color = RGB(.27, .33, .27),
#     ), c = RGB(0.9, 0.9, 0.9)
#     ))
# end

