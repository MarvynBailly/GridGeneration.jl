include("../../src/GridGeneration.jl")

include("data/GetAirfoilGrid.jl")
include("data/GetBoundary.jl")
include("metric/Metric.jl")
include("metric/CustomMetric.jl")

include("../../plotter/metric_grid_plotter.jl")
include("../../plotter/blocks_interfaces_boundaries.jl")

using MAT

function procPlot(p, saveFig, filePath)
    # if saveFig
    #     @info "saving fig at $filePath"
    #     savefig(p, filePath)
    # else
    display(p) 
    # end
end


problem = 2
name = ["Constant", "Leading Edge", "Trailing Edge", "Leading & Trailing Edge", "Custom", "Real"][problem]

# splitLocations = [
#     [ 300, 500 ],   # split along the x axis
#     [ 40, 80 ]      # split along the y axis
# ]

splitLocations = [
    [ ],   # split along the x axis
    [  ]      # split along the y axis
]

# set up real metric data
metricPath = "examples/airfoil/metric/A-airfoil_grid_data.mat"
# load metric
metricData = matread(metricPath)
# set up metric tree for fast nearest neighbor search
tree, refs = GridGeneration.setup_metric_tree(metricData)



# load in the initial grid - no trailing edge 
initialGrid = GetAirfoilGrid(airfoilPath="examples/airfoil/data/A-airfoil.txt", radius = 2)
# throw away trailing edge stuff
airfoilGrid = initialGrid[:, 101:end-100, :]
airfoil = airfoilGrid[:,:,1]
    
# reset
bndInfo = getBoundaryConditions(airfoilGrid)
interInfo = Any[]

##########################
# METRIC
##########################
scale = 0.001#80000
metricFunc = (x,y) -> Metric(x,y; problem = problem, scale = scale)

##########################
# SOLVE
##########################

# blocksSolve, bndInfoSolve, interInfoSolve = GridGeneration.SolveAllBlocks(metricFunc, blocksSplit, bndInfoSplit, interInfoSplit; solver="2ndorder")










##########################
# INFO
# ##########################
# verbose = true
# function printBlockInfo(blocks; verbose = false)
#     if verbose
#         println("===========================")
#         println("Number of blocks after splitting: ", length(blocksSolve))
#         println("The size of each block is: ", )
#         totalPoints = 0
#         for (i, block) in enumerate(blocksSolve)
#             println("- Block $i size: ", size(block))
#             totalPoints += size(block, 2) * size(block, 3)
#         end
#         println("Total number of points in all blocks: ", totalPoints)
#         println("===========================")
#     end
# end

# printBlockInfo(blocksSolve; verbose = verbose)





# @info "Continue to plot?"
# readline()




##########################
# # PLOTS
# ##########################
# bndInfo = getBoundaryConditions(airfoilGrid)
# interInfo = Any[]
# plt1 = plot_blocks_interfaces_boundaries([airfoilGrid], interInfo, bndInfo;
#     grid_stride= 1,
#     show_block_ids=true,
#     legend=false,
#     titleName="Before Splitting",
#     boundary_colors=Dict("BCWall"=>:purple, "BCInflow"=>:green, "BCOutflow"=>:blue) # optional
# )


# ########## METRIC
w = (x,y) -> first(metricFunc(x,y))
xl,xr = -0.5, 1.5
yl, yr = -1.0, 1.0
xs = range(xl, xr, length=450)
ys = range(yl, yr, length=300)

plt4, _ = plot_scalar_field(w, xs, ys; boundary=airfoil,
                        title="Isotropic $name Metric",
                        cb_label="M(x,y)", equal_aspect=true, colormap =cgrad([RGB(1,1,1), RGB(1,0,0)]))#colormap=cgrd(:imola, rev=true))

plot!(plt4, xlims=(xl, xr), ylims=(yl, yr))


# filename =  "airfoil_metric_$problem"
# filepath = "$path/$filename.svg"
# procPlot(plt4, saveFig, filepath)

display(plt4)

########## SPLIT INFORMATION

# filename =  "airfoil_split_$problem"
# filepath = "$path/$filename.svg"
# procPlot(plt2, saveFig, filepath)

# ########## ORIGINAL AND FINAL GRID
# plt5 = plot(title="Original Airfoil Grid",
#     xlabel="x", ylabel="y",
#     aspect_ratio=:equal)
# plot!(plt5, airfoilGrid[1, :, :], airfoilGrid[2, :, :], color=:black, lw=0.8, label=false, alpha=0.7)

# plt6 = plot(title="Final Grid", xlabel="x", ylabel="y",
#     aspect_ratio=:equal)

# for block in blocksSolve
#     X = block[1, :, :]
#     Y = block[2, :, :]

#     for j in 1:1:size(X, 1)
#         plot!(plt6, X[j, :], Y[j, :], color=:black, lw=0.05, label=false, alpha=1.0)
#         plot!(plt4, X[j, :], Y[j, :], color=:black, lw=0.6, label=false, alpha=0.7)
#     end

#     for i in 1:1:size(X, 2)
#         plot!(plt6, X[:, i], Y[:, i], color=:black, lw=0.05, label=false, alpha=1.0)
#         plot!(plt4, X[:, i], Y[:, i], color=:black, lw=0.6, label=false, alpha=0.7)
#     end
# end

# display(plt6)

# p4 = plot(plt5, plt4, plt6, layout=@layout([a ; b;  c]), size=(500, 1000),
#     xlabel="x", ylabel="y",
#     aspect_ratio=:equal)


# p3 = plot(plt1, plt2, plt3, layout=@layout([a ; b ; c]), size=(500, 1000),
#     xlabel="x", ylabel="y",
#     aspect_ratio=:equal)

# p6 = plot(p3, p4, layout=@layout([a b]), size=(1000, 1000),
#     xlabel="x", ylabel="y",
#     aspect_ratio=:equal)

# display(p6)