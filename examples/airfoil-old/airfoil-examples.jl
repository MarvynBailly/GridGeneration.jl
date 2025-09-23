include("../../src/GridGeneration.jl")

include("data/GetAirfoilGrid.jl")
include("data/GetBoundary.jl")
include("metric/Metric.jl")
include("metric/CustomMetric.jl")

include("../../plotter/metric_grid_plotter.jl")
include("../../plotter/blocks_interfaces_boundaries.jl")


# load in the initial grid - no trailing edge 
initialGrid = GetAirfoilGrid(airfoilPath="examples/airfoil/data/A-airfoil.txt", radius = 2)
# throw away trailing edge stuff
airfoilGrid = initialGrid[:, 101:end-100, :]
airfoil = airfoilGrid[:,:,1]


#################
# Make a custom metric
#################

function procPlot(p, saveFig, filePath)
    if saveFig
        @info "saving fig at $filePath"
        savefig(p, filePath)
    else
        display(p)
    end
end

# problem = 5
problems = [1, 2, 3, 4, 5,]
names = ["Constant", "Leading Edge", "Trailing Edge", "Leading & Trailing Edge", "Custom"]

splitLocationsArray = [ 
    [
        [ 300, 500 ],   # split along the x axis
        [ 40, 80 ]      # split along the y axis
    ],
    [
        [ 300, 500 ],   # split along the x axis
        [ 40, 80 ]      # split along the y axis
    ],
    [
        [ 300, 500 ],   # split along the x axis
        [ 40, 80 ]      # split along the y axis
    ],
    [
        [ 300, 500 ],   # split along the x axis
        [ 40, 80 ]      # split along the y axis
    ],
    [
        [ 300, 500 ],   # split along the x axis
        [ 40, 80 ]      # split along the y axis
    ] 
]


folder = "Examples/airfoil"
path = "docs/src/assets/images/$folder/"

saveFig = true

for (problem, name) in zip(problems, names)

    # reset
    bndInfo = getBoundaryConditions(airfoilGrid)
    interInfo = Any[]
    splitLocations = splitLocationsArray[problem]

    ##########################
    # METRIC
    ##########################
    metricFunc = (x,y) -> Metric(x,y; problem = problem, scale = 4000)

    ##########################
    # SPLIT
    ##########################

    blocksSplit, bndInfoSplit, interInfoSplit = GridGeneration.SplitBlock(airfoilGrid, splitLocations, bndInfo, interInfo)

    plt2 = plot_blocks_interfaces_boundaries(blocksSplit, interInfoSplit, bndInfoSplit;
        grid_stride= 1,
        show_block_ids=true,
        legend=false,
        titleName="Blocks After Splitting",
        boundary_colors=Dict("BCWall"=>:purple, "BCInflow"=>:green, "BCOutflow"=>:blue) # optional
    )


    ##########################
    # SOLVE
    ##########################

    blocksSolve, bndInfoSolve, interInfoSolve = GridGeneration.SolveAllBlocks(metricFunc, blocksSplit, bndInfoSplit, interInfoSplit)

    ##########################
    # PLOTS
    ##########################

    ########## METRIC
    ##########################
# PLOTS
##########################
bndInfo = getBoundaryConditions(airfoilGrid)
interInfo = Any[]
plt1 = plot_blocks_interfaces_boundaries([airfoilGrid], interInfo, bndInfo;
    grid_stride= 1,
    show_block_ids=true,
    legend=false,
    titleName="Before Splitting",
    boundary_colors=Dict("BCWall"=>:purple, "BCInflow"=>:green, "BCOutflow"=>:blue) # optional
)


########## METRIC
w = (x,y) -> first(metricFunc(x,y))
xl,xr = -0.5, 1.5
yl, yr = -1.0, 1.0
xs = range(xl, xr, length=450)
ys = range(yl, yr, length=300)

plt4, _ = plot_scalar_field(w, xs, ys; boundary=nothing,
                        title="Isotropic $name Metric",
                        cb_label="M(x,y)", equal_aspect=true, colormap =cgrad([RGB(1,1,1), RGB(1,0,0)]))#colormap=cgrd(:imola, rev=true))

plot!(plt4, xlims=(xl, xr), ylims=(yl, yr))

# filename =  "airfoil_metric_$problem"
# filepath = "$path/$filename.svg"
# procPlot(plt4, saveFig, filepath)

########## SPLIT INFORMATION

# filename =  "airfoil_split_$problem"
# filepath = "$path/$filename.svg"
# procPlot(plt2, saveFig, filepath)

########## SPLIT INFORMATION - post solver
plt3 = plot_blocks_interfaces_boundaries(blocksSolve, interInfoSolve, bndInfoSolve;
    grid_stride= 1,
    show_block_ids=true,
    legend=false,
    titleName="After Solving",
    boundary_colors=Dict("BCWall"=>:purple, "BCInflow"=>:green, "BCOutflow"=>:blue) # optional
)

########## ORIGINAL AND FINAL GRID
plt5 = plot(title="Original Airfoil Grid",
    xlabel="x", ylabel="y",
    aspect_ratio=:equal)
plot!(plt5, airfoilGrid[1, :, :], airfoilGrid[2, :, :], color=:black, lw=0.8, label=false, alpha=0.7)

plt6 = plot(title="Final Grid", xlabel="x", ylabel="y",
    aspect_ratio=:equal)

for block in blocksSolve
    X = block[1, :, :]
    Y = block[2, :, :]

    for j in 1:size(X, 1)
        plot!(plt6, X[j, :], Y[j, :], color=:black, lw=0.6, label=false, alpha=0.7)
        plot!(plt4, X[j, :], Y[j, :], color=:black, lw=0.6, label=false, alpha=0.7)
    end

    for i in 1:size(X, 2)
        plot!(plt6, X[:, i], Y[:, i], color=:black, lw=0.6, label=false, alpha=0.7)
        plot!(plt4, X[:, i], Y[:, i], color=:black, lw=0.6, label=false, alpha=0.7)
    end
end

p4 = plot(plt5, plt4, plt6, layout=@layout([a ; b;  c]), size=(500, 1000),
    xlabel="x", ylabel="y",
    aspect_ratio=:equal)


p3 = plot(plt1, plt2, plt3, layout=@layout([a ; b ; c]), size=(500, 1000),
    xlabel="x", ylabel="y",
    aspect_ratio=:equal)

p6 = plot(p3, p4, layout=@layout([a b]), size=(1000, 1000),
    xlabel="x", ylabel="y",
    aspect_ratio=:equal)

display(p6)

    filename =  "airfoil_all_$problem"
    filepath = "$path/$filename.svg"
    procPlot(p6, saveFig, filepath)
end





