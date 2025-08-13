include("../../src/GridGeneration.jl")

include("airfoil/GetAirfoilGrid.jl")
include("airfoil/GetBoundary.jl")

using Plots, MAT

##############
# functions
##############

function process_edge_pair(edgeA, edgeB, metricFunc, solver)
    metricA = GridGeneration.Get1DMetric(edgeA, metricFunc)
    metricB = GridGeneration.Get1DMetric(edgeB, metricFunc)

    xsA = GridGeneration.ProjectBoundary2Dto1D(edgeA)
    xsB = GridGeneration.ProjectBoundary2Dto1D(edgeB)

    solA = GridGeneration.SolveODE(metricA, metricA, length(xsA), xsA; method=solver)
    solB = GridGeneration.SolveODE(metricB, metricB, length(xsB), xsB; method=solver)

    optNA = GridGeneration.ComputeOptimalNumberofPoints(xsA, metricA, solA[1, :])
    optNB = GridGeneration.ComputeOptimalNumberofPoints(xsB, metricB, solB[1, :])

    optN = max(optNA, optNB)

    solOptA = GridGeneration.SolveODE(metricA, metricA, optN, xsA; method=solver)
    solOptB = GridGeneration.SolveODE(metricB, metricB, optN, xsB; method=solver)

    projectedA = GridGeneration.ProjectBoundary1Dto2D(edgeA, solOptA[1, :], xsA)
    projectedB = GridGeneration.ProjectBoundary1Dto2D(edgeB, solOptB[1, :], xsB)

    return projectedA, projectedB
end

function SolveBlock(block, bndInfo, interInfo, metricFunc; solver="analytic")
    left   = block[:, 1, :]
    right  = block[:, end, :]
    projectedLeft, projectedRight = process_edge_pair(left, right, metricFunc, solver)

    bottom = block[:, :, 1]
    top    = block[:, :, end]
    projectedBottom, projectedTop = process_edge_pair(bottom, top, metricFunc, solver)

    X,Y = GridGeneration.TFI_2D([projectedTop', projectedRight', projectedBottom', projectedLeft'])

    computedBlock = permutedims(cat(X, Y; dims=3), (3, 1, 2))

    return computedBlock, bndInfo, interInfo
end



#################
# Set up the input
#################

#################
# Real Metric Data
#################

# set up metric data
# metricPath = "examples/single_block_ns/airfoil/metric/A-airfoil_grid_data.mat"
# load metric
# metricData = matread(metricPath)
# set up metric function
# tree, refs = GridGeneration.setup_metric_tree(metricData)
# metricFunc = (x,y) -> GridGeneration.find_nearest_kd(metricData, tree, refs, x, y)

# load in the initial grid - no trailing edge 
initialGrid = GetAirfoilGrid(airfoilPath="examples/single_block_ns/airfoil/data/A-airfoil.txt", radius = 0.2)
# throw away trailing edge stuff
airfoilGrid = initialGrid[:, 101:end-100, :]

# define the boundary information
bndInfo = getBoundaryConditions(initialGrid)

# define interInfo
interInfo = Any[]

block, bndInfo, interInfo = SolveBlock(airfoilGrid, bndInfo, interInfo, metricFunc)

# plot the block
X = block[1, :, :]
Y = block[2, :, :]
p = plot()
for j in 1:size(X, 1)
    plot!(p, X[j, :], Y[j, :], color=:black, lw=0.5, label=false)
end

for i in 1:size(X, 2)
    plot!(p, X[:, i], Y[:, i], color=:black, lw=0.5, label=false)
end
display(p)
