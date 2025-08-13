include("../../src/GridGeneration.jl")

include("airfoil/GetAirfoilGrid.jl")
include("airfoil/GetBoundary.jl")
include("airfoil/metric/Metric.jl")
include("airfoil/metric/CustomMetric.jl")


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
    @info "Optimal number of points: $optN"

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
    # projectedLeft = left
    # projectedRight = right

    bottom = block[:, :, 1]
    top    = block[:, :, end]
    projectedBottom, projectedTop = process_edge_pair(bottom, top, metricFunc, solver)

    X,Y = GridGeneration.TFI_2D([projectedTop', projectedRight', projectedBottom', projectedLeft'])

    # p = scatter(projectedBottom[1, :], projectedBottom[2, :], label="Bottom", color=:blue)
    # scatter!(p, projectedTop[1, :], projectedTop[2, :], label="Top", color=:red)
    # display(p)

    computedBlock = permutedims(cat(X, Y; dims=3), (3, 1, 2))
    # computedBlock = block

    return computedBlock, bndInfo, interInfo
end



#################
# Set up the input
#################

#################
# Real Metric Data
#################

# # set up metric data
# metricPath = "examples/single_block_ns/airfoil/metric/A-airfoil_grid_data.mat"
# # load metric
# metricData = matread(metricPath)
# # set up metric function
# tree, refs = GridGeneration.setup_metric_tree(metricData)
# metricFunc = (x,y) -> 0.001 * GridGeneration.find_nearest_kd(metricData, tree, refs, x, y)


metricFunc = make_getMetric(airfoil;
    A_airfoil = 100.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   # tight decay from the airfoil
    A_origin  = 100.0,  ℓ_origin  = 0.10, p_origin  = 2,   # hotspot at (0,0)
    floor     = 1e-4,
profile   = :rational)  # or :gauss


# load in the initial grid - no trailing edge 
initialGrid = GetAirfoilGrid(airfoilPath="examples/single_block_ns/airfoil/data/A-airfoil.txt", radius = 10)
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
p = plot(aspect_ratio = 1 )

for j in 1:size(X, 1)
    plot!(p, X[j, :], Y[j, :], color=:red, lw=0.5, label=false)
end

for i in 1:size(X, 2)
    plot!(p, X[:, i], Y[:, i], color=:black, lw=0.5, label=false)
end
display(p)
