include("../../src/GridGeneration.jl")
include("airfoil/GetAirfoilGrid.jl")
include("airfoil/metric/CustomMetric.jl")

using Plots, MAT


# # set up metric data
# metricPath = "examples/single_block_ns/airfoil/metric/A-airfoil_grid_data.mat"
# # load metric
# metricData = matread(metricPath)
# # set up metric function
# tree, refs = GridGeneration.setup_metric_tree(metricData)
# metricFunc = (x,y) -> GridGeneration.find_nearest_kd(metricData, tree, refs, x, y)


initialGrid = GetAirfoilGrid(airfoilPath="examples/single_block_ns/airfoil/data/A-airfoil.txt", radius = 1)
airfoilGrid = initialGrid[:, 101:end-100, :]
airfoil = airfoilGrid[:,:,1]




getMetric = make_getMetric(airfoil;
    A_airfoil = 100.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   # tight decay from the airfoil
    A_origin  = 100.0,  ℓ_origin  = 0.10, p_origin  = 2,   # hotspot at (0,0)
    floor     = 1e-4,
profile   = :rational)  # or :gauss


metricValues = similar(airfoilGrid)

for i in 1:size(airfoilGrid, 2)
    for j in 1:size(airfoilGrid, 3)
        metricValues[1, i, j], metricValues[2, i, j] = getMetric(airfoilGrid[1, i, j], airfoilGrid[2, i, j])
    end
end



# topEdge = airfoilGrid[:, :, 1]

# topMetric = GridGeneration.Get1DMetric(topEdge, metricFunc)
# xs = GridGeneration.ProjectBoundary2Dto1D(topEdge)
# sol = GridGeneration.SolveODE(topMetric, topMetric, length(xs), xs)
# Nopt = GridGeneration.ComputeOptimalNumberofPoints(xs, topMetric, sol[1, :])

# color the points based on the metricFunc
p1 = scatter(airfoilGrid[1, :, :], airfoilGrid[2, :, :], title="M11", zcolor=metricValues[1, :, :], colorbar=true, label=false, markerstrokewidth=0)

p2 = scatter(airfoilGrid[1, :, :], airfoilGrid[2, :, :], title="M22", zcolor=metricValues[2, :, :], colorbar=true, label=false, markerstrokewidth=0)


# p2 = scatter(airfoilGrid[1, :, end], airfoilGrid[2, :, end], title="m", zcolor=topMetric, colorbar=true, label=false, markerstrokewidth=0)


p = plot(p1, p2, layout=(1, 2), size=(1000, 500))
display(p)