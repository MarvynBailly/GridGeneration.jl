include("../../src/GridGeneration.jl")
include("airfoil/GetAirfoilGrid.jl")

using Plots, MAT


# set up metric data
metricPath = "examples/single_block_ns/airfoil/metric/A-airfoil_grid_data.mat"
# load metric
metricData = matread(metricPath)
# set up metric function
tree, refs = GridGeneration.setup_metric_tree(metricData)
metricFunc = (x,y) -> GridGeneration.find_nearest_kd(metricData, tree, refs, x, y)


initialGrid = GetAirfoilGrid(airfoilPath="examples/single_block_ns/airfoil/data/A-airfoil.txt", radius = 1)
airfoilGrid = initialGrid[:, 101:end-100, :]


metricValues = similar(airfoilGrid)

for i in 1:size(airfoilGrid, 2)
    for j in 1:size(airfoilGrid, 3)
        metricValues[1, i, j], metricValues[2, i, j] = metricFunc(airfoilGrid[1, i, j], airfoilGrid[2, i, j])
    end
end





leftEdge = airfoilGrid[:, 1, :]


leftMetric = GridGeneration.Get1DMetric(leftEdge, metricFunc)
xs = GridGeneration.ProjectBoundary2Dto1D(leftEdge)
sol = GridGeneration.SolveODE(leftMetric, leftMetric, length(xs), xs)
Nopt = GridGeneration.ComputeOptimalNumberofPoints(xs, leftMetric, sol[1, :])
println(Nopt)

solOpt = GridGeneration.SolveODE(leftMetric, leftMetric, Nopt, xs)



d = diff(solOpt, dims=2)
d = [0, d...]

p = scatter(airfoilGrid[1, 1, :], airfoilGrid[2, 1, :], title="M22", zcolor=metricValues[2, 1, :], colorbar=true, label=false, markerstrokewidth=0)
p2 = scatter(xs, zeros(length(xs)), title="Projected Left Edge", zcolor=leftMetric, colorbar=true, label=false, markerstrokewidth=0)
p3 = scatter(solOpt[1, :], solOpt[2, :], title="Solution", label=false, markerstrokewidth=0)
plot!(p3, range(0,1, length=length(d)),d)

p4 = plot(p, p2, p3, layout=@layout([a b; c]), size=(1000, 500))