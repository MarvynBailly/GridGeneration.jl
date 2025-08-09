include("../../src/GridGeneration.jl")
include("GetAirfoilGrid.jl")
include("metric/Metric.jl")

using Plots, MAT, DelimitedFiles

function PlotProjectedPoints(projected_x, projected_y, projected_x_opt, projected_y_opt, x_bnd, y_bnd, metric, method, name)
    # 5. Plot the projected points on the airfoil
    p1 = plot(title = "Projected Points on Section (no ar)" )

    plot!(p1, x_bnd, y_bnd, 
        seriestype=:scatter, 
        label="Original Airfoil ($(length(x_bnd)))", 
        color=:black, 
        marker=:diamond,
        markersize=5, 
        markerstrokewidth=0)


    plot!(p1, x_bnd, y_bnd, 
        label="Original Airfoil", 
        color=:black, 

    )

    plot!(p1, projected_x, projected_y, 
        seriestype=:scatter, 
        label="Projected Points non-Optimal ($(length(projected_x)))", 
        color=:green, 
        markersize=3, 
        markerstrokewidth=0,
        xlabel="x",
        ylabel="y")

    plot!(p1, projected_x_opt, projected_y_opt, 
        seriestype=:scatter, 
        label="Projected Points Optimal ($(length(projected_x_opt)))", 
        color=:red, 
        marker=:rect,
        markersize=3, 
        markerstrokewidth=0,
        title="Projected Points on Airfoil Section",
        xlabel="x",
        ylabel="y")


    p2 = plot(projected_x_opt, projected_y_opt, 
        seriestype=:scatter, 
        label="Projected Points Optimal ($(length(projected_x_opt)))", 
        marker=:rect,
        color=:red, 
        markersize=2, 
        markerstrokewidth=0,
        title="Projected Points on Airfoil Section",
        xlabel="x",
        ylabel="y",
        aspect_ratio = 1)

    plot!(p2, x_bnd, y_bnd, 
        label="Original Airfoil", 
        color=:black, 
    )

    p3 = plot( title = "Discrete Boundary Metric Values", aspect_ratio = 1)
    plot!(p3, x_bnd, y_bnd, label="Boundary", color=:black, linewidth=1)
    
    scatter!(p3, x_bnd, y_bnd, label="Discrete Metric Values ($(length(x_bnd)))", marker_z = metric, markersize=4, markerstrokewidth=0, c = :brg)
    
    scatter!(p3, projected_x_opt, projected_y_opt, label="Projected Points Optimal ($(length(projected_x_opt)))", color = :red, marker=:rect, markersize=2, markerstrokewidth=0)


    titlepanel = plot(title = "1D Metric (method = $method) with $name clustering",
                framestyle = :none, grid = false, ticks = nothing)

    p5 = plot(titlepanel, p1, p2, p3, layout = @layout([A{0.001h}; b; c d]), size = (1000, 800))


    return p5
end


function ProjectBoundary1Dto2DGif(boundary, points, xs)
    projectedPoints = zeros(2, length(points))
    
    # make the plot have line and scatter points
    boundaryPlot = plot(boundary[1, :], boundary[2, :], 
        seriestype=:scatter, 
        label="Boundary", 
        color=:black, 
        markershape=:diamond,
        markersize=3, 
        markerstrokewidth=0,
        xlabel="x",
        ylabel="y",
        title="1D Boundary to 2D Projection")
    
    for (i, pnt) in enumerate(points)
        # clamp
        if i == 1
            intervalIndex = 1
            normalDist = 0
        elseif i == length(points)
            intervalIndex = length(xs) - 1
            normalDist = 1
        else
            intervalIndex = FindContainingIntervalIndex(pnt, xs)
            normalDist = ComputeNormalDistance(pnt, xs, intervalIndex)
        end
        
        projectPoint = ProjectPointOntoBoundary(normalDist, intervalIndex, boundary)
        projectedPoints[:, i] = projectPoint

        solPlot = plot(xs, zeros(length(points)), 
            seriestype=:scatter, 
            label="Original Points", 
            color=:black, 
            markersize=4, 
            markerstrokewidth=0,
            xlabel="x",
            ylabel="y",
            title="1D Point distribution")
        
        plot!(points, zeros(length(xs)), 
            seriestype=:scatter, 
            label="Solution Points", 
            color=:red, 
            markersize=3, 
            markerstrokewidth=0,
            xlabel="x",
            ylabel="y",
            title="1D Point distribution")

        # color the current point 
        scatter!(solPlot, [pnt], [0], 
            seriestype=:scatter, 
            label="Point $i", 
            color=:blue, 
            markersize=3, 
            markerstrokewidth=0) 

        distPlot = plot(xs[intervalIndex:intervalIndex+1], 
            zeros(2), 
            label="Interval $intervalIndex",
            color=:orange,
            markersize=3,
            markerstrokewidth=0,
            xlabel="x",
            ylabel="y",
            title="Interval containing point $i")
        scatter!(distPlot, [pnt], [0], 
            seriestype=:scatter, 
            label="Point $i", 
            color=:blue, 
            markersize=3, 
            markerstrokewidth=0)

        projDistPlot = scatter([boundary[1, intervalIndex]], [boundary[2, intervalIndex]],
            seriestype=:scatter,  
            label="Boundary Point i", 
            color=:green, 
            markersize=3, 
            markerstrokewidth=0,
            xlabel="x",
            ylabel="y",
            title="Boundary")

        scatter!(projDistPlot, [boundary[1, intervalIndex+1]], [boundary[2, intervalIndex + 1]],
            seriestype=:scatter,  
            label="Boundary Point i + 1", 
            color=:black, 
            markersize=3, 
            markerstrokewidth=0,
            xlabel="x",
            ylabel="y",
            title="Boundary")
        
        scatter!(projDistPlot, [projectPoint[1]], [projectPoint[2]], 
            seriestype=:scatter, 
            label="Projected Point", 
            color=:red, 
            markersize=5, 
            markerstrokewidth=0)
        
        scatter!(boundaryPlot, [projectPoint[1]], [projectPoint[2]], 
            seriestype=:scatter, 
            label="", 
            color=:red, 
            markersize=5, 
            markerstrokewidth=0)

        # println("Projecting point $i: $pnt onto boundary at index $intervalIndex with normalDist = $normalDist")    
        masterPlot = plot(solPlot, distPlot, projDistPlot, boundaryPlot, layout = @layout([a b; c d]), size = (800, 600))
        # display(masterPlot)
    end

    return projectedPoints
end

#####################
# SET UP DOMAIN
#####################

initialGrid = GetAirfoilGrid(airfoilPath = "examples/airfoil/data/A-airfoil.txt", radius = 3)

bottom = initialGrid[:,:,1]
airfoil = bottom[:, 101:end-100]  
sectionIndices = 100:300
# sectionIndices = 500:600

# uniform
# boundarySection = airfoil[:, sectionIndices]

# starting sparse and going denser
# Parameters
n = length(sectionIndices)
min_density = 0.1  # minimum probability of selecting a point (sparse)
max_density = 1.0  # maximum probability (dense)

# Create a function for density that increases as we move along the indices
density_function(i) = min_density + (max_density - min_density) * (i / n)

# Use the density function to decide which indices to keep
keep = [rand() < density_function(i) for i in 1:n]

# Apply the selection
newSectionIndices = sectionIndices[keep]
boundarySection = airfoil[:, newSectionIndices]

# gross
# keep = rand(length(sectionIndices)) .< 0.3  
# newSectionIndices = sectionIndices[keep]   
# boundarySection = airfoil[:, newSectionIndices]


N = length(boundarySection[1, :])

xs = GridGeneration.ProjectBoundary2Dto1D(boundarySection)

# build the metric
saveFig = false
method = "local"
numMethod = "2ndorder"
folder = "PointProjection"
path = "docs/src/assets/images/$folder/"

scale = 40000

problem = 2
name = "x=1"




M_func_test = (x,y) -> Metric(x, y, scale, problem)
M_u1_func_test = (x,y) -> MetricDerivative(x, y, scale, problem)

m = GridGeneration.Get1DMetric(boundarySection, M_func_test, method = method)
mx = GridGeneration.Get1DMetric(boundarySection, M_u1_func_test, method = method)

f = x -> mx_func(x) / (2 * m_func(x))

# get the direction of the boundary
dir = boundarySection[1, 1] > boundarySection[1, end] ? -1 : 1

sol_opt, sol = GridGeneration.GetOptimalSolution(m, mx, N, xs, dir; method= numMethod)

sol_x = sol[1, :]

m_func = GridGeneration.build_interps_linear(xs, m)
mx_func = GridGeneration.build_interps_linear(xs, mx)


p = plot(boundarySection[1, :], boundarySection[2, :], 
    seriestype=:scatter, 
    label="Boundary", 
    marker_z = m,
    markershape=:diamond,
    markersize=3, 
    markerstrokewidth=0,
    xlabel="x",
    ylabel="y",
    title="1D Boundary to 2D Projection")

p1 = scatter(xs, marker_z =m, zeros(N), label="1D projection")

p2 = scatter(xs, m, label="1d Metric values")
# scatter!(p2, xs, mx, label="1d Metric derivative values", color=:red)

p3 = scatter(xs, m_func.(xs), label="Interpolated 1d Metric values")
# scatter!(p3, xs, mx_func.(xs), label="Interpolated 1d Metric derivative values", color=:red)

p4 = scatter(xs, discreteForcing, label="Discrete forcing function")

p5 = scatter(sol_x, zeros(length(sol)), label="")

# p3 = plot(p, p1, p2, layout = @layout([a; b; c]), size = (800, 600))
p3 = plot(p, p1, p2, p3, p4, p5, layout = @layout([a b; c d; e f]), size = (800, 600))

