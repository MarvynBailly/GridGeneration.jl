include("../../../src/GridGeneration.jl")
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
# sectionIndices = 100:300

# 
sectionIndices = 1:length(airfoil[1, :]) 
# sectionIndices = 300:400

boundarySection = airfoil[:, sectionIndices]

# sectionIndices = 500:600
# sectionIndices = 1:length(airfoil[1, :])


# random subset of the points to make the distribution less uniform
# keep = rand(length(sectionIndices)) .< 0.1
# newSectionIndices = sectionIndices[keep]   
# boundarySection = airfoil[:, newSectionIndices]

N = length(boundarySection[1, :])

xs = GridGeneration.ProjectBoundary2Dto1D(boundarySection)

# build the metric
saveFig = false
method = "local"
numMethod = "analytic"
folder = "PointProjection"
path = "docs/src/assets/images/$folder/"

scale = 40000

problems = [2] # 1: x=0, 2: x=1, 3: uniform
names = ["x=0", "x=1", "uniform"]

for (name, problem) in zip(names, problems) 
# problem = 1
name = names[problem]




M_func_test = (x,y) -> Metric(x, y, scale, problem)
M_u1_func_test = (x,y) -> MetricDerivative(x, y, scale, problem)

# m = GridGeneration.Get1DMetric(boundarySection, M_func_test, method = method)
# mx = GridGeneration.Get1DMetric(boundarySection, M_u1_func_test, method = method)

m = GridGeneration.Get1DMetric(boundarySection, M_func_test, method = method)
mx = GridGeneration.Get1DMetric(boundarySection, M_u1_func_test, method = method)

dir = 1#boundarySection[1, 1] > boundarySection[1, end] ? 1 : -1

sol_opt, sol = GridGeneration.GetOptimalSolution(m, mx, N, xs; dir=dir, method= numMethod)

x_sol, x_sol_opt = sol[1, :], sol_opt[1, :]


projected_points = GridGeneration.ProjectBoundary1Dto2D(boundarySection, x_sol, xs)
projected_points_opt = GridGeneration.ProjectBoundary1Dto2D(boundarySection, x_sol_opt, xs)

projected_x, projected_y = projected_points[1, :], projected_points[2, :]
projected_x_opt, projected_y_opt = projected_points_opt[1, :], projected_points_opt[2, :]



p1 = plot(xs, m, title = "1D Metric with Boundary Spacing (method = $method, $name clustering)",
        xlabel = "s", ylabel = "m(x(s))", label = "m(x(s))",
        legend = :best, linewidth=2)
scatter!(p1, xs, zeros(length(xs)), markershape=:circle, markersize=2, markerstrokewidth=0, c = :black, label="1D Boundary Values")
scatter!(p1, xs, m, markershape=:diamond, markersize=2, markerstrokewidth=0, c = :black, label="2D Boundary Values")
# add red arrow to show the direction of xs
quiver!(p1, [xs[1]], [m[1]], quiver=([xs[5] - xs[1]], [m[5] - m[1]]), color=:red, label="s direction", lw=1)



p2 = plot(title = "Optimal and Non-Optimal 1D Distributions (method = $method, $name clustering)",
        xlabel = "s", ylabel = "m(x(s))", label = "m(x(s))",
        legend = :best, linewidth=2)
scatter!(p2, xs, zeros(length(xs)), markershape=:circle, markersize=3, markerstrokewidth=0, c = :black, label="1D Boundary Values")
scatter!(p2, sol[1,:], zeros(length(sol[1,:])), markershape=:circle, markersize=2.5, markerstrokewidth=0, c = :red, label="Non-Optimal Solution (N = $(length(sol[1,:])))")
scatter!(p2, sol_opt[1,:], zeros(length(sol_opt[1,:])), markershape=:circle, markersize=2, markerstrokewidth=0, c = :green, label="Optimal Solution (N = $(length(sol_opt[1,:])))")

p3 = plot( title = "Optimal Points Projected on Boundary", aspect_ratio = 1)

plot!(p3, boundarySection[1, :], boundarySection[2, :], label="Boundary", color=:black, linewidth=1)
scatter!(p3, boundarySection[1, :], boundarySection[2, :], label="1D Metric", marker_z = m, markersize=4, markerstrokewidth=0)
scatter!(p3, projected_x_opt, projected_y_opt, label="Projected Points Optimal ($(length(projected_x_opt)))", color = :green, marker=:circle, markersize=2, markerstrokewidth=0)
# scatter!(p3, boundarySection[1, :], boundarySection[2, :], label="Discrete Metric Values ($(length(boundarySection[1,:])))", marker_z = m_vals_section, markersize=2, markerstrokewidth=0, c = :brg)

p5 = plot()
# p5 = plot(aspect_ratio = 1)
plot!(p5, boundarySection[1, :], boundarySection[2, :], label="Boundary", color=:black, linewidth=2)
plot!(p5, projected_x_opt,projected_y_opt, label="Projected Solution", color=:red, linewidth=1)
scatter!(p5, projected_x_opt, projected_y_opt, label="Projected Points Optimal ($(length(projected_x_opt)))", color = :green, marker=:circle, markersize=2, markerstrokewidth=0)

# p6 = plot(xlims=[-0.05,0.05], ylims=[-0.1, 0.1])
# plot!(p6, boundarySection[1, :], boundarySection[2, :], label="Boundary", color=:black, linewidth=2, title="Projected Points")
# plot!(p6, projected_x_opt,projected_y_opt, label="Projected Solution", color=:red, linewidth=1, title="leading edge zoom in")

# p7 = plot(p5, p6, layout = @layout([a b]), size = (800, 400))
# display(p7)
# scatter!(p6, projected_x_opt, projected_y_opt, label="Projected Points Optimal ($(length(projected_x_opt)))", color = :green, marker=:circle, markersize=2, markerstrokewidth=0)


# p3 = plot(p1, p2, layout = @layout([A{0.001h}; b c]), size = (1000, 600))
p4 = plot(p1, p2, p3, p5, layout = @layout([a ; b ; c ; d]), size = (800, 1000))





if saveFig
    savefig(p4, path * "ode_solution_$(name)_$(method).svg")
else
    display(p4)
end
end