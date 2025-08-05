include("../../src/GridGeneration.jl")
include("GetAirfoilGrid.jl")
include("metric/Metric.jl")

using Plots


function plotBoundary(airfoil, boundarySection, m_vals, m_vals_section, method, name)
    p1 = plot(title = "Discrete Metric Values", aspect_ratio = 1)
    plot!(p1, airfoil[1, :], airfoil[2, :], label="Boundary", color=:black, linewidth=0.5)
    scatter!(p1, airfoil[1, :], airfoil[2, :], label="Airfoil", marker_z = m_vals, markersize=0.5, markerstrokewidth=0, c = :brg)

    p2 = plot( title = "Zoomed Leading Edge (no aspect_ratio)",
            xlims = (boundarySection[1, 1] - 0.1, boundarySection[1, end] + 0.1),
            ylims = (minimum(boundarySection[2, :]) - 0.1, maximum(boundarySection[2, :]) + 0.1))

    plot!(p2, airfoil[1, :], airfoil[2, :], label="Boundary", color=:black, linewidth=0.5)
    scatter!(p2, airfoil[1, :], airfoil[2, :], label="Airfoil", marker_z = m_vals, markersize=0.5, markerstrokewidth=0, c = :brg)


    p3 = plot( title = "Discrete Boundary Metric Values (no aspect_ratio)")
    plot!(p3, boundarySection[1, :], boundarySection[2, :], label="Boundary", color=:black, linewidth=1)
    scatter!(p3, boundarySection[1, :], boundarySection[2, :], label="Discrete Metric Values ($(length(boundarySection[1,:])))", marker_z = m_vals_section, markersize=2, markerstrokewidth=0, c = :brg)

    titlepanel = plot(title = "1D Metric (method = $method) with $name clustering",
                    framestyle = :none, grid = false, ticks = nothing)

    p4 = plot(titlepanel, p1, p2,  p3, layout = @layout([A{0.001h}; b c; d]), size = (1000, 600))

    return p4
end



# build the metric
saveFig = true
method = "local" # "local" or "nuclear"

folder = "Mapping2Dto1D"
path = "docs/src/assets/images/$folder/"

scale = 40000

problems = [1] # 1: x=0, 2: x=1, 3: uniform
names = ["x=0", "x=1", "uniform"]

initialGrid = GetAirfoilGrid(airfoilPath = "examples/airfoil/data/A-airfoil.txt", radius = 3)

bottom = initialGrid[:,:,1]
airfoil = bottom[:, 101:end-100]  
SectionIndices = 100:10:250
# SectionIndices = 1:length(airfoil[1, :])

boundarySection = airfoil[:, SectionIndices]

for (name, problem) in zip(names, problems)


M_func = (x,y) -> Metric(x, y, scale, problem)

m_vals = GridGeneration.Get1DMetric(airfoil, M_func, method = method)
m_vals_section = GridGeneration.Get1DMetric(boundarySection, M_func, method = method)


p = plotBoundary(airfoil, boundarySection, m_vals, m_vals_section, method, name)

display(p)

# imageName = "metricboundary_$(name)_$(method).svg"
# imagePath = "$path$imageName"


# if saveFig
#     savefig(p, imagePath)
# else
#     display(p)
# end

end