include("../../src/GridGeneration.jl")
include("GetAirfoilGrid.jl")
include("metric/Metric.jl")

using Plots


function plotBoundary(airfoil, boundarySection, m_vals, m_vals_section, xs, method, name)
    p1 = plot(title = "Discrete 1D Metric Values", aspect_ratio = 1)
    plot!(p1, boundarySection[1, :], boundarySection[2, :], label="Boundary Section", color=:yellow, linewidth=3)
    plot!(p1, airfoil[1, :], airfoil[2, :], label="Boundary", color=:black, linewidth=0.5)
    scatter!(p1, airfoil[1, :], airfoil[2, :], label="Discrete m(x)", marker_z = m_vals, markersize=0.5, markerstrokewidth=0, c = :brg)

    p2 = plot( title = "Zoomed Leading Edge (no aspect_ratio)",
            xlims = (boundarySection[1, 1] - 0.1, boundarySection[1, end] + 0.1),
            ylims = (minimum(boundarySection[2, :]) - 0.1, maximum(boundarySection[2, :]) + 0.1))

    plot!(p2, airfoil[1, :], airfoil[2, :], label="Boundary", color=:black, linewidth=0.5)
    scatter!(p2, airfoil[1, :], airfoil[2, :], label="Airfoil", marker_z = m_vals, markersize=0.5, markerstrokewidth=0, c = :brg)


    p3 = plot( title = "Discrete 1D Metric Boundary Metric Values (no aspect_ratio)")
    plot!(p3, boundarySection[1, :], boundarySection[2, :], label="Boundary", color=:black, linewidth=1)
    scatter!(p3, boundarySection[1, :], boundarySection[2, :], label="Discrete Metric Values ($(length(boundarySection[1,:])))", marker_z = m_vals_section, markersize=2, markerstrokewidth=0, c = :brg)

    titlepanel = plot(title = "1D Metric (method = $method) with $name clustering",
                    framestyle = :none, grid = false, ticks = nothing)

                    
    # Coordinates of the arrow base (e.g., near the leading edge)
    x0 = boundarySection[1, 1]
    y0 = boundarySection[2, 1]

    # Direction vector (e.g., pointing to the next boundary point)
    dx = boundarySection[1, 2] - boundarySection[1, 1]
    dy = boundarySection[2, 2] - boundarySection[2, 1]

    # Normalize and scale the vector to control arrow size
    arrow_scale = 0.05
    norm_factor = sqrt(dx^2 + dy^2)
    dx *= arrow_scale / norm_factor
    dy *= arrow_scale / norm_factor

    # Add arrow to p1, p2, and/or p3
    quiver!(p1, [x0], [y0], quiver=([dx], [dy]), color=:red, label="s direction", lw=1)
    quiver!(p3, [x0], [y0], quiver=([dx], [dy]), color=:red, label="s direction", lw=1)

    p5 = plot(xs, m_vals_section, title = "m(x(s)) along boundary section",
            xlabel = "s", ylabel = "m(x(s))", label = "m(x(s))",
            legend = :topright, linewidth=2, markershape=:circle)



    # Combine the plots
    # p4 = plot(titlepanel, p1, p2, p3, layout = @layout([A{0.001h}; b c; d]), size = (1000, 600))
    p4 = plot(titlepanel, p1, p3, p5, layout = @layout([A{0.001h}; b c; d]), size = (1000, 600))
    return p4
end



# build the metric
saveFig = false
method = "local" 

folder = "MetricReformulation/"
path = "docs/src/assets/images/$folder/"


scale = 40000

problems = [1] # 1: x=0, 2: x=1, 3: uniform
names = ["x=0", "x=1", "uniform"]

initialGrid = GetAirfoilGrid(airfoilPath = "examples/airfoil/data/A-airfoil.txt", radius = 3)

bottom = initialGrid[:,:,1]
airfoil = bottom[:, 101:end-100]  
SectionIndices = 10:10:250
# SectionIndices = 1:length(airfoil[1, :])

boundarySection = airfoil[:, SectionIndices]

x = boundarySection[1, :]
y = boundarySection[2, :]

# segment lengths (N-1)
Δx = diff(x)
Δy = diff(y)
Δs = sqrt.(Δx.^2 .+ Δy.^2)

xs = [0.0; cumsum(Δs)]  


for (name, problem) in zip(names, problems)
    M_func = (x,y) -> Metric(x, y, scale, problem)

    m_vals = GridGeneration.Get1DMetric(airfoil, M_func, method = method)
    m_vals_section = GridGeneration.Get1DMetric(boundarySection, M_func, method = method)



    p = plotBoundary(airfoil, boundarySection, m_vals, m_vals_section, xs, method, name)

    imageName = "metricboundary_$(name)_$(method).svg"
    imagePath = "$path$imageName"


    if saveFig
        savefig(p, imagePath)
    else
        display(p)
    end
end