include("../src/GridGeneration.jl")

using Plots, MAT, DelimitedFiles


function plotGrid!(p, X, Y, clr; sz = 1.0, skipV=1, skipH=1)
    for j in 1:skipV:size(X, 1)
        plot!(p, X[j, :], Y[j, :], color=clr, lw=sz, label=false)
    end

    for i in 1:skipH:size(X, 2)
        plot!(p, X[:, i], Y[:, i], color=clr, lw=sz, label=false)
    end
end














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



#####################
# SET UP DOMAIN
#####################

airfoilPath = "examples/airfoil/A-airfoil.txt"

# read the airfoil data
airfoilData = readdlm(airfoilPath, '\t', skipstart=1)

RADIUS = 3
vertN = 100
horzN = 100

# number of points in the airfoil
airfoilN = length(airfoilData[:, 1])

# set up a c grid around the provided inner boundary
boundary = GridGeneration.SetupDomain(
    airfoilData, 
    RADIUS, 
    vertN, 
    horzN;
    type = "cgrid"
)


initialGrid = GridGeneration.TFI_2D(boundary)
initialGrid = permutedims(cat(initialGrid...; dims=3), (3, 1, 2));

bottom = initialGrid[:,:,1]
airfoil = bottom[:, 101:end-100]  
SectionIndices = 100:250
# SectionIndices = 1:length(airfoil[1, :])
boundarySection = airfoil[:, SectionIndices]

arclength = sum(sqrt.(diff(boundarySection[1, :]).^2 .+ diff(boundarySection[2, :]).^2))
N = length(boundarySection[1, :])

# boundarySection is 2×N (rows: x,y; columns: points in order)
x = boundarySection[1, :]
y = boundarySection[2, :]

# segment lengths (N-1)
Δx = diff(x)
Δy = diff(y)
Δs = sqrt.(Δx.^2 .+ Δy.^2)

xs = [0.0; cumsum(Δs)]   # length N, xs[1]=0, xs[end]=arclength

# normalize 
xs = xs ./ xs[end]  # now xs is in [0, 1]


# build the metric
saveFig = true
method = "local" # "local" or "nuclear"


problems = [3] # 1: x=0, 2: x=1, 3: uniform
for problem in problems 
folder = "Mapping2Dto1D"
path = "docs/src/assets/images/$folder/"

if problem == 1
    name = "x=0"
elseif problem == 2
    name = "x=1"
elseif problem == 3
    name = "uniform"
end

scale = 40000

M_func = (x,y) -> Metric(x, y, scale, problem)
M_u1_func = (x,y) -> MetricDerivative(x, y, scale, problem)

M_func_values = M_func.(boundarySection[1, :], boundarySection[2, :])


m_vals = GridGeneration.Get1DMetric(airfoil, M_func, method = method)
m_vals_section = GridGeneration.Get1DMetric(boundarySection, M_func, method = method)





m = GridGeneration.Get1DMetric(boundarySection, M_func, method = method)
mx = GridGeneration.Get1DMetric(boundarySection, M_u1_func, method = method)
m_func = build_interps_linear(xs, m)
mx_func = build_interps_linear(xs, mx)

# pass the functions to the solver
sol = GridGeneration.SolveODE(m_func, mx_func, N, xs[1], xs[end])

sigma_opt = GridGeneration.ComputeOptimalSpacing(sol[1, :], m, xs)
N_opt = ceil(Int, 1 / sigma_opt)

@info("Optimal number of points: ", N_opt)

if N_opt < 2
    @info("Optimal number of points is less than 2, using N=3")
    N_opt = 10
end

sol_opt = GridGeneration.SolveODE(m_func, mx_func, N_opt, xs[1], xs[end])

x_sol = sol[1, :]
x_sol_opt = sol_opt[1, :]



# scale the solutions back to the correct size
x_sol = x_sol * xs[end]
x_sol_opt = x_sol_opt * xs[end]

projected_x, projected_y = InterpolateBoundaryManually(x_sol, boundarySection)
projected_x_opt, projected_y_opt = InterpolateBoundaryManually(x_sol_opt, boundarySection)




p = plotBoundary(airfoil, boundarySection, m_vals, m_vals_section, method, name)
imageName = "metricboundary_$(name)_$(method).svg"
imagePath = "$path$imageName"


if saveFig
    savefig(p, imagePath)
else
    display(p)
end



p = getPlots(m, x_sol, x_sol_opt, method, name)

imageName = "pointsmetric_$(name)_$(method).svg"
imagePath = "$path$imageName"


if saveFig
    savefig(p, imagePath)
else
    display(p)
end





p = PlotProjectedPoints(projected_x, projected_y, projected_x_opt, projected_y_opt, boundarySection[1, :], boundarySection[2, :], m_vals_section, method, name)

imageName = "result_$(name)_$(method).svg"
imagePath = "$path$imageName"


if saveFig
    savefig(p, imagePath)
else
    display(p)
end
end