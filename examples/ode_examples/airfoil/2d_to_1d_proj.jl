include("../../../src/GridGeneration.jl")

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









#####################
# SET UP DOMAIN
#####################

initialGrid = GetAirfoilGrid(airfoilPath = "examples/airfoil/data/A-airfoil.txt", radius = 3)

bottom = initialGrid[:,:,1]
airfoil = bottom[:, 101:end-100]  
SectionIndices = 100:250
# SectionIndices = 1:length(airfoil[1, :])
boundarySection = airfoil[:, SectionIndices]

arclength = sum(sqrt.(diff(boundarySection[1, :]).^2 .+ diff(boundarySection[2, :]).^2))
N = length(boundarySection[1, :])




xs = Boundary2Dto1D(boundarySection)

# build the metric
saveFig = true
method = "local" # "local" or "nuclear"
folder = "Mapping2Dto1D"
path = "docs/src/assets/images/$folder/"

scale = 40000

problems = [3] # 1: x=0, 2: x=1, 3: uniform
names = ["x=0", "x=1", "uniform"]

for (name, problem) in zip(names, problems) 

M_func = (x,y) -> Metric(x, y, scale, problem)
M_u1_func = (x,y) -> MetricDerivative(x, y, scale, problem)

m = GridGeneration.Get1DMetric(boundarySection, M_func, method = method)
mx = GridGeneration.Get1DMetric(boundarySection, M_u1_func, method = method)

m_func = build_interps_linear(xs, m)
mx_func = build_interps_linear(xs, mx)

# pass the functions to the solver
sol = GridGeneration.SolveODE(m_func, mx_func, N, xs[1], xs[end])

sigma_opt = GridGeneration.ComputeOptimalSpacing(sol[1, :], m, xs)

N_opt = ceil(Int, 1 / sigma_opt)
@info("Optimal number of points: ", N_opt)

sol_opt = GridGeneration.SolveODE(m_func, mx_func, N_opt, xs[1], xs[end])

x_sol = sol[1, :]
x_sol_opt = sol_opt[1, :]

# scale the solutions back to the correct size
x_sol = x_sol * xs[end]
x_sol_opt = x_sol_opt * xs[end]

projected_x, projected_y = InterpolateBoundaryManually(x_sol, boundarySection)
projected_x_opt, projected_y_opt = InterpolateBoundaryManually(x_sol_opt, boundarySection)

p = PlotProjectedPoints(projected_x, projected_y, projected_x_opt, projected_y_opt, boundarySection[1, :], boundarySection[2, :], m_vals_section, method, name)

imageName = "result_$(name)_$(method).svg"
imagePath = "$path$imageName"

if saveFig
    savefig(p, imagePath)
else
    display(p)
end
end