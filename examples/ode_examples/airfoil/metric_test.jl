include("../../../src/GridGeneration.jl")
include("GetAirfoilGrid.jl")
include("metric/Metric.jl")

using Plots

function getPlots(m_vals_section, x_sol, x_sol_opt, method, name)    
    diff = x_sol[2:end] - x_sol[1:end-1]
    diff_opt = x_sol_opt[2:end] - x_sol_opt[1:end-1]
    s = LinRange(0, 1, length(diff))
    s_opt = LinRange(0,1, length(diff_opt))

    p5 = scatter(x_sol, zeros(size(x_sol)) .+ 0.01, label="Non-opt Solution", color=:blue, linewidth=1.5,  title= "Grid Spacing with real metric")
    scatter!(p5, x_sol_opt, zeros(size(x_sol_opt)), label="opt Solution ($(length(x_sol_opt)))", color=:black, linewidth=1.5)
    scatter!(p5, s, diff, label="Non-opt Spacing", color=:red, markersize=2, markerstrokewidth=0)
    scatter!(p5, s_opt, diff_opt, label="opt Spacing", color=:green, markersize=2, markerstrokewidth=0)

    p6 = scatter(s, title = "M values", m_vals_section, label="M(s)", color=:black, markersize=2, markerstrokewidth=0)

        titlepanel = plot(title = "1D Metric (method = $method) with $name clustering",
                    framestyle = :none, grid = false, ticks = nothing)

    p7 = plot(titlepanel, p5, p6, layout = @layout([A{0.001h}; b c]), size = (1000, 600),)

    return p7
end


# build the metric
saveFig = true
method = "local" # "local" or "nuclear"

folder = "Mapping2Dto1D"
path = "docs/src/assets/images/$folder/"

scale = 40000

problems = [2] # 1: x=0, 2: x=1, 3: uniform
names = ["x=0", "x=1", "uniform"]


initialGrid = GetAirfoilGrid(airfoilPath = "examples/airfoil/data/A-airfoil.txt", radius = 3)

bottom = initialGrid[:,:,1]
airfoil = bottom[:, 101:end-100]  
SectionIndices = 100:10:250
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


for (problem, name) in zip(problems, names)
    @info "Processing problem: $problem with name: $name"
    M_func = (x,y) -> Metric(x, y, scale, problem)
    
    m = GridGeneration.Get1DMetric(boundarySection, M_func, method = method)
    
    m_func = GridGeneration.build_interps_linear(xs, m)
    
    pm = plot(xs, M_func(boundarySection[1, :], boundarySection[2, :])[1], markershape=:circle, label="M(x_i)", color=:red, linewidth=2, markerstrokewidth=0, xlabel="s", ylabel="M", title="Metric Values along Boundary (in s direction)")
    plot!(pm, xs, m, markershape=:circle, label="m(x_i)", color=:black, linewidth=1, markerstrokewidth=0, markersize=3.5, xlabel="s", ylabel="M")
    # xdense = LinRange(0, 1, 1000)
    # scatter!(pm, xdense, m_func.(xdense), label="m(s)", color=:blue, markershape=:circle, linewidth=1, markerstrokewidth=0, markersize=0.5)
    display(pm)
end