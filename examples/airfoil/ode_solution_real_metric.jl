include("../../src/GridGeneration.jl")
include("GetAirfoilGrid.jl")

using Plots, MAT, DelimitedFiles, NearestNeighbors


function setup_metric_tree(data, refs)
    coords = Float64[]
    for (b, (Xb,Yb)) in enumerate(zip(data["x"], data["y"]))
    Ny, Nz = size(Xb)
    for i in 1:Ny, j in 1:Nz
        push!(coords, Xb[i,j])
        push!(coords, Yb[i,j])
        push!(refs, (b,i,j))
    end
    end
    coords = reshape(coords, (2, length(refs)))  # now 2×N_points
    
    tree = KDTree(coords)
    return tree, refs
end


function find_nearest_kd(data, tree, refs, xq, yq)
    idxs, dists = knn(tree, [xq,yq], 1)   # 1‐NN
    ref = refs[idxs[1]]
    blk = ref.block
    i,j = ref.i, ref.j
    M11_val = data["M11"][blk][i, j]
    M22_val = data["M22"][blk][i, j]
    return M11_val, M22_val
end


#####################
# SET UP DOMAIN
#####################
initialGrid = GetAirfoilGrid(airfoilPath = "examples/airfoil/data/A-airfoil.txt", radius = 3)

bottom = initialGrid[:,:,1]
airfoil = bottom[:, 101:end-100]
sectionIndices = 1:length(airfoil[1, :])
boundarySection = airfoil[:, sectionIndices]


# set up real metric data
metric_path = "examples/airfoil/metric/A-airfoil_grid_data.mat"
# load metric
metric_data = matread(metric_path)
# set up metric tree for fast nearest neighbor search - this is used within metric
tree, refs = setup_metric_tree(metric_data, [])

# loop through points, use tree to get nearest neighbors, store neighbors, (interpolate?), and then do the real thing. Plot results and put up.

# start block stufff!



# N = length(boundarySection[1, :])

# xs = GridGeneration.ProjectBoundary2Dto1D(boundarySection)

# # build the metric
# saveFig = false
# folder = "PointProjection"
# path = "docs/src/assets/images/$folder/"


# M_func_test = (x,y) -> Metric(x, y, scale, problem)

# m = GridGeneration.Get1DMetric(boundarySection, M_func_test)

# sol_opt, sol = GridGeneration.GetOptimalSolution(m, m, N, xs)

# x_sol, x_sol_opt = sol[1, :], sol_opt[1, :]


# projected_points = GridGeneration.ProjectBoundary1Dto2D(boundarySection, x_sol, xs)
# projected_points_opt = GridGeneration.ProjectBoundary1Dto2D(boundarySection, x_sol_opt, xs)

# projected_x, projected_y = projected_points[1, :], projected_points[2, :]
# projected_x_opt, projected_y_opt = projected_points_opt[1, :], projected_points_opt[2, :]



# p1 = plot(xs, m, title = "1D Metric with Boundary Spacing (method = $method, $name clustering)",
#         xlabel = "s", ylabel = "m(x(s))", label = "m(x(s))",
#         legend = :best, linewidth=2)
# scatter!(p1, xs, zeros(length(xs)), markershape=:circle, markersize=2, markerstrokewidth=0, c = :black, label="1D Boundary Values")
# scatter!(p1, xs, m, markershape=:diamond, markersize=2, markerstrokewidth=0, c = :black, label="2D Boundary Values")
# # add red arrow to show the direction of xs
# quiver!(p1, [xs[1]], [m[1]], quiver=([xs[5] - xs[1]], [m[5] - m[1]]), color=:red, label="s direction", lw=1)



# p2 = plot(title = "Optimal and Non-Optimal 1D Distributions (method = $method, $name clustering)",
#         xlabel = "s", ylabel = "m(x(s))", label = "m(x(s))",
#         legend = :best, linewidth=2)
# scatter!(p2, xs, zeros(length(xs)), markershape=:circle, markersize=3, markerstrokewidth=0, c = :black, label="1D Boundary Values")
# scatter!(p2, sol[1,:], zeros(length(sol[1,:])), markershape=:circle, markersize=2.5, markerstrokewidth=0, c = :red, label="Non-Optimal Solution (N = $(length(sol[1,:])))")
# scatter!(p2, sol_opt[1,:], zeros(length(sol_opt[1,:])), markershape=:circle, markersize=2, markerstrokewidth=0, c = :green, label="Optimal Solution (N = $(length(sol_opt[1,:])))")

# p3 = plot( title = "Optimal Points Projected on Boundary", aspect_ratio = 1)

# plot!(p3, boundarySection[1, :], boundarySection[2, :], label="Boundary", color=:black, linewidth=1)
# scatter!(p3, boundarySection[1, :], boundarySection[2, :], label="1D Metric", marker_z = m, markersize=4, markerstrokewidth=0)
# scatter!(p3, projected_x_opt, projected_y_opt, label="Projected Points Optimal ($(length(projected_x_opt)))", color = :green, marker=:circle, markersize=2, markerstrokewidth=0)
# # scatter!(p3, boundarySection[1, :], boundarySection[2, :], label="Discrete Metric Values ($(length(boundarySection[1,:])))", marker_z = m_vals_section, markersize=2, markerstrokewidth=0, c = :brg)

# p5 = plot()
# plot!(p5, boundarySection[1, :], boundarySection[2, :], label="Boundary", color=:black, linewidth=2)
# plot!(p5, projected_x_opt,projected_y_opt, label="Projected Solution", color=:red, linewidth=1)
# scatter!(p5, projected_x_opt, projected_y_opt, label="Projected Points Optimal ($(length(projected_x_opt)))", color = :green, marker=:circle, markersize=2, markerstrokewidth=0)

# # p6 = plot(xlims=[-0.05,0.05], ylims=[-0.1, 0.1])
# # plot!(p6, boundarySection[1, :], boundarySection[2, :], label="Boundary", color=:black, linewidth=2, title="Projected Points")
# # plot!(p6, projected_x_opt,projected_y_opt, label="Projected Solution", color=:red, linewidth=1, title="leading edge zoom in")

# # p7 = plot(p5, p6, layout = @layout([a b]), size = (800, 400))
# # display(p7)
# # scatter!(p6, projected_x_opt, projected_y_opt, label="Projected Points Optimal ($(length(projected_x_opt)))", color = :green, marker=:circle, markersize=2, markerstrokewidth=0)


# p4 = plot(p1, p2, p3, p5, layout = @layout([a ; b ; c ; d]), size = (800, 1000))


# if saveFig
#     savefig(p4, path * "ode_solution_$(name)_$(method).svg")
# else
#     display(p4)
# end