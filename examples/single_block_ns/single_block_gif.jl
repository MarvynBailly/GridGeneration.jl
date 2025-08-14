include("../../src/GridGeneration.jl")

include("airfoil/GetAirfoilGrid.jl")
include("airfoil/GetBoundary.jl")
include("airfoil/metric/Metric.jl")
include("airfoil/metric/CustomMetric.jl")
include("../../plotter/metric_grid_plotter.jl")






function get_frame(origin, initialGrid)
    metricFunc = make_getMetric(airfoil;
    A_airfoil = 200.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
    A_origin  = 1000.0,  ℓ_origin  = 0.1, p_origin  = 10,   
    floor     = 1e-4,  origin_center=origin,
    profile   = :rational)  # or :gauss


    # define the boundary information
    bndInfo = getBoundaryConditions(initialGrid)

    # define interInfo
    interInfo = Any[]

    block, bndInfo, interInfo = GridGeneration.SolveBlock(initialGrid, bndInfo, interInfo, metricFunc; tfi_method="TFI")

    w = (x,y) -> first(metricFunc(x,y))   

    xl,xr = -2, 2 
    yl, yr = -2, 2
    xs = range(xl, xr, length=450)
    ys = range(yl, yr, length=300)

    p1, _ = plot_scalar_field(w, xs, ys; boundary=nothing,
                            title="Distance-based isotropic metric",
                            cb_label="M11(x,y)", equal_aspect=true, colormap =cgrad([RGB(1,1,1), RGB(1,0,0)]))#colormap=cgrd(:imola, rev=true))

    plot!(p1, xlims=(xl, xr), ylims=(yl, yr))

    X = block[1, :, :]
    Y = block[2, :, :]

    for j in 1:size(X, 1)
        plot!(p1, X[j, :], Y[j, :], color=:black, lw=0.8, label=false, alpha=0.7)
    end

    for i in 1:size(X, 2)
        plot!(p1, X[:, i], Y[:, i], color=:black, lw=0.8, label=false, alpha=0.7)
    end

    # add a zoomed in figure 
    p2 = deepcopy(p1)
    plot!(p2, xlims=(-0.15, 1.15), ylims=(-0.25, 0.25),
            title="Zoomed in view",
            xlabel="x", ylabel="y",
            cb_label="M11(x,y)", equal_aspect=true, colormap=:imola)

    for j in 1:size(X, 1)
        plot!(p2, X[j, :], Y[j, :], color=:black, lw=0.5, label=false, alpha=0.4)
    end

    for i in 1:size(X, 2)
        plot!(p2, X[:, i], Y[:, i], color=:black, lw=0.5, label=false, alpha=0.4)
    end

        
    scatter!(p2, X[:, :], Y[:, :], color=:blue, label="", markersize=1.5, 
            markerstrokewidth=0, shape=:diamond,)

    p3 = plot(p1, p2, layout=@layout([a b]), size=(1000, 500))
    return p3 
end









# load in the initial grid - no trailing edge 
initialGrid = GetAirfoilGrid(airfoilPath="examples/single_block_ns/airfoil/data/A-airfoil.txt", radius = 2)
# throw away trailing edge stuff
airfoilGrid = initialGrid[:, 101:end-100, :]
airfoil = airfoilGrid[:,:,1]


origins = [(airfoil[1,i], airfoil[2,i]) for i in 1:10:size(airfoil, 2)]

anim = @animate for (i, origin) in enumerate(origins)
    println("percentage: ", i / length(origins) * 100, "%")
    p = get_frame(origin, airfoilGrid)  # your function returns a Plots.Plot
    plot(p)                             # ensure p is the "current" plot for capture
end

# Write the GIF
gif(anim, "hotspot_along_airfoil_dense2.gif", fps=12)  # tweak fps to taste