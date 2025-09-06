include("../../src/GridGeneration.jl")

include("airfoil/GetAirfoilGrid.jl")
include("airfoil/GetBoundary.jl")
include("airfoil/metric/Metric.jl")
include("airfoil/metric/CustomMetric.jl")

include("../../plotter/metric_grid_plotter.jl")
include("../../plotter/blocks_interfaces_boundaries.jl")


using MAT


#################
# Load in initial single block grid
#################

# load in the initial grid - no trailing edge 
initialGrid = GetAirfoilGrid(airfoilPath="testing/single_block_ns/airfoil/data/A-airfoil.txt", radius = 2)
# throw away trailing edge stuff
airfoilGrid = initialGrid[:, 101:end-100, :]
airfoil = airfoilGrid[:,:,1]

# define the boundary information
bndInfo = getBoundaryConditions(airfoilGrid)

# define interInfo
interInfo = Any[]

#################
# Make a custom metric
#################

metricFunc1 = make_getMetric(airfoil;
    A_airfoil = 50.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
    A_origin  = 500.0,  ℓ_origin  = 0.1, p_origin  = 10,   
    floor     = 1e-4,  origin_center=(1, -0.6),
profile   = :rational)  # or :gauss

metricFunc2 = make_getMetric(airfoil;
    A_airfoil = 0.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
    A_origin  = 700.0,  ℓ_origin  = 0.1, p_origin  = 10,   
    floor     = 1e-4,  origin_center=(0.5, 0.1),
profile   = :rational)  # or :gauss

metricFunc = (x,y) -> metricFunc1(x,y) .+ metricFunc2(x,y)



#################
# Add split locations 
#################


splitLocations = [
    [  ],   # split along the x axis
    [  ]      # split along the y axis
]


plt1 = plot_blocks_interfaces_boundaries([airfoilGrid], interInfo, bndInfo;
    grid_stride= 1,
    show_block_ids=true,
    legend=false,
    titleName="Before Splitting",
    boundary_colors=Dict("BCWall"=>:purple, "BCInflow"=>:green, "BCOutflow"=>:blue) # optional
)

#################
## Split the blocks
#################

blocks, bndInfo, interInfo = GridGeneration.SplitBlock(airfoilGrid, splitLocations, bndInfo, interInfo)

plt2 = plot_blocks_interfaces_boundaries(blocks, interInfo, bndInfo;
    grid_stride= 1,
    show_block_ids=true,
    legend=false,
    titleName="After Splitting",
    boundary_colors=Dict("BCWall"=>:purple, "BCInflow"=>:green, "BCOutflow"=>:blue) # optional
)

#################
##### Run the solver on them
#################

# blocks, bndInfo, interInfo = GridGeneration.SolveAllBlocks(metricFunc, blocks, bndInfo, interInfo; solver="analytic")
# blocks, bndInfo, interInfo = GridGeneration.SolveAllBlocks(metricFunc, blocks, bndInfo, interInfo; solver="2ndorder")


plt3 = plot_blocks_interfaces_boundaries(blocks, interInfo, bndInfo;
    grid_stride= 1,
    show_block_ids=true,
    legend=false,
    titleName="After Solving",
    boundary_colors=Dict("BCWall"=>:purple, "BCInflow"=>:green, "BCOutflow"=>:blue) # optional
)


#################
##### Plotting stuff
#################

p3 = plot(plt1, plt2, plt3, layout=@layout([a ; b ; c]), size=(500, 1000),
    xlabel="x", ylabel="y",
    aspect_ratio=:equal)



w = (x,y) -> first(metricFunc(x,y))   

xl,xr = -2, 2 
yl, yr = -2, 2
xs = range(xl, xr, length=450)
ys = range(yl, yr, length=300)

plt4, _ = plot_scalar_field(w, xs, ys; boundary=nothing,
                         title="Distance-based symmetric isotropic metric",
                         cb_label="M(x,y)", equal_aspect=true, colormap =cgrad([RGB(1,1,1), RGB(1,0,0)]))#colormap=cgrd(:imola, rev=true))

plot!(plt4, xlims=(xl, xr), ylims=(yl, yr))

for block in blocks
    X = block[1, :, :]
    Y = block[2, :, :]

    for j in 1:size(X, 1)
        plot!(plt4, X[j, :], Y[j, :], color=:black, lw=0.8, label=false, alpha=0.7)
    end

    for i in 1:size(X, 2)
        plot!(plt4, X[:, i], Y[:, i], color=:black, lw=0.8, label=false, alpha=0.7)
    end
end


plt5 = plot(title="Original Airfoil Grid",
    xlabel="x", ylabel="y",
    aspect_ratio=:equal)
plot!(plt5, airfoilGrid[1, :, :], airfoilGrid[2, :, :], color=:black, lw=0.8, label=false, alpha=0.7)

plt6 = plot(title="Final Grid", xlabel="x", ylabel="y",
    aspect_ratio=:equal)

for block in blocks
    X = block[1, :, :]
    Y = block[2, :, :]

    for j in 1:size(X, 1)
        plot!(plt6, X[j, :], Y[j, :], color=:black, lw=0.6, label=false, alpha=0.7)
    end

    for i in 1:size(X, 2)
        plot!(plt6, X[:, i], Y[:, i], color=:black, lw=0.6, label=false, alpha=0.7)
    end
end

p4 = plot(plt5, plt4, plt6, layout=@layout([a ; b;  c]), size=(500, 1000),
    xlabel="x", ylabel="y",
    aspect_ratio=:equal)


p5 = plot(p3, p4, layout=@layout([a b]), size=(1000, 1000))

display(p5)