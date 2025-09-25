include("airfoil/AirfoilExample.jl")
include("plotter/metric_grid_plotter.jl")
include("plotter/angle_deviation.jl")


##############################################
################ PARAM STRUCT ################
##############################################

struct SimParams
    initialGrid::Array{Float64,3}
    bndInfo::Array{Any,1}
    interInfo::Array{Any,1}
    splitLocations::Vector{Vector{Int}}
    M::Function
    max_iter::Int
    tol::Float64
    ω::Float64
    use_top_wall::Bool
    s_top::Float64
    a_decay_top::Float64
    b_decay_top::Float64
    use_left_wall::Bool
    s_left::Float64
    a_decay_left::Float64
    b_decay_left::Float64
    use_right_wall::Bool
    s_right::Float64
    a_decay_right::Float64
    b_decay_right::Float64
    use_bottom_wall::Bool
    s_bottom::Float64
    a_decay_bottom::Float64
    b_decay_bottom::Float64
    EllipticSolver_verbose::Bool
    showPlots::Bool
end


##############################################
################## BOUNDARY ##################
##############################################

airfoilOuterRadius::Float64 = 0.8
initialGrid = GetAirfoilGrid(airfoilPath="examples/airfoil/data/A-airfoil.txt", radius = airfoilOuterRadius)

# throw away trailing edge stuff 
airfoilGrid = initialGrid[:, 101:end-100, :]
airfoil = airfoilGrid[:,:,1]

# define the boundary information
bndInfo = getBoundaryConditions(airfoilGrid)

# define interInfo
interInfo = Any[]

##############################################
################### METRIC ###################
##############################################
# define the metric using either a custom function or load from file
function GetMetric(problem)
    if problem == 1
        metricFunc1 = GridGeneration.make_getMetric(airfoil;
        A_airfoil = 10000.0,  ℓ_airfoil = 0.1, p_airfoil = 5,   
        A_origin  = 1000.0,  ℓ_origin  = 0.25, p_origin  = 6,   
        floor     = 1000,  origin_center=(-0.25, 0.0),
        profile   = :rational)
        
        M = (x,y) -> metricFunc1(x, y)
    elseif problem == 2
        metricFunc1 = GridGeneration.make_getMetric(airfoil;
        A_airfoil = 10000.0,  ℓ_airfoil = 0.1, p_airfoil = 5,   
        A_origin  = 5000.0,  ℓ_origin  = 0.25, p_origin  = 6,   
        floor     = 1e-4,  origin_center=(-0.25, 0.0),
        profile   = :rational)
        
        M = (x,y) -> metricFunc1(x, y)
    elseif problem == 3
        metricFunc1 = GridGeneration.make_getMetric(airfoil;
        A_airfoil = 500.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
        A_origin  = 5000.0,  ℓ_origin  = 0.2, p_origin  = 10,   
        floor     = 1e-4,  origin_center=(0.0, 0.0),
        profile   = :rational)

        metricFunc2 = GridGeneration.make_getMetric(airfoil;
        A_airfoil = 500.0,  ℓ_airfoil = 0.5, p_airfoil = 2,   
        A_origin  = 5000.0,  ℓ_origin  = 0.2, p_origin  = 10,   
        floor     = 1e-4,  origin_center=(1.0, 0.0),
        profile   = :rational)

        M = (x,y) -> metricFunc1(x, y) .+ metricFunc2(x, y)
    elseif problem == 4
        # set up metric data
        metricPath = "examples/airfoil/metric/A-airfoil_grid_data.mat"
        # load metric
        metricData = matread(metricPath)
        # set up metric function
        tree, refs = GridGeneration.setup_metric_tree(metricData)
        M = (x,y) -> 0.05 * GridGeneration.find_nearest_kd(metricData, tree, refs, x, y)
    end
    return M
end

problem::Int = 4
M::Function = GetMetric(problem)  # M(x,y)


##############################################
######### BLOCK SPLITTING PARAMETERS #########
##############################################
splitLocations::Vector{Vector{Int}} = [
    [ 300 , 400],     # split along the x axis
    [ 30]      # split along the y axis
]



##############################################
############ Elliptic Parameters #############
##############################################
max_iter::Int=5000
tol::Float64=1e-6
ω::Float64=0.2

# note user defined `s_side` is being ignored in favor of automatic calculation based on grid spacing

use_top_wall::Bool=true 
s_top::Float64=-1
a_decay_top::Float64=0.4
b_decay_top::Float64=0.4

use_left_wall::Bool=true 
s_left::Float64=-1
a_decay_left::Float64=0.4
b_decay_left::Float64=0.4


use_right_wall::Bool=true
s_right::Float64=-1
a_decay_right::Float64=0.4
b_decay_right::Float64=0.4

use_bottom_wall::Bool=true 
s_bottom::Float64=-1
a_decay_bottom::Float64=0.4
b_decay_bottom::Float64=0.4

EllipticSolver_verbose::Bool=false

showPlots::Bool=false


##############################################
################### PARAMS ###################
##############################################

SimParams(;
    initialGrid = airfoilGrid, 
    bndInfo = bndInfo, 
    interInfo = interInfo, 
    splitLocations = splitLocations, 
    M = M,
    max_iter = max_iter,
    tol = tol,
    ω = ω,
    use_top_wall = use_top_wall,
    s_top = s_top,
    a_decay_top = a_decay_top,
    b_decay_top = b_decay_top,
    use_left_wall = use_left_wall,
    s_left = s_left,
    a_decay_left = a_decay_left,
    b_decay_left = b_decay_left,
    use_right_wall = use_right_wall,
    s_right = s_right,
    a_decay_right = a_decay_right,
    b_decay_right = b_decay_right,
    use_bottom_wall = use_bottom_wall,
    s_bottom = s_bottom,
    a_decay_bottom = a_decay_bottom,
    b_decay_bottom = b_decay_bottom,
    EllipticSolver_verbose = EllipticSolver_verbose,
    showPlots = showPlots
) = SimParams(initialGrid, bndInfo, interInfo, splitLocations, M, max_iter, tol, ω, use_top_wall, s_top, a_decay_top, b_decay_top,
    use_left_wall, s_left, a_decay_left, b_decay_left,
    use_right_wall, s_right, a_decay_right, b_decay_right,
    use_bottom_wall, s_bottom, a_decay_bottom, b_decay_bottom,
    EllipticSolver_verbose, showPlots)




############## POST PROCESSING ##############


#### Post Processing
function ComputeAngleDeviation(blocks)
    angleDeviations = []
    for block in blocks
        x = block[1, :, :]
        y = block[2, :, :]
        nrows, ncols = size(x)
        
        # Initialize a matrix to store the deviation for each cell corner.
        # The result will be (nrows-1, ncols-1), so the last row/col will be 0.
        deviations = zeros(nrows, ncols)

        # --- Corrected Vector Calculation ---
        
        # Vector 1: Along the xi-direction (horizontal grid lines)
        # For each cell (i,j), this is the vector from point (i,j) to (i,j+1)
        v1_x = x[1:end-1, 2:end]   - x[1:end-1, 1:end-1]
        v1_y = y[1:end-1, 2:end]   - y[1:end-1, 1:end-1]

        # Vector 2: Along the eta-direction (vertical grid lines)
        # For each cell (i,j), this is the vector from point (i,j) to (i+1,j)
        v2_x = x[2:end,   1:end-1] - x[1:end-1, 1:end-1]
        v2_y = y[2:end,   1:end-1] - y[1:end-1, 1:end-1]

        # --- Vectorized Angle Calculation ---
        
        # Dot product: v1 ⋅ v2
        dot_product = v1_x .* v2_x .+ v1_y .* v2_y

        # Magnitude of v1: ||v1||
        norm1 = sqrt.(v1_x.^2 .+ v1_y.^2)
        
        # Magnitude of v2: ||v2||
        norm2 = sqrt.(v2_x.^2 .+ v2_y.^2)

        # Angle from dot product formula: θ = acos((v1 ⋅ v2) / (||v1|| * ||v2||))
        epsilon = 1e-12
        theta = acosd.(dot_product ./ (norm1 .* norm2))

        # The deviation is the absolute difference from 90 degrees.
        deviations[1:end-1, 1:end-1] = abs.(theta .- 90.0)
        
        push!(angleDeviations, deviations)
    end
    return angleDeviations
end










##############################################
#################### Main ####################
##############################################
runSolver = true
saveplots = false
case = 4

@info "runSolver = $runSolver"

if case == 1
    params = SimParams(use_left_wall = false,
                use_right_wall = false,
                use_bottom_wall = false,
                use_top_wall = false)
elseif case == 2
    params = SimParams(use_left_wall = false,
                use_right_wall = false,
                use_bottom_wall = true,
                use_top_wall = false)
elseif case == 3
    params = SimParams(use_left_wall = false,
                use_right_wall = false,
                use_bottom_wall = true,
                use_top_wall = true)
elseif case == 4
    params = SimParams(use_left_wall = true,
                use_right_wall = true,
                use_bottom_wall = true,
                use_top_wall = true)
end



if runSolver
    blocks_smooth, blocks = AirfoilSolver(params)
end

##############################################
################### PLOTS ####################
##############################################
function plot_grid(x, y, title_str; plt = nothing, c = nothing)
    p = plt == nothing ? plot() : plt
    plot!(p, title=title_str, aspect_ratio=:equal, legend=false, framestyle=:box, grid = false)
    # Plot η-lines (lines of constant j)
    for j in 1:size(x, 2)
        plot!(p, x[:, j], y[:, j], lw=0.1, color=c)
    end
    # Plot ξ-lines (lines of constant i)
    for i in 1:2:size(x, 1)
        plot!(p, x[i, :], y[i, :], lw=0.1, color=c)
    end
    return p
end

pInit = plot_grid(airfoilGrid[1, :, :], airfoilGrid[2, :, :], "Initial Grid", c = RGB(0.0, 0.0, 0.0))

pSmooth = plot()

for i in 1:length(blocks)
    plot_grid(blocks_smooth[i][1, :, :], blocks_smooth[i][2, :, :], "Block Elliptic Grid Case $case", plt=pSmooth, c = RGB(0.0, 0.0, 0.0))
end

finalPlot = plot(pInit, pSmooth, layout=(1,2), size=(1200,600))
display(finalPlot)

if saveplots
    savefig(finalPlot, "images/bs_el_$case.pdf")
end


# savefig(finalPlot, "bs_el_$case.svg")


w = (x,y) -> first(M(x,y))   
xl,xr = -1, 1 
yl, yr = -1, 1
xs = range(xl, xr, length=450)
ys = range(yl, yr, length=300)

plt4, _ = plot_scalar_field(w, xs, ys; boundary=airfoil,
                         title="0.05*M(x,y)",
                         cb_label="M(x,y)", equal_aspect=true, colormap =cgrad([RGB(1,1,1), RGB(1,0,0)]))#colormap=cgrd(:imola, rev=true))

plot!(plt4, xlims=(xl, xr), ylims=(yl, yr))
if saveplots
    savefig(plt4, "images/metric_05.pdf")
end



postProcessing = true

@info "postProcessing = $postProcessing"
if postProcessing
    angleDeviations = ComputeAngleDeviation(blocks_smooth)
    qualityPlot = PlotGridAngleDeviation(blocks_smooth, angleDeviations)
    display(qualityPlot)
end