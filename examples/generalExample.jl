#=
General example: build an initial grid, split it into blocks, redistribute the block
edges according to a metric, and smooth the result.

Run from the repository root with

    julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'   # once
    julia --project=examples examples/generalExample.jl

The final grid is plotted to examples/output/.
=#
using GridGeneration
using Plots

include(joinpath(@__DIR__, "airfoil", "airfoil.jl"))
include(joinpath(@__DIR__, "rectangle", "rectangle.jl"))
include(joinpath(@__DIR__, "..", "plotters", "plot_grid.jl"))


case  = :airfoil  # :rectangle

##########################################
##########################################
#                 Domain                 #
##########################################
##########################################
"""
Define the initial grid here:
- initial grid
- boundary conditions
- interfaces

Note that an initial grid can be generated using TFI(boundary)
where boundary is saved in [top, right, bottom, left] format.
"""

###################
#   Airfoil Grid  #
###################
if case == :airfoil
    initialGrid, bndInfo, interInfo = GetAirfoilSetup(radius = 3, type =:cgrid)
elseif case == :rectangle
    initDomain = GetRectangleDomain()
    initialGrid = TFI(initDomain)
    bndInfo = Any[]
    interInfo = Any[]
end


##########################################
##########################################
#                 Metric                 #
##########################################
##########################################
"""
Define the metric tensor field M(x, y) -> (M11, M22)

Can either load a metric from file or use a custom metric function.
"""

###################
#  Airfoil Metric #
###################
if case == :airfoil
    problem = 6
    M = GetAirfoilMetric(problem; scale = 0.05)
elseif case == :rectangle
    problem = 1
    M = GetRectangleMetric(problem; scale = 10000)
end


##############################################
##############################################
#              Parameters                    #
##############################################
##############################################

"""
Define parameters
"""
# define split locations using the indices of the initial grid
splitLocations::Vector{Vector{Int}} = [
    [ 300 , 400],                           # split along the x axis
    [ 30 ]                                  # split along the y axis
]


# set up the parameters
params = SimParams(
    useSplitting = true,
    splitLocations = splitLocations,
    useEdgeSolver = true,
    boundarySolver = :analytic,     # :numeric
    useSmoothing = true,
    smoothMethod = :ellipticSS,
    elliptic = EllipticParams(
        max_iter = 5000,
        tol = 1e-6,
        ω = 0.2,
        useTopWall = true,
        useBottomWall = true,
        useLeftWall = true,
        useRightWall = true,
        a_decay_left = 0.4, b_decay_left = 0.4,
        a_decay_right = 0.4, b_decay_right = 0.4,
        a_decay_top = 0.4, b_decay_top = 0.4,
        a_decay_bottom = 0.4, b_decay_bottom = 0.4,
        verbose = false
    )
)


##############################################
##############################################
#              Run the Method                #
##############################################
##############################################

smoothBlocks, blocks, bndInfo, interInfo, finalErrors, finalIterations = GenerateGrid(initialGrid, bndInfo, interInfo, M, params=params)


##############################################
##############################################
#              Plot the Result               #
##############################################
##############################################

outdir = mkpath(joinpath(@__DIR__, "output"))
p = plot(
    plot_blocks(blocks, "After edge solve"),
    plot_blocks(smoothBlocks, "After smoothing"),
    layout = (1, 2), size = (1400, 600)
)
savefig(p, joinpath(outdir, "generalExample_$(case).png"))
println("Saved plot to ", joinpath(outdir, "generalExample_$(case).png"))
