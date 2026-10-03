#=
Template: regenerate an existing Tortuga ("Turtle") grid using a metric field file.

Fill in `metricFile` and `gridFile`, choose split locations, then run from the repository root:

    julia --project=examples examples/generalExample_blank.jl

The output grid is written as a .grid file next to this script (examples/output/).
=#
using GridGeneration

# Turtle grid/metric loading helpers (pure Julia; shared with the GUI)
include(joinpath(@__DIR__, "..", "gui", "setup", "domain_setup.jl"))

##############################################
##############################################
#                  Set Up                    #
##############################################
##############################################

metricFile = ""   # path to a Turtle .metric field file
gridFile = ""     # path to a Turtle .grid file
(isfile(metricFile) && isfile(gridFile)) || error("Set metricFile and gridFile at the top of $(@__FILE__)")

initialGrid, bndInfo, interInfo, M = setup_turtle_grid_domain(metricFile, gridFile)


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
#              Save the Grid                 #
##############################################
##############################################

# extrude the smoothed 2D grid into a thin 3D slab and write it as a Turtle grid file
extrusion_length = 0.1
k_layers = 20
mesh3D, bndInfo3D, interInfo3D = GridGeneration.convert_2D_to_3D(smoothBlocks, bndInfo, interInfo, extrusion_length, k_layers)

outfile = joinpath(mkpath(joinpath(@__DIR__, "output")), "regenerated.grid")
GridGeneration.write_turtle_grid(mesh3D, interInfo3D, bndInfo3D, outfile)
println("Wrote ", outfile)
