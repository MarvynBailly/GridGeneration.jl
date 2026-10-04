#=
Template: regenerate an existing Tortuga ("Turtle") grid using a metric field file.

Fill in `metricFile` and `gridFile`, choose split locations, then run from the repository root:

    julia --project=examples examples/generalExample_blank.jl

The output grid is written as a .grid file next to this script (examples/output/).
=#
using GridGeneration

##############################################
##############################################
#                  Set Up                    #
##############################################
##############################################

metricFile = ""   # path to a Turtle .metric field file
gridFile = ""     # path to the Turtle .grid file to regenerate
# Grid the metric was computed on. `nothing` uses the grid named in the metric file's header
# (if it sits next to metricFile), otherwise gridFile.
metricGridFile = nothing
(isfile(metricFile) && isfile(gridFile)) || error("Set metricFile and gridFile at the top of $(@__FILE__)")

initialGrid, bndInfo, interInfo, M = setup_turtle_grid_domain(metricFile, gridFile; metricGridFile = metricGridFile)


# split requests: (blockId, [[i_splits...], [j_splits...]]) using the indices of each block;
# splits propagate across interfaces into neighbouring blocks automatically
splitRequests = [
    (1, [[20], [30]]),
]

# set up the parameters
params = SimParams(
    useSplitting = true,
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

smoothBlocks, blocks, bndInfo, interInfo, finalErrors, finalIterations = GenerateGrid(initialGrid, bndInfo, interInfo, M; params = params, splitRequests = splitRequests)


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
