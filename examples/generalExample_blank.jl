using GridGeneration

##############################################
##############################################
#                  Set Up                    #
##############################################
##############################################

metricFile = "" 
gridFile = "" 
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

# save the final grid to a turtle grid file
# extrusion
extrusion_length = 0.1
k_layers = 20
mesh3D, bndInfo3D, interfaceInfo3D = GridGeneration.convert_2D_to_3D(blocks, bndInfo, interfaceInfo, extrusion_length, k_layers)

# write .grid file
GridGeneration.write_turtle_grid(mesh3D, interfaceInfo3D, bndInfo3D, filename)