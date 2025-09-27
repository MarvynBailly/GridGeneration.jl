using Plots
using DelimitedFiles
using MAT: matread

include("../src/GridGeneration.jl")
include("airfoil/data/GetAirfoilSetup.jl")
include("airfoil/metric/GetMetric.jl")
include("plotter/plot_grid.jl")

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

Note that an initial grid can be generated using GridGeneration.TFI(boundary) 
where boundary is saved in [top, right, bottom, left] format. 
"""

###################
#   Airfoil Grid  #
###################
initialGrid, bndInfo, interInfo = GetAirfoilSetup(airfoilPath = "examples/airfoil/data/A-airfoil.txt", radius = 3, type =:cgrid)


##########################################
##########################################
#                 Metric                 #
##########################################
##########################################
"""
Define the metric tensor field M 

Can either provide the filepath or use custom metric function.
"""

###################
#  Airfoil Metric #
###################
problem = 6
M = GetMetric(problem; scale = 0.05)



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
params = GridGeneration.SimParams(
    useSplitting = true, 
    splitLocations = splitLocations, 
    useEdgeSolver = true, 
    boundarySolver = :analytic,     # :numeric
    useSmoothing = true,
    smoothMethod = :ellipticSS,
    elliptic = GridGeneration.EllipticParams(
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

smoothBlocks, blocks, bndInfo, interInfo, finalErrors, finalIterations = GridGeneration.GenerateGrid(initialGrid, bndInfo, interInfo, M, params=params)




# pSmooth = plot()

# for i in 1:length(blocks)
#     plot_grid(smoothBlocks[i][1, :, :], smoothBlocks[i][2, :, :], "Block Elliptic Grid Case", plt=pSmooth, c = RGB(0.0, 0.0, 0.0))
# end

# pZoom = deepcopy(pSmooth)
# plot!(pZoom, xlims=(0.5, 1.5), ylims=(-0.3, 0.3))
# p = plot(pSmooth, pZoom, layout = (1,2), size=(1200, 600))
# savefig(p, "examples/generalExample.pdf")