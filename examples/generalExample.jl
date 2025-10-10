using Plots
using DelimitedFiles
using MAT: matread

include("../src/GridGeneration.jl")
using .GridGeneration

include("airfoil/airfoil.jl")
include("rectangle/rectangle.jl")


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
    initialGrid, bndInfo, interInfo = GetAirfoilSetup(airfoilPath = "examples/airfoil/data/A-airfoil.txt", radius = 3, type =:cgrid)
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
Define the metric tensor field M 

Can either provide the filepath or use custom metric function.
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
