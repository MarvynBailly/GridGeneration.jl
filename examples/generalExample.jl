using Plots
using DelimitedFiles
using MAT: matread

include("../src/GridGeneration.jl")
include("airfoil/data/GetAirfoilSetup.jl")
include("airfoil/metric/GetMetric.jl")

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
###################
#  Split Params  #
###################

# should the method split
useSplitting = true
# define split locations using the indices of the initial grid 
splitLocations::Vector{Vector{Int}} = [
    [ 300 , 400],                           # split along the x axis
    [ 30 ]                                  # split along the y axis
]

###################
#  Edge Solver Params  #
###################
# should the method solve along the edges
useEdgeSolver = true
# which edge solver to use
boundarySolver = :analytic  # :analytic or :numeric



###################
# Smooth Params  #
###################
# should the method smooth
useSmoothing = true

# which smoothing method to use
smoothMethod = :ellipticSS 

# solver parameters 
max_iter::Int=5000
tol::Float64=1e-6
ω::Float64=0.2

# verbose
ellipticSolverVerbose::Bool=false

#### SS parameters
# SS forcing terms
useLeftWall, useRightWall, useTopWall, useBottomWall = false, false, true, true

# SS decay parameters
a_decay_top::Float64=0.4
b_decay_top::Float64=0.4

a_decay_left::Float64=0.4
b_decay_left::Float64=0.4


a_decay_right::Float64=0.4
b_decay_right::Float64=0.4

a_decay_bottom::Float64=0.4
b_decay_bottom::Float64=0.4


# create parameter struct
# find a way to have this be automatic and just allow the user to change what they want
params = SimParams(
    useSplitting = useSplitting, 
    splitLocations = splitLocations, 
    useEdgeSolver = useEdgeSolver, 
    boundarySolver = boundarySolver,
    useSmoothing = useSmoothing, 
    smoothMethod = smoothMethod,
    elliptic = EllipticParams(
        max_iter = max_iter, tol = tol, ω = ω,
        useTopWall = useTopWall, useBottomWall = useBottomWall, useLeftWall = useLeftWall, useRightWall = useRightWall,
        a_decay_left = a_decay_left, b_decay_left = b_decay_left,
        a_decay_right = a_decay_right, b_decay_right = b_decay_right,
        a_decay_top = a_decay_top, b_decay_top = b_decay_top,
        a_decay_bottom = a_decay_bottom, b_decay_bottom = b_decay_bottom,
        verbose = ellipticSolverVerbose
    )
)






##############################################
##############################################
#              Run the Method                 #
##############################################
##############################################




function MAIN(initialGrid, bndInfo, interInfo, M; params=params)
    if params.useSplitting
        blocks, bndInfo, interInfo = GridGeneration.SplitBlock(initialGrid, params.splitLocations, bndInfo, interInfo)
    end
    if params.useEdgeSolver
        blocks, bndInfo, interInfo = GridGeneration.SolveAllBlocks(M, blocks, bndInfo, interInfo; solver=params.boundarySolver)
    end

    if params.useSmoothing
        smoothBlocks, finalErrors, finalIterations = GridGeneration.SmoothBlocks(blocks; solver=params.smoothMethod, params=params.elliptic)
    end

    return smoothBlocks, blocks, bndInfo, interInfo
end


MAIN(initialGrid, bndInfo, interInfo, M, params=params)