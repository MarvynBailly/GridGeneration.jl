using DelimitedFiles

include("SetupDomain.jl")
include("GetBoundary.jl")

"""
    GetAirfoilSetup(; airfoilPath, radius = 3, cutN = 100, type = :cgrid)

Build the initial single-block grid around the airfoil, plus its boundary conditions.
Returns `(airfoilGrid, bndInfo, interInfo)`.
"""
function GetAirfoilSetup(; airfoilPath = joinpath(@__DIR__, "A-airfoil.txt"), radius = 3, cutN = 100, type =:cgrid)

    # read the airfoil data
    airfoilData = readdlm(airfoilPath, '\t', skipstart=1)

    # set up a c grid around the provided inner boundary
    boundary = SetupDomain(
        airfoilData, 
        radius, 
        cutN, 
        cutN;
        type = type
    )

    initialGrid = GridGeneration.TFI(boundary)

    airfoilGrid = initialGrid[:, 101:end-100,:]
    
    bndInfo = getBoundaryConditions(airfoilGrid)

    interInfo = Any[]

    return airfoilGrid, bndInfo, interInfo
end