using DelimitedFiles

include("SetupDomain.jl")
include("GetBoundary.jl")

function GetAirfoilSetup(; airfoilPath = "examples/airfoil/A-airfoil.txt", radius = 3, cutN = 100, type =:cgrid)
    """
    Get the initial grid for the airfoil.
    """

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