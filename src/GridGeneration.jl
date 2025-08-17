module GridGeneration

include("numerics/NumSolver.jl")
include("SetupDomain.jl")
include("TFI.jl")
include("Get1DMetric.jl")
include("Interpolators.jl")
include("1D2DFuctions.jl")
include("Metric.jl")
include("BlockFunctions.jl")
include("BoundaryUtils.jl")
include("InterfaceUtils.jl")
include("SplittingFunctions.jl")
include("SolveAllBlocks.jl")


end