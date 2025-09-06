module GridGeneration

include("numerics/NumSolver.jl")
include("numerics/CentralDiff.jl")

include("SetupDomain.jl")
include("TFI.jl")
include("Get1DMetric.jl")
include("Interpolators.jl")
include("1D2DFunctions.jl")
include("Metric.jl")
include("BlockFunctions.jl")
include("BoundaryUtils.jl")
include("InterfaceUtils.jl")
include("SplittingFunctions.jl")
include("SolveAllBlocks.jl")


end