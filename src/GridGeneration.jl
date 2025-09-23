module GridGeneration

include("numerics/NumSolver.jl")

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
include("CustomMetric.jl")


end