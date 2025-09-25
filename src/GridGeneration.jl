module GridGeneration

include("numerics/NumSolver.jl")

include("metric/Metric.jl")
include("metric/CustomMetric.jl")

include("Interpolators/TFI.jl")
include("Interpolators/LinearInterpolator.jl")

include("projections/Get1DMetric.jl")
include("projections/1D2DFunctions.jl")

include("blocksplitting/GetNeighbors.jl")
include("blocksplitting/BlockFunctions.jl")
include("blocksplitting/BoundaryUtils.jl")
include("blocksplitting/InterfaceUtils.jl")
include("blocksplitting/SplittingFunctions.jl")
include("blocksplitting/SolveAllBlocks.jl")
end