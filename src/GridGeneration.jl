module GridGeneration


include("types.jl")
export SimParams, EllipticParams


include("numerics/NumSolver.jl")
include("numerics/EllipticSolver.jl")

include("metric/Metric.jl")
include("metric/CustomMetric.jl")

include("interpolators/TFI.jl")
include("interpolators/LinearInterpolator.jl")

include("projections/Get1DMetric.jl")
include("projections/1D2DFunctions.jl")

include("blocksplitting/GetNeighbors.jl")
include("blocksplitting/BlockFunctions.jl")
include("blocksplitting/BoundaryUtils.jl")
include("blocksplitting/InterfaceUtils.jl")
include("blocksplitting/SplittingFunctions.jl")
include("blocksplitting/SolveAllBlocks.jl")

include("smoothing/SmoothBlocks.jl")



end