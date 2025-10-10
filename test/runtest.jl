using Test
using GridGeneration

@testset "GridGeneration.jl" begin
    
    # Test modules organized by functionality
    include("test_exports.jl")
    include("test_parameters.jl")
    include("test_interpolators.jl")
    include("test_projections.jl")
    # include("test_numerics.jl")
    include("test_blocks.jl")
    include("test_multiblock_splitting.jl")
    # include("test_integration.jl")
end
