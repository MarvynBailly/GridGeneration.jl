using Test
using GridGeneration

@testset "Module Exports" begin
    @test isdefined(GridGeneration, :GenerateGrid)
    @test isdefined(GridGeneration, :SimParams)
    @test isdefined(GridGeneration, :EllipticParams)
    @test isdefined(GridGeneration, :TFI)
    @test isdefined(GridGeneration, :SplitMultiBlock)
    for name in (:SolveAllBlocks, :SmoothBlocks, :ImportTurtleGrid, :readTurtleFields,
                 :convert_2D_to_3D, :write_turtle_grid)
        @test name in names(GridGeneration)
    end
    @test isdefined(GridGeneration, :make_getMetric)
    @test isdefined(GridGeneration, :setup_metric_tree)
    @test isdefined(GridGeneration, :find_nearest_kd)
end
