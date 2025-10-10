using Test
using GridGeneration

@testset "Module Exports" begin
    @test isdefined(GridGeneration, :GenerateGrid)
    @test isdefined(GridGeneration, :SimParams)
    @test isdefined(GridGeneration, :EllipticParams)
    @test isdefined(GridGeneration, :TFI)
    @test isdefined(GridGeneration, :make_getMetric)
    @test isdefined(GridGeneration, :setup_metric_tree)
    @test isdefined(GridGeneration, :find_nearest_kd)
end
