using Test
using GridGeneration
using LinearAlgebra

@testset "GridGeneration.jl" begin
    
    @testset "Module Structure" begin
        @testset "Exports" begin
            @test isdefined(GridGeneration, :GenerateGrid)
            @test isdefined(GridGeneration, :SimParams)
            @test isdefined(GridGeneration, :EllipticParams)
            @test isdefined(GridGeneration, :TFI)
            @test isdefined(GridGeneration, :TFI_2D_Hermite)
            @test isdefined(GridGeneration, :make_getMetric)
            @test isdefined(GridGeneration, :setup_metric_tree)
            @test isdefined(GridGeneration, :find_nearest_kd)
        end
    end
    
    include("test_interpolators.jl")
    include("test_projections.jl")
    include("test_parameters.jl")
    include("test_metrics.jl")
    include("test_numerics.jl")
    
end
