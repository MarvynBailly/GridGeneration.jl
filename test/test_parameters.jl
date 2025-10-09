using Test
using GridGeneration

@testset "Configuration Types" begin
    @testset "SimParams Construction" begin
        # Test default construction
        params = SimParams()
        @test params isa SimParams
        @test params.useSplitting == true
        @test params.boundarySolver == :none
        
        # Test with custom values
        params2 = SimParams(
            useSplitting = false,
            boundarySolver = :analytic,
            useSmoothing = true,
            smoothMethod = :ellipticSS
        )
        @test params2.useSplitting == false
        @test params2.boundarySolver == :analytic
        @test params2.useSmoothing == true
        @test params2.smoothMethod == :ellipticSS
    end

    @testset "EllipticParams Construction" begin
        # Test default construction
        params = EllipticParams()
        @test params isa EllipticParams
        @test params.max_iter == 5000
        @test params.tol == 1e-6
        @test params.ω == 0.2
        
        # Test with custom values
        params2 = EllipticParams(
            max_iter = 1000,
            tol = 1e-5,
            ω = 0.3,
            verbose = true
        )
        @test params2.max_iter == 1000
        @test params2.tol == 1e-5
        @test params2.ω == 0.3
        @test params2.verbose == true
    end
end
