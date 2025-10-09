using Test
using GridGeneration

@testset "Transfinite Interpolation (TFI)" begin
    @testset "Unit Square - Basic TFI" begin
        # Create a unit square: [top, right, bottom, left]
        N = 5
        top = [range(0, 1, length=N) ones(N)]
        right = [ones(N) range(0, 1, length=N)]
        bottom = [range(0, 1, length=N) zeros(N)]
        left = [zeros(N) range(0, 1, length=N)]
        
        boundary = [top, right, bottom, left]
        grid = TFI(boundary)
        
        # Check grid dimensions: [coordinate, i, j]
        @test size(grid) == (2, N, N)
        
        # Check corner points
        @test grid[1, 1, 1] ≈ 0.0 atol=1e-10  # bottom-left x
        @test grid[2, 1, 1] ≈ 0.0 atol=1e-10  # bottom-left y
        @test grid[1, N, 1] ≈ 1.0 atol=1e-10  # bottom-right x
        @test grid[2, N, 1] ≈ 0.0 atol=1e-10  # bottom-right y
        @test grid[1, 1, N] ≈ 0.0 atol=1e-10  # top-left x
        @test grid[2, 1, N] ≈ 1.0 atol=1e-10  # top-left y
        @test grid[1, N, N] ≈ 1.0 atol=1e-10  # top-right x
        @test grid[2, N, N] ≈ 1.0 atol=1e-10  # top-right y
        
        # Check that interior points are reasonable
        @test all(0 .<= grid[1, :, :] .<= 1)  # x in [0,1]
        @test all(0 .<= grid[2, :, :] .<= 1)  # y in [0,1]
    end

    @testset "Rectangular Domain" begin
        # Create a 2×1 rectangle
        L, H = 2.0, 1.0
        N = 7
        top = [range(0, L, length=N) fill(H, N)]
        right = [fill(L, N) range(0, H, length=N)]
        bottom = [range(0, L, length=N) zeros(N)]
        left = [zeros(N) range(0, H, length=N)]
        
        boundary = [top, right, bottom, left]
        grid = TFI(boundary)
        
        @test size(grid) == (2, N, N)
        
        # Check corners
        @test grid[1, 1, 1] ≈ 0.0 atol=1e-10
        @test grid[2, 1, 1] ≈ 0.0 atol=1e-10
        @test grid[1, N, 1] ≈ L atol=1e-10
        @test grid[2, N, 1] ≈ 0.0 atol=1e-10
        @test grid[2, N, N] ≈ H atol=1e-10
    end

    @testset "TFI Boundary Preservation" begin
        # TFI should exactly preserve boundary points
        N = 6
        top = [range(0, 1, length=N) ones(N)]
        right = [ones(N) range(0, 1, length=N)]
        bottom = [range(0, 1, length=N) zeros(N)]
        left = [zeros(N) range(0, 1, length=N)]
        
        grid = TFI([top, right, bottom, left])
        
        # Check bottom boundary
        for i in 1:N
            @test grid[1, i, 1] ≈ bottom[i, 1] atol=1e-10
            @test grid[2, i, 1] ≈ bottom[i, 2] atol=1e-10
        end
        
        # Check top boundary
        for i in 1:N
            @test grid[1, i, N] ≈ top[i, 1] atol=1e-10
            @test grid[2, i, N] ≈ top[i, 2] atol=1e-10
        end
        
        # Check left boundary
        for j in 1:N
            @test grid[1, 1, j] ≈ left[j, 1] atol=1e-10
            @test grid[2, 1, j] ≈ left[j, 2] atol=1e-10
        end
        
        # Check right boundary
        for j in 1:N
            @test grid[1, N, j] ≈ right[j, 1] atol=1e-10
            @test grid[2, N, j] ≈ right[j, 2] atol=1e-10
        end
    end
end
