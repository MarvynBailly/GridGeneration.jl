using Test
using GridGeneration

@testset "Solvers" begin
    G = GridGeneration

    @testset "1D ODE - numeric agrees with analytic" begin
        xs = collect(range(0, 1, length=101))
        for m in (x -> 1 + 4x^2, x -> 1 + 50exp(-200(x - 0.5)^2))
            mf = G.LinearInterpolate(xs, m.(xs))
            ref = G.SolveODE(mf, xs; solver=:analytic)
            sol = G.SolveODE(mf, xs; solver=:numeric)

            # Both solvers return a plain vector of point locations
            @test sol isa AbstractVector
            @test length(sol) == length(ref)
            @test sol[1] ≈ xs[1]
            @test sol[end] ≈ xs[end]
            @test all(diff(sol) .> 0)
            @test maximum(abs.(sol .- ref)) < 1e-2
        end

        @test_throws ArgumentError G.SolveODE(x -> 1.0, xs; solver=:numerical)
        @test_throws ArgumentError G.SolveODEFixedN(x -> 1.0, xs, 10; solver=:numerical)
    end

    @testset "Optimal number of points" begin
        # Constant metric m on a unit-length edge: equidistributed spacing 1/sqrt(m)
        for m in (25.0, 100.0, 400.0)
            x = collect(range(0, 1, length=201))
            N = G.ComputeOptimalNumberofPoints(x, s -> m)
            @test N == floor(Int, sqrt(m))
        end
        # Never fewer than 3 points
        x = collect(range(0, 1, length=11))
        @test G.ComputeOptimalNumberofPoints(x, s -> 1e-4) == 3
    end

    @testset "Elliptic smoothing" begin
        # Wall-clustered grid whose first ξ-spacing at the left wall varies along η (corners = 0.02)
        Ni, Nj, h = 21, 21, 0.02
        x = zeros(Ni, Nj); y = zeros(Ni, Nj)
        for j in 1:Nj, i in 1:Ni
            η = (j - 1) / (Nj - 1); ξ = (i - 1) / (Ni - 1)
            a = h * (1 + 2sin(pi * η)) * (Ni - 1)
            x[i, j] = a * ξ + (1 - a) * ξ^2
            y[i, j] = η + 0.05sin(pi * ξ) * sin(pi * η)
        end

        p = EllipticParams(max_iter=20000, tol=1e-9, useBottomWall=false, useLeftWall=true)
        xr, yr, err, iters = G.EllipticSolver(copy(x), copy(y); params=p)

        @test err < 1e-9
        @test iters < 20000
        # Boundary points are fixed
        @test xr[1, :] == x[1, :] && xr[end, :] == x[end, :]
        @test yr[:, 1] == y[:, 1] && yr[:, end] == y[:, end]
        # Left-wall forcing pulls the first-cell spacing towards the corner spacing h
        # (the ξ-wall metric coefficients were previously swapped, giving spacings up to ~0.033)
        spacing = [hypot(xr[2, j] - xr[1, j], yr[2, j] - yr[1, j]) for j in 2:Nj-1]
        @test maximum(spacing) < 0.028
        # Grid stays unfolded
        J = [(xr[i+1, j] - xr[i, j]) * (yr[i, j+1] - yr[i, j]) - (xr[i, j+1] - xr[i, j]) * (yr[i+1, j] - yr[i, j])
             for i in 1:Ni-1, j in 1:Nj-1]
        @test minimum(J) > 0
    end

    @testset "SmoothBlocks parameter forms" begin
        N = 6
        grid = TFI([[range(0, 1, length=N) ones(N)], [ones(N) range(0, 1, length=N)],
                    [range(0, 1, length=N) zeros(N)], [zeros(N) range(0, 1, length=N)]])
        blocks = [grid, grid]
        p = EllipticParams(useBottomWall=false)

        single = G.SmoothBlocks(blocks; params=p)
        perblock = G.SmoothBlocks(blocks; params=[p, EllipticParams(skipBlock=true)])
        @test length(single[1]) == 2
        @test perblock[3][2] == 0   # skipped block reports zero iterations
        @test_throws ArgumentError G.SmoothBlocks(blocks; params=[p])
        @test_throws ArgumentError G.SmoothBlocks(blocks; solver=:none, params=p)
    end
end
