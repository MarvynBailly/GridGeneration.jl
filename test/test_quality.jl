using Test
using GridGeneration

@testset "Grid quality" begin
    N = 6
    square = TFI([[range(0, 1, length=N) ones(N)], [ones(N) range(0, 1, length=N)],
                  [range(0, 1, length=N) zeros(N)], [zeros(N) range(0, 1, length=N)]])

    # sheared grid whose grid lines meet at 60°
    sheared = similar(square)
    for j in 1:N, i in 1:N
        ξ, η = (i - 1) / (N - 1), (j - 1) / (N - 1)
        sheared[1, i, j] = ξ + η * cosd(60)
        sheared[2, i, j] = η * sind(60)
    end

    devs, blk, maxdev, mi, mj = ComputeAngleDeviation([square, sheared])
    @test length(devs) == 2 && size(devs[1]) == (N, N)
    @test maximum(devs[1]) ≈ 0 atol=1e-10
    @test all(isapprox.(devs[2][1:end-1, 1:end-1], 30.0; atol=1e-10))
    @test blk == 2
    @test maxdev ≈ 30.0
    @test devs[blk][mi, mj] == maxdev

    # a collapsed edge gives NaN rather than an error
    degenerate = copy(square)
    degenerate[:, 2, 1] = degenerate[:, 1, 1]
    d, _, m, _, _ = ComputeAngleDeviation([degenerate])
    @test isnan(d[1][1, 1])
    @test isfinite(m)
end
