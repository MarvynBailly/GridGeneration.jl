# Wavy channel with a staggered row of hotspots, like the cores of a vortex street.
# The metric reaches the grid through the block edges, so the split lines are placed
# through the cores.

wave(x) = 0.15 * sin(2π * x / 2)

function wavy_channel_case()
    L = 6.0
    initial = block_from_curves(
        t -> (L * t, wave(L * t)),
        t -> (L * t, 1.0 + wave(L * t)),
        line((0.0, 0.0), (0.0, 1.0)),
        line((L, wave(L)), (L, 1.0 + wave(L)));
        Ni = 361, Nj = 121)

    # cores at x = 1, ..., 5, alternating above and below the centreline
    cores = [(1.0 + k, 0.5 + (isodd(k) ? 0.18 : -0.18)) for k in 0:4]
    M = metric_sum(background(200.0), (hotspot(c, 20_000.0, 0.16) for c in cores)...)

    # i-splits at every core (x = 1..5) and half-way between; j-splits through both rows of cores
    params = SimParams(splitLocations = [collect(31:30:331), [39, 83]], useSmoothing = false)
    return (; input = initial, bndInfo = Any[], interInfo = Any[], M, params,
              zoom = ((1.5, 3.5), (-0.2, 1.2)))
end
