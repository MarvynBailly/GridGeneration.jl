# Channel with a Gaussian bump on the lower wall.
# Anisotropic wall layers (only M22) cluster cells across both walls without
# over-refining along them; a hotspot refines the flow over the bump crest.

bump(x) = 0.12 * exp(-((x - 2.0) / 0.35)^2)

function bump_channel_case()
    L, H = 4.0, 1.0
    # The initial grid only describes the geometry and where the metric is sampled,
    # so make it fine (the generated grid chooses its own resolution).
    initial = block_from_curves(
        t -> (L * t, bump(L * t)),      # bottom wall with the bump
        t -> (L * t, H),                # top wall
        line((0.0, 0.0), (0.0, H)),     # inlet
        line((L, 0.0), (L, H));         # outlet
        Ni = 321, Nj = 201)

    # distance to the nearest wall (vertical distance is accurate enough for a gentle bump)
    wall_dist(x, y) = min(y - bump(x), H - y)

    M = metric_sum(
        background(400.0),                                    # spacing ≈ 0.05 away from walls
        wall_layer(wall_dist, 400_000.0, 0.012; normal = :y), # wall-normal clustering only
        hotspot((2.0, 0.22), 3_000.0, 0.3),                   # refine over the bump
    )

    params = SimParams(
        splitLocations = [[81, 161, 241], [16, 186]],   # thin wall blocks hold the layers
        boundarySolver = :analytic,
        useSmoothing = false,      # keep the metric-driven clustering from the edge solve + TFI
    )
    return (; input = initial, bndInfo = Any[], interInfo = Any[], M, params,
              zoom = ((1.5, 2.5), (0.0, 0.3)))
end
