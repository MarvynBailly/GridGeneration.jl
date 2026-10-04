# Quarter annulus 1 ≤ r ≤ 3: an isotropic layer on the inner wall plus a circular
# "front" that cuts diagonally across the blocks.

function annulus_case()
    arc(R) = t -> (R * cos(π / 2 * (1 - t)), R * sin(π / 2 * (1 - t)))   # from (0, R) to (R, 0)
    initial = block_from_curves(
        arc(1.0),                        # inner wall
        arc(3.0),                        # outer boundary
        line((0.0, 1.0), (0.0, 3.0)),
        line((1.0, 0.0), (3.0, 0.0));
        Ni = 241, Nj = 161)

    M = metric_sum(
        background(300.0),
        wall_layer((x, y) -> hypot(x, y) - 1.0, 60_000.0, 0.04),   # inner-wall layer
        front((3.2, 3.2), 2.6, 4_000.0, 0.08),                     # curved front
    )

    params = SimParams(splitLocations = [[61, 121, 181], [9, 61, 111]], useSmoothing = false)
    return (; input = initial, bndInfo = Any[], interInfo = Any[], M, params,
              zoom = ((0.6, 1.6), (0.6, 1.6)))
end
