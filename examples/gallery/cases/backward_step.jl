# Backward-facing step as a genuine multi-block input (3 blocks with interfaces).
# Splits requested on two blocks propagate across the interfaces; the metric refines
# the shear layer leaving the step corner, the corner itself, and the walls.

function backward_step_case()
    h = 0.5
    rect(x0, x1, y0, y1, ni, nj) = block_from_curves(
        line((x0, y0), (x1, y0)), line((x0, y1), (x1, y1)),
        line((x0, y0), (x0, y1)), line((x1, y0), (x1, y1)); Ni = ni, Nj = nj)
    inlet = rect(-1.0, 0.0, h, 1.0, 41, 41)     # block 1: channel upstream of the step
    upper = rect(0.0, 4.0, h, 1.0, 161, 41)     # block 2: downstream, above the step height
    lower = rect(0.0, 4.0, 0.0, h, 161, 41)     # block 3: downstream, below (behind the step)
    blocks = [inlet, upper, lower]

    face(b, s, e) = Dict{String,Any}("block" => b, "start" => s, "end" => e)
    interInfo = Any[
        Dict{String,Any}("blockA" => 1, "blockB" => 2, "start_blkA" => [41, 1], "end_blkA" => [41, 41],
                         "start_blkB" => [1, 1], "end_blkB" => [1, 41], "offset" => [0.0, 0.0, 0.0], "angle" => 0.0),
        Dict{String,Any}("blockA" => 2, "blockB" => 3, "start_blkA" => [1, 1], "end_blkA" => [161, 1],
                         "start_blkB" => [1, 41], "end_blkB" => [161, 41], "offset" => [0.0, 0.0, 0.0], "angle" => 0.0),
    ]
    bndInfo = Any[
        Dict{String,Any}("name" => "inflow",  "faces" => Any[face(1, [1, 1], [1, 41])]),
        Dict{String,Any}("name" => "outflow", "faces" => Any[face(2, [161, 1], [161, 41]), face(3, [161, 1], [161, 41])]),
        Dict{String,Any}("name" => "wall",    "faces" => Any[face(1, [1, 1], [41, 1]), face(1, [1, 41], [41, 41]),
                                                             face(2, [1, 41], [161, 41]), face(3, [1, 1], [161, 1]),
                                                             face(3, [1, 1], [1, 41])]),
    ]

    # distances to the horizontal walls (top, upstream floor at y = h, downstream floor at 0)
    # and to the vertical step face (x = 0, y < h)
    floor_dist(x, y) = min(1.0 - y, x < 0 ? y - h : y)
    face_dist(x, y) = (x >= 0 && y < h) ? x : Inf

    # shear layer leaving the corner along y = h, spreading downstream (only across it: M22)
    shear_layer(x, y) = x <= 0 ? (0.0, 0.0) :
        (0.0, 30_000.0 * exp(-((y - h) / (0.015 + 0.03x))^2) / (1 + 2x))

    M = metric_sum(
        background(250.0),
        wall_layer(floor_dist, 15_000.0, 0.02; normal = :y),   # across the horizontal walls only
        wall_layer(face_dist, 15_000.0, 0.02; normal = :x),    # across the step face only
        shear_layer,
        hotspot((0.0, h), 15_000.0, 0.06),                      # step corner
        hotspot((2.6, 0.0), 4_000.0, 0.35),                     # reattachment region
    )

    # split block 2 in i (propagates into block 3) and in j (propagates into block 1),
    # and block 3 near its floor
    splitRequests = [(2, [[21, 61, 101], [21]]), (3, [Int[], [9]])]
    params = SimParams(useSmoothing = false)
    return (; input = blocks, bndInfo, interInfo, M, params, splitRequests,
              zoom = ((-0.3, 0.9), (0.2, 0.8)))
end
