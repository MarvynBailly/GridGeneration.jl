# Elliptic smoothing on a distorted domain with curved, skewed sides.
# TFI alone carries the skew of the boundaries into the interior; elliptic smoothing
# straightens and evens out the grid lines, and wall forcing on the curved bottom wall
# additionally makes the grid meet that wall at right angles.

function smoothing_cases()
    initial = block_from_curves(
        t -> (t, 0.25 * sin(π * t)),                   # curved bottom wall
        t -> (-0.3 + 1.5t, 1.0 + 0.15 * sin(2π * t)),  # wavy top
        line((0.0, 0.0), (-0.3, 1.0)),                  # skewed left side
        line((1.0, 0.0), (1.2, 1.0));                   # skewed right side
        Ni = 121, Nj = 121)
    M = metric_sum(background(900.0),
                   wall_layer((x, y) -> y - 0.25 * sin(π * clamp(x, 0, 1)), 20_000.0, 0.06))
    base = (; splitLocations = Vector{Vector{Int}}(), useSplitting = false)
    variants = [
        "Edge solve + TFI" => SimParams(; base..., useSmoothing = false),
        "Elliptic smoothing" => SimParams(; base...,
            elliptic = EllipticParams(sweep = :line, ω = 1.3, useBottomWall = false, max_iter = 20_000, tol = 1e-8)),
        "Smoothing + bottom-wall forcing" => SimParams(; base...,
            elliptic = EllipticParams(useBottomWall = true, max_iter = 20_000, tol = 1e-8)),
    ]
    return initial, M, variants
end

# angle between the grid line leaving the bottom wall and the wall itself, averaged
function wall_angle_deviation(b)
    x, y = b[1, :, :], b[2, :, :]
    devs = map(2:size(x, 1)-1) do i
        t = (x[i+1, 1] - x[i-1, 1], y[i+1, 1] - y[i-1, 1])
        n = (x[i, 2] - x[i, 1], y[i, 2] - y[i, 1])
        abs(acosd(clamp((t[1] * n[1] + t[2] * n[2]) / (hypot(t...) * hypot(n...)), -1, 1)) - 90)
    end
    return sum(devs) / length(devs)
end

function smoothing_figure(path)
    initial, M, variants = smoothing_cases()
    kw = (aspect_ratio = :equal, framestyle = :box, grid = false, legend = false,
          xlims = (-0.35, 1.25), ylims = (-0.05, 1.2), titlefontsize = 11, tickfontsize = 7)
    panels = Any[]; stats = []
    for (name, params) in variants
        result, s = run_case(initial, Any[], Any[], M; params)
        wall = wall_angle_deviation(result.smoothBlocks[1])
        push!(stats, name => merge(s, (; wall)))
        p = plot(; title = "$name\nmean angle dev $(round(s.meandev, digits=1))°, at wall $(round(wall, digits=1))°", kw...)
        draw_blocks!(p, result.smoothBlocks; lw = 0.4, outline_lw = 1.2)
        push!(panels, p)
    end
    fig = plot(panels...; layout = (1, 3), size = (1500, 560), left_margin = 3Plots.mm, bottom_margin = 3Plots.mm)
    savefig(fig, path)
    return stats
end
