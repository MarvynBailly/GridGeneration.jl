# One domain, one block layout, four metrics: the unit square split into 8×8 blocks,
# adapted to a uniform metric, two hotspots, an oblique "shock" and a circular front.
# Clustering on the block edges is joined across each block by TFI, so straight
# features crossing a block are followed even when they cut it diagonally.

function square_metrics_cases()
    square = block_from_curves(line((0, 0), (1, 0)), line((0, 1), (1, 1)),
                               line((0, 0), (0, 1)), line((1, 0), (1, 1)); Ni = 161, Nj = 161)
    params = SimParams(splitLocations = [collect(21:20:141), collect(21:20:141)], useSmoothing = false)
    metrics = [
        "Uniform"        => background(900.0),
        "Two hotspots"   => metric_sum(background(400.0), hotspot((0.3, 0.3), 30_000.0, 0.1),
                                       hotspot((0.7, 0.68), 30_000.0, 0.1)),
        "Oblique shock"  => metric_sum(background(400.0), feature_line((0.5, 0.46), (1.0, 0.6), 25_000.0, 0.03)),
        "Circular front" => metric_sum(background(400.0), front((0.5, 0.5), 0.3, 25_000.0, 0.03)),
    ]
    return square, params, metrics
end

function square_metrics_figure(path)
    square, params, metrics = square_metrics_cases()
    kw = (aspect_ratio = :equal, framestyle = :box, grid = false, legend = false,
          xlims = (-0.02, 1.02), ylims = (-0.02, 1.02), titlefontsize = 12, tickfontsize = 7)
    tops = Any[]; bottoms = Any[]; stats = []
    for (name, M) in metrics
        result, s = run_case(square, Any[], Any[], M; params)
        push!(stats, name => s)
        p = plot(; title = name, kw...)
        draw_metric!(p, result.smoothBlocks, M; npix = 300)
        plot!(p; colorbar = false)
        push!(tops, p)
        q = plot(; title = "$(s.cells) cells", kw...)
        draw_blocks!(q, result.smoothBlocks; lw = 0.3, outline_lw = 1.0)
        push!(bottoms, q)
    end
    fig = plot(tops..., bottoms...; layout = (2, 4), size = (1600, 820),
               left_margin = 2Plots.mm, bottom_margin = 2Plots.mm)
    savefig(fig, path)
    return stats
end
